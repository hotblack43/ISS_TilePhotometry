#!/usr/bin/env python3
# Earthshine MVP v0.9.12 — adds --no-dcraw to force rawpy-only; retains skip/overwrite flags.
#
# This code creates aan archive of FITS files made from NEF and CR2 files
#
# Example usage: Example usage:
# uv run   es_mvp_v0_9_12.py --filelist all_iss_files.txt   --root /media/pth/LenovoPS8/ISS_ARCHIVE_TEST   --no-dcraw
#
import argparse, hashlib, json, os, re, shutil, sqlite3, subprocess, sys, fnmatch
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
from tifffile import imread, imwrite
from astropy.io import fits
import exifread

import tarfile, gzip

ARCHIVE_EXTS = ('.gz', '.tgz', '.tar.gz', '.tar')

RAW_EXTS_ALL = {".NEF", ".CR2", ".RAF", ".ARW", ".ORF", ".RW2", ".CR3",
                ".nef", ".cr2", ".raf", ".arw", ".orf", ".rw2", ".cr3"}

def _is_archive(p: Path) -> bool:
    n = p.name.lower()
    return n.endswith('.gz') or n.endswith('.tgz') or n.endswith('.tar.gz') or n.endswith('.tar')

def _decompress_gz_to(raw_gz: Path, out_dir: Path):
    name = raw_gz.name
    if not name.endswith('.gz'):
        return []
    inner = name[:-3]
    if not any(inner.endswith(ext) for ext in RAW_EXTS_ALL):
        return []
    out_path = out_dir / inner
    ensure_dirs(out_dir)
    with gzip.open(raw_gz, 'rb') as fin, open(out_path, 'wb') as fout:
        shutil.copyfileobj(fin, fout)
    ARCHIVE_MAP[str(out_path)] = (str(raw_gz), raw_gz.name)
    return [out_path]

def _extract_tar_members(tar_path: Path, out_dir: Path):
    ensure_dirs(out_dir)
    out = []
    mode = 'r:*'
    with tarfile.open(tar_path, mode) as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            base = Path(m.name).name
            if not any(base.endswith(ext) for ext in RAW_EXTS_ALL):
                continue
            target = out_dir / base
            if target.exists():
                import hashlib as _hashlib
                h = _hashlib.sha256(m.name.encode('utf-8')).hexdigest()[:8]
                target = out_dir / f"{h}_{base}"
            with tf.extractfile(m) as fin, open(target, 'wb') as fout:
                shutil.copyfileobj(fin, fout)
            ARCHIVE_MAP[str(target)] = (str(tar_path), m.name)
            out.append(target)
    return out

ARCHIVE_MAP = {}


try:
    import rawpy
    RAWPY_AVAILABLE = True
except Exception:
    RAWPY_AVAILABLE = False

DEFAULT_TEMPLATE = "{date}_{time}_{orig}_{sha8}_L1_rawcube.fits"

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1<<20), b''):
            h.update(chunk)
    return h.hexdigest()

def run(cmd, check=True):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and p.returncode != 0:
        raise RuntimeError(f"CMD failed: {' '.join(cmd)}\nstdout:\n{p.stdout}\nstderr:\n{p.stderr}")
    return p

def ensure_dirs(*paths):
    for p in paths:
        Path(p).mkdir(parents=True, exist_ok=True)

def has_tool(name: str) -> bool:
    return shutil.which(name) is not None

def now_utc_iso():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

def sanitize(s: str, maxlen=180):
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    if len(s) > maxlen:
        s = s[:maxlen]
    s = s.strip("._")
    return s or "file"

def read_exif(raw_path: Path) -> dict:
    with open(raw_path, 'rb') as f:
        # Parse ALL tags; include MakerNote/unknowns; be lenient
        tags = exifread.process_file(f, stop_tag=None, details=True, strict=False)
    exif = {}
    for k, v in tags.items():
        ks = str(k)  # original key, e.g., "EXIF DateTimeOriginal" or "MakerNote SerialNumber"
        vs = str(v)
        # 1) original
        exif[ks] = vs
        # 2) space->colon normalisation used elsewhere in this script
        exif[ks.replace(' ', ':')] = vs
        # 3) lowercased variants to enable case-insensitive lookups
        exif[ks.lower()] = vs
        exif[ks.lower().replace(' ', ':')] = vs
    return exif

def parse_exif_dt_parts(exif: dict):
    candidates = [
        "Composite:SubSecDateTimeOriginal",
        "EXIF:SubSecDateTimeOriginal",
        "EXIF:DateTimeOriginal",
        "EXIF:CreateDate",
        "QuickTime:CreateDate",
        "DateTimeOriginal",
        "CreateDate",
    ]
    s = None
    for k in candidates:
        v = exif.get(k)
        if v:
            s = str(v).strip()
            break
    if not s:
        return "00000000", "000000"
    s = s.replace('T', ' ').replace('-', ':').replace('Z', '').strip()
    parts = s.split()
    if len(parts) == 1:
        date_part, time_part = parts[0], "00:00:00"
    else:
        date_part, time_part = parts[0], parts[1]
    d = date_part.split(':') + ["01","01"]
    t = time_part.split(':') + ["00","00","00"]
    y, m, da = d[0], d[1], d[2]
    hh, mm, ss = t[0], t[1], t[2].split('.')[0]
    return f"{y}{m}{da}", f"{hh}{mm}{ss}"

def dateobs_from_exif_or_mtime(exif: dict, raw_path: Path) -> str:
    date, time = parse_exif_dt_parts(exif)
    if date != "00000000":
        return f"{date[0:4]}-{date[4:6]}-{date[6:8]}T{time[0:2]}:{time[2:4]}:{time[4:6]}.000"
    dt = datetime.fromtimestamp(raw_path.stat().st_mtime, tz=timezone.utc)
    return dt.strftime("%Y-%m-%dT%H:%M:%S.000")

def _pick_first(exif: dict, keys):
    for k in keys:
        v = exif.get(k)
        if v and str(v).strip():
            return str(v).strip()
    return ""

def _floatish(s):
    s = str(s).strip()
    try:
        if "/" in s:
            a,b = s.split("/",1)
            return float(a)/float(b)
        return float(s)
    except Exception:
        m = re.search(r"[-+]?\d+(\.\d+)?", s)
        if m:
            try: return float(m.group(0))
            except Exception: pass
    return None

def extract_lens_info(exif: dict):
    out = {}
    lens_model = _pick_first(exif, [
        "EXIF:LensModel",
        "EXIF:LensMake",
        "MakerNotes:LensModel",
        "MakerNotes:Lens",
        "MakerNotes:LensType",
        "Nikon:Lens",
        "Nikon3:Lens",
        "Composite:LensID",
        "EXIF:LensType",
        "Image:LensModel",
    ])
    if lens_model:
        out["TELESCOP"] = lens_model
    lens_id = _pick_first(exif, ["Composite:LensID","MakerNotes:LensID","EXIF:LensID"])
    if lens_id:
        out["LENSID"] = lens_id
    fl = _pick_first(exif, ["EXIF:FocalLength","Image:FocalLength"])
    flv = _floatish(fl) if fl else None
    if flv is not None:
        out["FOCALLEN"] = float(flv)
    spec = _pick_first(exif, ["EXIF:LensSpecification","Image:LensSpecification"])
    if spec:
        nums = re.findall(r"[-+]?\d+(?:\.\d+)?", spec.replace(",", " "))
        try:
            nums_f = [float(n) for n in nums]
            if len(nums_f) >= 2:
                out["LNSMIN"] = nums_f[0]
                out["LNSMAX"] = nums_f[1]
            if len(nums_f) >= 4:
                out["APMIN"]  = nums_f[2]
                out["APMAX"]  = nums_f[3]
        except Exception:
            pass
    return out


def extract_body_serial(exif: dict) -> str:
    """
    Extract camera body serial number from exifread-derived tags only (no exiftool).
    We check many plausible keys (case-insensitive) and then fall back to a generic scan
    of any tag name containing 'serial' that yields a plausible serial-like value.
    """
    import re as _re
    candidates = [
        # Common EXIF/IFD names (space->colon and lowercase forms too)
        "EXIF:BodySerialNumber", "EXIF:CameraSerialNumber", "EXIF:SerialNumber",
        "Image:BodySerialNumber", "Image:SerialNumber",
        "ExifIFD:BodySerialNumber", "IFD0:SerialNumber", "TIFF:SerialNumber",
        "exif:bodyserialnumber", "exif:cameraserialnumber", "exif:serialnumber",
        "image:bodyserialnumber", "image:serialnumber",
        "exififd:bodyserialnumber", "ifd0:serialnumber", "tiff:serialnumber",
        # MakerNote variants
        "MakerNote:SerialNumber", "MakerNote:CameraSerialNumber", "MakerNote:InternalSerialNumber",
        "makernote:serialnumber", "makernote:cameraserialnumber", "makernote:internalserialnumber",
        # Brand hints (rare with exifread but harmless)
        "Nikon:SerialNumber", "Nikon3:SerialNumber",
        "Canon:SerialNumber", "Sony:SerialNumber", "Fujifilm:SerialNumber",
        "Pentax:SerialNumber", "Olympus:SerialNumber",
        "nikon:serialnumber", "nikon3:serialnumber",
        "canon:serialnumber", "sony:serialnumber", "fujifilm:serialnumber",
        "pentax:serialnumber", "olympus:serialnumber",
        # Very generic
        "Serial:Number", "SerialNumber", "serial:number", "serialnumber",
    ]
    # Direct lookups first (honour our exif dict which includes original/lower/colonised keys)
    for key in candidates:
        if key in exif and str(exif[key]).strip():
            val = str(exif[key]).strip()
            m = _re.search(r'(?:^|[^0-9])([0-9]{5,12})(?:[^0-9]|$)', val)
            return m.group(1) if m else val
    # Case-insensitive sweep across all keys for anything with 'serial' in the name
    for k, v in exif.items():
        if 'serial' in str(k).lower():
            s = str(v).strip()
            if not s:
                continue
            m = _re.search(r'(?:^|[^0-9])([0-9]{5,12})(?:[^0-9]|$)', s)
            if m:
                return m.group(1)
            if 3 <= len(s) <= 32:
                return s
    return ""


def render_name(template: str, raw_path: Path, exif: dict, sha: str):
    orig = raw_path.stem
    ext  = raw_path.suffix.lstrip(".")
    date, time = parse_exif_dt_parts(exif)
    model = str(exif.get("EXIF:Model") or exif.get("Image:Model") or "").strip().replace(" ", "_") or "Unknown"
    subs = {"orig":orig,"ext":ext,"sha8":sha[:8],"sha16":sha[:16],"sha256":sha,"date":date,"time":time,"model":model}
    try:
        name = template.format(**subs)
    except KeyError as e:
        raise RuntimeError(f"Unknown key {e} in --name-template. Allowed: {', '.join(subs.keys())}")
    name = sanitize(name)
    if not name.lower().endswith(".fits"):
        name += ".fits"
    return name

def pick_dcraw_binary():
    for cand in ("dcraw_emu", "dcraw"):
        p = shutil.which(cand)
        if p: return cand
    return None

def decode_raw_to_tiff_dcraw(raw_path: Path, out_tiff: Path, binary: str):
    cmd = [binary, "-4", "-D", "-W", "-o", "0", "-g", "1", "1", "-j", "-T", str(raw_path)]
    run(cmd)
    cand = raw_path.with_suffix(".tiff")
    if not cand.exists():
        cand = raw_path.with_suffix(".tif")
    if not cand.exists():
        raise RuntimeError(f"{binary} did not produce TIFF")
    shutil.move(str(cand), str(out_tiff))

def cfa_pattern_from_rawpy(rp) -> str:
    try:
        pat = rp.raw_pattern
        desc = rp.color_desc.decode() if isinstance(rp.color_desc, bytes) else rp.color_desc
        s = ''.join([desc[pat[0,0]], desc[pat[0,1]], desc[pat[1,0]], desc[pat[1,1]]]).upper()
        return s if s in ("RGGB","BGGR","GBRG","GRBG") else "RGGB"
    except Exception:
        return "RGGB"

def decode_raw_to_tiff_rawpy(raw_path: Path, out_tiff: Path):
    if not RAWPY_AVAILABLE:
        raise RuntimeError("rawpy not available (pip install rawpy)")
    with rawpy.imread(str(raw_path)) as rp:
        arr = rp.raw_image_visible.astype(np.uint16) if rp.raw_image_visible is not None else rp.raw_image.astype(np.uint16)
        imwrite(out_tiff, arr, photometric='minisblack')
        white_level = None; black_levels = None; bayerpat_hint = None
        try: white_level = int(rp.white_level)
        except Exception: pass
        try:
            bl = rp.black_level_per_channel
            if bl is not None:
                bl = list(bl)
                if len(bl)==1: bl = bl*4
                elif len(bl)==2: bl = [bl[0], bl[1], bl[1], bl[0]]
                black_levels = [int(b) for b in bl[:4]]
        except Exception: pass
        try: bayerpat_hint = cfa_pattern_from_rawpy(rp)
        except Exception: pass
    return white_level, black_levels, bayerpat_hint

def guess_cfa_pattern_from_exif(exif: dict, default="RGGB"):
    for k in ("MakerNotes:CFAPattern", "EXIF:CFAPattern", "CFAPattern", "CFAPattern2"):
        v = exif.get(k)
        if not v: continue
        s = str(v).upper().replace(" ", "").replace(",", "")
        if "RGGB" in s: return "RGGB"
        if "BGGR" in s: return "BGGR"
        if "GBRG" in s: return "GBRG"
        if "GRBG" in s: return "GRBG"
    return default

def split_cfa_planes(data_f32, pattern: str):
    pat = pattern.upper()
    if pat == "RGGB":
        r  = data_f32[0::2,0::2]; g1 = data_f32[0::2,1::2]
        g2 = data_f32[1::2,0::2]; b  = data_f32[1::2,1::2]
    elif pat == "BGGR":
        b  = data_f32[0::2,0::2]; g1 = data_f32[0::2,1::2]
        g2 = data_f32[1::2,0::2]; r  = data_f32[1::2,1::2]
    elif pat == "GBRG":
        g1 = data_f32[0::2,0::2]; b  = data_f32[0::2,1::2]
        r  = data_f32[1::2,0::2]; g2 = data_f32[1::2,1::2]
    elif pat == "GRBG":
        g1 = data_f32[0::2,0::2]; r  = data_f32[0::2,1::2]
        b  = data_f32[1::2,0::2]; g2 = data_f32[1::2,1::2]
    else:
        r  = data_f32[0::2,0::2]; g1 = data_f32[0::2,1::2]
        g2 = data_f32[1::2,0::2]; b  = data_f32[1::2,1::2]; pat = "RGGB"
    return r, g1, g2, b, pat

def write_fits_cfa_cube(tiff_path: Path, exif: dict, raw_path: Path, out_fits: Path, procver: str, extras: dict,
                        adc_bits_hint=None, black_levels=None, white_level=None, bayerpat_hint=None):
    arr_u16 = imread(tiff_path)
    if arr_u16.dtype != 'uint16':
        raise RuntimeError(f"Expected uint16 TIFF, got {arr_u16.dtype}")

    adc_bits = adc_bits_hint or 14
    for k in ("EXIF:BitsPerSample","Image:BitsPerSample","MakerNotes:BitsPerSample","BitsPerSample"):
        if k in exif:
            try:
                adc_bits = int(str(exif[k]).split()[0]); break
            except Exception:
                pass
    adc_limit = 4095 if adc_bits == 12 else 16383
    observed_max = int(np.max(arr_u16))
    if observed_max > adc_limit:
        raise RuntimeError(f"Decoder scaled data (max={observed_max} > limit for {adc_bits}-bit = {adc_limit}).")

    data_f32 = arr_u16.astype(np.float32)
    pattern = bayerpat_hint or guess_cfa_pattern_from_exif(exif, default=extras.get("BAYERPAT","RGGB"))
    r, g1, g2, b, pat = split_cfa_planes(data_f32, pattern)
    cube = np.stack([r, g1, g2, b], axis=0)
    cube_u16 = np.clip(np.round(cube), 0, 65535).astype(np.uint16)

    hdr = fits.Header()
    hdr["SIMPLE"]   = True
    hdr["BITPIX"]   = 16
    hdr["NAXIS"]    = 3
    hdr["BUNIT"]    = "adu"
    hdr["CHANNELS"] = "R,G1,G2,B"
    hdr["BAYERPAT"] = pat
    hdr["ADC_BITS"] = adc_bits
    hdr["DATAMAX"]  = observed_max
    hdr["DATAMIN"]  = int(np.min(arr_u16))
    hdr["BSCAL"]    = 1
    hdr["BZERO"]    = 0

    if black_levels is not None and len(black_levels)==4:
        hdr["BLACK_R"]  = int(black_levels[0])
        hdr["BLACK_G1"] = int(black_levels[1])
        hdr["BLACK_G2"] = int(black_levels[2])
        hdr["BLACK_B"]  = int(black_levels[3])
    if white_level is not None:
        hdr["WHITELEV"] = int(white_level)

    hdr["DATE-OBS"] = extras["DATE-OBS"]
    hdr["TIMESYS"]  = "UTC"

    def _floatish_local(s):
        return _floatish(s)
    et = exif.get("EXIF:ExposureTime") or exif.get("EXIF:ExposureBiasValue")
    if et:
        try: hdr["EXPTIME"] = _floatish_local(et)
        except Exception: pass
    fnum = exif.get("EXIF:FNumber")
    if fnum:
        try: hdr["FNUMBER"] = _floatish_local(fnum)
        except Exception: pass
    iso = exif.get("EXIF:ISOSpeedRatings") or exif.get("EXIF:ISO")
    if iso:
        try: hdr["ISOSPEED"] = int(str(iso).split()[0])
        except Exception: pass

    cam = (exif.get("EXIF:Model") or exif.get("Image:Model") or "").strip()
    if cam:
        hdr["INSTRUME"] = cam
    body_ser = extract_body_serial(exif)
    if body_ser:
        hdr["BODYSER"] = body_ser

    lens_meta = extract_lens_info(exif)
    if "TELESCOP" in lens_meta and lens_meta["TELESCOP"]:
        hdr["TELESCOP"] = lens_meta["TELESCOP"]
    if "LENSID" in lens_meta and lens_meta["LENSID"]:
        hdr["LENSID"]   = lens_meta["LENSID"]
    if "FOCALLEN" in lens_meta and lens_meta["FOCALLEN"] is not None:
        try: hdr["FOCALLEN"] = float(lens_meta["FOCALLEN"])
        except Exception: pass
    if "LNSMIN" in lens_meta:
        try: hdr["LNSMIN"]   = float(lens_meta["LNSMIN"])
        except Exception: pass
    if "LNSMAX" in lens_meta:
        try: hdr["LNSMAX"]   = float(lens_meta["LNSMAX"])
        except Exception: pass
    if "APMIN" in lens_meta:
        try: hdr["APMIN"]    = float(lens_meta["APMIN"])
        except Exception: pass
    if "APMAX" in lens_meta:
        try: hdr["APMAX"]    = float(lens_meta["APMAX"])
        except Exception: pass

    filt = extras.get("FILTER")
    if filt:
        hdr["FILTER"] = filt
    plat = extras.get("PLATFORM")
    if plat:
        hdr["PLATFORM"] = plat
    route = extras.get("ROUTE")
    if route:
        hdr["ROUTE"] = route

    hdr["ORIGNAME"] = raw_path.name
    hdr["ORIGHASH"] = "sha256:" + sha256(raw_path)
    hdr["PROCVER"]  = procver
    try:
        from pathlib import Path as _P
        if str(raw_path).startswith(str(extract_root)):
            src_mem = ARCHIVE_MAP.get(str(raw_path))
            if src_mem:
                src, mem = src_mem
                hdr["ARCHIVE"] = _P(src).name[:68]
                hdr["MEMBER"]  = _P(mem).name[:68]
                hdr["ORIGSRC"] = "archive"
    except Exception:
        pass

    for k, v in (extras or {}).items():
        if k in ("DATE-OBS","MJD-OBS","BAYERPAT","INSTRUME","TELESCOP","FILTER","PLATFORM","ROUTE","HISTORY",
                 "LENSID","FOCALLEN","LNSMIN","LNSMAX","APMIN","APMAX", "_OVERWRITE_FITS_"):
            continue
        hdr[k[:8]] = v

    for hline in (extras.get("HISTORY") or []):
        try: hdr.add_history(hline)
        except Exception: pass

    fits.HDUList([fits.PrimaryHDU(data=cube_u16, header=hdr)]).writeto(out_fits, overwrite=extras.get("_OVERWRITE_FITS_", False))

SCHEMA = """
CREATE TABLE IF NOT EXISTS items(
  raw_sha256 TEXT PRIMARY KEY,
  raw_path   TEXT,
  format     TEXT,
  camera     TEXT,
  route      TEXT,
  procver    TEXT,
  l1_path    TEXT,
  status     TEXT,
  first_seen_utc TEXT,
  last_processed_utc TEXT,
  notes      TEXT
);
CREATE INDEX IF NOT EXISTS idx_status ON items(status);
"""

def db_open(path: Path):
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=WAL;")
    for stmt in SCHEMA.strip().split(";"):
        s = stmt.strip()
        if s:
            conn.execute(s)
    return conn

def db_upsert(conn, row: dict):
    conn.execute("""
    INSERT INTO items(raw_sha256,raw_path,format,camera,route,procver,l1_path,status,first_seen_utc,last_processed_utc,notes)
    VALUES(:raw_sha256,:raw_path,:format,:camera,:route,:procver,:l1_path,:status,:first_seen_utc,:last_processed_utc,:notes)
    ON CONFLICT(raw_sha256) DO UPDATE SET
      raw_path=excluded.raw_path,
      format=excluded.format,
      camera=excluded.camera,
      route=excluded.route,
      procver=excluded.procver,
      l1_path=excluded.l1_path,
      status=excluded.status,
      last_processed_utc=excluded.last_processed_utc,
      notes=excluded.notes;
    """, row)
    conn.commit()

def already_done(conn, h: str, procver: str) -> bool:
    cur = conn.execute("SELECT status, procver FROM items WHERE raw_sha256=?", (h,))
    r = cur.fetchone()
    if not r: return False
    status, pv = r
    return status == "OK" and pv == procver

def size_of(p: Path) -> int:
    return p.stat().st_size

def pack_files(l1_root: Path, packs_root: Path, max_bytes=49*1024**3, max_files=100):
    files = sorted([p for p in l1_root.glob("*.fits")], key=lambda p: p.stat().st_mtime)
    if not files:
        return
    pack_idx, cur_bytes, cur_count = 1, 0, 0
    staged = []
    def flush_pack(idx, items):
        if not items: return
        pdir = packs_root / f"pack_{idx:03d}"
        ensure_dirs(pdir)
        for f in items:
            target = pdir / f.name
            if not target.exists():
                os.link(f, target)
    for f in files:
        b = size_of(f)
        if cur_bytes + b > max_bytes or cur_count + 1 > max_files:
            flush_pack(pack_idx, staged)
            pack_idx += 1
            cur_bytes, cur_count = 0, 0
            staged = []
        staged.append(f); cur_bytes += b; cur_count += 1
    flush_pack(pack_idx, staged)

DEFAULT_EXTS = {".NEF", ".CR2", ".RAF", ".ARW", ".ORF", ".RW2", ".CR3"}
CAN_USE_DCRAW = {".NEF", ".CR2", ".RAF", ".ARW", ".ORF", ".RW2"}
NEEDS_VENDOR  = {".CR3"}

def looks_like_nikon_he(raw_path: Path) -> bool:
    return raw_path.suffix.upper()==".NEF" and ("HE" in raw_path.stem.upper() or "Z9" in raw_path.stem.upper())


def load_filelist(path: Path, extract_dir: Path | None = None, accept_archives: bool = True):
    base = path.parent
    items, missing, extracted = [], [], []
    for i, line in enumerate(path.read_text().splitlines(), start=1):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        s = os.path.expanduser(s)
        p = Path(s)
        if not p.is_absolute():
            p = (base / s).resolve()
        if not p.exists() or not p.is_file():
            missing.append((i, s))
            continue
        if accept_archives and _is_archive(p):
            if extract_dir is None:
                extract_dir = (base / "EXTRACT").resolve()
            ensure_dirs(extract_dir)
            outs = []
            name = p.name.lower()
            if name.endswith('.gz') and not name.endswith('.tar.gz') and not name.endswith('.tgz'):
                outs = _decompress_gz_to(p, extract_dir)
            elif name.endswith('.tar.gz') or name.endswith('.tgz') or name.endswith('.tar'):
                outs = _extract_tar_members(p, extract_dir)
            extracted.extend(outs)
            items.extend(outs)
        else:
            items.append(p)
    return items, missing, extracted

def discover_files(roots, globs, exts, follow_symlinks=False):
    seen = set()
    out = []
    def want(p: Path):
        if exts and p.suffix.upper() not in exts:
            return False
        if globs:
            return any(fnmatch.fnmatch(p.name, pat) for pat in globs)
        return True
    for r in roots:
        r = Path(r)
        if not r.exists(): continue
        for dirpath, dirnames, filenames in os.walk(r, followlinks=follow_symlinks):
            for name in filenames:
                p = Path(dirpath) / name
                if want(p):
                    key = (p.stat().st_ino, p.stat().st_dev) if hasattr(os, "stat") else (str(p),)
                    if key in seen: continue
                    seen.add(key); out.append(p)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", action="append")
    ap.add_argument("--filelist")
    ap.add_argument("--glob", action="append")
    ap.add_argument("--exts")
    ap.add_argument("--follow-symlinks", action="store_true")
    ap.add_argument("--root",  required=True)
    ap.add_argument("--procver", default=os.environ.get("PROCVER","earthshine-pipeline v0.9.12"))
    ap.add_argument("--max-pack-gb", type=float, default=49.0)
    ap.add_argument("--max-pack-files", type=int, default=100)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--accept-archives", action="store_true", default=True)
    ap.add_argument("--keep-extracted", action="store_true")
    ap.add_argument("--extract-dir", default=None)
    ap.add_argument("--skip-if-fits-exists", action="store_true",
                    help="If the target L1 FITS already exists, skip processing regardless of manifest DB.")
    ap.add_argument("--overwrite-fits", action="store_true",
                    help="Allow overwriting existing L1 FITS files (default is to fail if exists).")
    ap.add_argument("--no-dcraw", action="store_true",
                    help="Force RAW decoding via rawpy only (ignore dcraw_emu even if present).")
    ap.add_argument("--name-template", default=DEFAULT_TEMPLATE,
                    help="Filename pattern; keys: {orig},{ext},{sha8},{sha16},{sha256},{date},{time},{model}")
    args = ap.parse_args()

    root   = Path(args.root)
    l1     = root/"L1"
    packs  = root/"PACKS"
    man    = root/"MANIFEST"/"earthshine_manifest.sqlite"
    logs   = root/"LOGS"
    ensure_dirs(l1, packs, man.parent, logs)

    dcraw_bin = pick_dcraw_binary()
    have_dcraw = dcraw_bin is not None
    if args.no_dcraw:
        have_dcraw = False
        dcraw_bin = None

    candidates = []; missing_report = []

    extract_root = Path(args.extract_dir).resolve() if args.extract_dir else (root / "EXTRACT")
    extracted_paths_all = []
    if args.filelist:
        fl = Path(args.filelist)
        if not fl.exists():
            print(f"ERROR: filelist not found: {fl}", file=sys.stderr); sys.exit(2)
        got, missing, extracted = load_filelist(fl, extract_dir=extract_root, accept_archives=args.accept_archives)
        extracted_paths_all.extend(extracted)
        candidates.extend(got); missing_report = missing

    exts = set(e.strip().upper() for e in args.exts.split(",")) if args.exts else DEFAULT_EXTS
    globs = args.glob or []
    if args.input:
        candidates.extend(discover_files(args.input, globs, exts, follow_symlinks=args.follow_symlinks))

    if not candidates and missing_report:
        print("ERROR: filelist contained no existing files. First few missing entries:", file=sys.stderr)
        for i,(lineno, s) in enumerate(missing_report[:10], start=1):
            print(f"  line {lineno}: {s}", file=sys.stderr)
        sys.exit(2)

    if not candidates:
        print("ERROR: No inputs. Provide --input and/or --filelist.", file=sys.stderr); sys.exit(2)

    candidates = sorted(set(Path(p) for p in candidates if Path(p).is_file()), key=lambda p: p.name)

    if args.dry_run:
        n_total = len(candidates)
        n_dcraw = sum(1 for p in candidates if have_dcraw and p.suffix.upper() in CAN_USE_DCRAW and not looks_like_nikon_he(p))
        n_rawpy = n_total - n_dcraw
        print(f"[DRY-RUN] Candidates: {n_total} (dcraw route: {n_dcraw}, rawpy route: {n_rawpy})")
        if missing_report:
            print(f"[DRY-RUN] Missing entries in filelist: {len(missing_report)} (showing up to 10):")
            for i,(lineno, s) in enumerate(missing_report[:10], start=1):
                print(f"  line {lineno}: {s}")
        return

    conn = db_open(man)
    procver = args.procver
    max_bytes = int(args.max_pack_gb * (1024**3))
    max_files = args.max_pack_files

    for raw in candidates:
        ext = raw.suffix.upper()
        h = sha256(raw)
        if already_done(conn, h, procver):
            print(f"SKIP (cached {procver}): {raw}"); continue

        fmt = ext.lstrip(".")
        camera = ""; route = ""; status = "ERROR"; notes = ""; l1_path = ""

        try:
            exif = read_exif(raw)
            camera = str(exif.get("EXIF:Model") or exif.get("Image:Model") or "")

            interm_tiff = (root/"LOGS"/f"{h}.mosaic.tiff")
            ensure_dirs(interm_tiff.parent)

            extras = {}
            _arc = ARCHIVE_MAP.get(str(raw))
            if _arc:
                extras["_ARCHIVE_SRC_"] = _arc[0]
                extras["_ARCHIVE_MEMBER_"] = _arc[1]
            extras["DATE-OBS"] = dateobs_from_exif_or_mtime(exif, raw)
            if raw.name.lower().startswith("iss"): extras["PLATFORM"] = "ISS"
            extras["FILTER"] = "RGB"

            adc_bits_hint = None; black_levels=None; white_level=None; bayerpat_hint=None

            if have_dcraw and ext in CAN_USE_DCRAW and not looks_like_nikon_he(raw):
                decode_raw_to_tiff_dcraw(raw, raw.with_suffix(".mosaic.tiff"), dcraw_bin)
                src = raw.with_suffix(".mosaic.tiff")
                shutil.move(str(src), str(interm_tiff))
                route = "dcraw"; extras["ROUTE"] = "dcraw"
                extras["HISTORY"] = [f"RAW->TIFF via {dcraw_bin} (-4 -D -W -o 0 -g 1 1 -j -T)."]
            else:
                wl, bl, bp = decode_raw_to_tiff_rawpy(raw, interm_tiff)
                white_level, black_levels, bayerpat_hint = wl, bl, bp
                route = "rawpy"; extras["ROUTE"] = "rawpy"
                extras["HISTORY"] = ["RAW read via rawpy (LibRaw); mosaic written to TIFF (no demosaic/scale)."]

            out_name = render_name(args.name_template, raw, exif, h)
            out_fits = (Path(args.root)/"L1") / out_name
            if args.skip_if_fits_exists and out_fits.exists():
                print(f"SKIP (exists): {raw} -> {out_fits}")
                try: interm_tiff.unlink(missing_ok=True)
                except Exception: pass
                continue

            extras["_OVERWRITE_FITS_"] = bool(args.overwrite_fits)

            write_fits_cfa_cube(interm_tiff, exif, raw, out_fits, procver, extras,
                                adc_bits_hint=adc_bits_hint,
                                black_levels=black_levels,
                                white_level=white_level,
                                bayerpat_hint=bayerpat_hint)
            l1_path = str(out_fits); status = "OK"; notes = ""

            try: interm_tiff.unlink(missing_ok=True)
            except Exception: pass

        except Exception as e:
            status = "ERROR"; notes = f"{type(e).__name__}: {e}"

        row = {"raw_sha256": h, "raw_path": str(raw), "format": fmt, "camera": camera,
               "route": route, "procver": procver, "l1_path": l1_path,
               "status": status, "first_seen_utc": now_utc_iso(),
               "last_processed_utc": now_utc_iso(), "notes": notes}
        db_upsert(conn, row)

        msg = f"{status}: {raw}"
        if status != "OK" and notes: msg += f" -> {notes}"
        print(msg)

    pack_files(Path(args.root)/"L1", Path(args.root)/"PACKS",
               max_bytes=int(args.max_pack_gb*(1024**3)), max_files=args.max_pack_files)
    print("Done.")
    try:
        if not args.keep_extracted and extract_root.exists():
            shutil.rmtree(extract_root)
    except Exception:
        pass

import sqlite3
def db_open(path: Path):
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=WAL;")
    for stmt in SCHEMA.strip().split(";"):
        s = stmt.strip()
        if s: conn.execute(s)
    return conn

if __name__ == "__main__":
    main()
