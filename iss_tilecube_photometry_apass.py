#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Peter Thejll, DMI November 2025

iss_tilecube_photometry_apass.py

Pipeline step:
- Read WCS'd 4-plane tile cubes: *_tNN_cube_rg1g2b.fits  (shape: 4 x ny x nx)
  Planes: [0:R, 1:G1, 2:G2, 3:B]
- DETECT stars on (G1 + G2)/2 to avoid R saturation and low B S/N.
- Measure aperture photometry on R, G1, G2, B at the same centroids.
- Convert (x,y) -> (RA,Dec) using 2-D celestial WCS (TAN/TAN-SIP).
- Query APASS DR9 (VizieR table II/336/apass9) around the tile center
  with radius = min(half-diagonal, 2.0 deg) to keep queries stable.
- Cross-match detections to APASS within --match arcsec.
- Write per-tile CSVs + merged field CSV.

Requires: numpy, pandas, astropy, sep, astroquery

Example usage:
    uv run iss_tilecube_photometry_apass.py TESTOUT4 TESTOUT4 --maxrad 2.5 --match 60
"""

import os, re, glob, math, argparse, logging, time
import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

import sep
from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = 50000

# --------------------
# Logging
# --------------------
logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")
LOG = logging.getLogger("iss_tilecube_phot_apass")

# --------------------
# Defaults (CLI overridable)
# --------------------
DEF_APRADIUS       = 4.0     # aperture radius [pix]
DEF_BKG_INNER      = 12.0    # annulus inner radius [pix]
DEF_BKG_OUTER      = 20.0    # annulus outer radius [pix]
DEF_DETECT_THRESH  = 5.0     # detection threshold [sigma]
DEF_DETECT_MINAREA = 5       # min connected pixels > thresh
DEF_DEBLEND_CONT   = 0.005   # SEP deblend contrast
DEF_MATCH_ARCSEC   = 2.5     # match radius [arcsec]
DEF_MAX_APASS_RAD  = 2.0     # cap per-tile search radius [deg]
DEF_ROWS_LIMIT     = 200000  # Vizier row limit
DEF_TIMEOUT        = 120     # Vizier timeout seconds
DEF_RETRIES        = 3       # Vizier retries

# APASS via VizieR
VIZIER_CAT  = "II/336/apass9"
VIZIER_COLS = ["RAJ2000","DEJ2000","Vmag","e_Vmag","Bmag","e_Bmag"]

CUBE_GLOB   = "*_t??_cube_rg1g2b.fits"

# --------------------
# Helpers
# --------------------
def diagnostics_offsets(ra_det, dec_det, apass_df, ap_idx, ap_sep, tag, p=95):
    """
    Print separation and signed offset stats for matches.
    Offsets are in arcsec, with dRA already multiplied by cos(dec).
    """
    import numpy as np
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    m = ~np.isnan(ap_idx)
    if not np.any(m):
        LOG.info(f"{tag}: diagnostics: no matches; cannot compute offsets")
        return

    # Build matched pairs
    det = SkyCoord(ra=np.array(ra_det[m])*u.deg, dec=np.array(dec_det[m])*u.deg, frame="icrs")
    ap  = SkyCoord(ra=apass_df["RAJ2000"].to_numpy(float)[ap_idx[m].astype(int)]*u.deg,
                   dec=apass_df["DEJ2000"].to_numpy(float)[ap_idx[m].astype(int)]*u.deg, frame="icrs")

    # Vector offsets (small-angle)
    dra  = (ap.ra  - det.ra ).to(u.arcsec)
    ddec = (ap.dec - det.dec).to(u.arcsec)
    # Convert dRA to true on-sky arcsec (cos(dec) factor)
    dra_on_sky = (dra * np.cos(det.dec)).value
    ddec_val   = ddec.value
    sep_val    = ap_sep[m]

    def pct(a, p): return float(np.nanpercentile(a, p))

    LOG.info(
        f"{tag}: diagnostics:"
        f" Nmatch={m.sum()}, sep[med/{p}%]={np.nanmedian(sep_val):.2f}/{pct(sep_val,p):.2f}\""
        f", dRA*cos(dec) [med/σ]={np.nanmedian(dra_on_sky):.2f}/{np.nanstd(dra_on_sky):.2f}\""
        f", dDec [med/σ]={np.nanmedian(ddec_val):.2f}/{np.nanstd(ddec_val):.2f}\""
    )


def load_cube(path):
    with fits.open(path) as hdul:
        data = hdul[0].data
        hdr  = hdul[0].header.copy()
    if data is None or data.ndim != 3 or data.shape[0] != 4:
        raise ValueError(f"Expected cube shape (4, ny, nx); got {None if data is None else data.shape} for {path}")
    return data.astype(np.float32, copy=False), hdr

def tile_tag(path):
    m = re.search(r"_t(\d{2})_cube_rg1g2b\.fits$", path)
    return f"t{m.group(1)}" if m else "t??"

def stem_from_cube(path):
    return re.sub(r"_t\d{2}_cube_rg1g2b\.fits$", "", path)

def tile_center_radius_deg(wcs2, ny, nx):
    # Center in 1-based FITS convention:
    x = nx / 2.0; y = ny / 2.0
    ra_c, dec_c = wcs2.wcs_pix2world([[x+1, y+1]], 1)[0]
    # Estimate width/height from pixel scale matrix (deg/pix)
    cd = wcs2.pixel_scale_matrix
    sx = float(np.hypot(cd[0,0], cd[1,0]))
    sy = float(np.hypot(cd[0,1], cd[1,1]))
    width_deg, height_deg = sx * nx, sy * ny
    half_diag = 0.5 * float(np.hypot(width_deg, height_deg))
    return float(ra_c), float(dec_c), float(half_diag)

def detect_sources(detect_img, thresh, minarea, deblend):
    arr = np.ascontiguousarray(detect_img)
    bkg = sep.Background(arr)
    arr_sub = arr - bkg.back()
    objs = sep.extract(arr_sub, thresh=thresh, err=bkg.rms(),
                       minarea=minarea, deblend_cont=deblend)
    return objs, bkg

def aperture_phot(img, x, y, r, rin, rout):
    a = np.ascontiguousarray(img)
    flux, fluxerr, flag = sep.sum_circle(a, x, y, r, subpix=1)

    # Background from annulus:
    try:
        bkg_ann, bkgerr, bkgflag = sep.sum_circann(a, x, y, rin, rout, subpix=1)
    except AttributeError:
        # Older SEP fallback: annulus = diff of two circles
        b_outer, _, _ = sep.sum_circle(a, x, y, rout, subpix=1)
        b_inner, _, _ = sep.sum_circle(a, x, y, rin,  subpix=1)
        bkg_ann = b_outer - b_inner

    area_ap  = math.pi * (r**2)
    area_ann = math.pi * (rout**2 - rin**2)
    bkg_per_pix = bkg_ann / area_ann
    flux_net = flux - bkg_per_pix * area_ap

    flux_net = np.where(flux_net > 0, flux_net, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        mag_inst = -2.5 * np.log10(flux_net)
    return flux_net, mag_inst, flag

def query_apass_cone(ra_deg, dec_deg, radius_deg, retries=3):
    from astroquery.vizier import Vizier
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import pandas as pd
    import numpy as np

    Vizier.ROW_LIMIT = -1
    cols = [
        "RAJ2000","DEJ2000",
        "Vmag","e_Vmag","Bmag","e_Bmag",
        "g'mag","e_g'mag","r'mag","e_r'mag","i'mag","e_i'mag"
    ]
    v = Vizier(columns=cols, row_limit=-1)
    center = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame="icrs")

    last_exc = None
    for _ in range(max(1, retries)):
        try:
            res = v.query_region(center, radius=radius_deg*u.deg, catalog="II/336/apass9")
            if len(res) == 0:
                return pd.DataFrame(columns=cols)
            df = res[0].to_pandas()

            # Coerce numeric cols (important for SkyCoord + comparisons)
            for c in ["RAJ2000","DEJ2000","Vmag","e_Vmag","Bmag","e_Bmag",
                      "g'mag","e_g'mag","r'mag","e_r'mag","i'mag","e_i'mag"]:
                if c in df.columns:
                    df[c] = pd.to_numeric(df[c], errors="coerce")

            # Optional: trim to V < 14 (adjust if you like)
            if "Vmag" in df.columns:
                df = df[df["Vmag"].notna() & (df["Vmag"] < 14.0)].reset_index(drop=True)

            return df
        except Exception as e:
            last_exc = e
    # If we fall through retries
    LOG.warning(f"APASS query failed after {retries} tries: {last_exc}")
    return pd.DataFrame(columns=cols)

def badly_query_apass_cone(ra_deg, dec_deg, radius_deg, retries=DEF_RETRIES):
    Vizier.ROW_LIMIT = DEF_ROWS_LIMIT
    Vizier.TIMEOUT = DEF_TIMEOUT
    v = Vizier(columns=VIZIER_COLS)
    last_err = None
    for attempt in range(1, retries+1):
        try:
            res = v.query_region(
                SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame="icrs"),
                radius=radius_deg * u.deg,
                catalog=VIZIER_CAT
            )
            if len(res) == 0:
                LOG.warning(f"APASS query returned no results at RA={ra_deg}, Dec={dec_deg}")
                return pd.DataFrame(columns=VIZIER_COLS)

            df = res[0].to_pandas()
            if len(df) == 0:
                LOG.warning(f"No matches found for RA={ra_deg}, Dec={dec_deg} in APASS.")
            
            # Debug print for the APASS query results
            LOG.info(f"APASS query result sample:\n{df.head()}")
            
            for c in ["Vmag", "Bmag", "e_Vmag", "e_Bmag", "RAJ2000", "DEJ2000"]:
                if c not in df.columns:
                    df[c] = np.nan
            return df
        except Exception as e:
            last_err = e
            LOG.warning(f"APASS query attempt {attempt}/{retries} failed: {e}")
            time.sleep(1.5 * attempt)
    
    LOG.error(f"APASS/VizieR permanently failed after {retries} attempts: {last_err}")
    return pd.DataFrame(columns=VIZIER_COLS)

def old_query_apass_cone(ra_deg, dec_deg, radius_deg, retries=DEF_RETRIES):
    Vizier.ROW_LIMIT = DEF_ROWS_LIMIT
    Vizier.TIMEOUT = DEF_TIMEOUT
    v = Vizier(columns=VIZIER_COLS)
    last_err = None
    for attempt in range(1, retries+1):
        try:
            res = v.query_region(
                SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame="icrs"),
                radius=radius_deg * u.deg,
                catalog=VIZIER_CAT
            )
            if len(res) == 0:
                return pd.DataFrame(columns=VIZIER_COLS)
            df = res[0].to_pandas()
            for c in ["Vmag","Bmag","e_Vmag","e_Bmag","RAJ2000","DEJ2000"]:
                if c not in df.columns:
                    df[c] = np.nan
            return df
        except Exception as e:
            last_err = e
            LOG.warning(f"APASS query attempt {attempt}/{retries} failed: {e}")
            time.sleep(1.5 * attempt)
    LOG.error(f"APASS/VizieR permanently failed after {retries} attempts: {last_err}")
    return pd.DataFrame(columns=VIZIER_COLS)

def crossmatch_apass(ra_deg, dec_deg, apass_df, match_arcsec):
    """
    ra_deg, dec_deg: numpy arrays of detections (degrees)
    apass_df: DataFrame containing RAJ2000, DEJ2000 (degrees)
    match_arcsec: maximum separation to accept (arcseconds)

    Returns:
      idx: float array same length as detections; index into apass_df or NaN if no match
      sep_arcsec: float array of separations for nearest APASS neighbour
    """
    import numpy as np
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    if apass_df is None or len(apass_df) == 0:
        return np.full(len(ra_deg), np.nan, dtype=float), np.full(len(ra_deg), np.nan, dtype=float)

    # Ensure numeric (in case upstream didn’t coerce)
    ra_ap = apass_df["RAJ2000"].to_numpy(dtype=float)
    de_ap = apass_df["DEJ2000"].to_numpy(dtype=float)

    # Build SkyCoord sets
    det = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame="icrs")
    cat = SkyCoord(ra=ra_ap*u.deg,  dec=de_ap*u.deg,  frame="icrs")

    # Nearest-neighbour match
    idx_cat, sep2d, _ = det.match_to_catalog_sky(cat)
    sep_arcsec = sep2d.to(u.arcsec).value

    # Accept only within threshold
    ok = sep_arcsec <= float(match_arcsec)

    # Return index (float with NaN for no match) + separations
    idx_float = idx_cat.astype(float)
    idx_float[~ok] = np.nan
    return idx_float, sep_arcsec

def no_crossmatch_apass(ra, dec, apass_df, max_sep_arcsec):
    if len(apass_df) > 0:
        import pandas as _pd
        apass_df["RAJ2000"] = _pd.to_numeric(apass_df["RAJ2000"], errors="coerce")
        apass_df["DEJ2000"] = _pd.to_numeric(apass_df["DEJ2000"], errors="coerce")
    
    # Handle the case where there is no valid data
    if len(apass_df.dropna(subset=["RAJ2000", "DEJ2000"])) == 0:
        return np.full(ra.shape, np.nan), np.full(ra.shape, np.nan)
    
    det = SkyCoord(ra*u.deg, dec*u.deg)
    cat = SkyCoord(apass_df["RAJ2000"].values*u.deg, apass_df["DEJ2000"].values*u.deg)

    # Debug print for RA/Dec values of detections and catalog
    LOG.info(f"Detected RA/Dec: {ra}, {dec}")
    LOG.info(f"APASS RA/Dec: {apass_df[['RAJ2000', 'DEJ2000']].head()}")

    # Cross-match
    idx, sep2d, _ = det.match_to_catalog_sky(cat)
    sep_arcsec = sep2d.arcsec
    
    # Check for NaN values in separation
    if np.any(np.isnan(sep_arcsec)):
        LOG.warning(f"NaN values found in separation distances: {sep_arcsec}")

    # Mask the values that exceed max_sep_arcsec
    idx_out = idx.astype(float)
    idx_out[sep_arcsec > max_sep_arcsec] = np.nan
    sep_out = np.where(sep_arcsec > max_sep_arcsec, np.nan, sep_arcsec)

    # Debugging the cross-matching result
    LOG.info(f"Matched {len(sep_out)} sources from APASS with separation <= {max_sep_arcsec} arcsec.")
    if len(sep_out) > 0:
        LOG.info(f"Sample of cross-matching distances: {sep_out[:10]}")
    
    return idx_out, sep_out

def vile_crossmatch_apass(ra, dec, apass_df, max_sep_arcsec):
    if len(apass_df) > 0:
        import pandas as _pd
        apass_df["RAJ2000"] = _pd.to_numeric(apass_df["RAJ2000"], errors="coerce")
        apass_df["DEJ2000"] = _pd.to_numeric(apass_df["DEJ2000"], errors="coerce")
    if len(apass_df) == 0:
        return np.full(ra.shape, np.nan), np.full(ra.shape, np.nan)
    det = SkyCoord(ra*u.deg, dec*u.deg)
    cat = SkyCoord(apass_df["RAJ2000"].values*u.deg, apass_df["DEJ2000"].values*u.deg)
    
    # Debug print for RA/Dec values of detections and catalog
    LOG.info(f"Detected RA/Dec: {ra}, {dec}")
    LOG.info(f"APASS RA/Dec: {apass_df[['RAJ2000', 'DEJ2000']].head()}")

    # Cross-match
    idx, sep2d, _ = det.match_to_catalog_sky(cat)
    sep_arcsec = sep2d.arcsec
    idx_out = idx.astype(float)
    idx_out[sep_arcsec > max_sep_arcsec] = np.nan
    sep_out = np.where(sep_arcsec > max_sep_arcsec, np.nan, sep_arcsec)
    
    # Debugging the cross-matching result
    LOG.info(f"Matched {len(sep_out)} sources from APASS with separation <= {max_sep_arcsec} arcsec.")
    if len(sep_out) > 0:
        LOG.info(f"Sample of cross-matching distances: {sep_out[:10]}")
    
    return idx_out, sep_out

def really_offensive_crossmatch_apass(ra, dec, apass_df, max_sep_arcsec):
    if len(apass_df) > 0:
        import pandas as _pd
        apass_df["RAJ2000"] = _pd.to_numeric(apass_df["RAJ2000"], errors="coerce")
        apass_df["DEJ2000"] = _pd.to_numeric(apass_df["DEJ2000"], errors="coerce")
    if len(apass_df) == 0:
        return np.full(ra.shape, np.nan), np.full(ra.shape, np.nan)
    det = SkyCoord(ra*u.deg, dec*u.deg)
    cat = SkyCoord(apass_df["RAJ2000"].values*u.deg, apass_df["DEJ2000"].values*u.deg)
    idx, sep2d, _ = det.match_to_catalog_sky(cat)
    sep_arcsec = sep2d.arcsec
    idx_out = idx.astype(float)
    idx_out[sep_arcsec > max_sep_arcsec] = np.nan
    sep_out = np.where(sep_arcsec > max_sep_arcsec, np.nan, sep_arcsec)
    
    # Debugging the cross-matching result
    LOG.info(f"Matched {len(sep_out)} sources from APASS with separation <= {max_sep_arcsec} arcsec.")
    if len(sep_out) > 0:
        LOG.info(f"Sample of cross-matching distances: {sep_out[:10]}")
    
    return idx_out, sep_out

def old_crossmatch_apass(ra, dec, apass_df, max_sep_arcsec):
    if len(apass_df) > 0:
        import pandas as _pd
        apass_df["RAJ2000"] = _pd.to_numeric(apass_df["RAJ2000"], errors="coerce")
        apass_df["DEJ2000"] = _pd.to_numeric(apass_df["DEJ2000"], errors="coerce")
    if len(apass_df) == 0:
        return np.full(ra.shape, np.nan), np.full(ra.shape, np.nan)
    det = SkyCoord(ra*u.deg, dec*u.deg)
    cat = SkyCoord(apass_df["RAJ2000"].values*u.deg, apass_df["DEJ2000"].values*u.deg)
    idx, sep2d, _ = det.match_to_catalog_sky(cat)
    sep_arcsec = sep2d.arcsec
    idx_out = idx.astype(float)
    idx_out[sep_arcsec > max_sep_arcsec] = np.nan
    sep_out = np.where(sep_arcsec > max_sep_arcsec, np.nan, sep_arcsec)
    return idx_out, sep_out


def process_tile_cube(cube_path, output_dir, s):
    cube, hdr = load_cube(cube_path)
    # Force 2-D celestial WCS (valid for all planes)
    w = WCS(hdr, naxis=2)

    R  = cube[0]
    G1 = cube[1]
    G2 = cube[2]
    B  = cube[3]

    ny, nx = R.shape
    tag = tile_tag(cube_path)
    LOG.info(f"{tag}: {os.path.basename(cube_path)}  shape={cube.shape}")

    detect_img = 0.5 * (G1 + G2)
    objs, bkg = detect_sources(detect_img, s["detect_thresh"], s["minarea"], s["deblend"])
    if len(objs) == 0:
        LOG.warning(f"{tag}: no detections on G1+G2; skipping tile")
        return None

    x = objs["x"].astype(np.float64)
    y = objs["y"].astype(np.float64)

    flux_R,  mag_R,  flag_R  = aperture_phot(R,  x, y, s["apr"], s["rin"], s["rout"])
    flux_G1, mag_G1, flag_G1 = aperture_phot(G1, x, y, s["apr"], s["rin"], s["rout"])
    flux_G2, mag_G2, flag_G2 = aperture_phot(G2, x, y, s["apr"], s["rin"], s["rout"])
    flux_B,  mag_B,  flag_B  = aperture_phot(B,  x, y, s["apr"], s["rin"], s["rout"])

    xy = np.vstack([x, y]).T
    sky = w.wcs_pix2world(xy, 0)
    ra  = sky[:,0].astype(np.float64)
    dec = sky[:,1].astype(np.float64)

    # Compute cone centre and half-diagonal (deg), then clamp to max radius
    ra_c, dec_c, r_half = tile_center_radius_deg(w, ny, nx)
    radius_deg = min(r_half, s["max_apass_rad"])
    LOG.info(f"{tag}: APASS cone: RA={ra_c:.6f}, Dec={dec_c:.6f}, r={radius_deg:.3f} deg (half-diag={r_half:.3f})")

    apass = query_apass_cone(ra_c, dec_c, radius_deg, retries=s["retries"])
    LOG.info(f"{tag}: APASS query @ RA={ra_c:.3f}, Dec={dec_c:.3f}, r={radius_deg:.2f}°")
    LOG.info(f"{tag}: APASS rows returned = {len(apass)}")

    ap_idx, ap_sep = crossmatch_apass(ra, dec, apass, s["match_arcsec"])
    diagnostics_offsets(ra, dec, apass, ap_idx, ap_sep, tag, p=95)


    out = pd.DataFrame({
        "tile": tag,
        "x": x, "y": y, "ra": ra, "dec": dec,
        "flux_R": flux_R,   "mag_inst_R": mag_R,   "flag_R": flag_R,
        "flux_G1": flux_G1, "mag_inst_G1": mag_G1, "flag_G1": flag_G1,
        "flux_G2": flux_G2, "mag_inst_G2": mag_G2, "flag_G2": flag_G2,
        "flux_B": flux_B,   "mag_inst_B": mag_B,   "flag_B": flag_B,
        "match_sep_arcsec": ap_sep
    })

    for c in ["Vmag","e_Vmag","Bmag","e_Bmag","RAJ2000","DEJ2000"]:
        out[c] = np.nan

    if len(apass) > 0:
        for i in range(len(out)):
            j = ap_idx[i]
            if not np.isnan(j):
                row = apass.iloc[int(j)]
                out.at[i, "Vmag"]    = row.get("Vmag", np.nan)
                out.at[i, "e_Vmag"]  = row.get("e_Vmag", np.nan)
                out.at[i, "Bmag"]    = row.get("Bmag", np.nan)
                out.at[i, "e_Bmag"]  = row.get("e_Bmag", np.nan)
                out.at[i, "RAJ2000"] = row.get("RAJ2000", np.nan)
                out.at[i, "DEJ2000"] = row.get("DEJ2000", np.nan)

    out["matched"] = ~np.isnan(out["Vmag"])

    os.makedirs(output_dir, exist_ok=True)
    stem = stem_from_cube(cube_path)
    tile_csv = os.path.join(output_dir, f"{os.path.basename(stem)}_{tag}_photometry.csv")
    out.to_csv(tile_csv, index=False)
    LOG.info(f"{tag}: wrote {tile_csv}  (detections={len(out)}, matched={int(out['matched'].sum())})")

    return out, stem

def old_process_tile_cube(cube_path, output_dir, s):
    cube, hdr = load_cube(cube_path)
    # Force 2-D celestial WCS (valid for all planes)
    w = WCS(hdr, naxis=2)

    R  = cube[0]
    G1 = cube[1]
    G2 = cube[2]
    B  = cube[3]

    ny, nx = R.shape
    tag = tile_tag(cube_path)
    LOG.info(f"{tag}: {os.path.basename(cube_path)}  shape={cube.shape}")

    detect_img = 0.5 * (G1 + G2)
    objs, bkg = detect_sources(detect_img, s["detect_thresh"], s["minarea"], s["deblend"])
    if len(objs) == 0:
        LOG.warning(f"{tag}: no detections on G1+G2; skipping tile")
        return None

    x = objs["x"].astype(np.float64)
    y = objs["y"].astype(np.float64)

    flux_R,  mag_R,  flag_R  = aperture_phot(R,  x, y, s["apr"], s["rin"], s["rout"])
    flux_G1, mag_G1, flag_G1 = aperture_phot(G1, x, y, s["apr"], s["rin"], s["rout"])
    flux_G2, mag_G2, flag_G2 = aperture_phot(G2, x, y, s["apr"], s["rin"], s["rout"])
    flux_B,  mag_B,  flag_B  = aperture_phot(B,  x, y, s["apr"], s["rin"], s["rout"])

    xy = np.vstack([x, y]).T
    sky = w.wcs_pix2world(xy, 0)
    ra  = sky[:,0].astype(np.float64)
    dec = sky[:,1].astype(np.float64)

    ra_c, dec_c, r_half = tile_center_radius_deg(w, ny, nx)
    LOG.info(f"RA: {ra_c}, Dec: {dec_c}, Radius: {radius_deg}")

    radius_deg = min(r_half, s["max_apass_rad"])

    apass = query_apass_cone(ra_c, dec_c, radius_deg, retries=s["retries"])
    LOG.info(f"{tag}: APASS query @ RA={ra_c:.3f}, Dec={dec_c:.3f}, r={radius_deg:.2f}°")
    LOG.info(f"{tag}: APASS rows returned = {len(apass)}")

    ap_idx, ap_sep = crossmatch_apass(ra, dec, apass, s["match_arcsec"])

    out = pd.DataFrame({
        "tile": tag,
        "x": x, "y": y, "ra": ra, "dec": dec,
        "flux_R": flux_R,   "mag_inst_R": mag_R,   "flag_R": flag_R,
        "flux_G1": flux_G1, "mag_inst_G1": mag_G1, "flag_G1": flag_G1,
        "flux_G2": flux_G2, "mag_inst_G2": mag_G2, "flag_G2": flag_G2,
        "flux_B": flux_B,   "mag_inst_B": mag_B,   "flag_B": flag_B,
        "match_sep_arcsec": ap_sep
    })

    for c in ["Vmag","e_Vmag","Bmag","e_Bmag","RAJ2000","DEJ2000"]:
        out[c] = np.nan

    if len(apass) > 0:
        for i in range(len(out)):
            j = ap_idx[i]
            if not np.isnan(j):
                row = apass.iloc[int(j)]
                out.at[i, "Vmag"]    = row.get("Vmag", np.nan)
                out.at[i, "e_Vmag"]  = row.get("e_Vmag", np.nan)
                out.at[i, "Bmag"]    = row.get("Bmag", np.nan)
                out.at[i, "e_Bmag"]  = row.get("e_Bmag", np.nan)
                out.at[i, "RAJ2000"] = row.get("RAJ2000", np.nan)
                out.at[i, "DEJ2000"] = row.get("DEJ2000", np.nan)

    out["matched"] = ~np.isnan(out["Vmag"])

    os.makedirs(output_dir, exist_ok=True)
    stem = stem_from_cube(cube_path)
    tile_csv = os.path.join(output_dir, f"{os.path.basename(stem)}_{tag}_photometry.csv")
    out.to_csv(tile_csv, index=False)
    LOG.info(f"{tag}: wrote {tile_csv}  (detections={len(out)}, matched={int(out['matched'].sum())})")

    return out, stem

def main():
    ap = argparse.ArgumentParser(description="Photometry + APASS match on WCS'd 4-plane tile cubes.")
    ap.add_argument("tile_dir",   help="Directory with *_t??_cube_rg1g2b.fits (e.g., ISS_ARCHIVE/)")
    ap.add_argument("output_dir", help="Directory to write CSVs")
    ap.add_argument("--apr",      type=float, default=DEF_APRADIUS,       help="Aperture radius [pix]")
    ap.add_argument("--rin",      type=float, default=DEF_BKG_INNER,      help="Annulus inner radius [pix]")
    ap.add_argument("--rout",     type=float, default=DEF_BKG_OUTER,      help="Annulus outer radius [pix]")
    ap.add_argument("--thresh",   type=float, default=DEF_DETECT_THRESH,  help="Detection threshold [sigma]")
    ap.add_argument("--minarea",  type=int,   default=DEF_DETECT_MINAREA, help="Detection min area [pix]")
    ap.add_argument("--deblend",  type=float, default=DEF_DEBLEND_CONT,   help="SEP deblend contrast")
    ap.add_argument("--match",    type=float, default=DEF_MATCH_ARCSEC,   help="Match radius [arcsec]")
    ap.add_argument("--maxrad",   type=float, default=DEF_MAX_APASS_RAD,  help="Max APASS cone radius per tile [deg]")
    ap.add_argument("--retries",  type=int,   default=DEF_RETRIES,        help="APASS query retries")
    args = ap.parse_args()

    settings = dict(
        apr=args.apr, rin=args.rin, rout=args.rout,
        detect_thresh=args.thresh, minarea=args.minarea, deblend=args.deblend,
        match_arcsec=args.match, max_apass_rad=args.maxrad, retries=args.retries
    )

    cubes = sorted(glob.glob(os.path.join(args.tile_dir, CUBE_GLOB)))
    if not cubes:
        LOG.error(f"No cubes matching {CUBE_GLOB} in {args.tile_dir}")
        return

    merged = []
    out_stem = None
    for cp in cubes:
        try:
            res = process_tile_cube(cp, args.output_dir, settings)
            if res is not None:
                df, stem = res
                merged.append(df)
                out_stem = stem
        except Exception as e:
            LOG.error(f"Tile failed: {cp} -> {e}")

    if merged and out_stem is not None:
        field = pd.concat(merged, ignore_index=True)
        field_csv = os.path.join(args.output_dir, f"{os.path.basename(out_stem)}_field_catalog.csv")
        field.to_csv(field_csv, index=False)
        LOG.info(f"Wrote merged field catalog: {field_csv}  (rows={len(field)}, matched={int(field['matched'].sum())})")
    else:
        LOG.warning("No per-tile outputs to merge; nothing written.")

if __name__ == "__main__":
    main()

