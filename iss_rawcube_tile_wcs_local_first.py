#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Peter Thejll, DMI November 2025

iss_rawcube_tile_wcs_local_first.py

Tile a 4-plane (R, G1, G2, B) rawcube FITS, solve WCS per tile.
Preference order:
  1) Local astrometry.net `solve-field` (fast, offline)
  2) nova.astrometry.net API (online fallback)

If solved, the WCS header is copied onto each colour-plane tile and a 4-plane cube tile.

Prereqs: astropy, numpy, requests (only needed for nova fallback)

Example usage:

    uv run iss_rawcube_tile_wcs_local_first.py /home/pth/pCloudDrive/ACRYLIC/IMAGEFILES/ISS066E0321/20211031_162721_iss066e032140_5c2cc127_L1_rawcube.fits TESTOUT6
"""
import os
import sys
import time
import json
import argparse
import logging
import shutil
import subprocess
from io import BytesIO
from typing import Tuple, List, Optional

import numpy as np
from astropy.io import fits

# ----------------------------
# Logging
# ----------------------------
logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")
LOG = logging.getLogger("iss_rawcube_tile_wcs")

# ----------------------------
# NOVA endpoints
# ----------------------------
NOVA_BASE = "https://nova.astrometry.net"
API_LOGIN = f"{NOVA_BASE}/api/login"
API_UPLOAD = f"{NOVA_BASE}/api/upload"
API_SUB_STATUS = f"{NOVA_BASE}/api/submissions"
API_JOB_STATUS = f"{NOVA_BASE}/api/jobs"
WCS_FITS_URL = f"{NOVA_BASE}/wcs_file"

# ----------------------------
# Utils
# ----------------------------
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def read_fits_cube(path: str) -> np.ndarray:
    with fits.open(path) as hdul:
        data = hdul[0].data
    if data is None or data.ndim != 3 or data.shape[0] != 4:
        raise ValueError(f"Expected 4-plane cube (R,G1,G2,B). Got {None if data is None else data.shape} from {path}")
    return data

def write_fits(path: str, data: np.ndarray, header: Optional[fits.Header] = None) -> None:
    fits.PrimaryHDU(data=data, header=header).writeto(path, overwrite=True)

def tile_bounds(ny: int, nx: int, grid: int, overlap: float) -> List[Tuple[int, int, int, int]]:
    if not (0.0 <= overlap < 0.49):
        raise ValueError("overlap must be in [0, 0.49)")
    step_y = int(round(ny / grid))
    step_x = int(round(nx / grid))
    oy = int(round(step_y * overlap))
    ox = int(round(step_x * overlap))
    tiles = []
    for gy in range(grid):
        for gx in range(grid):
            y0 = max(0, gy * step_y - (oy if gy > 0 else 0))
            x0 = max(0, gx * step_x - (ox if gx > 0 else 0))
            y1 = min(ny, (gy + 1) * step_y + (oy if gy < grid - 1 else 0))
            x1 = min(nx, (gx + 1) * step_x + (ox if gx < grid - 1 else 0))
            if (y1 - y0) > 10 and (x1 - x0) > 10:
                tiles.append((y0, y1, x0, x1))
    return tiles

# ----------------------------
# Local solve-field backend
# ----------------------------
def have_solve_field() -> bool:
    return shutil.which("solve-field") is not None

def build_solve_field_cmd(
    inp: str,
    out_wcs: str,
    scale_lower: Optional[float],
    scale_upper: Optional[float],
    center_ra: Optional[float],
    center_dec: Optional[float],
    radius: Optional[float],
    downsample: Optional[int],
    cpulimit: Optional[int],
) -> List[str]:
    cmd = [
        "solve-field",
        "--overwrite",
        "--no-plots",
        "--wcs", out_wcs,
    ]
    # Hints: scale in degrees across (degwidth)
    if scale_lower is not None and scale_upper is not None:
        cmd += ["--scale-units", "degwidth", "--scale-low", str(scale_lower), "--scale-high", str(scale_upper)]
    # Hints: center + radius (degrees)
    if center_ra is not None and center_dec is not None and radius is not None:
        cmd += ["--ra", str(center_ra), "--dec", str(center_dec), "--radius", str(radius)]
    # Downsample factor
    if downsample is not None and downsample > 1:
        cmd += ["--downsample", str(downsample)]
    # CPU time limit (sec)
    if cpulimit is not None and cpulimit > 0:
        cmd += ["--cpulimit", str(cpulimit)]
    # Input file last
    cmd.append(inp)
    return cmd

def try_local_solve(
    summed_path: str,
    hints: dict,
) -> Optional[fits.Header]:
    out_wcs = os.path.splitext(summed_path)[0] + ".wcs"
    cmd = build_solve_field_cmd(
        inp=summed_path,
        out_wcs=out_wcs,
        scale_lower=hints.get("scale_lower"),
        scale_upper=hints.get("scale_upper"),
        center_ra=hints.get("center_ra"),
        center_dec=hints.get("center_dec"),
        radius=hints.get("radius"),
        downsample=hints.get("downsample"),
        cpulimit=hints.get("cpulimit"),
    )
    LOG.info("Local solve-field: " + " ".join(cmd))
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=hints.get("timeout_local", 900))
        LOG.info("solve-field output:\n" + res.stdout.strip())
        if res.returncode != 0:
            LOG.warning(f"solve-field exited with code {res.returncode}")
            return None
        if not os.path.exists(out_wcs):
            LOG.warning("solve-field claimed success but .wcs not found")
            return None
        with fits.open(out_wcs) as hdul:
            return hdul[0].header.copy()
    except subprocess.TimeoutExpired:
        LOG.error("solve-field timed out.")
        return None
    except Exception as e:
        LOG.error(f"solve-field error: {e}")
        return None

# ----------------------------
# Nova backend
# ----------------------------
import requests

class NovaClient:
    def __init__(self, apikey: str, timeout_s: int = 3600, poll_s: int = 10, max_upload_retries: int = 5):
        self.apikey = apikey
        self.timeout_s = timeout_s
        self.poll_s = poll_s
        self.session = None
        self.max_upload_retries = max_upload_retries

    def login(self) -> None:
        LOG.info("Logging in to nova.astrometry.net ...")
        resp = requests.post(API_LOGIN, data={"request-json": json.dumps({"apikey": self.apikey})}, timeout=30)
        resp.raise_for_status()
        js = resp.json()
        if js.get("status") != "success":
            raise RuntimeError(f"Login failed: {js}")
        self.session = js.get("session")
        if not self.session:
            raise RuntimeError("Login did not return a session key.")
        LOG.info("Login OK; session acquired.")

    def upload_image(self, image_path: str) -> int:
        if not self.session:
            raise RuntimeError("Not logged in.")
        for attempt in range(1, self.max_upload_retries + 1):
            try:
                LOG.info(f"Uploading {os.path.basename(image_path)} (attempt {attempt}/{self.max_upload_retries}) ...")
                with open(image_path, "rb") as f:
                    files = {"file": f}
                    data = {"request-json": json.dumps({"session": self.session})}
                    resp = requests.post(API_UPLOAD, files=files, data=data, timeout=180)
                resp.raise_for_status()
                js = resp.json()
                if js.get("status") != "success":
                    raise RuntimeError(f"Upload failed: {js}")
                subid = js.get("subid")
                if subid is None:
                    raise RuntimeError(f"No submission id returned: {js}")
                LOG.info(f"Upload OK; submission id = {subid}")
                return int(subid)
            except Exception as e:
                LOG.warning(f"Upload attempt {attempt} failed: {e}")
                if attempt < self.max_upload_retries:
                    time.sleep(min(60, 2 ** attempt))  # backoff
                else:
                    raise

    def wait_for_job(self, subid: int) -> int:
        start = time.time()
        jobid = None
        while time.time() - start < self.timeout_s:
            r = requests.get(f"{API_SUB_STATUS}/{subid}", timeout=30)
            r.raise_for_status()
            js = r.json()
            jobs = js.get("jobs") or []
            if jobs:
                jobid = jobs[0]
                LOG.info(f"Submission {subid} → job {jobid}")
                break
            time.sleep(self.poll_s)
        if not jobid:
            raise TimeoutError(f"Timed out waiting for job from submission {subid}")
        while time.time() - start < self.timeout_s:
            r = requests.get(f"{API_JOB_STATUS}/{jobid}", timeout=30)
            r.raise_for_status()
            js = r.json()
            status = js.get("status")
            if status == "success":
                LOG.info(f"Job {jobid} solved.")
                return int(jobid)
            elif status in ("failure", "aborted"):
                raise RuntimeError(f"Job {jobid} failed: {js}")
            time.sleep(self.poll_s)
        raise TimeoutError(f"Timed out waiting for job {jobid} to solve")

    def fetch_wcs_header(self, jobid: int) -> fits.Header:
        url = f"{WCS_FITS_URL}/{jobid}"
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        with fits.open(BytesIO(r.content)) as hdul:
            return hdul[0].header.copy()

# ----------------------------
# Core
# ----------------------------
def process_file(
    in_fits: str,
    out_dir: str,
    grid: int,
    overlap: float,
    hints: dict,
    apikey: Optional[str],
    timeout_nova: int,
    poll_nova: int,
) -> None:
    stem = os.path.splitext(os.path.basename(in_fits))[0]
    ensure_dir(out_dir)

    cube = read_fits_cube(in_fits)   # (4, ny, nx)
    ny, nx = cube.shape[1], cube.shape[2]
    tiles = tile_bounds(ny, nx, grid, overlap)
    LOG.info(f"Image {stem}: shape (ny={ny}, nx={nx}), grid={grid} → {len(tiles)} tiles")

    use_local = have_solve_field()
    if use_local:
        LOG.info("Local solve-field found — will try local solving first.")
    else:
        LOG.info("Local solve-field NOT found — will use nova fallback.")

    nova = None
    if not use_local:
        if not apikey:
            LOG.error("No local solve-field and no nova API key; cannot solve.")
            return
        nova = NovaClient(apikey=apikey, timeout_s=timeout_nova, poll_s=poll_nova)
        nova.login()

    for t_idx, (y0, y1, x0, x1) in enumerate(tiles, start=1):
        tag = f"t{t_idx:02d}"
        LOG.info(f"Processing tile {tag}: y[{y0}:{y1}) x[{x0}:{x1})")

        # Extract per-plane tiles
        r_tile  = cube[0, y0:y1, x0:x1]
        g1_tile = cube[1, y0:y1, x0:x1]
        g2_tile = cube[2, y0:y1, x0:x1]
        b_tile  = cube[3, y0:y1, x0:x1]

        # Summed tile for solving
        summed = r_tile.astype(np.float64) + g1_tile + g2_tile + b_tile
        summed_path = os.path.join(out_dir, f"{stem}_{tag}_summed.fits")
        write_fits(summed_path, summed)

        wcs_hdr = None

        # Try local first
        if use_local:
            wcs_hdr = try_local_solve(summed_path, hints=hints)

        # Fallback to nova if needed
        if wcs_hdr is None and apikey:
            if nova is None:
                nova = NovaClient(apikey=apikey, timeout_s=timeout_nova, poll_s=poll_nova)
                nova.login()
            try:
                subid = nova.upload_image(summed_path)
                jobid = nova.wait_for_job(subid)
                wcs_hdr = nova.fetch_wcs_header(jobid)
            except Exception as e:
                LOG.error(f"{tag}: nova solving failed: {e}")

        # Write outputs (with or without WCS)
        if wcs_hdr is not None:
            hdr = wcs_hdr.copy()
            hdr["HISTORY"] = "WCS from local solve-field" if use_local and wcs_hdr else "WCS from nova.astrometry.net"
            hdr["HISTORY"] = f"tile: y0={y0}, y1={y1}, x0={x0}, x1={x1}"
            hdr["HISTORY"] = "planes: 1=R, 2=G1, 3=G2, 4=B"
            write_fits(os.path.join(out_dir, f"{stem}_{tag}_plane_R.fits"),  r_tile,  header=hdr)
            write_fits(os.path.join(out_dir, f"{stem}_{tag}_plane_G1.fits"), g1_tile, header=hdr)
            write_fits(os.path.join(out_dir, f"{stem}_{tag}_plane_G2.fits"), g2_tile, header=hdr)
            write_fits(os.path.join(out_dir, f"{stem}_{tag}_plane_B.fits"),  b_tile,  header=hdr)
            tile_cube = np.stack([r_tile, g1_tile, g2_tile, b_tile], axis=0)
            write_fits(os.path.join(out_dir, f"{stem}_{tag}_cube_rg1g2b.fits"), tile_cube, header=hdr)
        else:
            LOG.error(f"{tag}: No WCS obtained (local nor nova). Writing tiles without WCS for inspection.")
            write_fits(os.path.join(out_dir, f"{stem}_{tag}_plane_R.fits"),  r_tile)
            write_fits(os.path.join(out_dir, f"{stem}_{tag}_plane_G1.fits"), g1_tile)
            write_fits(os.path.join(out_dir, f"{stem}_{tag}_plane_G2.fits"), g2_tile)
            write_fits(os.path.join(out_dir, f"{stem}_{tag}_plane_B.fits"),  b_tile)
            tile_cube = np.stack([r_tile, g1_tile, g2_tile, b_tile], axis=0)
            write_fits(os.path.join(out_dir, f"{stem}_{tag}_cube_rg1g2b.fits"), tile_cube)

    LOG.info("Done.")

# ----------------------------
# CLI
# ----------------------------
def main():
    p = argparse.ArgumentParser(description="Tile a 4-plane rawcube and solve WCS per tile (local solve-field first, nova fallback).")
    p.add_argument("input_fits", help="Path to 4-plane rawcube FITS (R,G1,G2,B).")
    p.add_argument("output_dir", help="Directory for tiles and outputs.")
    p.add_argument("--grid", type=int, default=3, help="Grid size (grid×grid). Default 3.")
    p.add_argument("--overlap", type=float, default=0.10, help="Tile overlap fraction. Default 0.10.")
    # Hints (applied to both backends)
    p.add_argument("--scale-lower", type=float, default=None, help="Lower bound of field width (degrees).")
    p.add_argument("--scale-upper", type=float, default=None, help="Upper bound of field width (degrees).")
    p.add_argument("--center-ra", type=float, default=None, help="Approx RA (deg).")
    p.add_argument("--center-dec", type=float, default=None, help="Approx Dec (deg).")
    p.add_argument("--radius", type=float, default=None, help="Search radius (deg).")
    p.add_argument("--downsample", type=int, default=None, help="Downsample factor for solving (int > 1).")
    # Local solve options
    p.add_argument("--cpulimit", type=int, default=None, help="solve-field CPU time limit (seconds).")
    p.add_argument("--timeout-local", type=int, default=900, help="Timeout for local solve-field call (seconds).")
    # Nova fallback
    p.add_argument("--apikey", default=os.environ.get("ASTROMETRY_API_KEY", ""), help="nova API key (or env ASTROMETRY_API_KEY).")
    p.add_argument("--timeout-nova", type=int, default=3600, help="Per-tile timeout for nova (seconds).")
    p.add_argument("--poll-nova", type=int, default=10, help="Nova poll interval (seconds).")

    args = p.parse_args()

    hints = dict(
        scale_lower=args.scale_lower,
        scale_upper=args.scale_upper,
        center_ra=args.center_ra,
        center_dec=args.center_dec,
        radius=args.radius,
        downsample=args.downsample,
        cpulimit=args.cpulimit,
        timeout_local=args.timeout_local,
    )

    try:
        process_file(
            in_fits=args.input_fits,
            out_dir=args.output_dir,
            grid=args.grid,
            overlap=args.overlap,
            hints=hints,
            apikey=args.apikey,
            timeout_nova=args.timeout_nova,
            poll_nova=args.poll_nova,
        )
    except Exception as e:
        LOG.error(f"Fatal: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

