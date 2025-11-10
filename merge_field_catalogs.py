#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Peter Thejll, DMI, November 2025

Merge overlapping ISS tile field catalogs by sky position and average machine magnitudes.

Inputs:
  - One or more glob patterns pointing to *field_catalog.csv files.
    These are the outputs you listed, e.g. TESTOUT*/..._field_catalog.csv

Output:
  - A merged CSV with one row per clustered star, containing:
    * median RA/Dec (deg)
    * counts of detections
    * per-band machine magnitudes: median, mean, std, N_used
    * APASS catalog photometry (median of Vmag, e_Vmag, Bmag, e_Bmag)
    * cluster diagnostics

Requires: numpy, pandas, astropy


Example usage:
    uv run merge_field_catalogs.py \
  --glob './TESTOUT*/**/*_field_catalog.csv' \
  --match-arcsec 8 \
  --only-matched \
  --outfile merged_field_catalog_averaged.csv \
  --members-outfile merged_members_listing.csv

"""
import argparse
import glob
import os
import sys
from typing import List, Dict
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

BANDS = [
    ("R",  "mag_inst_R"),
    ("G1", "mag_inst_G1"),
    ("G2", "mag_inst_G2"),
    ("B",  "mag_inst_B"),
]
APASS_COLS = ["Vmag","e_Vmag","Bmag","e_Bmag"]

def load_catalogs(paths: List[str], only_matched: bool) -> pd.DataFrame:
    rows = []
    for p in paths:
        try:
            df = pd.read_csv(p)
        except Exception as e:
            print(f"WARNING: failed to read {p}: {e}", file=sys.stderr)
            continue
        df["source_file"] = p
        # prefer APASS RA/Dec if present & finite; else fall back to WCS RA/Dec
        ra_ap = pd.to_numeric(df.get("RAJ2000"), errors="coerce")
        de_ap = pd.to_numeric(df.get("DEJ2000"), errors="coerce")
        ra_im = pd.to_numeric(df.get("ra"), errors="coerce")
        de_im = pd.to_numeric(df.get("dec"), errors="coerce")
        ra_use = ra_ap.where(ra_ap.notna(), ra_im)
        de_use = de_ap.where(de_ap.notna(), de_im)
        df["ra_use"] = ra_use
        df["dec_use"] = de_use

        if only_matched and "matched" in df.columns:
            df = df[df["matched"] == True].copy()

        # keep only rows with usable coords
        df = df[np.isfinite(df["ra_use"]) & np.isfinite(df["dec_use"])].copy()
        if len(df):
            rows.append(df)

    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True, sort=False)

# Union-Find (Disjoint Set) for clustering
def uf_find(parent, i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i

def uf_union(parent, rank, i, j):
    ri, rj = uf_find(parent, i), uf_find(parent, j)
    if ri == rj: return
    if rank[ri] < rank[rj]:
        parent[ri] = rj
    elif rank[ri] > rank[rj]:
        parent[rj] = ri
    else:
        parent[rj] = ri
        rank[ri] += 1

def cluster_by_sky(ra_deg: np.ndarray, dec_deg: np.ndarray, match_arcsec: float) -> Dict[int, List[int]]:
    """Friends-of-friends clustering via search_around_sky within match_arcsec."""
    coords = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame="icrs")
    # Find all neighbors within the match radius
    idx1, idx2, sep2d, _ = coords.search_around_sky(coords, (match_arcsec * u.arcsec))
    n = len(coords)
    parent = list(range(n))
    rank = [0]*n
    for i, j in zip(idx1, idx2):
        if i == j:  # skip self-pairs
            continue
        # unify both ways (search_around_sky returns both orders; union handles it)
        uf_union(parent, rank, int(i), int(j))
    # Build components
    comp: Dict[int, List[int]] = {}
    for i in range(n):
        r = uf_find(parent, i)
        comp.setdefault(r, []).append(i)
    return comp

def robust_stats(x: pd.Series):
    x = pd.to_numeric(x, errors="coerce")
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan, 0, np.nan
    return float(np.nanmedian(x)), float(np.nanmean(x)), float(np.nanstd(x, ddof=1) if len(x) > 1 else 0.0), int(len(x))

def aggregate_cluster(df: pd.DataFrame, idxs: List[int], cluster_id: int) -> dict:
    sub = df.iloc[idxs].copy()
    # Sky position: robust central value (median)
    ra_med  = float(np.nanmedian(sub["ra_use"]))
    dec_med = float(np.nanmedian(sub["dec_use"]))
    n_obs = int(len(sub))

    row = {
        "cluster_id": cluster_id,
        "n_obs": n_obs,
        "ra_deg": ra_med,
        "dec_deg": dec_med,
    }

    # Per-band machine magnitudes
    for band, col in BANDS:
        med, mean, std, n = robust_stats(sub.get(col))
        row[f"{band}_mag_med"]  = med
        row[f"{band}_mag_mean"] = mean
        row[f"{band}_mag_std"]  = std
        row[f"{band}_mag_n"]    = n

    # APASS photometry (median across members)
    for c in APASS_COLS:
        if c in sub.columns:
            row[c] = float(np.nanmedian(pd.to_numeric(sub[c], errors="coerce")))
        else:
            row[c] = np.nan

    # How tight is the cluster? (median/max separation to the cluster median)
    try:
        cc = SkyCoord(sub["ra_use"].to_numpy()*u.deg, sub["dec_use"].to_numpy()*u.deg)
        c0 = SkyCoord(ra_med*u.deg, dec_med*u.deg)
        sep = cc.separation(c0).arcsec
        row["sep_median_arcsec"] = float(np.nanmedian(sep))
        row["sep_max_arcsec"]    = float(np.nanmax(sep))
    except Exception:
        row["sep_median_arcsec"] = np.nan
        row["sep_max_arcsec"]    = np.nan

    # Optional diagnostics
    row["source_files_unique"] = len(set(sub["source_file"]))
    row["source_files"] = ";".join(sorted(set(sub["source_file"])))
    return row

def main():
    ap = argparse.ArgumentParser(description="Merge overlapping field catalogs and average machine magnitudes.")
    ap.add_argument("--glob", dest="globs", action="append",
                    default=["./**/*_field_catalog.csv"],
                    help="Glob pattern(s) for input field catalogs (repeatable). Default: ./**/*_field_catalog.csv")
    ap.add_argument("--match-arcsec", type=float, default=8.0,
                    help="Clustering radius in arcseconds (friends-of-friends). Default 8.0")
    ap.add_argument("--only-matched", action="store_true",
                    help="If set, include only rows where 'matched' == True (APASS-matched).")
    ap.add_argument("--outfile", type=str, default="merged_field_catalog_averaged.csv",
                    help="Output CSV (one row per clustered star).")
    ap.add_argument("--members-outfile", type=str, default="merged_members_listing.csv",
                    help="Optional CSV listing members for each cluster.")
    args = ap.parse_args()

    # Resolve globs
    paths = []
    for g in args.globs:
        paths.extend(glob.glob(g, recursive=True))
    paths = sorted(set(p for p in paths if p.lower().endswith("_field_catalog.csv")))

    if not paths:
        print("ERROR: No input catalogs found. Check your --glob pattern(s).", file=sys.stderr)
        sys.exit(2)

    print(f"Found {len(paths)} input catalogs.")
    df = load_catalogs(paths, only_matched=args.only_matched)
    if df.empty:
        print("ERROR: After loading/filtering, no rows remain.", file=sys.stderr)
        sys.exit(3)

    # Build clusters
    ra = df["ra_use"].to_numpy(dtype=float)
    dec = df["dec_use"].to_numpy(dtype=float)
    comps = cluster_by_sky(ra, dec, match_arcsec=args.match_arcsec)
    print(f"Formed {len(comps)} clusters with match radius {args.match_arcsec:.2f}\" from {len(df)} detections.")

    # Aggregate per cluster
    rows = []
    for k, idxs in comps.items():
        rows.append(aggregate_cluster(df, idxs, cluster_id=len(rows)+1))
    out = pd.DataFrame(rows)
    # Order columns: basic, bands, APASS, diagnostics
    base_cols = ["cluster_id","n_obs","ra_deg","dec_deg"]
    band_cols = []
    for b,_ in BANDS:
        band_cols += [f"{b}_mag_med", f"{b}_mag_mean", f"{b}_mag_std", f"{b}_mag_n"]
    apass_cols = APASS_COLS
    diag_cols  = ["sep_median_arcsec","sep_max_arcsec","source_files_unique","source_files"]
    ordered = base_cols + band_cols + apass_cols + diag_cols
    other = [c for c in out.columns if c not in ordered]
    out = out[ordered + other]

    out.to_csv(args.outfile, index=False)
    print(f"Wrote merged catalog with {len(out)} stars → {args.outfile}")

    # Optional members dump
    members = []
    for cid, idxs in zip(out["cluster_id"], comps.values()):
        sub = df.iloc[idxs].copy()
        sub.insert(0, "cluster_id", cid)
        members.append(sub)
    if members:
        members_df = pd.concat(members, ignore_index=True, sort=False)
        members_df.to_csv(args.members_outfile, index=False)
        print(f"Wrote cluster membership listing → {args.members_outfile}")

if __name__ == "__main__":
    main()

