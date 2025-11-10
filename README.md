BEGIN README.md
ISS Rawcube Tiling, WCS Solving, and APASS Photometry Pipeline

This repository provides a complete pipeline for processing ISS
Earth–Moon rawcube FITS images into WCS-solved tiles, performing
photometry, matching sources to the APASS DR9 catalogue, and finally
merging overlapping tile catalogues into a single per-field photometric
dataset.

The input data (the 'rawcube' files) are generated elsewhere. See
code es_mvp_v0_9_12.py in earthshine_mvp_v0_3/ read the README there.

The system is designed for high-precision Earthshine photometry but
may be useful for other applications involving 4-plane (R, G1, G2, B)
Bayer rawcube images.

Overview of the Pipeline

The processing consists of four major stages:

Tile the rawcube and perform WCS solving
(iss_rawcube_tile_wcs_local_first.py)

Perform photometry and APASS crossmatching
(iss_tilecube_photometry_apass.py)

Merge overlapping tile catalogues
(merge_field_catalogs.py)

Plot and analyse the merged results
(R notebooks: plot_merged_averaged_catalogs.Rmd, plot_results_from_tilesolver_1.Rmd)

A shorter description is available in HOWTO.txt.

1. Tiling and WCS Solving

Script: iss_rawcube_tile_wcs_local_first.py

This script:

Loads a 4-plane FITS rawcube containing R, G1, G2 and B channels.

Splits the rawcube into an overlapping grid of tiles (typically 3×3).

Creates a white-light tile for WCS solving.

Attempts local astrometry.net solution first; if that fails, falls back
to nova.astrometry.net.

Writes tile products:

*_tNN_plane_R.fits, *_tNN_plane_G1.fits, *_tNN_plane_G2.fits, *_tNN_plane_B.fits

*_tNN_cube_rg1g2b.fits (4-plane tile with WCS headers)

Example
uv run iss_rawcube_tile_wcs_local_first.py RAWCUBE.fits TESTOUT1


Batch processing:

./go.sh

2. Photometry and APASS Crossmatching

Script: iss_tilecube_photometry_apass.py

This step:

Reads WCS-solved tiles (*_tNN_cube_rg1g2b.fits)

Detects sources (using combined G1/G2)

Performs aperture photometry in all four bands

Converts (x,y) → (RA,Dec) via WCS

Queries APASS DR9

Matches detections to APASS

Produces per-tile and per-rawcube field catalogues

Example
uv run iss_tilecube_photometry_apass.py TESTOUT1 TESTOUT1 --maxrad 2.5 --match 60


Batch version:

./go2.sh

3. Merge Field Catalogues

Script: merge_field_catalogs.py

Features:

Loads all *_field_catalog.csv files

Clusters detections using friends-of-friends in sky coordinates

Computes median RA/Dec per cluster

Computes robust per-band machine magnitudes (median/mean/std/count)

Averages APASS V/B magnitudes

Outputs:

merged_field_catalog_averaged.csv

merged_members_listing.csv (optional)

Example
uv run merge_field_catalogs.py \
    --glob './TESTOUT*/**/*_field_catalog.csv' \
    --match-arcsec 8 \
    --only-matched \
    --outfile merged_field_catalog_averaged.csv \
    --members-outfile merged_members_listing.csv


Convenience command:

./go_merge_catalogs.sh

4. Plotting and Analysis

Two R notebooks are provided:

plot_merged_averaged_catalogs.Rmd

plot_results_from_tilesolver_1.Rmd

Render using RStudio or:

Rscript -e "rmarkdown::render('plot_merged_averaged_catalogs.Rmd')"

Dependencies
Python

numpy

pandas

astropy

sep

astroquery

requests

tqdm (optional)

uv (optional, recommended)

External Tools

astrometry.net (solve-field + index files)

Or: nova.astrometry.net API (set ASTROMETRY_API_KEY)

R

tidyverse

ggplot2

dplyr

rmarkdown

Example End-to-End Workflow
# 1. Tile and WCS-solve
uv run iss_rawcube_tile_wcs_local_first.py RAWCUBE.fits TESTOUT1

# 2. Photometry + APASS matching
uv run iss_tilecube_photometry_apass.py TESTOUT1 TESTOUT1 --maxrad 2.5 --match 60

# 3. Merge all tile catalogs
uv run merge_field_catalogs.py \
    --glob './TESTOUT1/**/*_field_catalog.csv' \
    --match-arcsec 8 \
    --outfile merged_catalog.csv

# 4. Plot in R
Rscript -e "rmarkdown::render('plot_merged_averaged_catalogs.Rmd')"

Repository Structure
iss_rawcube_tile_wcs_local_first.py     # Step 1: tiling + WCS
iss_tilecube_photometry_apass.py        # Step 2: photometry + APASS match
merge_field_catalogs.py                 # Step 3: merge catalogs
go.sh                                   # Batch WCS solving
go2.sh                                  # Batch photometry
go_merge_catalogs.sh                    # Batch merge
plot_merged_averaged_catalogs.Rmd       # R analysis
plot_results_from_tilesolver_1.Rmd      # Diagnostics
HOWTO.txt                               # Short user guide

Notes

Scripts are designed to run independently but are usually executed in the order described.

If local astrometry.net is installed with index files, WCS solving is fast and reproducible.

Rawcubes must be 4-plane monochrome FITS (R,G1,G2,B).

END README.md
