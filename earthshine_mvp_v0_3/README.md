Earthshine MVP v0.9.12 — RAW → 4-plane FITS archiver

This tool ingests RAW camera files (NEF/CR2/… or RAWs inside .gz/.tgz/.tar) and produces Level-1 4-plane FITS cubes (R, G1, G2, B) with rich headers suitable for downstream Earthshine/photometry work. It can decode via dcraw_emu/dcraw or rawpy, keeps a SQLite manifest, and optionally groups outputs into size-bounded PACKS/ folders. 

es_mvp_v0_9_12

Features (brief)

Inputs: RAW files (.NEF .CR2 .CR3 .ARW .RAF .ORF .RW2, case-insensitive) and/or archives containing RAWs (.gz .tgz .tar.gz .tar). Archives are auto-extracted to EXTRACT/.  Note that the high-efficiency HE* NEF files from Nikon Z9 camera cannot be archived with this software. Use of Adobe DNGConverter is required for that.

es_mvp_v0_9_12

Decoding routes:

dcraw_emu/dcraw (fully linear, no demosaic) when available, or

rawpy (LibRaw) with --no-dcraw to force this route. 

es_mvp_v0_9_12

Output: 4-plane FITS (R,G1,G2,B) with key headers: BAYERPAT, ADC_BITS, INSTRUME, BODYSER, lens fields, DATE-OBS, FILTER, PLATFORM (ISS if filename starts with iss…). File name is templated and includes SHA tags. 

es_mvp_v0_9_12

Manifest DB: MANIFEST/earthshine_manifest.sqlite tracks per-RAW status/hash/outputs (idempotent processing). 

es_mvp_v0_9_12

Packing: hard-link packs under PACKS/pack_### with size/file count limits (defaults: 49 GB, 100 files). 

es_mvp_v0_9_12

Safety: refuses scaled data (e.g., if decoder expanded beyond ADC limit); supports --skip-if-fits-exists and --overwrite-fits. 

es_mvp_v0_9_12

Quick start
# Example (recommended pattern)
uv run es_mvp_v0_9_12.py \
  --filelist all_iss_files.txt \
  --root /path/to/ISS_ARCHIVE_TEST \
  --no-dcraw


--filelist is a text file with one path per line (RAWs and/or archives).

Outputs go to: /path/to/ISS_ARCHIVE_TEST/L1/*.fits, plus PACKS/ and MANIFEST/. 

es_mvp_v0_9_12

CLI (most useful flags)

Inputs: --filelist FILE, --input DIR (repeatable), --glob PATTERN (repeatable), --exts ".NEF,.CR2,...", --follow-symlinks, --accept-archives/--extract-dir. 

es_mvp_v0_9_12

Output root (required): --root /path/to/archive_root (creates L1/ PACKS/ MANIFEST/ LOGS/). 

es_mvp_v0_9_12

Decoding: --no-dcraw (force rawpy), otherwise auto-choose dcraw_emu/dcraw when suitable. 

es_mvp_v0_9_12

Idempotency & overwrite: --skip-if-fits-exists, --overwrite-fits. 

es_mvp_v0_9_12

Naming: --name-template (keys: {orig},{ext},{sha8},{sha16},{sha256},{date},{time},{model}; default ends with _L1_rawcube.fits). 

es_mvp_v0_9_12

Packs: --max-pack-gb 49.0, --max-pack-files 100. 

es_mvp_v0_9_12

Dry run: --dry-run (summarises candidates and missing entries). 

es_mvp_v0_9_12

FITS specifics

Cube order: [R, G1, G2, B]; BITPIX=16, BUNIT=adu.

Header includes ADC/black/white levels when available; camera body serial; lens model/spec; exposure, f-number, ISO; SHA of original file (ORIGHASH). Archive provenance stored when decoded from inside .gz/.tar. 

es_mvp_v0_9_12

Requirements

Python 3.10+; packages: numpy, tifffile, astropy, exifread, optional rawpy. (SQLite is stdlib.)

Optional external: dcraw_emu or dcraw on PATH (for the dcraw route). 

es_mvp_v0_9_12

Troubleshooting

Scaled data error (e.g., “max > limit for 14-bit”): ensure decoder is truly linear; try --no-dcraw (rawpy route) or correct dcraw_emu flags. 

es_mvp_v0_9_12

No inputs: check --filelist paths; the tool prints the first missing lines on error. 

es_mvp_v0_9_12

Re-runs: processed files are skipped via the manifest status unless you change --procver or force overwrite. 

es_mvp_v0_9_12

License

Choose a permissive license (e.g., MIT) for broad reuse in research/software.
