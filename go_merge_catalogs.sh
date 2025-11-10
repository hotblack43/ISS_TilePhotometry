uv run merge_field_catalogs.py \
  --glob './TESTOUT*/**/*_field_catalog.csv' \
  --match-arcsec 8 \
  --only-matched \
  --outfile merged_field_catalog_averaged.csv \
  --members-outfile merged_members_listing.csv

