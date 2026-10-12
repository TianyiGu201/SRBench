# R2 payload compression assessment

Measured on the current `web_data` snapshot on 2026-10-12. Sizes below are
decimal GB (1 GB = 1,000,000,000 bytes). No source data was changed.

## What the static website actually reads

The baseline upload script sends 20,969 files totaling 28.255 GB. It already
excludes 36 `.bak` files (0.268 GB). A path-by-path comparison with the Vite
frontend and generated static manifests found another 237 files (0.325 GB)
that the deployed site never requests: 36 original out-of-sample CSVs (0.315
GB) already represented by `site/oos/.../*.json`, plus 201 extra raw-resolution
CSVs (0.010 GB) absent from the published manifests. Keep these sources on the
server for future regeneration; omit them from the R2 payload. The remaining
20,732 files total 27.930 GB.

The 13 JPEGs are used for overview previews. The 13 TIFFs are used for the
full-resolution HE-image download and must remain available. The NPYs are
fetched by the browser to draw expression images. The site JSON and the
manifest-selected raw-resolution CSVs drive the remaining views.

## Lossless compression measurements

| Required file group | Current GB | Proposed GB | Method |
| --- | ---: | ---: | --- |
| 3,108 NPY arrays | 16.813 | 4.289 | gzip level 6, restore original NPY bytes in the browser |
| 13 TIFF images | 9.527 | 4.520 | TIFF Deflate with horizontal predictor (`tiffcp -c zip:2:p6`) |
| 14,470 CSV files | 1.148 | 0.550 | gzip level 6, explicitly decoded by browser |
| 3,128 JSON files | 0.439 | 0.108 | gzip level 6, explicitly decoded by browser |
| 13 JPEG previews | 0.003 | 0.003 | Keep as-is |
| **Total** | **27.930** | **9.470** | **66.5% below the baseline upload size** |

Gzip sizes were calculated for every required NPY, CSV, and JSON file, without
replacing any source file. All 13 TIFFs were independently re-encoded in a
temporary directory to measure their resulting sizes. Pixel equality was
verified for a small TIFF and the 1.93 GB pancreas TIFF; the production build
must check all 13 images and preserve required metadata. SHA-256 duplicate
scanning of NPY, CSV, and JSON found only 0.049 GB of duplicate copies, too
little to justify a special cross-object deduplication scheme.

The compression plan places gzipped NPY, CSV, and JSON at distinct R2 keys
ending in `.gz`, without an HTTP `Content-Encoding` header. This avoids a
Cloudflare/R2 behavior reported by rclone in which gzip-encoded objects may be
decompressed during transfer. With `VITE_DATA_COMPRESSION=gzip`, the frontend
uses the browser's `DecompressionStream` to recover the original NPY/CSV/JSON
bytes. The TIFFs remain ordinary lossless TIFF files at the same keys. Merely
enabling CDN/network compression does not reduce R2 storage.

## Deployment checks and limits

1. Build a separate staged payload; keep the current `web_data` untouched.
2. Validate gzip round trips against the source and all TIFF pixels/important
   tags. Compare the staged inventory with the static manifest paths.
3. Upload only the staged bytes to a fresh bucket and compare remote sizes
   with the staged files. Do not put both original and compressed copies in
   the same bucket; do not set `Content-Encoding: gzip` on these `.gz` objects.
4. Confirm public `Content-Type` and CORS headers and the absence of
   `Content-Encoding: gzip` for `.gz` objects. Run the
   browser smoke flow against the real R2 URL: catalog, NPY image, OOS view,
   raw-resolution view, CSV download, and TIFF download.

At 9.470 GB, the staged payload is about 0.530 GB below R2 Standard's monthly
10 GB-month free allowance if the account has no other R2 storage. This is a
measured **payload estimate**, not a verified R2 bill; check actual bucket
usage after upload. Large compressed TIFFs still exceed Cloudflare Free's
512 MB CDN cacheable-file limit, so their full-resolution downloads will be
uncached. Avoid converting arrays to float32 or TIFFs to lossy JPEG if exact
numeric/pixel preservation is required.
