# Private local Pages + R2 simulation

This preview uses two loopback origins. The frontend is a static Vite build under
`/spatial-transcriptomics-viz/`; the read-only data origin serves the original
`web_data/data_root` and generated `web_data/site` files. No public bucket or
GitHub Pages site is involved.

## Generate static data

Run from this repository root with the Python environment installed for the backend:

```bash
.venv-server/bin/python scripts/export_static_data.py \
  --data-root /mnt/disk30/gty/project/super_resolution/web_data/data_root \
  --output /mnt/disk30/gty/project/super_resolution/web_data/site
```

The exporter creates the dataset catalog, out-of-sample region tables, and
different-resolution plot summaries. Original NPY, CSV, JPG, TIFF, and metadata
JSON files are served directly from `data_root`.

## Build and preview

Run from `front/`:

```bash
export PATH=/mnt/disk10/gty/tools/node-v24.17.0-linux-x64/bin:$PATH
VITE_DATA_BASE_URL=http://127.0.0.1:7801 \
VITE_BASE_PATH=/spatial-transcriptomics-viz/ npm run build
```

Create the Pages-style path and start its server:

```bash
mkdir -p /mnt/disk30/gty/project/super_resolution/web_data/pages_preview
ln -sfn /mnt/disk30/gty/project/super_resolution/spatial-transcriptomics-viz/front/dist \
  /mnt/disk30/gty/project/super_resolution/web_data/pages_preview/spatial-transcriptomics-viz
python3 -m http.server 7800 --bind 127.0.0.1 \
  --directory /mnt/disk30/gty/project/super_resolution/web_data/pages_preview
```

In a second terminal, start the read-only data server from this repository root:

```bash
python3 scripts/serve_static_data.py \
  --root /mnt/disk30/gty/project/super_resolution/web_data \
  --port 7801 --origin http://127.0.0.1:7800
```

From a workstation, forward both ports over SSH, replacing `<server-alias>`
with the SSH host name you use for this server:

```bash
ssh -L 7800:127.0.0.1:7800 -L 7801:127.0.0.1:7801 <server-alias>
```

Open `http://127.0.0.1:7800/spatial-transcriptomics-viz/` on that workstation.

The frontend also supports the original FastAPI mode when `VITE_DATA_BASE_URL`
is unset. For a real R2 deployment, upload `data_root` and `site` as bucket
prefixes, use an HTTPS custom domain as `VITE_DATA_BASE_URL`, and configure R2
CORS for the Pages origin. The local preview does not measure public network
speed or validate Cloudflare account configuration.
