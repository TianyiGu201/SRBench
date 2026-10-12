# GitHub Pages frontend + Cloudflare R2 data

`TianyiGu201/SRBench` can publish its frontend at
`https://tianyigu201.github.io/SRBench/` once the repository is public and
GitHub Pages is enabled. The existing Cloudflare Pages site at
`https://srbench.pages.dev/` can remain available. Both sites use the same
public R2 bucket; no research data is committed to GitHub.

The static adapter already maps website requests to two bucket prefixes:

| Website request | R2 object key |
| --- | --- |
| Dataset catalog and generated summaries | `site/...*.json.gz` |
| NPY and CSV used by the browser | `data_root/...*.npy.gz`, `data_root/...*.csv.gz` |
| JPEG previews and lossless TIFF downloads | `data_root/...*.jpg`, `data_root/...*.tiff` |

The data bucket is separate from GitHub. No R2 access key belongs in this
repository or in the frontend build: `VITE_DATA_BASE_URL` is a public URL
embedded in browser JavaScript.

## 1. Validate the local payload

Run from the repository root:

```bash
python3 scripts/audit_static_payload.py \
  --web-data /mnt/disk30/gty/project/super_resolution/web_data
```

The current export has 13 datasets. All required references exist; 44 optional
repeat CSVs are absent in the source data and appear as empty results, matching
the original API. If source data changes, rerun `scripts/export_static_data.py`
before the audit and upload.

## 2. Build the optimized payload on this server

The frontend uses `VITE_DATA_COMPRESSION=gzip` to fetch `.gz` objects and
decompress them in the browser. The `.gz` objects are uploaded as ordinary
files **without** HTTP `Content-Encoding`; Cloudflare may transparently
expand objects that use that header. The source `web_data` is not modified.

```bash
python3 scripts/build_optimized_payload.py build \
  --web-data /mnt/disk30/gty/project/super_resolution/web_data \
  --output /mnt/disk30/gty/project/super_resolution/web_data/r2_payload_gzip_v1
python3 scripts/build_optimized_payload.py verify \
  --web-data /mnt/disk30/gty/project/super_resolution/web_data \
  --output /mnt/disk30/gty/project/super_resolution/web_data/r2_payload_gzip_v1
```

The builder refuses to write into an existing output directory. The verifier
checks every gzip round trip, every TIFF's pixel values and important tags,
all JPEG copies, and the full published file inventory. The measured target
is 20,732 files and about 9.47 GB. See
[`compression-assessment.md`](compression-assessment.md).

## 3. Create and fill an R2 bucket

First activate the R2 subscription in Cloudflare if the dashboard requests it;
the account can remain on the Free plan, but storage beyond the monthly free
allowance is billed. Create a **Standard** storage bucket such as
`srbench-data`. Leave public access disabled while uploading and checking
files. In **R2 > Overview > API Tokens > Manage**, create a bucket-scoped
**Object Read & Write** API token. Configure the server's installed `rclone`
with S3 provider `Cloudflare`, region `auto`, the R2 account endpoint
`https://<account-id>.r2.cloudflarestorage.com`, and the API key. Keep the
credentials in the server's private rclone configuration; do not paste them
into GitHub or browser environment variables.

Run `rclone config` on the server, create a new remote named `r2`, choose
**Amazon S3 Compliant Storage Providers** and **Cloudflare R2**, then enter
the Access Key ID, Secret Access Key, and endpoint shown by Cloudflare. Use
`rclone listremotes` to confirm that `r2:` appears. If you used an object-only
token and rclone reports a bucket-check permission error, set
`no_check_bucket = true` for that remote.

Check the **optimized** upload plan, then apply it from this repository root:

```bash
bash scripts/upload_r2_optimized.sh \
  /mnt/disk30/gty/project/super_resolution/web_data \
  /mnt/disk30/gty/project/super_resolution/web_data/r2_payload_gzip_v1 \
  r2:srbench-data
bash scripts/upload_r2_optimized.sh \
  /mnt/disk30/gty/project/super_resolution/web_data \
  /mnt/disk30/gty/project/super_resolution/web_data/r2_payload_gzip_v1 \
  r2:srbench-data --apply
```

The script verifies local staged bytes, uploads `data_root` then `site`, never
deletes remote objects, and compares staged and remote file sizes. Use an
empty bucket to avoid paying for duplicate original and compressed objects.
The older `upload_r2.sh` is retained only as an uncompressed baseline.

## 4. Choose access mode

**Public website and data:** For the first end-to-end test, enable the bucket's
Public Development URL in **R2 > bucket > Settings**, and copy its
`https://...r2.dev` URL. It is rate limited and has no CDN caching, so use a
custom domain such as `https://data.example.org` for production. A custom
domain requires a domain you control in Cloudflare. Set the bucket CORS policy
from `r2-cors.example.json`. It lists both current frontend origins:
`https://srbench.pages.dev` and `https://tianyigu201.github.io`. The GitHub
origin does not include `/SRBench/`; CORS origins contain only the scheme and
hostname. In the bucket's **Settings > CORS Policy > Add CORS policy > JSON**
editor, paste the policy. Check public access and CORS, for example:

```bash
curl -I https://<your-r2-public-host>/site/datasets.json.gz
curl -I -H 'Origin: https://tianyigu201.github.io' \
  https://<your-r2-public-host>/site/datasets.json.gz
```

The second response must include
`Access-Control-Allow-Origin: https://tianyigu201.github.io`.
Configure caching for the data custom domain after correctness checks. The
`.json.gz` objects are ordinary gzip files, so set their cache rules based on
the public data URL and verify the response headers. On the Free plan, files
above 512 MB cannot be cached by Cloudflare, so the largest TIFFs will remain
cold reads.

**Private website or data:** Ordinary GitHub Pages and a public R2 domain do
not meet this requirement. Use a frontend host with access control and a
private R2 bucket behind an authenticated Worker or equivalent proxy. The
current static frontend would then need an authenticated data transport.

## 5. Publish the frontend with GitHub Pages

After checking the source and Git history for credentials, change the
`TianyiGu201/SRBench` repository visibility to **Public** in **Settings >
General > Danger Zone > Change repository visibility**. GitHub Free cannot
publish Pages directly from a private repository. The repository includes
website source and deployment scripts, but no R2 credentials or research data.

In **Settings > Pages > Build and deployment**, set **Source** to **GitHub
Actions**. The committed `.github/workflows/pages.yml` installs frontend
dependencies, checks TypeScript, builds with `VITE_BASE_PATH=/SRBench/`, the
public R2 URL, and `VITE_DATA_COMPRESSION=gzip`, then uploads only `front/dist`
to GitHub Pages. Trigger **Actions > Deploy SRBench to GitHub Pages > Run
workflow** once if changing visibility does not start a deployment. The
expected URL is `https://tianyigu201.github.io/SRBench/`.

Before the browser test, add `https://tianyigu201.github.io` to R2 CORS while
keeping the existing `https://srbench.pages.dev` origin. The frontend uses
hash routing, so links under `/#/` require no rewrite rule. Verify the dataset
list, an NPY image, out-of-sample regions, CSV download, and TIFF download.
The R2 `r2.dev` hostname is suitable for testing; use an R2 custom domain for
production caching and to avoid development-endpoint rate limits.

## Existing Cloudflare Pages deployment

`https://srbench.pages.dev/` is already connected to the same GitHub repository
and R2 bucket. Its build root is `front`, command is `npm run build`, output is
`dist`, and build variables are `VITE_BASE_PATH=/`, the same public R2 URL,
and `VITE_DATA_COMPRESSION=gzip`. Keep its CORS origin while both sites run.

## Cost estimate for this payload

The baseline upload scope is 28.26 GB, while the verified optimized target is
about 9.47 GB. If this bucket has no other objects, the optimized payload fits
R2 Standard's 10 GB-month free storage allowance with about 0.53 GB of room.
The first upload is far below the included one million Class A operations;
10 million Class B reads are included monthly. Standard retrieval and outbound
transfer have no R2 charge. GitHub Pages is available from this public
repository on GitHub Free, subject to GitHub Pages usage limits. Cloudflare
Pages Free includes 500 builds/month. A custom domain, if purchased, has a
separate registration cost. These are estimates, not a measured bill.
