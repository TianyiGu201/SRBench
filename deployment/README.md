# Private GitHub source + Cloudflare Pages frontend + R2 data

This is the deployment path for a private `TianyiGu201/SRBench` repository on
GitHub Free. Cloudflare Pages can build from a private GitHub repository and
publish the resulting frontend publicly. GitHub Pages from that private
repository would require GitHub Pro or another eligible paid plan.

The static adapter already maps website requests to two bucket prefixes:

| Website request | R2 object key |
| --- | --- |
| Dataset catalog and generated summaries | `site/...` |
| JPG, TIFF, NPY, CSV and original metadata | `data_root/...` |

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

## 2. Create and fill an R2 bucket

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

Check the upload plan, then apply it from this repository root:

```bash
bash scripts/upload_r2.sh \
  /mnt/disk30/gty/project/super_resolution/web_data r2:srbench-data
bash scripts/upload_r2.sh \
  /mnt/disk30/gty/project/super_resolution/web_data r2:srbench-data --apply
```

The script uploads `data_root` first and `site` second, excludes `.bak` backup
files, never deletes remote objects, and compares uploaded file sizes. The
current payload is 20,969 objects and 28.26 GB before upload.

## 3. Choose access mode

**Public website and data:** For the first end-to-end test, enable the bucket's
Public Development URL in **R2 > bucket > Settings**, and copy its
`https://...r2.dev` URL. It is rate limited and has no CDN caching, so use a
custom domain such as `https://data.example.org` for production. A custom
domain requires a domain you control in Cloudflare. Set the bucket CORS policy
from `r2-cors.example.json`, replacing the Pages origin with the actual
`https://<your-pages-project>.pages.dev` hostname. The origin contains only
the scheme and hostname; add a second origin if the frontend later uses a
custom domain. In the bucket's **Settings > CORS Policy > Add CORS policy >
JSON** editor, paste the edited policy. The example file uses the dashboard
JSON format. Check both public access and CORS, for example:

```bash
curl -I https://<your-r2-public-host>/site/datasets.json
curl -I -H 'Origin: https://<your-pages-project>.pages.dev' \
  https://<your-r2-public-host>/site/datasets.json
```

The second response must include
`Access-Control-Allow-Origin: https://<your-pages-project>.pages.dev`.
Configure caching for the data custom domain after correctness checks. JSON
is not cached by default. On the Free plan, files above 512 MB cannot be
cached by Cloudflare, so the largest TIFFs will remain cold reads.

**Private website or data:** Ordinary GitHub Pages and a public R2 domain do
not meet this requirement. Use a frontend host with access control and a
private R2 bucket behind an authenticated Worker or equivalent proxy. The
current static frontend would then need an authenticated data transport.

## 4. Publish the frontend from the private GitHub repository

In Cloudflare, select **Workers & Pages > Create application > Pages > Connect
to Git**, grant the Cloudflare Pages GitHub app access to `TianyiGu201/SRBench`,
and select that private repository. Configure:

| Setting | Value |
| --- | --- |
| Production branch | `main` |
| Root directory | `front` |
| Build command | `npm run build` |
| Build output directory | `dist` |
| Environment variable `VITE_BASE_PATH` | `/` |
| Environment variable `VITE_DATA_BASE_URL` | Public R2 URL without trailing slash |

The frontend uses hash routing, so URLs under `/#/` work on Pages without a
rewrite rule. Set the build variables before the first deployment. The R2 URL
is public browser configuration, not a secret. Cloudflare Pages will give the
site a `*.pages.dev` address; use that exact origin in the R2 CORS policy.

If the project name changes during setup, update the CORS origin. After
deployment, verify the dataset list, an NPY image, an out-of-sample region, a
CSV download, and a raw TIFF link from a browser outside the server. Measure
one cold load and one repeat load; local loopback timings do not represent
public network speed.

## Optional: GitHub Pages if the GitHub plan changes

The deploy workflow is staged at `pages-public.yml.example` and is not active.
If the account is upgraded to GitHub Pro or another eligible plan, move it to
`.github/workflows/pages-public.yml`, create the GitHub Actions repository
variable `R2_PUBLIC_BASE_URL` with the HTTPS data origin (no trailing slash),
and set **Settings > Pages > Build and deployment > Source** to
**GitHub Actions**. The workflow builds with `/SRBench/` as its base path and
deploys only `front/dist`.

Alternatively, a separate public repository containing only built assets can
serve GitHub Pages on GitHub Free while `SRBench` stays private.

## Cost estimate for this payload

The current upload scope is 28.26 GB. At the current R2 Standard price of
US$0.015 per GB-month with 10 GB-month included, full-month storage is about
US$0.29/month before billing rounding or other usage. The first upload is far
below the included one million Class A operations; 10 million Class B reads
are included monthly. Standard retrieval and outbound transfer have no R2
charge. Pages Free includes 500 builds/month. A custom domain, if purchased,
has a separate registration cost. These are estimates, not a measured bill.
