# SRBench private deployment test

This repository contains the spatial transcriptomics comparison website source.
The 27 GiB research dataset and the generated static data catalog stay outside
GitHub. The frontend can read those files from a separate data origin when
`VITE_DATA_BASE_URL` is set at build time.

## Private GitHub test

The `Private static build` GitHub Actions workflow checks the TypeScript source,
builds the frontend with `/SRBench/` as its base path, and stores the static
output as a private workflow artifact for seven days. The artifact points at
`http://127.0.0.1:7801` so it can be tested against the local data server via
SSH port forwarding. It is not a public deployment. GitHub Pages is not enabled.

For the original local two-origin browser test and data export instructions,
see [README.pages-r2-local.md](README.pages-r2-local.md).

## Production configuration

For an eventual public Pages plus R2 deployment, build `front/` with:

```bash
npm ci
VITE_BASE_PATH=/SRBench/ VITE_DATA_BASE_URL=https://<data-domain> npm run build
```

The data origin must serve matching `/data_root/` and `/site/` paths and allow
the website origin through CORS. A public static deployment will make the site
and any publicly served data URLs accessible to others; a private hosted site
requires an access-control layer.
