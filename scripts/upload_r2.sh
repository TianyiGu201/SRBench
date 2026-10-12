#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 WEB_DATA_DIR RCLONE_REMOTE:BUCKET [--apply]" >&2
  exit 2
fi

web_data_dir=${1%/}
destination=${2%/}
mode=${3:---plan}
if [[ "$destination" != *:* || "$destination" == *: ]]; then
  echo "Destination must be an rclone remote and bucket, such as r2:srbench-data" >&2
  exit 2
fi
if [[ "$mode" != --plan && "$mode" != --apply ]]; then
  echo "Third argument must be --apply or omitted for a plan" >&2
  exit 2
fi

python3 "$(dirname "$0")/audit_static_payload.py" --web-data "$web_data_dir"

echo "Data target: $destination/data_root"
echo "Derived JSON target: $destination/site"
echo "Backup .bak files will be excluded; no remote objects will be deleted."
if [[ "$mode" == --plan ]]; then
  echo "Plan only. Run again with --apply after configuring the R2 rclone remote."
  exit 0
fi

command -v rclone >/dev/null
rclone copy "$web_data_dir/data_root" "$destination/data_root" \
  --exclude '*.bak' --transfers 8 --checkers 16 --s3-chunk-size 64M --progress
rclone copy "$web_data_dir/site" "$destination/site" \
  --transfers 8 --checkers 16 --progress
rclone check "$web_data_dir/data_root" "$destination/data_root" \
  --exclude '*.bak' --one-way --size-only
rclone check "$web_data_dir/site" "$destination/site" \
  --one-way --size-only
echo "R2 upload and size checks complete."
