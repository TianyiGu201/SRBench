#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 WEB_DATA_DIR STAGED_PAYLOAD_DIR RCLONE_REMOTE:BUCKET [--apply]" >&2
  exit 2
fi

web_data_dir=${1%/}
stage_dir=${2%/}
destination=${3%/}
mode=${4:---plan}
if [[ "$destination" != *:* || "$destination" == *: ]]; then
  echo "Destination must include an rclone remote and bucket, such as r2:srbench-data" >&2
  exit 2
fi
if [[ "$mode" != --plan && "$mode" != --apply ]]; then
  echo "Fourth argument must be --apply or omitted for a plan" >&2
  exit 2
fi
if [[ ! -d "$stage_dir/data_root" || ! -d "$stage_dir/site" ]]; then
  echo "Missing staged data_root or site: $stage_dir" >&2
  exit 2
fi

echo "Target: $destination; staged bytes will be uploaded without Content-Encoding."
echo "No remote objects will be deleted. Use an empty bucket to stay within the free storage allowance."
if [[ "$mode" == --plan ]]; then
  du -sh "$stage_dir/data_root" "$stage_dir/site"
  echo "Plan only. Run again with --apply after configuring the R2 rclone remote."
  exit 0
fi

command -v rclone >/dev/null
python3 "$(dirname "$0")/build_optimized_payload.py" verify \
  --web-data "$web_data_dir" --output "$stage_dir"
rclone copy "$stage_dir/data_root" "$destination/data_root" \
  --transfers 8 --checkers 16 --s3-chunk-size 64M --progress
rclone copy "$stage_dir/site" "$destination/site" \
  --transfers 8 --checkers 16 --progress
rclone check "$stage_dir/data_root" "$destination/data_root" --one-way --size-only
rclone check "$stage_dir/site" "$destination/site" --one-way --size-only
rclone size "$destination"
echo "Optimized R2 upload and size checks complete."
