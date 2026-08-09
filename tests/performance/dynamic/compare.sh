#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 CPU_LOG GPU_LOG CPU_STATOUT GPU_STATOUT" >&2
  exit 2
fi

script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
tmp_dir=$(mktemp -d /tmp/jetspin-dynamic-compare.XXXXXX)
trap 'rm -rf "$tmp_dir"' EXIT
rg '^Topology event:' "$1" > "$tmp_dir/cpu.events"
rg '^Topology event:' "$2" > "$tmp_dir/gpu.events"
diff -u "$script_dir/topology-events.txt" "$tmp_dir/cpu.events"
diff -u "$script_dir/topology-events.txt" "$tmp_dir/gpu.events"
python3 "$script_dir/../../regression/compare_statout.py" \
  --rtol "${RTOL:-3e-4}" --atol "${ATOL:-1e-6}" "$3" "$4"
