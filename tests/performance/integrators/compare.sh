#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 {euler|rk2} CPU_STATOUT [GPU_STATOUT]" >&2
  exit 2
fi

case "$1" in
  euler) prefix=euler ;;
  rk2) prefix=rk2 ;;
  *) echo "Unknown integrator: $1" >&2; exit 2 ;;
esac

script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
compare="$script_dir/../../regression/compare_statout.py"
rtol=${RTOL:-1e-6}
atol=${ATOL:-1e-9}

python3 "$compare" --rtol "$rtol" --atol "$atol" \
  "$script_dir/${prefix}-nvfortran-cpu.statout" "$2"

if [[ $# -eq 3 ]]; then
  python3 "$compare" --rtol "$rtol" --atol "$atol" \
    "$script_dir/${prefix}-openacc-a30.statout" "$3"
  python3 "$compare" --rtol "$rtol" --atol "$atol" "$2" "$3"
fi
