#!/bin/sh

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
record_dir="$repo_root/tests/performance/test9"
compare="$repo_root/tests/regression/compare_statout.py"
rtol=${JETSPIN_TEST9_RTOL:-1e-6}
atol=${JETSPIN_TEST9_ATOL:-1e-9}

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 CPU_STATOUT [GPU_STATOUT]" >&2
    exit 2
fi

cpu_result=$1
python3 "$compare" --rtol "$rtol" --atol "$atol" \
    "$record_dir/baseline-nvfortran-cpu.statout" "$cpu_result"

if [ "$#" -eq 2 ]; then
    gpu_result=$2
    python3 "$compare" --rtol "$rtol" --atol "$atol" \
        "$record_dir/baseline-openacc-a30.statout" "$gpu_result"
    python3 "$compare" --rtol "$rtol" --atol "$atol" \
        "$cpu_result" "$gpu_result"
fi

echo "Test Case 9 numerical comparison passed (rtol=$rtol, atol=$atol)"
