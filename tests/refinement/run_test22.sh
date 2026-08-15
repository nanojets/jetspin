#!/bin/sh

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
backend=${1:-nvfortran}
case "$backend" in
    nvfortran)
        build_target=nvfortran
        build_label="NVFORTRAN CPU reference"
        ;;
    openacc)
        build_target=nvfortran-openacc
        build_label="NVFORTRAN OpenACC"
        ;;
    force-oracle)
        build_target=nvfortran-openacc-force-oracle
        build_label="NVFORTRAN OpenACC host-force oracle"
        ;;
    *)
        echo "Usage: $0 [nvfortran|openacc|force-oracle]" >&2
        exit 2
        ;;
esac

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/jetspin-refinement-stress.XXXXXX")
if [ "${JETSPIN_REFINEMENT_KEEP:-0}" = 1 ]; then
    echo "Keeping Test 22 files in $work_dir"
else
    trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
fi

mkdir -p "$work_dir/source" "$work_dir/execute"
cp "$repo_root"/source/*.f90 "$work_dir/source/"
cp "$repo_root/build/Makefile" "$work_dir/source/Makefile"

echo "Building the Test 22 $build_label executable"
make -C "$work_dir/source" "$build_target" \
    GPUCC="${GPUCC:-80}" BINROOT="$work_dir/execute"

cp "$repo_root/examples/input-22/input.dat" "$work_dir/execute/input.dat"
echo "Running Test 22"
(
    cd "$work_dir/execute"
    JETSPIN_REFINEMENT_INITIAL_RESERVE=20 \
      JETSPIN_REFINEMENT_GROWTH_INCREMENT=20 \
      timeout "${JETSPIN_REFINEMENT_TIMEOUT:-180}" ./main.x > run.log 2>&1
)

python3 "$repo_root/tests/refinement/check_test22.py" \
    --run-log "$work_dir/execute/run.log" \
    --statout "$work_dir/execute/statout.dat"

if [ "$backend" = openacc ] || [ "$backend" = force-oracle ]; then
    rebind_count=$(grep -c "OpenACC refinement capacity rebind:" \
        "$work_dir/execute/run.log" || true)
    if [ "$rebind_count" -ne 3 ]; then
        echo "Test 22 validation failed: expected 3 OpenACC capacity rebinds, got $rebind_count" >&2
        exit 1
    fi
fi
