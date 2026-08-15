#!/bin/sh

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
backend=${1:-gfortran}
mode=${2:-standard}
case "$backend" in
    gfortran)
        build_target=gfortran
        build_label="GFortran reference"
        ;;
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
        echo "Usage: $0 [gfortran|nvfortran|openacc|force-oracle] [standard|capacity-growth]" >&2
        exit 2
        ;;
esac
case "$mode" in
    standard)
        initial_reserve=100
        checker_mode=
        ;;
    capacity-growth)
        # Keep enough room for nozzle insertion before the event, but force
        # the accepted Akima mesh to exceed the original device capacity.
        initial_reserve=50
        checker_mode=--require-capacity-growth
        ;;
    *)
        echo "Usage: $0 [gfortran|nvfortran|openacc|force-oracle] [standard|capacity-growth]" >&2
        exit 2
        ;;
esac
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/jetspin-refinement.XXXXXX")
if [ "${JETSPIN_REFINEMENT_KEEP:-0}" = 1 ]; then
    echo "Keeping Test 21 files in $work_dir"
else
    trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
fi

mkdir -p "$work_dir/source" "$work_dir/execute"
cp "$repo_root"/source/*.f90 "$work_dir/source/"
cp "$repo_root/build/Makefile" "$work_dir/source/Makefile"

echo "Building the Test 21 $build_label executable"
make -C "$work_dir/source" "$build_target" \
    GPUCC="${GPUCC:-80}" BINROOT="$work_dir/execute"

cp "$repo_root/examples/input-21/input.dat" "$work_dir/execute/input.dat"
echo "Running Test 21 ($mode)"
(
    cd "$work_dir/execute"
    JETSPIN_REFINEMENT_INITIAL_RESERVE="$initial_reserve" \
      timeout "${JETSPIN_REFINEMENT_TIMEOUT:-90}" ./main.x > run.log 2>&1
)

python3 "$repo_root/tests/refinement/check_test21.py" \
    --run-log "$work_dir/execute/run.log" \
    --statout "$work_dir/execute/statout.dat" $checker_mode

if [ "$mode" = capacity-growth ] && \
   { [ "$backend" = openacc ] || [ "$backend" = force-oracle ]; }; then
    if ! grep -q "OpenACC refinement capacity rebind:" \
        "$work_dir/execute/run.log"; then
        echo "Test 21 validation failed: OpenACC capacity rebind is missing" >&2
        exit 1
    fi
fi
