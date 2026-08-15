#!/bin/sh

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
backend=${1:-nvfortran}
accelerator_backend=0
case "$backend" in
    nvfortran)
        build_target=nvfortran
        build_label="NVFORTRAN CPU reference"
        ;;
    openacc)
        build_target=nvfortran-openacc
        build_label="NVFORTRAN OpenACC"
        accelerator_backend=1
        ;;
    force-oracle)
        build_target=nvfortran-openacc-force-oracle
        build_label="NVFORTRAN OpenACC host-force oracle"
        accelerator_backend=1
        ;;
    host-akima)
        build_target=nvfortran-openacc-host-akima
        build_label="NVFORTRAN OpenACC host-Akima oracle"
        accelerator_backend=1
        ;;
    akima-compare)
        build_target=nvfortran-openacc-compare-akima
        build_label="NVFORTRAN OpenACC Akima A/B comparison"
        accelerator_backend=1
        ;;
    refinement-compare)
        build_target=nvfortran-openacc-compare-refinement
        build_label="NVFORTRAN OpenACC refinement-assembly A/B comparison"
        accelerator_backend=1
        ;;
    *)
        echo "Usage: $0 [nvfortran|openacc|force-oracle|host-akima|akima-compare|refinement-compare]" >&2
        exit 2
        ;;
esac

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/jetspin-refinement-removal.XXXXXX")
if [ "${JETSPIN_REFINEMENT_KEEP:-0}" = 1 ]; then
    echo "Keeping Test 23 files in $work_dir"
else
    trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
fi

mkdir -p "$work_dir/source" "$work_dir/execute"
cp "$repo_root"/source/*.f90 "$work_dir/source/"
cp "$repo_root/build/Makefile" "$work_dir/source/Makefile"

echo "Building the Test 23 $build_label executable"
make -C "$work_dir/source" "$build_target" \
    GPUCC="${GPUCC:-80}" BINROOT="$work_dir/execute"

cp "$repo_root/examples/input-23/input.dat" "$work_dir/execute/input.dat"
echo "Running Test 23"
(
    cd "$work_dir/execute"
    JETSPIN_REFINEMENT_INITIAL_RESERVE=20 \
      JETSPIN_REFINEMENT_GROWTH_INCREMENT=20 \
      timeout "${JETSPIN_REFINEMENT_TIMEOUT:-240}" ./main.x > run.log 2>&1
)

akima_check_arg=
if [ "$backend" = akima-compare ]; then
    akima_check_arg=--require-akima-comparison
fi
python3 "$repo_root/tests/refinement/check_test23.py" \
    --run-log "$work_dir/execute/run.log" \
    --statout "$work_dir/execute/statout.dat" $akima_check_arg

if [ "$accelerator_backend" -eq 1 ]; then
    rebind_count=$(grep -c "OpenACC refinement capacity rebind:" \
        "$work_dir/execute/run.log" || true)
    if [ "$rebind_count" -ne 3 ]; then
        echo "Test 23 validation failed: expected 3 OpenACC capacity rebinds, got $rebind_count" >&2
        exit 1
    fi
fi

if [ "$backend" = akima-compare ]; then
    comparison_count=$(grep -c "Akima device comparison:" \
        "$work_dir/execute/run.log" || true)
    if [ "$comparison_count" -ne 33 ]; then
        echo "Test 23 validation failed: expected 33 Akima A/B comparisons, got $comparison_count" >&2
        exit 1
    fi
fi

if [ "$backend" = refinement-compare ]; then
    comparison_count=$(grep -c "Refinement assembly device comparison:" \
        "$work_dir/execute/run.log" || true)
    if [ "$comparison_count" -ne 3 ]; then
        echo "Test 23 validation failed: expected 3 refinement-assembly A/B comparisons, got $comparison_count" >&2
        exit 1
    fi
    evaporation_comparison_count=$(grep -c \
        "Refinement assembly evaporation comparison:" \
        "$work_dir/execute/run.log" || true)
    if [ "$evaporation_comparison_count" -ne 3 ]; then
        echo "Test 23 validation failed: expected 3 refinement-assembly evaporation comparisons, got $evaporation_comparison_count" >&2
        exit 1
    fi
fi
