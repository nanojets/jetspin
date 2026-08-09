#!/bin/sh

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/jetspin-regression.XXXXXX")
if [ "${JETSPIN_REGRESSION_KEEP:-0}" = 1 ]; then
    echo "Keeping regression files in $work_dir"
else
    trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
fi

backend=${1:-gfortran}
action=${2:-check}
case "$backend" in
    check|update)
        # Preserve the original `run.sh update` interface.
        action=$backend
        backend=gfortran
        ;;
esac
case "$backend" in
    gfortran)
        baseline_dir="$repo_root/tests/regression/baselines"
        serial_target=gfortran
        comparison_target=gfortran-mpi
        comparison_mode=mpi
        comparison_label="two-rank MPI"
        ;;
    nvfortran)
        baseline_dir="$repo_root/tests/regression/baselines/nvfortran"
        serial_target=nvfortran
        comparison_target=nvfortran-mpi
        comparison_mode=mpi
        comparison_label="two-rank MPI/NVFORTRAN"
        ;;
    openacc)
        baseline_dir="$repo_root/tests/regression/baselines/nvfortran"
        serial_target=nvfortran
        comparison_target=nvfortran-openacc
        comparison_mode=serial
        comparison_label="NVFORTRAN OpenACC"
        if [ "$action" = update ]; then
            echo "Use '$0 nvfortran update' to update NVFORTRAN CPU baselines" >&2
            exit 2
        fi
        ;;
    *)
        echo "Usage: $0 [gfortran|nvfortran|openacc] [check|update]" >&2
        exit 2
        ;;
esac
case "$action" in
    check|update) ;;
    *) echo "Usage: $0 [gfortran|nvfortran|openacc] [check|update]" >&2; exit 2 ;;
esac

rtol=${JETSPIN_REGRESSION_RTOL:-1e-7}
atol=${JETSPIN_REGRESSION_ATOL:-1e-10}
mpi_rtol=${JETSPIN_MPI_REGRESSION_RTOL:-1e-6}
mpi_atol=${JETSPIN_MPI_REGRESSION_ATOL:-1e-9}
openacc_rtol=${JETSPIN_OPENACC_REGRESSION_RTOL:-1e-6}
openacc_atol=${JETSPIN_OPENACC_REGRESSION_ATOL:-1e-9}
timeout_seconds=${JETSPIN_REGRESSION_TIMEOUT:-30}
mpiexec_command=${MPIEXEC:-mpirun}
mpi_fc=${JETSPIN_MPIFC:-mpif90}
first_case=${JETSPIN_REGRESSION_FIRST_CASE:-1}
last_case=${JETSPIN_REGRESSION_LAST_CASE:-8}
case "$first_case:$last_case" in
    *[!0-9:]*|:*|*:) echo "Invalid regression case range" >&2; exit 2 ;;
esac
if [ "$first_case" -lt 1 ] || [ "$last_case" -gt 8 ] || \
   [ "$first_case" -gt "$last_case" ]; then
    echo "Regression case range must be within 1..8" >&2
    exit 2
fi

mkdir -p "$work_dir/serial-source" "$work_dir/comparison-source"
cp "$repo_root"/source/*.f90 "$work_dir/serial-source/"
cp "$repo_root"/source/*.f90 "$work_dir/comparison-source/"
cp "$repo_root/build/Makefile" "$work_dir/serial-source/Makefile"
cp "$repo_root/build/Makefile" "$work_dir/comparison-source/Makefile"
mkdir -p "$work_dir/serial-bin" "$work_dir/comparison-bin"

echo "Building $backend serial CPU regression executable"
make -C "$work_dir/serial-source" "$serial_target" \
    BINROOT="$work_dir/serial-bin"
mpi_compat_flag=
if [ "$backend" = gfortran ]; then
    if printf 'end\n' | "$mpi_fc" -fallow-argument-mismatch -x f95 \
        -c -o "$work_dir/mpi-flag-test.o" - >/dev/null 2>&1; then
        mpi_compat_flag=-fallow-argument-mismatch
    else
        mpi_compat_flag=-Wno-argument-mismatch
    fi
fi
echo "Building $comparison_label regression executable"
make -C "$work_dir/comparison-source" "$comparison_target" \
    BINROOT="$work_dir/comparison-bin" MPI_COMPAT_FLAG="$mpi_compat_flag" \
    MPIFC="$mpi_fc" GPUCC="${GPUCC:-80}" \
    CUDA_VERSION="${CUDA_VERSION:-12.3}"

prepare_case() {
    case_number=$1
    destination=$2
    mkdir -p "$destination"
    cp "$repo_root/examples/input-$case_number/input.dat" "$destination/input.dat"
    sed -i 's/\r$//' "$destination/input.dat"
    timestep=$(awk '
        tolower($1) == "timestep" {
            value = $2
            gsub(/[dD]/, "e", value)
            printf "%.17e", value + 0
            exit
        }
    ' "$destination/input.dat")
    if [ -z "$timestep" ]; then
        echo "No timestep found for case $case_number" >&2
        exit 1
    fi
    final_time=$(awk -v value="$timestep" 'BEGIN { printf "%.17e", value * 1000 }')
    print_time=$(awk -v value="$timestep" 'BEGIN { printf "%.17e", value * 200 }')
    sed -i \
        -e "s/^[[:space:]]*final time[[:space:]].*/ final time $final_time/" \
        -e "s/^[[:space:]]*print time[[:space:]].*/ print time $print_time/" \
        -e 's/^[[:space:]]*printstat list[[:space:]].*/ printstat list t x y z vx vy vz st n rc curn curc nref/' \
        -e '/^[[:space:]]*printstat binary/d' \
        "$destination/input.dat"
    if [ "$case_number" -eq 4 ]; then
        # Force the stochastic Platen case above the 10-bead minimum chunk,
        # so both ranks consume their own part of the shared Gaussian block.
        sed -i 's/^[[:space:]]*points[[:space:]].*/ points 25/' \
            "$destination/input.dat"
        sed -i 's/^[[:space:]]*initial length[[:space:]].*/ initial length 0.2d0/' \
            "$destination/input.dat"
    fi
}

run_case() {
    executable=$1
    case_dir=$2
    mode=$3
    cp "$executable" "$case_dir/main.x"
    if [ "$mode" = mpi ]; then
        (cd "$case_dir" && timeout "$timeout_seconds" \
            "$mpiexec_command" -np 2 ./main.x > run.log 2>&1)
    else
        (cd "$case_dir" && timeout "$timeout_seconds" \
            ./main.x > run.log 2>&1)
    fi
    grep -q 'Program closed correctly' "$case_dir/run.log"
    test -s "$case_dir/statout.dat"
}

case_number=$first_case
while [ "$case_number" -le "$last_case" ]; do
    serial_dir="$work_dir/serial-$case_number"
    comparison_dir="$work_dir/comparison-$case_number"
    prepare_case "$case_number" "$serial_dir"
    prepare_case "$case_number" "$comparison_dir"
    echo "Running $backend serial CPU regression case $case_number"
    run_case "$work_dir/serial-bin/main.x" "$serial_dir" serial

    baseline="$baseline_dir/case-$case_number.statout"
    if [ "$action" = update ]; then
        mkdir -p "$baseline_dir"
        cp "$serial_dir/statout.dat" "$baseline"
        sed -i 's/[[:space:]]*$//' "$baseline"
        echo "Updated $baseline"
    else
        python3 "$repo_root/tests/regression/compare_statout.py" \
            --rtol "$rtol" --atol "$atol" "$baseline" "$serial_dir/statout.dat"
    fi

    echo "Running $comparison_label regression case $case_number"
    run_case "$work_dir/comparison-bin/main.x" "$comparison_dir" \
        "$comparison_mode"
    if [ "$backend" = openacc ]; then
        comparison_rtol=$openacc_rtol
        comparison_atol=$openacc_atol
    else
        comparison_rtol=$mpi_rtol
        comparison_atol=$mpi_atol
    fi
    python3 "$repo_root/tests/regression/compare_statout.py" \
        --rtol "$comparison_rtol" --atol "$comparison_atol" \
        "$serial_dir/statout.dat" "$comparison_dir/statout.dat"
    case_number=$((case_number + 1))
done

echo "Numerical regression suite passed for $backend cases $first_case..$last_case"
