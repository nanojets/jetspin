#!/bin/sh

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/jetspin-smoke.XXXXXX")
if [ "${JETSPIN_SMOKE_KEEP:-0}" = 1 ]; then
    echo "Keeping smoke-test files in $work_dir"
else
    trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
fi

mode=${1:-serial}
case "$mode" in
    serial)
        build_target=gfortran
        last_case=7
        ;;
    debug)
        build_target=gfortran-debugger
        last_case=7
        ;;
    mpi)
        build_target=gfortran-mpi
        last_case=1
        ;;
    *)
        echo "Usage: $0 [serial|debug|mpi]" >&2
        exit 2
        ;;
esac

mkdir -p "$work_dir/source" "$work_dir/execute"
cp "$repo_root"/source/*.f90 "$work_dir/source/"
cp "$repo_root/build/Makefile" "$work_dir/source/Makefile"

echo "Building JETSPIN executable ($mode)"
if [ "$mode" = mpi ]; then
    if printf 'end\n' | mpif90 -fallow-argument-mismatch -x f95 \
        -c -o "$work_dir/mpi-flag-test.o" - >/dev/null 2>&1; then
        mpi_compat_flag=-fallow-argument-mismatch
    else
        mpi_compat_flag=-Wno-argument-mismatch
    fi
    make -C "$work_dir/source" "$build_target" \
        BINROOT="$work_dir/execute" MPI_COMPAT_FLAG="$mpi_compat_flag"
else
    make -C "$work_dir/source" "$build_target" BINROOT="$work_dir/execute"
fi

case_number=1
while [ "$case_number" -le "$last_case" ]; do
    case_dir="$work_dir/case-$case_number"
    mkdir -p "$case_dir"
    cp "$work_dir/execute/main.x" "$case_dir/main.x"
    cp "$repo_root/examples/input-$case_number/input.dat" "$case_dir/input.dat"

    # Normalize historical DOS line endings and shorten the run to 1,000
    # integration steps while retaining every case-specific directive.
    sed -i 's/\r$//' "$case_dir/input.dat"
    sed -i \
        -e 's/^[[:space:]]*timestep[[:space:]].*/ timestep 1.d-6/' \
        -e 's/^[[:space:]]*final time[[:space:]].*/ final time 1.d-3/' \
        -e 's/^[[:space:]]*print time[[:space:]].*/ print time 2.d-4/' \
        "$case_dir/input.dat"

    echo "Running smoke case $case_number"
    (
        cd "$case_dir"
        if [ "$mode" = mpi ]; then
            timeout "${JETSPIN_SMOKE_TIMEOUT:-30}" \
                mpirun -np 2 ./main.x > run.log 2>&1
        else
            timeout "${JETSPIN_SMOKE_TIMEOUT:-30}" \
                ./main.x > run.log 2>&1
        fi
    )

    if ! grep -q 'Program closed correctly' "$case_dir/run.log"; then
        echo "Case $case_number did not close correctly" >&2
        tail -40 "$case_dir/run.log" >&2
        exit 1
    fi

    if grep -Eiq '(^|[^[:alpha:]])(error|nan|[-+]?inf(inity)?)([^[:alpha:]]|$)' \
        "$case_dir/run.log" "$case_dir/statout.dat"; then
        echo "Case $case_number produced an error or non-finite value" >&2
        tail -40 "$case_dir/run.log" >&2
        exit 1
    fi

    data_rows=$(awk '!/^#/ && NF { rows++ } END { print rows + 0 }' \
        "$case_dir/statout.dat")
    if [ "$data_rows" -lt 2 ]; then
        echo "Case $case_number produced only $data_rows numerical rows" >&2
        exit 1
    fi

    if ! awk '
        /^#/ || !NF { next }
        !seen { for (i = 2; i <= NF; i++) first[i] = $i; seen = 1; next }
        { for (i = 2; i <= NF; i++) if ($i != first[i]) changed = 1 }
        END { exit changed ? 0 : 1 }
    ' "$case_dir/statout.dat"; then
        echo "Case $case_number produced no evolving numerical observable" >&2
        exit 1
    fi

    case "$case_number" in
        2)
            grep -q 'inserting yes - inserting mode' "$case_dir/run.log"
            ;;
        5)
            grep -q 'dynamic refinement yes' "$case_dir/run.log"
            ;;
        7)
            grep -q 'external potential type.*3' "$case_dir/run.log"
            grep -q 'evaporation yes' "$case_dir/run.log"
            ;;
    esac

    echo "Case $case_number passed ($data_rows numerical rows)"
    case_number=$((case_number + 1))
done

echo "Smoke mode $mode passed ($last_case case(s))"
