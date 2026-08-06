#!/bin/sh

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/jetspin-smoke.XXXXXX")
trap 'rm -rf "$work_dir"' EXIT HUP INT TERM

mkdir -p "$work_dir/source" "$work_dir/execute"
cp "$repo_root"/source/*.f90 "$work_dir/source/"
cp "$repo_root/build/Makefile" "$work_dir/source/Makefile"

echo "Building serial JETSPIN executable"
make -C "$work_dir/source" gfortran BINROOT="$work_dir/execute"

case_number=1
while [ "$case_number" -le 7 ]; do
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
        timeout "${JETSPIN_SMOKE_TIMEOUT:-30}" ./main.x > run.log 2>&1
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

    echo "Case $case_number passed ($data_rows numerical rows)"
    case_number=$((case_number + 1))
done

echo "All 7 smoke cases passed"
