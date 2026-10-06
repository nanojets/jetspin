#!/bin/sh
# Exact restart of the dynamic stochastic Platen runs (Tests 24 and 25,
# 2026-10-06): for each test, run A stops at step 100,000 and writes
# save.dat, run B restarts from it and stops at step 200,000, and run C runs
# the same 200,000 steps without interruption.  Every output row and frame
# of B after step 100,000 must equal C's.
#
# Usage: tests/restart/run.sh [nvfortran|openacc]

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
backend=${1:-nvfortran}
case "$backend" in
    nvfortran) build_target=nvfortran ;;
    openacc) build_target=nvfortran-openacc ;;
    *)
        echo "Usage: $0 [nvfortran|openacc]" >&2
        exit 2
        ;;
esac

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/jetspin-restart.XXXXXX")
if [ "${JETSPIN_RESTART_KEEP:-0}" = 1 ]; then
    echo "Keeping restart test files in $work_dir"
else
    trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
fi

mkdir -p "$work_dir/source" "$work_dir/bin"
cp "$repo_root"/source/*.f90 "$work_dir/source/"
cp "$repo_root/build/Makefile" "$work_dir/source/Makefile"
echo "Building the $build_target executable"
make -C "$work_dir/source" "$build_target" GPUCC="${GPUCC:-80}" \
    BINROOT="$work_dir/bin" > "$work_dir/build.log" 2>&1 || {
    tail -20 "$work_dir/build.log" >&2
    exit 1
}

# make_input TEST FINALTIME RESTART: the test's input with 100,000-step
# restart dumps, a 20,000-step print interval and the given final time;
# RESTART=yes also reads restart.dat.
make_input() {
    restart_line=""
    if [ "$3" = yes ]; then
        restart_line="\\
 restart yes"
    fi
    sed -e "s/^ final time .*/ final time $2/" \
        -e "s/^ print time .*/ print time 1.d-4/" \
        -e "s/^ seed \(.*\)/ seed \1\\
 restart dump 100000$restart_line/" "$repo_root/examples/input-$1/input.dat"
}

run() {
    (
        cd "$1"
        timeout "${JETSPIN_RESTART_TIMEOUT:-600}" "$work_dir/bin/main.x" \
            > run.log 2>&1
    )
}

status=0
for test in 24 25; do
    for name in a b c; do
        mkdir -p "$work_dir/t$test-$name"
    done
    make_input $test 5.d-4 no > "$work_dir/t$test-a/input.dat"
    make_input $test 1.d-3 yes > "$work_dir/t$test-b/input.dat"
    make_input $test 1.d-3 no > "$work_dir/t$test-c/input.dat"
    echo "Test $test: run A (100,000 steps)"
    run "$work_dir/t$test-a"
    cp "$work_dir/t$test-a/save.dat" "$work_dir/t$test-b/restart.dat"
    echo "Test $test: run B (restart to 200,000 steps) and run C (uninterrupted)"
    run "$work_dir/t$test-b" &
    run "$work_dir/t$test-c"
    wait
    if ! python3 "$repo_root/tests/restart/check_restart.py" \
        "$work_dir/t$test-b" "$work_dir/t$test-c" 100000; then
        status=1
    fi
done
exit $status
