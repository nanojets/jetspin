#!/bin/sh
# Exact restart: for each case, run A stops at the restart step and writes
# save.dat, run B restarts from it, and run C runs the same steps without
# interruption.  Every output row and frame of B after the restart step
# must equal C's.
#
# Cases (2026-10-06): the dynamic stochastic Platen runs (Tests 24 and 25),
# the Platen scheme on a fixed jet (Test 12), RK4 with insertion and
# removal (Test 13), with evaporation (Test 16) and with the Kelvin-Voigt
# model and evaporation (Test 17), the Platen scheme without refinement
# (Example 4) and RK4 with evaporation restarted just after a compaction
# (Example 8).  In the OpenACC build Tests 12, 13, 16 and 17 run on the
# device step from the first step, and Example 8 restarts from a state the
# device step had engaged with fewer than 100 beads in the arrays.  The
# MPI backends run every copy on two ranks ($MPIEXEC -np 2, the executable
# built with $JETSPIN_MPIFC), as tests/regression/run.sh does.
#
# Usage: tests/restart/run.sh
#        [nvfortran|openacc|gfortran|gfortran-mpi|nvfortran-mpi] [CASE...]

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
backend=${1:-nvfortran}
mpi=no
case "$backend" in
    nvfortran|gfortran) build_target=$backend ;;
    openacc) build_target=nvfortran-openacc ;;
    gfortran-mpi|nvfortran-mpi) build_target=$backend mpi=yes ;;
    *)
        echo "Usage: $0 [nvfortran|openacc|gfortran|gfortran-mpi|nvfortran-mpi] [CASE...]" >&2
        exit 2
        ;;
esac
[ $# -gt 0 ] && shift
mpiexec_command=${MPIEXEC:-mpirun}
mpi_fc=${JETSPIN_MPIFC:-mpif90}

# NAME INPUT RESTART_STEP FINAL_TIME_A FINAL_TIME_BC PRINT_TIME XYZ_TIME:
# final times half a step short of the step count (Tests 24 and 25: exact
# multiples), so that rounding could not add a step before the step count
# became final time / timestep rounded up (2026-10-07); the restart step is
# a multiple of the print interval.
case_table='
t24 24 100000 5.d-4 1.d-3 1.d-4 1.d-4
t25 25 100000 5.d-4 1.d-3 1.d-4 1.d-4
t12 12 500 4.995d-6 9.995d-6 1.d-6 1.d-6
t13 13 500 4.995d-5 9.995d-5 1.d-5 1.d-5
t16 16 500 4.995d-5 9.995d-5 1.d-5 1.d-5
t17 17 500 4.995d-5 9.995d-5 1.d-5 1.d-5
ex4 4 10000 9.9995d-5 1.99995d-4 1.d-6 1.d-5
ex8 8 14000 1.39995d-4 1.99995d-4 1.d-6 1.d-5
'
selected=${*:-t24 t25 t12 t13 t16 t17 ex4 ex8}
for name in $selected; do
    if ! printf '%s\n' "$case_table" | grep -q "^$name "; then
        echo "Unknown case $name" >&2
        exit 2
    fi
done

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/jetspin-restart.XXXXXX")
if [ "${JETSPIN_RESTART_KEEP:-0}" = 1 ]; then
    echo "Keeping restart test files in $work_dir"
else
    trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
fi

mkdir -p "$work_dir/source" "$work_dir/bin"
if [ -n "${JETSPIN_RESTART_EXE:-}" ]; then
    echo "Using $JETSPIN_RESTART_EXE"
    cp "$JETSPIN_RESTART_EXE" "$work_dir/bin/main.x"
else
    cp "$repo_root"/source/*.f90 "$work_dir/source/"
    cp "$repo_root/build/Makefile" "$work_dir/source/Makefile"
    mpi_compat_flag=
    if [ "$backend" = gfortran-mpi ]; then
        if printf 'end\n' | "$mpi_fc" -fallow-argument-mismatch -x f95 \
            -c -o "$work_dir/mpi-flag-test.o" - >/dev/null 2>&1; then
            mpi_compat_flag=-fallow-argument-mismatch
        else
            mpi_compat_flag=-Wno-argument-mismatch
        fi
    fi
    echo "Building the $build_target executable"
    make -C "$work_dir/source" "$build_target" GPUCC="${GPUCC:-80}" \
        MPIFC="$mpi_fc" MPI_COMPAT_FLAG="$mpi_compat_flag" \
        BINROOT="$work_dir/bin" > "$work_dir/build.log" 2>&1 || {
        tail -20 "$work_dir/build.log" >&2
        exit 1
    }
fi

# make_input INPUT FINALTIME PRINTTIME XYZTIME DUMPSTEP RESTART: the
# input with the given final time, print and trajectory intervals and
# restart dumps; RESTART=yes also reads restart.dat.  Line ends are made
# Unix ones (the input of Example 4 has DOS ones, which the parser does not
# accept mixed with the added lines).
make_input() {
    restart_line=""
    if [ "$6" = yes ]; then
        restart_line="\\
 restart yes"
    fi
    sed -e 's/\r$//' \
        -e "/^ *print xyz *[0-9.]/d" \
        -e "s/^ final time .*/ final time $2/" \
        -e "s/^ print time .*/ print time $3/" \
        -e "s/^ seed \(.*\)/ seed \1\\
 print xyz $4\\
 restart dump $5$restart_line/" "$repo_root/examples/input-$1/input.dat"
}

run() {
    (
        cd "$1"
        if [ "$mpi" = yes ]; then
            timeout "${JETSPIN_RESTART_TIMEOUT:-600}" "$mpiexec_command" \
                -np 2 "$work_dir/bin/main.x" > run.log 2>&1
        else
            timeout "${JETSPIN_RESTART_TIMEOUT:-600}" "$work_dir/bin/main.x" \
                > run.log 2>&1
        fi
    )
}

status=0
for name in $selected; do
    set -- $(printf '%s\n' "$case_table" | grep "^$name ")
    input=$2 step=$3 final_a=$4 final_bc=$5 print_time=$6 xyz_time=$7
    for run_name in a b c; do
        mkdir -p "$work_dir/$name-$run_name"
    done
    make_input "$input" "$final_a" "$print_time" "$xyz_time" "$step" no \
        > "$work_dir/$name-a/input.dat"
    make_input "$input" "$final_bc" "$print_time" "$xyz_time" "$step" yes \
        > "$work_dir/$name-b/input.dat"
    make_input "$input" "$final_bc" "$print_time" "$xyz_time" "$step" no \
        > "$work_dir/$name-c/input.dat"
    echo "Case $name: run A ($step steps)"
    if ! run "$work_dir/$name-a" || [ ! -f "$work_dir/$name-a/save.dat" ]; then
        echo "Case $name: run A failed" >&2
        tail -5 "$work_dir/$name-a/run.log" >&2
        status=1
        continue
    fi
    cp "$work_dir/$name-a/save.dat" "$work_dir/$name-b/restart.dat"
    echo "Case $name: run B (restart) and run C (uninterrupted)"
    if [ "$mpi" = yes ]; then
        # One two-rank launch at a time, within the ranks allotted.
        run "$work_dir/$name-b" || true
        run "$work_dir/$name-c" || true
    else
        run "$work_dir/$name-b" &
        run_b=$!
        run "$work_dir/$name-c" || true
        wait $run_b || true
    fi
    for run_name in a b; do
        grep 'OpenACC device step' "$work_dir/$name-$run_name/run.log" |
            sed "s/^ */  run $run_name: /" || true
    done
    if ! python3 "$repo_root/tests/restart/check_restart.py" \
        "$work_dir/$name-b" "$work_dir/$name-c" "$step"; then
        status=1
    fi
done
exit $status
