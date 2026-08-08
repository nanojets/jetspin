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
        last_case=8
        ;;
    debug)
        build_target=gfortran-debugger
        last_case=8
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
        8)
            grep -q 'evaporation yes' "$case_dir/run.log"
            python3 "$repo_root/tests/evaporation/check_yarin2001.py" \
                --input "$repo_root/examples/input-8/input.dat" \
                --run-log "$case_dir/run.log" \
                --nanojet-source "$repo_root/source/nanojet_mod.f90" \
                --eom-source "$repo_root/source/eom_ev_mod.f90"
            ;;
    esac

    echo "Case $case_number passed ($data_rows numerical rows)"
    case_number=$((case_number + 1))
done

if [ "$mode" != mpi ]; then
    echo "Running Kelvin-Voigt evaporation regression"
    python3 "$repo_root/tests/evaporation/check_kv_evaporation.py" \
        --eom-source "$repo_root/source/eom_ev_mod.f90" \
        --integrator-source "$repo_root/source/integrator_kv_ev_mod.f90"

    kv_integrator=1
    while [ "$kv_integrator" -le 3 ]; do
        kv_dir="$work_dir/kv-evap-$kv_integrator"
        mkdir -p "$kv_dir"
        cp "$work_dir/execute/main.x" "$kv_dir/main.x"
        cp "$repo_root/examples/input-8/input.dat" "$kv_dir/input.dat"
        sed -i 's/\r$//' "$kv_dir/input.dat"
        sed -i \
            -e "s/^[[:space:]]*integrator[[:space:]].*/ integrator $kv_integrator/" \
            -e 's/^[[:space:]]*timestep[[:space:]].*/ timestep 1.d-8/' \
            -e 's/^[[:space:]]*final time[[:space:]].*/ final time 1.d-6/' \
            -e 's/^[[:space:]]*print time[[:space:]].*/ print time 2.d-7/' \
            "$kv_dir/input.dat"
        sed -i '/^[[:space:]]*Finish/i\ kvfluid yes' "$kv_dir/input.dat"

        echo "Running Kelvin-Voigt evaporation with integrator $kv_integrator"
        (
            cd "$kv_dir"
            timeout "${JETSPIN_SMOKE_TIMEOUT:-30}" ./main.x > run.log 2>&1
        )
        grep -q 'Program closed correctly' "$kv_dir/run.log"
        grep -q 'Kelvin-Voigt evaporation integrator active' "$kv_dir/run.log"
        if grep -Eiq '(^|[^[:alpha:]])(error|nan|[-+]?inf(inity)?)([^[:alpha:]]|$)' \
            "$kv_dir/run.log" "$kv_dir/statout.dat"; then
            echo "Kelvin-Voigt evaporation integrator $kv_integrator failed" >&2
            tail -60 "$kv_dir/run.log" >&2
            exit 1
        fi
        kv_integrator=$((kv_integrator + 1))
    done

    # Exercise the three-way coupling explicitly.  Test Case 5 is the
    # historical dynamic-refinement example; a lower allowed refinement
    # threshold and a short accelerated timestep make at least one Akima
    # remeshing event occur during the smoke run.  Evaporation and the
    # Kelvin-Voigt extension are then enabled on the same trajectory.
    refine_dir="$work_dir/kv-evap-refine"
    mkdir -p "$refine_dir"
    cp "$work_dir/execute/main.x" "$refine_dir/main.x"
    cp "$repo_root/examples/input-5/input.dat" "$refine_dir/input.dat"
    sed -i 's/\r$//' "$refine_dir/input.dat"
    sed -i \
        -e 's/^[[:space:]]*integrator[[:space:]].*/ integrator 3/' \
        -e 's/^[[:space:]]*timestep[[:space:]].*/ timestep 1.d-6/' \
        -e 's/^[[:space:]]*final time[[:space:]].*/ final time 6.d-3/' \
        -e 's/^[[:space:]]*print time[[:space:]].*/ print time 5.d-4/' \
        -e 's/^[[:space:]]*dynamic refinement every[[:space:]].*/ dynamic refinement every 5.d-4/' \
        -e 's/^[[:space:]]*dynamic refinement threshold[[:space:]].*/ dynamic refinement threshold 0.1d0/' \
        -e 's/^[[:space:]]*print list[[:space:]].*/ print list t n nref/' \
        "$refine_dir/input.dat"
    sed -i '/^[[:space:]]*Finish/i\ kvfluid yes\
 evaporation yes\
 evaporation polymer frac 0.06d0\
 evaporation airviscosity 0.15d0\
 evaporation temperature 293.15d0\
 evaporation umidity 0.165d0\
 evaporation bconstant 7.d0\
 evaporation mconstant 0.1d0\
 evaporation tconstant 1.d0\
 evaporation diffusivity 0.242d0' "$refine_dir/input.dat"

    echo "Running Kelvin-Voigt evaporation with dynamic refinement"
    (
        cd "$refine_dir"
        timeout "${JETSPIN_SMOKE_TIMEOUT:-30}" ./main.x > run.log 2>&1
    )
    grep -q 'Program closed correctly' "$refine_dir/run.log"
    grep -q 'dynamic refinement yes' "$refine_dir/run.log"
    grep -q 'evaporation yes' "$refine_dir/run.log"
    grep -q 'Kelvin-Voigt evaporation integrator active' "$refine_dir/run.log"
    if ! awk '
        $1 ~ /^[0-9]+$/ && NF == 4 {
            if (($4 + 0) > 0) refined = 1
            if (($3 + 0) > maxbeads) maxbeads = $3 + 0
        }
        END { exit (refined && maxbeads > 100) ? 0 : 1 }
    ' "$refine_dir/run.log"; then
        echo "Dynamic-refinement evaporation run did not force refinement/capacity growth" >&2
        tail -80 "$refine_dir/run.log" >&2
        exit 1
    fi
    if grep -Eiq '(^|[^[:alpha:]])(error|nan|[-+]?inf(inity)?)([^[:alpha:]]|$)' \
        "$refine_dir/run.log" "$refine_dir/statout.dat"; then
        echo "Kelvin-Voigt evaporation dynamic-refinement run failed" >&2
        tail -80 "$refine_dir/run.log" >&2
        exit 1
    fi
fi

echo "Smoke mode $mode passed ($last_case case(s))"
