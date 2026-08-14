#!/bin/sh

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/jetspin-dynamic-evaporation.XXXXXX")
if [ "${JETSPIN_DYNAMIC_EVAP_KEEP:-0}" = 1 ]; then
    echo "Keeping validation files in $work_dir"
else
    trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
fi

gpu_cc=${GPUCC:-80}
cuda_version=${CUDA_VERSION:-12.3}
timeout_seconds=${JETSPIN_DYNAMIC_EVAP_TIMEOUT:-60}

if ! command -v nvfortran >/dev/null 2>&1; then
    echo "nvfortran is not available; load NVIDIA HPC SDK first" >&2
    exit 2
fi

mkdir -p "$work_dir/cpu-source" "$work_dir/gpu-source"
mkdir -p "$work_dir/cpu-bin" "$work_dir/gpu-bin"
cp "$repo_root"/source/*.f90 "$work_dir/cpu-source/"
cp "$repo_root"/source/*.f90 "$work_dir/gpu-source/"
cp "$repo_root/build/Makefile" "$work_dir/cpu-source/Makefile"
cp "$repo_root/build/Makefile" "$work_dir/gpu-source/Makefile"

nvfortran --version > "$work_dir/toolchain.txt"
printf '%s\n' \
    "GPUCC=$gpu_cc" \
    "CUDA_VERSION=$cuda_version" \
    "OpenACC flags=-O3 -acc=gpu -gpu=cc${gpu_cc},cuda${cuda_version},nofma" \
    >> "$work_dir/toolchain.txt"

echo "Building the NVFORTRAN CPU reference"
if ! make -C "$work_dir/cpu-source" nvfortran \
    BINROOT="$work_dir/cpu-bin" FPPFLAGS_EXTRA= \
    > "$work_dir/cpu-build.log" 2>&1; then
    tail -80 "$work_dir/cpu-build.log" >&2
    exit 1
fi

echo "Building the standard OpenACC executable (no diagnostic macros)"
if ! make -C "$work_dir/gpu-source" nvfortran-openacc \
    BINROOT="$work_dir/gpu-bin" GPUCC="$gpu_cc" \
    CUDA_VERSION="$cuda_version" FPPFLAGS_EXTRA= \
    > "$work_dir/gpu-build.log" 2>&1; then
    tail -80 "$work_dir/gpu-build.log" >&2
    exit 1
fi

prepare_input() {
    case_number=$1
    integrator_number=$2
    duration=$3
    destination=$4

    cp "$repo_root/examples/input-$case_number/input.dat" "$destination/input.dat"
    sed -i 's/\r$//' "$destination/input.dat"
    sed -i \
        -e "s/^[[:space:]]*integrator[[:space:]].*/ integrator $integrator_number/" \
        "$destination/input.dat"

    if [ "$duration" = short ]; then
        sed -i \
            -e 's/^[[:space:]]*final time[[:space:]].*/ final time 3.d-7/' \
            -e 's/^[[:space:]]*print time[[:space:]].*/ print time 1.d-7/' \
            "$destination/input.dat"
        if [ "$case_number" -eq 16 ]; then
            sed -i '/^[[:space:]]*Finish/i\ print xyz 1.d-7\
 print xyz maxnumber 1000' "$destination/input.dat"
        fi
    fi
}

run_case() {
    executable=$1
    destination=$2

    cp "$executable" "$destination/main.x"
    (
        cd "$destination"
        timeout "$timeout_seconds" ./main.x > run.log 2>&1
    )
    grep -q 'Program closed correctly' "$destination/run.log"
    test -s "$destination/statout.dat"
    if grep -Eiq '(^|[^[:alpha:]])(error|nan|[-+]?inf(inity)?)([^[:alpha:]]|$)' \
        "$destination/run.log" "$destination/statout.dat"; then
        echo "Non-finite value or runtime error in $destination" >&2
        tail -80 "$destination/run.log" >&2
        exit 1
    fi
}

check_topology() {
    run_log=$1
    grep -q '^Topology additions: 111$' "$run_log"
    grep -q '^Topology removals: 122$' "$run_log"
    grep -q '^Array reallocations: 2$' "$run_log"
    grep -q '^Topology active beads: 89$' "$run_log"
}

case_number=16
while [ "$case_number" -le 17 ]; do
    integrator_number=1
    while [ "$integrator_number" -le 3 ]; do
        label="test-${case_number}-integrator-${integrator_number}"
        cpu_full="$work_dir/$label-cpu-full"
        gpu_full="$work_dir/$label-gpu-full"
        cpu_short="$work_dir/$label-cpu-short"
        gpu_short="$work_dir/$label-gpu-short"
        mkdir -p "$cpu_full" "$gpu_full" "$cpu_short" "$gpu_short"

        prepare_input "$case_number" "$integrator_number" full "$cpu_full"
        cp "$cpu_full/input.dat" "$gpu_full/input.dat"
        prepare_input "$case_number" "$integrator_number" short "$cpu_short"
        cp "$cpu_short/input.dat" "$gpu_short/input.dat"

        echo "Running Test $case_number integrator $integrator_number (CPU/GPU, 1,000 steps)"
        run_case "$work_dir/cpu-bin/main.x" "$cpu_full"
        run_case "$work_dir/gpu-bin/main.x" "$gpu_full"
        check_topology "$cpu_full/run.log"
        check_topology "$gpu_full/run.log"

        echo "Comparing Test $case_number integrator $integrator_number before topology events"
        run_case "$work_dir/cpu-bin/main.x" "$cpu_short"
        run_case "$work_dir/gpu-bin/main.x" "$gpu_short"
        python3 "$repo_root/tests/regression/compare_statout.py" \
            --rtol 1e-12 --atol 1e-13 \
            "$cpu_short/statout.dat" "$gpu_short/statout.dat"

        if [ "$case_number" -eq 16 ]; then
            if ! cmp -s "$cpu_short/traj.xyz" "$gpu_short/traj.xyz"; then
                echo "Maxwell pre-event XYZ geometry differs for integrator $integrator_number" >&2
                exit 1
            fi
            echo "PASS Test 16 integrator $integrator_number XYZ geometry: byte-identical"
        fi

        integrator_number=$((integrator_number + 1))
    done
    case_number=$((case_number + 1))
done

echo "Dynamic evaporation validation passed for Tests 16 and 17, integrators 1--3"
