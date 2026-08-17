#!/usr/bin/env python3
"""Validate repeated refinement interleaved with collector removal."""

import argparse
from collections import Counter
import math
import re
from pathlib import Path

from check_test21 import (
    AMOUNT_RE,
    ANCHOR_FIELD_RE,
    CAPACITY_RE,
    EVENT_RE,
    INITIAL_RE,
    INVARIANT_RE,
    numerical_rows,
)


EXPECTED_EVENTS = 3
TOLERANCE = 1.0e-12
TOPOLOGY_RE = re.compile(
    r"Topology event: step=(\d+) add=([TF]) remove=(\d+) active=(\d+)"
)
AKIMA_RE = re.compile(
    r"Akima device comparison: field=(\S+) "
    r"coefficient_max_abs=\s*([+\-0-9.Ee]+) "
    r"coefficient_max_rel=\s*([+\-0-9.Ee]+) "
    r"value_max_abs=\s*([+\-0-9.Ee]+) "
    r"value_max_rel=\s*([+\-0-9.Ee]+)"
)
AKIMA_FIELDS = {
    "x", "y", "z", "vx", "vy", "vz", "stress", "mass_density",
    "charge_density", "radius_area", "evap_radius_area",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"Test 23 validation failed: {message}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-log", type=Path, required=True)
    parser.add_argument("--statout", type=Path, required=True)
    parser.add_argument("--require-akima-comparison", action="store_true")
    args = parser.parse_args()

    log = args.run_log.read_text(encoding="utf-8")
    require("Program closed correctly" in log, "run did not close correctly")
    require("ERROR - numerical instability" not in log, "numerical instability")

    initial = INITIAL_RE.search(log)
    require(initial is not None, "missing initial-anchor diagnostic")
    initial_active, initial_anchors = map(int, initial.group(1, 2))
    spacing = float(initial.group(3))
    require(initial_active == 400,
            f"expected 400 initial elements, got {initial_active}")
    require(initial_anchors == 79,
            f"expected 79 interior anchors, got {initial_anchors}")
    require(math.isclose(spacing, 0.10, rel_tol=0.0, abs_tol=TOLERANCE),
            f"expected 0.10 cm anchor spacing, got {spacing}")

    events = [tuple(map(int, match)) for match in EVENT_RE.findall(log)]
    require(len(events) == EXPECTED_EVENTS,
            f"expected {EXPECTED_EVENTS} refinement events, got {len(events)}")
    previous_step = 0
    for index, event in enumerate(events, 1):
        step, active_before, active_after, anchors_before, anchors_after = event
        require(step > previous_step,
                f"event {index} step is not strictly increasing")
        require(active_after > active_before,
                f"event {index} did not add resolution")
        require(anchors_after == anchors_before,
                f"event {index} changed the active anchor count")
        require(anchors_before >= initial_anchors,
                f"event {index} lost an initial anchor prematurely")
        previous_step = step

    topology = [
        (int(step), added == "T", int(removed), int(active))
        for step, added, removed, active in TOPOLOGY_RE.findall(log)
    ]
    removals = [event for event in topology if event[2] > 0]
    require(removals, "collector removal was never exercised")
    require(all(not added and removed == 1 for _, added, removed, _ in removals),
            "unexpected combined or multi-bead collector event")
    require(any(step < events[-1][0] for step, _, _, _ in removals),
            "no collector removal occurred before the final remesh")
    require(any(step > events[-1][0] for step, _, _, _ in removals),
            "no collector removal occurred after the final remesh")

    invariants = INVARIANT_RE.findall(log)
    require(len(invariants) == EXPECTED_EVENTS,
            "missing per-event refinement invariants")
    for index, values in enumerate(invariants, 1):
        anchor_error, volume_error, evaporation_error = map(float, values[:3])
        require(values[3] == "T", f"event {index} path is not ordered")
        require(anchor_error <= TOLERANCE,
                f"event {index} anchor displacement exceeds tolerance")
        require(volume_error <= TOLERANCE,
                f"event {index} reference-volume error exceeds tolerance")
        require(evaporation_error <= TOLERANCE,
                f"event {index} evaporation-volume error exceeds tolerance")

    anchor_fields = ANCHOR_FIELD_RE.findall(log)
    require(len(anchor_fields) == EXPECTED_EVENTS,
            "missing per-event anchor-field invariants")
    for index, values in enumerate(anchor_fields, 1):
        require(all(float(value) <= TOLERANCE for value in values),
                f"event {index} changed an anchor field")

    amounts = AMOUNT_RE.findall(log)
    require(len(amounts) == EXPECTED_EVENTS,
            "missing per-event mass/charge invariants")
    for index, values in enumerate(amounts, 1):
        require(all(float(value) <= TOLERANCE for value in values),
                f"event {index} failed mass/charge conservation")

    capacities = [tuple(map(int, match)) for match in CAPACITY_RE.findall(log)]
    require(len(capacities) == EXPECTED_EVENTS,
            f"expected {EXPECTED_EVENTS} capacity increases, got {len(capacities)}")
    for index, (old, new, retained_steps, values) in enumerate(capacities, 1):
        require(new > old, f"growth {index} did not increase capacity")
        if index > 1:
            require(old == capacities[index - 2][1],
                    f"growth {index} does not continue the previous capacity")
        require(retained_steps > 0 and values == (new + 1) * 6 * retained_steps,
                f"growth {index} has an inconsistent Gaussian history")

    rows = numerical_rows(args.statout)
    require(rows, "statout.dat contains no numerical rows")
    require(all(math.isfinite(value) for row in rows for value in row),
            "statout.dat contains a non-finite value")
    require(int(round(rows[-1][-1])) == EXPECTED_EVENTS,
            f"final nref is not {EXPECTED_EVENTS}")

    akima_summary = ""
    if args.require_akima_comparison:
        comparisons = AKIMA_RE.findall(log)
        require(len(comparisons) == EXPECTED_EVENTS * len(AKIMA_FIELDS),
                f"expected 33 Akima comparisons, got {len(comparisons)}")
        field_counts = Counter(values[0] for values in comparisons)
        require(set(field_counts) == AKIMA_FIELDS,
                "Akima comparison field set is incomplete")
        require(all(count == EXPECTED_EVENTS for count in field_counts.values()),
                "each Akima field was not compared at every event")
        parsed_metrics = []
        for field, coefficient_abs, coefficient_rel, value_abs, value_rel in comparisons:
            metrics = tuple(map(float, (
                coefficient_abs, coefficient_rel, value_abs, value_rel,
            )))
            parsed_metrics.append(metrics)
            require(all(math.isfinite(value) for value in metrics),
                    f"non-finite Akima comparison for {field}")
            require(metrics[1] <= 1.0e-12,
                    f"Akima coefficient relative error for {field} exceeds 1e-12")
            require(metrics[2] <= 1.0e-9,
                    f"Akima interpolated absolute error for {field} exceeds 1e-9")
            require(metrics[3] <= 1.0e-12,
                    f"Akima interpolated relative error for {field} exceeds 1e-12")
        akima_summary = (
            ", Akima max coefficient rel="
            f"{max(values[1] for values in parsed_metrics):.4e}, "
            "max value abs="
            f"{max(values[2] for values in parsed_metrics):.4e}, "
            "max value rel="
            f"{max(values[3] for values in parsed_metrics):.4e}"
        )

    event_summary = ", ".join(
        f"{step}:{before}->{after}"
        for step, before, after, _, _ in events
    )
    removal_summary = ", ".join(str(step) for step, _, _, _ in removals)
    print(
        "Test 23 passed: "
        f"events=[{event_summary}], removals=[{removal_summary}], "
        f"final_active={int(round(rows[-1][9]))}{akima_summary}"
    )


if __name__ == "__main__":
    main()
