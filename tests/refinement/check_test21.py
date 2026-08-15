#!/usr/bin/env python3
"""Validate the event-level invariants of JETSPIN Test Case 21."""

import argparse
import math
import re
from pathlib import Path
from typing import List


INITIAL_RE = re.compile(
    r"Initial dynamic-refinement anchors: active=(\d+) anchors=(\d+) "
    r"spacing_cm=\s*([+\-0-9.Ee]+)"
)
EVENT_RE = re.compile(
    r"Dynamic refinement event: step=(\d+) active_before=(\d+) "
    r"active_after=(\d+) anchors_before=(\d+) anchors_after=(\d+)"
)
INVARIANT_RE = re.compile(
    r"Dynamic refinement invariants: "
    r"anchor_position_max_displacement_cm=\s*([+\-0-9.Ee]+) "
    r"reference_volume_relative_difference=\s*([+\-0-9.Ee]+) "
    r"evaporation_volume_relative_difference=\s*([+\-0-9.Ee]+) "
    r"ordered_path=([TF])"
)
ANCHOR_FIELD_RE = re.compile(
    r"Dynamic refinement anchor fields: "
    r"velocity_max_difference_internal=\s*([+\-0-9.Ee]+) "
    r"stress_max_difference_internal=\s*([+\-0-9.Ee]+) "
    r"radius_max_difference_cm=\s*([+\-0-9.Ee]+) "
    r"evaporation_radius_max_difference_cm=\s*([+\-0-9.Ee]+)"
)
AMOUNT_RE = re.compile(
    r"Dynamic refinement conserved amounts: "
    r"mass_relative_difference=\s*([+\-0-9.Ee]+) "
    r"charge_relative_difference=\s*([+\-0-9.Ee]+)"
)
CAPACITY_RE = re.compile(
    r"Gaussian history capacity: old=(\d+) new=(\d+) "
    r"retained_steps=(\d+) values=(\d+)"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"Test 21 validation failed: {message}")


def numerical_rows(path: Path) -> List[List[float]]:
    rows = []  # type: List[List[float]]
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        try:
            rows.append([float(value) for value in stripped.split()])
        except ValueError as exc:
            raise SystemExit(f"Invalid numerical row in {path}: {line}") from exc
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-log", type=Path, required=True)
    parser.add_argument("--statout", type=Path, required=True)
    parser.add_argument("--require-capacity-growth", action="store_true")
    args = parser.parse_args()

    log = args.run_log.read_text(encoding="utf-8")
    require("Program closed correctly" in log, "run did not close correctly")
    require("ERROR - numerical instability" not in log, "numerical instability")

    initial = INITIAL_RE.search(log)
    require(initial is not None, "missing initial-anchor diagnostic")
    initial_active, initial_anchors = map(int, initial.group(1, 2))
    spacing = float(initial.group(3))
    require(initial_active == 400, f"expected 400 initial elements, got {initial_active}")
    require(initial_anchors == 79, f"expected 79 interior anchors, got {initial_anchors}")
    require(math.isclose(spacing, 0.10, rel_tol=0.0, abs_tol=1.0e-12),
            f"expected 0.10 cm anchor spacing, got {spacing}")

    events = EVENT_RE.findall(log)
    require(len(events) == 1,
            f"expected one accepted refinement event, got {len(events)}")
    step, active_before, active_after, anchors_before, anchors_after = map(
        int, events[0]
    )
    require(step > 0, "invalid refinement step")
    require(active_after > active_before, "refinement did not add resolution")
    require(anchors_after >= anchors_before, "an existing anchor was lost")

    invariants = INVARIANT_RE.search(log)
    require(invariants is not None, "missing refinement-invariant diagnostic")
    anchor_error, volume_error, evaporation_error = map(
        float, invariants.group(1, 2, 3)
    )
    require(invariants.group(4) == "T", "refined path is not strictly ordered")
    require(anchor_error <= 1.0e-12,
            f"anchor displacement {anchor_error:.3e} cm exceeds tolerance")
    require(volume_error <= 1.0e-12,
            f"reference-volume error {volume_error:.3e} exceeds tolerance")
    require(evaporation_error <= 1.0e-12,
            f"evaporation-volume error {evaporation_error:.3e} exceeds tolerance")

    anchor_fields = ANCHOR_FIELD_RE.search(log)
    require(anchor_fields is not None, "missing anchor-field diagnostic")
    velocity_difference, stress_difference, radius_difference, evradius_difference = map(
        float, anchor_fields.groups()
    )
    require(velocity_difference <= 1.0e-12,
            f"anchor velocity difference {velocity_difference:.3e} exceeds tolerance")
    require(stress_difference <= 1.0e-12,
            f"anchor stress difference {stress_difference:.3e} exceeds tolerance")
    require(radius_difference <= 1.0e-12,
            f"anchor radius difference {radius_difference:.3e} cm exceeds tolerance")
    require(evradius_difference <= 1.0e-12,
            f"anchor evaporation-radius difference {evradius_difference:.3e} cm exceeds tolerance")

    amounts = AMOUNT_RE.search(log)
    require(amounts is not None, "missing conserved-amount diagnostic")
    mass_difference, charge_difference = map(float, amounts.groups())
    require(mass_difference <= 1.0e-12,
            f"mass difference {mass_difference:.3e} exceeds tolerance")
    require(charge_difference <= 1.0e-12,
            f"charge difference {charge_difference:.3e} exceeds tolerance")

    rows = numerical_rows(args.statout)
    require(rows, "statout.dat contains no numerical rows")
    require(all(math.isfinite(value) for row in rows for value in row),
            "statout.dat contains a non-finite value")
    require(rows[-1][-1] >= 1.0, "final nref is zero")

    capacity = CAPACITY_RE.search(log)
    if args.require_capacity_growth:
        require(capacity is not None, "capacity-growth diagnostic is missing")
        old_capacity, new_capacity, retained_steps, history_values = map(
            int, capacity.groups()
        )
        require(new_capacity > old_capacity,
                "capacity-growth run did not increase capacity")
        require(new_capacity >= active_after,
                "new capacity is smaller than the refined topology")
        require(retained_steps > 0 and history_values > 0,
                "resized Gaussian history is empty")

    print(
        "Test 21 passed: "
        f"step={step}, active={active_before}->{active_after}, "
        f"anchors={anchors_before}->{anchors_after}, "
        f"anchor_displacement={anchor_error:.3e} cm, "
        f"volume_differences=({volume_error:.3e}, {evaporation_error:.3e})"
    )


if __name__ == "__main__":
    main()
