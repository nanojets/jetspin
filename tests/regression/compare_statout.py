#!/usr/bin/env python3
"""Compare two JETSPIN statout.dat files with mixed exact/tolerant fields."""

import argparse
import math
from pathlib import Path
from typing import List, Optional, Tuple

EXACT_COLUMNS = {"nstep", "n", "nref"}


def read_table(path: Path) -> Tuple[List[str], List[List[float]]]:
    columns = None  # type: Optional[List[str]]
    rows = []  # type: List[List[float]]
    for line_number, raw in enumerate(
        path.read_text(encoding="utf-8", errors="strict").splitlines(), 1
    ):
        fields = raw.split()
        if not fields:
            continue
        if raw.startswith("#    nstep"):
            # io_mod writes nstep in a 10-character field followed by one
            # 20-character field per observable. Units inside those fields
            # contain spaces, so ordinary whitespace splitting is ambiguous.
            columns = ["nstep"]
            offset = 10
            while offset < len(raw):
                label_field = raw[offset:offset + 20].strip()
                if label_field:
                    columns.append(label_field.split()[0])
                offset += 20
            continue
        if fields[0].startswith("#"):
            continue
        try:
            rows.append([float(value.replace("D", "E").replace("d", "e"))
                         for value in fields])
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: invalid numerical row") from error

    if columns is None:
        raise ValueError(f"{path}: statout column header not found")
    if not rows:
        raise ValueError(f"{path}: no numerical rows found")
    for index, row in enumerate(rows, 1):
        if len(row) != len(columns):
            raise ValueError(
                f"{path}: data row {index} has {len(row)} fields; "
                f"expected {len(columns)}"
            )
        if not all(math.isfinite(value) for value in row):
            raise ValueError(f"{path}: data row {index} contains NaN or infinity")
    return columns, rows


def compare(reference: Path, actual: Path, rtol: float, atol: float) -> None:
    ref_columns, ref_rows = read_table(reference)
    got_columns, got_rows = read_table(actual)
    if got_columns != ref_columns:
        raise AssertionError(
            f"column mismatch:\nreference {ref_columns}\nactual    {got_columns}"
        )
    if len(got_rows) != len(ref_rows):
        raise AssertionError(
            f"row-count mismatch: reference {len(ref_rows)}, actual {len(got_rows)}"
        )

    worst_ratio = 0.0
    worst_detail = ""
    failures = []  # type: List[str]
    for row_index, (ref_row, got_row) in enumerate(zip(ref_rows, got_rows), 1):
        for column, expected, actual_value in zip(ref_columns, ref_row, got_row):
            if column in EXACT_COLUMNS:
                if actual_value != expected:
                    failures.append(
                        f"row {row_index}, {column}: {actual_value:g} != {expected:g}"
                    )
                continue
            allowed = atol + rtol * abs(expected)
            difference = abs(actual_value - expected)
            ratio = difference / allowed if allowed else float("inf")
            detail = (
                f"row {row_index}, {column}: actual={actual_value:.12g}, "
                f"reference={expected:.12g}, abs_diff={difference:.3g}, "
                f"allowed={allowed:.3g}"
            )
            if ratio > worst_ratio:
                worst_ratio = ratio
                worst_detail = detail
            if difference > allowed:
                failures.append(detail)

    if failures:
        preview = "\n".join(f"  {item}" for item in failures[:12])
        suffix = "" if len(failures) <= 12 else f"\n  ... {len(failures) - 12} more"
        raise AssertionError(
            f"{len(failures)} numerical mismatch(es) with "
            f"rtol={rtol:g}, atol={atol:g}:\n{preview}{suffix}"
        )

    print(
        f"PASS {actual}: {len(got_rows)} rows, {len(got_columns)} columns; "
        f"worst normalized difference={worst_ratio:.3g}"
    )
    if worst_detail:
        print(f"     {worst_detail}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("actual", type=Path)
    parser.add_argument("--rtol", type=float, default=1.0e-7)
    parser.add_argument("--atol", type=float, default=1.0e-10)
    args = parser.parse_args()
    if args.rtol < 0.0 or args.atol < 0.0:
        parser.error("tolerances must be non-negative")
    compare(args.reference, args.actual, args.rtol, args.atol)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
