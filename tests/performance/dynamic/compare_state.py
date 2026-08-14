#!/usr/bin/env python3
"""Compare full bead states written by JETSPIN_TOPOLOGY_SNAPSHOT=1."""

import argparse
from pathlib import Path


COLUMNS = ("x", "y", "z", "stress", "vx", "vy", "vz", "mass", "charge", "volume", "ev_volume", "ev_radius")


def read(path):
    events = []
    current = None
    for lineno, line in enumerate(Path(path).read_text().splitlines(), 1):
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "event":
            current = {"meta": tuple(fields[1:]), "beads": []}
            events.append(current)
        elif current is None:
            raise ValueError(f"{path}:{lineno}: bead before event")
        else:
            current["beads"].append(
                (int(fields[0]), fields[1] == "T", *map(float, fields[2:]))
            )
    return events


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("actual")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument("--allow-event-divergence", action="store_true",
                        help="compare the common event prefix and report topology drift")
    args = parser.parse_args()
    ref, actual = read(args.reference), read(args.actual)
    if len(ref) != len(actual) and not args.allow_event_divergence:
        raise AssertionError(f"event count differs: {len(ref)} != {len(actual)}")
    event_limit = min(len(ref), len(actual))
    if len(ref) != len(actual):
        print(f"WARNING: event count differs: CPU={len(ref)} GPU={len(actual)}; "
              f"comparing first {event_limit} events")
    worst = None
    first = None
    mismatches = 0
    for event_index, (revent, aevent) in enumerate(zip(ref[:event_limit], actual[:event_limit]), 1):
        if revent["meta"] != aevent["meta"]:
            if not args.allow_event_divergence:
                raise AssertionError(f"event {event_index} metadata differs")
            print(f"WARNING: event {event_index} metadata differs: "
                  f"CPU={revent['meta']} GPU={aevent['meta']}")
            continue
        if len(revent["beads"]) != len(aevent["beads"]):
            raise AssertionError(f"event {event_index} bead count differs")
        step = int(revent["meta"][0])
        for rbead, abead in zip(revent["beads"], aevent["beads"]):
            if rbead[:2] != abead[:2]:
                raise AssertionError(f"step {step}: bead identity/flag differs")
            bead = rbead[0]
            for column, reference, value in zip(COLUMNS, rbead[2:], abead[2:]):
                difference = abs(value - reference)
                allowed = args.atol + args.rtol * abs(reference)
                normalized = difference / allowed if allowed else float("inf")
                record = (normalized, step, bead, column, reference, value, difference, allowed)
                if worst is None or record[0] > worst[0]:
                    worst = record
                if difference > allowed:
                    mismatches += 1
                    if first is None:
                        first = record
    if first:
        _, step, bead, column, reference, value, difference, allowed = first
        raise AssertionError(
            f"{mismatches} mismatches; first at step={step}, bead={bead}, "
            f"column={column}: cpu={reference:.17g}, gpu={value:.17g}, "
            f"abs_diff={difference:.3g}, allowed={allowed:.3g}"
        )
    if worst is None:
        print("PASS: no common events to compare")
        return
    normalized, step, bead, column, reference, value, difference, allowed = worst
    print(
        f"PASS: {event_limit} common events; worst={normalized:.3g} at step={step}, "
        f"bead={bead}, column={column}, abs_diff={difference:.3g}, allowed={allowed:.3g}"
    )


if __name__ == "__main__":
    main()
