#!/usr/bin/env python3
"""Compare a restarted run with an uninterrupted one after the restart step.

Usage: check_restart.py RESTARTED_DIR CONTINUOUS_DIR RESTART_STEP

Both directories hold statout.dat, traj.xyz and run.log.  Every statout.dat
row, terminal row and traj.xyz frame written after RESTART_STEP by the
restarted run must equal the uninterrupted run's, character for character;
the exit status is 1 otherwise, or if fewer than two rows were compared.
"""
import os
import sys


def statout_rows(path):
    rows = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split()
            rows[int(fields[0])] = fields[1:]
    return rows


def terminal_rows(path):
    rows = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3 or not fields[0].isdigit():
                continue
            try:
                [float(value) for value in fields[1:]]
            except ValueError:
                continue
            rows[int(fields[0])] = fields[1:]
    return rows


def xyz_frames(path):
    frames = {}
    with open(path) as handle:
        lines = handle.read().split('\n')
    index = 0
    while index < len(lines):
        if not lines[index].strip():
            index += 1
            continue
        count = int(lines[index])
        header = lines[index + 1]
        frames[int(header.split()[-1])] = lines[index + 2:index + 2 + count]
        index += 2 + count
    return frames


def compare(label, restarted, continuous, first_step):
    steps = [step for step in sorted(restarted)
             if step > first_step and step in continuous]
    differing = [step for step in steps if restarted[step] != continuous[step]]
    print('%s: %d after step %d, %d differ%s' % (
        label, len(steps), first_step, len(differing),
        (' (first at step %d)' % differing[0]) if differing else ''))
    return len(steps), len(differing)


def main():
    restarted, continuous, first_step = sys.argv[1], sys.argv[2], int(sys.argv[3])
    failed = False
    for label, reader, name in (('statout.dat rows', statout_rows, 'statout.dat'),
                                ('terminal rows', terminal_rows, 'run.log'),
                                ('traj.xyz frames', xyz_frames, 'traj.xyz')):
        paths = [os.path.join(restarted, name), os.path.join(continuous, name)]
        missing = [path for path in paths if not os.path.exists(path)]
        if missing:
            print('%s: missing %s' % (label, ', '.join(missing)))
            failed = True
            continue
        count, differing = compare(label, reader(paths[0]), reader(paths[1]),
                                   first_step)
        if count < 2 or differing:
            failed = True
    if failed:
        print('Restart check FAILED', file=sys.stderr)
        return 1
    print('Restart check passed: the restarted run continues the uninterrupted one exactly')
    return 0


if __name__ == '__main__':
    sys.exit(main())
