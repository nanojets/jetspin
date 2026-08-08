#!/usr/bin/env python3
"""Regression checks for the Yarin et al. (2001) evaporation model.

The checker reads the dimensional parameters of Test Case 8 and recomputes
fixed reference quantities.  When source files are supplied it also guards
the production expressions used for the cutoff and rheological update.  A
JETSPIN run log can be supplied to verify the actual run-time fallback value
of the solvent diffusivity.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


CS_CUTOFF = 0.1
DA_PREFAC = 0.211
DA_TREF = 273.15
DA_EXP = 1.94

# Fixed Test Case 8 regression targets, independently recorded here.
REF_THETA0 = 1.0e-2
REF_DA = 0.2420017580740141
REF_EV_LIM = 0.06666666666666667
REF_MU_RATIO = 43.978335002057946
REF_THETA_RATIO = 15.0
REF_G_RATIO = 2.931889000137196


def fortran_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def read_directive(path: Path, directive: str) -> float:
    wanted = directive.lower().split()
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith(("!", "#")):
            continue
        line = line.split("!", 1)[0].split("#", 1)[0].strip()
        fields = line.split()
        if [x.lower() for x in fields[: len(wanted)]] == wanted:
            if len(fields) <= len(wanted):
                raise ValueError(f"No value after directive {directive!r}")
            return fortran_float(fields[len(wanted)])
    raise KeyError(f"Directive {directive!r} not found in {path}")


def assert_close(name: str, actual: float, expected: float,
                 rel: float = 5.0e-9, abs_: float = 5.0e-12) -> None:
    if not math.isclose(actual, expected, rel_tol=rel, abs_tol=abs_):
        raise AssertionError(
            f"{name}: got {actual:.15g}, expected {expected:.15g}"
        )
    print(f"PASS {name:24s} {actual:.12g}")


def parse_runtime_value(log: str, label: str) -> float:
    pattern = re.compile(re.escape(label) + r"\s*([-+0-9.eEdD]+)", re.I)
    match = pattern.search(log)
    if not match:
        raise AssertionError(f"Run-time diagnostic not found: {label!r}")
    return fortran_float(match.group(1))


def compact_fortran(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace").lower()
    # Strip comments and whitespace so line continuations/indentation do not
    # make the regression check sensitive to formatting.
    text = "\n".join(line.split("!", 1)[0] for line in text.splitlines())
    return re.sub(r"\s+", "", text)


def assert_source_expression(path: Path, expression: str, name: str) -> None:
    compact = compact_fortran(path)
    expected = re.sub(r"\s+", "", expression.lower())
    if expected not in compact:
        raise AssertionError(f"{name}: production expression not found in {path}")
    print(f"PASS {name}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True,
                        help="JETSPIN Test Case 8 input.dat")
    parser.add_argument("--run-log", type=Path,
                        help="optional JETSPIN run.log for run-time checks")
    parser.add_argument("--nanojet-source", type=Path,
                        help="optional source/nanojet_mod.f90")
    parser.add_argument("--eom-source", type=Path,
                        help="optional source/eom_ev_mod.f90")
    args = parser.parse_args()

    cp0 = read_directive(args.input, "evaporation polymer frac")
    temperature = read_directive(args.input, "evaporation temperature")
    bconst = read_directive(args.input, "evaporation bconstant")
    mconst = read_directive(args.input, "evaporation mconstant")
    tconst = read_directive(args.input, "evaporation tconstant")
    mu0 = read_directive(args.input, "viscosity")
    g0 = read_directive(args.input, "elastic modulus")

    theta0 = mu0 / g0
    da = DA_PREFAC * (temperature / DA_TREF) ** DA_EXP
    evlim = cp0 / (1.0 - CS_CUTOFF)
    cp_cut = cp0 / evlim
    mu_ratio = 10.0 ** (bconst * (cp_cut ** mconst - cp0 ** mconst))
    theta_ratio = (cp_cut / cp0) ** tconst
    g_ratio = mu_ratio / theta_ratio

    assert_close("theta0 [s]", theta0, REF_THETA0)
    assert_close("fallback D_a [cm2/s]", da, REF_DA)
    assert_close("cutoff V/V0", evlim, REF_EV_LIM)
    assert_close("cutoff cp", cp_cut, 0.9)
    assert_close("mu/mu0 at cutoff", mu_ratio, REF_MU_RATIO)
    assert_close("theta/theta0 cutoff", theta_ratio, REF_THETA_RATIO)
    assert_close("G/G0 at cutoff", g_ratio, REF_G_RATIO)
    assert_close("Yarin tconstant", tconst, 1.0)

    if args.nanojet_source is not None:
        assert_source_expression(
            args.nanojet_source,
            "evlim=min(1.d0,cp0/(1.d0-evsolvlim))",
            "production Yarin cutoff",
        )
        assert_source_expression(
            args.nanojet_source,
            "evmasscoeff=0.211d0*((evtemp/273.15d0)**1.94d0)",
            "production diffusivity fallback",
        )

    if args.eom_source is not None:
        assert_source_expression(
            args.eom_source,
            "cp=cp0*yvl(ipoint)/yve(ipoint)",
            "production polymer concentration",
        )
        assert_source_expression(
            args.eom_source,
            "rattao=(cp/cp0)**tev",
            "production relaxation-time ratio",
        )
        assert_source_expression(
            args.eom_source,
            "ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))",
            "production viscosity ratio",
        )

    if args.run_log is not None:
        log = args.run_log.read_text(encoding="utf-8", errors="replace")
        runtime_da = parse_runtime_value(
            log, "mass diffusivity of solvent automatically set equal to"
        )
        # error_mod prints this value with g20.10, so allow the precision of
        # the legacy diagnostic format rather than full Python precision.
        assert_close("runtime fallback D_a", runtime_da, da,
                     rel=5.0e-9, abs_=5.0e-10)

    print("Yarin-2001 evaporation regression checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
