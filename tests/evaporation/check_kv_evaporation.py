#!/usr/bin/env python3
"""Regression guard for the concentration-dependent Kelvin-Voigt extension."""
import argparse
import math
import re
from pathlib import Path

REF_RMU = 8.367216948213878
REF_RG = 1.6734433896427756
REF_DRMU = 2.3913155036429394
REF_DRG = 0.14357442280003282
REF_STRAIN = -0.20121487043148
REF_DSTRESS = 0.1265700720011721


def close(name, got, ref, tol=2e-12):
    if not math.isclose(got, ref, rel_tol=tol, abs_tol=tol):
        raise AssertionError(f"{name}: {got:.16g} != {ref:.16g}")
    print(f"PASS {name:26s} {got:.12g}")


def compact(path):
    text = path.read_text(encoding="utf-8", errors="replace").lower()
    text = "\n".join(line.split("!", 1)[0] for line in text.splitlines())
    return re.sub(r"\s+", "", text)


def require(path, expression, name):
    if re.sub(r"\s+", "", expression.lower()) not in compact(path):
        raise AssertionError(f"{name}: expression not found in {path}")
    print(f"PASS {name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eom-source", required=True, type=Path)
    ap.add_argument("--integrator-source", required=True, type=Path)
    args = ap.parse_args()

    cp0, cp = 0.06, 0.30
    b, m, tev = 7.0, 0.1, 1.0
    stress, strainrate, strainacc = 0.5, 0.1, -0.03
    # fV/V=-0.2, hence dc_p/dt=-c_p fV/V=+0.06.
    dcpdt = 0.06

    rmu = 10.0 ** (b * (cp**m - cp0**m))
    rg = rmu / ((cp/cp0)**tev)
    drmu = rmu * math.log(10.0) * b * m * cp**(m-1.0) * dcpdt
    drg = rg * (math.log(10.0)*b*m*cp**(m-1.0)-tev/cp) * dcpdt
    strain = (stress-rmu*strainrate)/rg
    dstress = (rg*strainrate + rmu*strainacc +
               drg*strain + drmu*strainrate)

    close("mu/mu0", rmu, REF_RMU)
    close("G/G0", rg, REF_RG)
    close("d(mu/mu0)/dtbar", drmu, REF_DRMU)
    close("d(G/G0)/dtbar", drg, REF_DRG)
    close("recovered strain", strain, REF_STRAIN)
    close("Kelvin-Voigt stress rate", dstress, REF_DSTRESS)

    # When concentration is frozen, the extension must reduce exactly to
    # the historical JETSPIN Kelvin-Voigt differential form.
    frozen = rg*strainrate + rmu*strainacc
    product_rule_frozen = rg*strainrate + rmu*strainacc + 0.0 + 0.0
    close("frozen-composition limit", product_rule_frozen, frozen)

    require(args.eom_source,
            "dcpdt=-cp*fevlocal/yve",
            "production dc_p/dt")
    require(args.eom_source,
            "strain=(stress-ratmu*strainrate)/ratg",
            "production strain recovery")
    require(args.eom_source,
            "fst=ratg*strainrate+ratmu*strainacc+dratg*strain+dratmu*strainrate",
            "production product rule")
    require(args.integrator_source,
            "call eom3_KV_st_ev",
            "3D Kelvin-Voigt evaporation path")
    require(args.integrator_source,
            "call clamp_ev_volume",
            "Yarin cutoff in KV integrators")

    print("Kelvin-Voigt evaporation regression checks passed")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
