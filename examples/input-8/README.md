# Test Case 8 -- Yarin et al. (2001) evaporation reference case

This example is a JETSPIN reference/regression case for the solvent-evaporation and rheological-solidification model of A. L. Yarin, S. Koombhongse and D. H. Reneker, *Journal of Applied Physics* **89**, 3018--3026 (2001), DOI 10.1063/1.1333035.

The input follows the dimensional parameter set used in Sec. VI of that paper as closely as possible within the current JETSPIN discretization and input conventions:

| quantity | Test 8 value | source value |
| --- | ---: | ---: |
| initial polymer mass fraction | 0.06 | 6 wt% PEO |
| initial radius `a0` | 0.015 cm | 150 micrometers |
| density | 1 g/cm^3 | 1000 kg/m^3 |
| surface tension | 70 dyn/cm | 0.07 N/m |
| initial viscosity `mu0` | 1.0e4 P | 1000 kg/(m s) |
| initial relaxation time `theta0` | 0.01 s | 10 ms |
| elastic modulus `G0=mu0/theta0` | 1.0e6 Ba | derived |
| charge density | 2.99792458e6 statC/cm^3 | 1 C/L |
| collector distance | 20 cm | 20 cm |
| external field | 0.05 statV/cm | 1.5 kV/m |
| applied potential in JETSPIN | 1 statV | derived from `E*h` |
| relative humidity | 0.165 | 16.5% |
| temperature | 293.15 K | 20 C |
| `B` | 7 | 7 |
| `m` | 0.1 | 0.1 |
| relaxation exponent | 1 | `theta/theta0 = cp/cp0` |

The paper reports the dimensionless perturbation frequency `K_s=100`. With `theta0=mu0/G0=0.01 s`, Test 8 therefore uses `omega=1.0e4 s^-1`. The perturbation amplitude `1.0e-3 cm` is only a small JETSPIN seed for the bending instability; it is not one of the dimensional values reported in the parameter list of Yarin et al.

The original 2001 calculation and the current JETSPIN bead-insertion algorithm are not identical discretizations. Test 8 therefore uses `lstep=0.0325 cm`, the electrostatic cutoff length estimated in the same paper, which also satisfies the JETSPIN requirement that the discretization length exceed the initial jet radius. Consequently this case should be used to verify the evaporation/rheology implementation and parameter mapping, not as a point-by-point reproduction of the original Fig. 2 trajectory.

## Expected evaporation/rheology checks

At 293.15 K, with no explicit `evaporation diffusivity` directive, JETSPIN uses the atmospheric-pressure water-vapour correlation

`D_a = 0.211*(T/273.15)^1.94 = 0.242001758 cm^2/s`.

The Yarin solidification cutoff is `cs=0.1`. Since

`cp = cp0*V0/V`,

Test 8 must stop evaporation at

`V/V0 = cp0/(1-cs) = 0.06/0.9 = 0.0666666667`.

At that cutoff `cp=0.9`, so the Yarin concentration laws give approximately

- `mu/mu0 = 43.9783350`,
- `theta/theta0 = 15`,
- `G/G0 = 2.9318890`.

These values are useful regression targets for checking the `evrc`, `visc` and `gc` observables and for diagnosing future changes to the evaporation implementation.
