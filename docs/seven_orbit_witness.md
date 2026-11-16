# A common seven-orbit mixture

The inputs in `data/seven_orbit_witness.json` make the mixture in article Sec. V A and Supplemental Sec. S5 directly inspectable. This is a retrospective fitted feasibility witness. Its weights do not identify a microscopic preparation, and the response away from the specified observables remains undetermined.

## Reproduce the scalar check

From the collection root, run:

```sh
python3 scripts/check_seven_orbit_witness.py --output outputs/witness/seven_orbit_witness_check.json
```

Only the Python standard library is required. The script locates its default input relative to its own location, so the input JSON and script can also be copied alone into `data/` and `scripts/` in a separate directory. No original repository, optimizer, state vector, or network connection is needed. `--input` selects another input JSON; omitting `--output` prints the result. Invalid input returns a nonzero exit status.

The deterministic generated JSON matches `results/witness/seven_orbit_witness_check.json` byte for byte. It contains exact reconstructed endpoints, joint populations, moments, and containment slacks. The checker recomputes these values from the seven rows rather than accepting a saved success flag.

## Common weights and observables

The occupied orbit IDs are **92, 3614, 4909, 5596, 6711, 7952, and 8163**. Their positive rational weights sum to one exactly. All other IDs in the 8,192-orbit universe have zero weight. These fitted weights are independent of the fixed disorder-law probabilities supplied elsewhere in this collection; no catalog-probability factors enter this check.

The finite system has 24 sites and torus translation columns `(2,2)` and `(-2,4)`. The field is along unit `[111]`, with `b/J = 0.50`. The probe is the local site-0 z spin, `S = sigma/2`. In working energy units, `J = 1/2`, `b = 1/4`, and `J_sigma = 1/8`.

One common weight vector is used for:

- Three raw face-dephased intensities at `c/J_sigma = 1, 2, 3`, with Lorentzian half-width `gamma/J_sigma = 0.05` (`1/160` in working units).
- The means of face 0, face 11, and their product, in that order.
- The four ordered joint populations `++`, `+-`, `-+`, and `--` of those pre-insertion face signs.

The filter is an unsubtracted, peak-one Lorentzian. The band enclosures retain the elastic contribution and use absolute spin intensity per initial orbit, whose integrated local strength is `1/4`. For row `i`, signs `(s_i,t_i)` give moments `(s_i,t_i,s_i*t_i)` and one unit joint-population indicator in the corresponding ordered cell. The checker weights those quantities and both endpoints of every band with the same `w_i`.

The downward observable is the unbroadened full-spin mass at signed transfer `E_final - E_initial < 0`, with equal additive block constants. Every occupied row has interval `[0,0]`, so the mixture has exactly zero downward mass. Nonnegative spectral weights and this attained zero establish the minimum for this local probe; they do not establish zero downward response for every observable.

## Exact input representation

Every numerical endpoint and weight is an explicit rational string: an integer or `numerator/denominator`. There is no endpoint rounding in the export.

| JSON field | Interpretation |
|---|---|
| `orbits[].weight` | Exact positive common-mixture weight |
| `orbits[].face_signs` | Ordered pre-insertion signs of faces 0 and 11; sufficient to reconstruct the three moments and four cell indicators |
| `orbits[].band_intervals` | Three full-response `[lower, upper]` pairs, in increasing probe-center order |
| `orbits[].downward_interval` | Unbroadened negative-transfer mass enclosure |
| `orbits[].archive_qualification`, `initial_ground_rank` | Retained enclosure status and occupied initial ground-space rank; each supplied row has rank one |
| `targets.bands` | Three physical raw face-dephased target intervals |
| `targets.moments`, `targets.joint_populations` | Physical stored-ray moment and ordered-cell target intervals |
| `saved_witness` | Saved numeric reference values and slacks, compared after exact-number normalization |

Orbit band/downward endpoints and physical band targets preserve the exact rational values of their binary64 representations, equivalent to `Fraction.from_float(float(endpoint))` in the archived interpretation. Physical population and moment target strings instead preserve exact decimal rationals, equivalent to `Fraction(decimal_string)`. Converting all displayed decimals by one rule would change the input problem.

For each predicted interval `[a,b]` and target `[L,U]`, the checker calculates `[a-L,U-b]` and requires both slacks to be nonnegative. This verifies whole-interval containment, not merely intersection. Moments and cells are exact points represented as degenerate intervals. Every reconstructed value and slack must also equal its saved reference after exact rational normalization.

The physical targets refer to a normalized stored ray. Their interpretation as lowest-state observations inherits the rank-one ground-space completeness and true-separation conditions. This supplied-scalar replay retains those conditions; it does not regenerate the spectral enclosures or supply the missing physical state vectors, elastic matrices, or complete production runtime.
