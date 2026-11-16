# Scalar inputs and conventions

The data are saved numerical enclosures and exact probability metadata used in the paper. CSV endpoint strings are retained as supplied. The reducers recover the original binary64 endpoints exactly with `Fraction.from_float(float(cell))` or `Decimal.from_float(float(cell))`; they then use rational or directed-decimal arithmetic. Reinterpreting those cells as exact decimal measurements would change the numerical convention. Exact orbit numerators and denominators are integers.

## Shared scientific definitions

`P` is the intensity-matched one-particle control; `S` is the full response. Both use the same frozen background and local spin quench. The filtered intensity uses a peak-height-one Lorentzian. `g005` and `g010` denote filter half-widths `gamma/J_sigma = 0.05` and `0.10`. The `width` columns in these tables refer to this filter parameter unless a filename or field explicitly identifies an enclosure or integration window.

For a triplet at `c/J_sigma = 1, 2, 3`, `d(U) = [U(1)+U(3)]/2-U(2)` and `D = d(P)-d(S)`. Thus positive `D` is shallowing and negative `D` is deepening. `X(c)=S(c)-P(c)`. Normalized total local intensity is one; absolute spin intensity is one quarter of that normalization. Output fields explicitly state the factor of `10^3` where applied.

The finite-cohort signed `D` kernel evaluates `U(2)-[U(1)+U(3)]/2 = -d(U)`. Consequently its saved `X = S-P` equals the manuscript's paired `D`. Do not negate that saved correction a second time.

## Finite cohorts

`data/scalars/fresh_endpoint_boxes.csv` contains 8,759 rows covering all 1,920 source IDs, including shared even sources. The pairing key is `(geometry, id)`; the kernel adds a filter setting and observable. No sources are discarded on numerical precision grounds.

| Geometry ID | Sites | Shape | Eligible odd sources | Shared even sources |
|---|---:|---|---:|---:|
| C3n02 | 96 | Tilted | 164 | 156 |
| C3n04 | 384 | Tilted | 156 | 164 |
| C3n08 | 1,536 | Tilted | 175 | 145 |
| C3n12 | 3,456 | Tilted | 177 | 143 |
| C3n16 | 6,144 | Tilted | 170 | 150 |
| R24x32 | 1,536 | Rectangular | 135 | 185 |

Each geometry has IDs 0 through 319. Eligible rows have `eligibility=known`, `parity=odd`, and eight kernel rows: two filter settings times `c180`, `c200`, `c220`, and `D`. The first three centers are `c/J_sigma=1.8,2.0,2.2`. Excluded even rows have an empty kernel and a shared response with zero correction; their blank endpoints are not missing odd-source observations.

| Column family | Meaning |
|---|---|
| `geometry`, `N`, `id` | Geometry, site count, and stable within-geometry source ID |
| `eligibility`, `parity`, `core_status`, `disposition`, `issues` | Source eligibility, parity, and recorded numerical status |
| `offset_lo`, `offset_hi` | Saved energy-transfer offset enclosure |
| `kernel` | Filter setting and central, neighboring, or signed-depth kernel |
| `P_lo/hi`, `S_lo/hi` | Matched-control and full-response kernel enclosures |
| `X_lo/hi` | Paired full-minus-control kernel enclosure; for the signed kernel this is `D` |
| `C1_lo/hi` | Unrenormalized one-particle contribution, with its original mass |
| `w1_lo/hi` | One-particle mass enclosure |
| `integration_tail_halfwidth`, `integration_tail_scope`, `error_inventory` | Recorded integration-tail allowance and numerical-error inventory |

`protected_medians.csv` contains the accepted source-level rank reductions. It retains geometry, endpoint, eligible/ambiguous counts, lower and upper bounds, fixed means/medians, error allocation, coverage, and rank-completion metadata. For the population-median checks the central and depth endpoints use error allocation `1/800`, and neighboring-center endpoints use `1/2400`. The reducer computes the corresponding binomial ranks exactly; it does not choose new ranks.

## Exact 24-site probability distribution

The six `data/scalars/n24/` files retain all 8,192 gauge orbits and the 3,892 eligible odd sources. The reference torus has translation columns `(2,2)` and `(-2,4)`. The negative-link probability is the exact rational encoded in `n24_exact_law.json`, equal to `66871911161826185187153288884631757060/2^128`. Its gauge-orbit probabilities include the Wilson sectors.

| File | Contents |
|---|---|
| `n24_exact_law.json` | Common denominator, probability digest, cycle/face/Wilson masks, edge ordering, and gauge convention |
| `n24_probability_classes.csv` | Exact integer numerator for each probability class |
| `n24_orbit_catalog.csv` | Stable orbit `id`, `classification`, allowed initial-ground-space rank, probability class, and saved probability bounds |
| `n24_odd_scalar_boxes.csv` | One row for each eligible orbit and filter setting; `D`, `X18`, `X2`, `X22`, `w1`, and the `P1..3`, `S1..3` enclosures |
| `n24_numerical_error_metadata.csv` | Saved per-orbit/filter numerical-error contributions, including energy, frame, mass, and projector allowances |
| `n24_published_reference_summary.json` | Accepted conditional-odd mean and median targets used for comparison |

The catalog uses `odd_modified` for eligible odd sources and `shared_even` for shared even responses. Exact class numerators divided by the common denominator give orbit probabilities; all are positive and sum to one. The reducer verifies the digest of the ordered fixed-width numerators. Weighted median lower endpoints cross cumulative weight at `>=1/2`; upper endpoints cross at `>1/2`, retaining exact-half median sets. Odd conditioning divides by the exact odd probability. The unconditional median correction is zero because the shared even class has probability greater than one half.

## Relative mean contrast

`data/relative_contrast/original_mean_inputs.csv` has four rows: N24 and N6144 at both filter settings. `M_lo/hi` enclose the odd-source mean paired correction, and `F_lo/hi` enclose the full odd-source mean depth. `eligible` and `population` retain the conditioning: exact probability-weighted odd sources at N24 and the realized eligible-odd cohort at N6144. These endpoints are in `10^-3` absolute spin intensity.

`printed_input_rows.csv` is a separate calculation using rounded paper inputs. Its `S1_lo/hi`, `S2_lo/hi`, and `S3_lo/hi` give the printed full-response ordinates; `M_lo/hi` give the printed paired mean. The two routes must not be conflated: original-record bounds can be tighter.

For each route the reducer forms `R=M/F`, with `F>0`, taking all four interval-corner quotients. It reports `100R` with outward rounding to 0.1 percentage points. The `R_percent_lo/hi` columns already contain percentage numbers; no second factor of 100 belongs in their display. The independent checker uses directed 70-digit Decimal division rather than importing the rational reducer.

## Diagnostics and plot inputs

The separate `data/seven_orbit_witness.json` schema is documented in [seven_orbit_witness.md](seven_orbit_witness.md). Its explicit rational strings preserve binary64 orbit/band endpoints and exact-decimal physical population/moment targets without applying a single decimal interpretation to both. Pre-insertion face signs determine each orbit's moments and ordered joint-population indicators. One common exact weight vector reconstructs all observables; the checker requires full interval containment and zero unbroadened downward mass for the specified local probe. Raw intensities retain elastic weight, the signed transfer convention, and the documented stored-ray/ground-space qualifications.

All scalar units and geometry/source labels are retained in the file headers. The following compact tables support saved-value checks and local plots.

| File | Meaning |
|---|---|
| `absolute_model_ordinates.csv` | P/S absolute-intensity enclosures at the three probes, for each geometry and estimand |
| `N6144_sign_frequencies.csv` | Definite/possible counts for central gain, deepening, and their conjunction |
| `endpoint_mean_table.csv` | Higher-sector, removed one-particle, and net mean contributions at the center and in contrast |
| `odd_depth_means.csv` | Fixed-cohort mean enclosures and inherited pointwise uncertainty intervals |
| `odd_depth_intercepts.csv` | Saved intercept intervals for specified finite-size forms and size windows |
| `physical_050_comparison.csv` | Fixed benchmark model/control and physical-response enclosures at `b/J=0.50` |
| `physical_field_curves.csv` | 401 raw-response probe points per field at `b/J=0.50,0.55`, from `c/J_sigma=1.5` to `3.5` in increments of `0.005` |
| `physical_minimum_candidates.csv` | Saved minimizing-grid candidates by field and raw/inelastic layer |
| `window_fraction_bounds_N6144.csv` | All 170 eligible source IDs, constituent window/response bounds, and the difference of higher-sector and one-particle fractions |
| `window_fraction_counts.json` | Accepted window-classification counts and parameter settings |
| `fixed_class_bound.json` | Existing four-class restricted-reweighting bound and its inputs |
| `plot_summary.json` | Reviewed summary intervals and counts consumed by the diagnostic plotting script |

In the window table, `P_n/P_b` and `h_n/h_b` denote narrow/broad filtered one-particle and higher-sector inputs; `R_P` and `R_h` enclose broad/narrow filter ratios, while `f_P` and `f_h` enclose the corresponding within-window fractions defined in the paper. The window radius is a hard integration-window half-width, distinct from the Lorentzian filter settings. `f_h_minus_f_P_lo/hi` determine `forced_order`: higher-sector when strictly positive, one-particle when strictly negative, and unresolved otherwise.

The four `paper_tables/*.txt` files are numeric table excerpts in TeX-style row syntax used only as acceptance targets, not manuscript sources. They preserve the printed bounds and source-ID lists exactly. The diagnostic checker confirms that all 96 Table S1 boxes contain the saved ordinates, classifies all 170 window rows, checks the forced-ID list, and reproduces the discrete field-displacement sets. It does not treat the Cartesian set of allowed grid displacements as a continuous interval.
