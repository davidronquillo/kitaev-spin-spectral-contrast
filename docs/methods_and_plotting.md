# Methods, saved validation, and plotting

## Saved validation reductions

Run from the repository root using Python 3.9 or newer:

```sh
python3 scripts/summarize_validation.py
```

The command uses the standard library, reads only `data/validation/`, and writes
`outputs/validation/validation_summary.json` and `outputs/validation/direct_summary.csv`.
The expected results are in `results/validation/`.

| Input | Contents and checks |
| --- | --- |
| `fixtures_24_site.json` | Eight selected 24-site orbit fixtures: occupation probabilities, component errors and allowances, intersecting Fock/Gaussian observable boxes, and first moments. |
| `high_sector_controls.json` | 620 saved checks of the sector estimator, including coherent amplitudes for five- and seven-particle sectors. |
| `transform_controls.json` | 406 saved independent-transform checks at the two filter settings. |
| `core_controls.json` | Low-sector precision and Gaussian-chart controls, plus an assembly control whose status remains unresolved. |
| `direct_precision.csv` | Eight kernel comparisons for each of six selected geometry/source pairs, including stochastic, direct, and combined numerical-enclosure half-widths. |
| `direct_selection.csv` | The recorded source-selection qualifications for those six pairs. |
| `geometry_metadata.json` | Torus periods, site counts, canonical labeling, source edge, link encoding, and available sampling specifications. |

The summary checks the recorded statuses, interval intersections, errors versus
allowances, and comparison availability, then reports maxima and counts. The
unresolved assembly control is retained explicitly; a passing reduction does not
change its status. Scientific geometry/source identifiers, seeds, decimal
endpoints, diagnostic values, and selection qualifications are preserved.

These commands reproduce reductions of saved scalar records. Recreating the
underlying spectral calculations or validation experiments requires the full
solver runtime, state and transform records, and physical-response evidence not
contained here. The method selections below do not constitute that complete
runtime.

## Selected implementation modules

| Module in `scripts/methods/` | Purpose |
| --- | --- |
| `filter_kernels.py` | Ten prescribed rational Lorentzians and 24 paired linear kernels, including coefficient ranges and floating-point allowances. |
| `numerical_certificates.py` | Complete-frame source and energy-operator certificates, provided a represented eigensystem and the necessary source records. |
| `gaussian_chart.py` | Phase-consistent Gaussian-chart construction, source contractions, occupation spectrum, and error diagnostics. |
| `sector_importance.py` | Reference importance estimator for an odd pure Gaussian chart in the five- and seven-particle sectors. |
| `geometry.py` | Canonical quotient-honeycomb labeling, oriented links and faces, parity conventions, and hopping construction for a supplied link assignment. |

These selected modules preserve the numerical algorithms. NumPy and SciPy are
needed to use the numerical modules. The sampling-domain identifier and external
reference-labeling map are not supplied; the geometry selection does not
regenerate the original background assignments.
The estimator's synthetic-test entry point is retained for inspection; it is not
part of the saved-record reproduction commands. No historical solver or
synthetic validation experiment was run during export.

## Optional local drawing

The tested drawing environment is Python 3.9.18, NumPy 1.24.4, SciPy 1.11.2, and
Matplotlib 3.7.2. Install `requirements-plotting.txt` in a Python 3.9–3.11
environment, then run from the repository root:

```sh
python3 scripts/plotting/draw_local_quench_schematic.py
python3 scripts/plotting/plot_saved_diagnostics.py
```

Both scripts accept `--output-dir PATH`; the default is `outputs/plots/`.
Generated graphics are local outputs and are excluded from version control.

`draw_local_quench_schematic.py` produces PDF, SVG, and PNG versions of the
conceptual local-quench diagram, with the unchanged approved panel geometry.
Its bar partition and depth triplets are illustrative layout coordinates, not
numerical results, representative spectra, or pointwise medians.

`plot_saved_diagnostics.py` produces five PDF diagnostic views of the supplied
scalar tables: paired depth reversal, absolute spin intensities, filter ratios
and hard-window fractions, the prespecified physical benchmark, and the
physical field displacement. Its plotting operations use saved ordinates,
interval midpoints, differences, and classifications. The code is an earlier
plotting implementation and does not exactly regenerate the approved manuscript
artwork. The manuscript's approved vector figures remain in the separately
delivered typesetting package.

Both drawing commands were exercised locally using only supplied inputs; all
five diagnostic files and the conceptual schematic were produced. OpenAI
ChatGPT and Codex assisted with plot and figure preparation. No license grant
is added by this export; existing attribution is preserved.
