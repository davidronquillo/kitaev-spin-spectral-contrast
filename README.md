# Kitaev spin spectral contrast

Data and analysis scripts accompanying **“Finite-size reversal of a multiparticle correction to local spin spectra in a flux-disordered Kitaev model”**, by **David C. Ronquillo**.

The companion is currently private during author review. Public access is planned to coincide with the preprint. Licensing remains subject to the author's release decision; this preparation does not grant a new license. No DOI, arXiv identifier, or publication status is assigned here.

## Reproduce the numerical reductions

Use Python 3 and run from this repository's root. These commands require only the Python standard library; no credentials, original project installation, or network access are needed.

```sh
python3 scripts/scalars/reproduce_scalar_tables.py
python3 scripts/relative_contrast/extract_original_means.py
python3 scripts/relative_contrast/reduce_relative_contrast.py
python3 scripts/relative_contrast/verify_relative_contrast.py
python3 scripts/check_reported_diagnostics.py
python3 scripts/summarize_validation.py
```

Expected results:

- All 56 Table S2 intervals, all 8 Table S6 intervals, and both complete-N24 Figure 3 annotations agree with the supplied targets.
- Four original-record mean-input rows are reconstructed exactly. The relative-contrast reducer reproduces eight original/printed-input intervals, and a separate directed-Decimal calculation agrees with each ratio.
- All 96 Table S1 ordinate boxes are compatible. The 170 saved window bounds yield 37 strictly higher-sector, 0 strictly one-particle, and 133 unresolved classifications. The physical-field grid and discrete displacement sets agree.
- Saved validation records summarize 8 N24 fixtures, 620 higher-sector checks, 406 transform checks, and 6 direct-comparison sources.

Outputs are written under `outputs/`. Accepted expected files are retained under `results/`; the data, exact weights, endpoint conventions, and source IDs remain unchanged by these commands. To compare every generated report with its expected file:

```sh
python3 - <<'PY'
from pathlib import Path
expected = sorted(Path('results').rglob('*'))
count = 0
for file in expected:
    if file.is_file():
        generated = Path('outputs') / file.relative_to('results')
        assert generated.read_bytes() == file.read_bytes(), str(file)
        count += 1
print(f'PASS: {count} expected result files match byte-for-byte.')
PY
```

The complete file-to-paper mapping, exact outputs, and reproduction scope are in [docs/reproduction_map.md](docs/reproduction_map.md). Input columns, normalization, interval arithmetic, conditioning, and probability conventions are documented in [docs/scalar_schema.md](docs/scalar_schema.md). [docs/file_manifest.json](docs/file_manifest.json) lists the exported scientific files and their SHA-256 hashes.

## Inspect methods and make local plots

Selected method code, recorded validation results, graph conventions, the editable vector schematic, and diagnostic plotting code are included. Their dependencies and specific capabilities are described in [docs/methods_and_plotting.md](docs/methods_and_plotting.md).

Optional plotting dependencies are listed in `requirements-plotting.txt`. With those installed:

```sh
python3 scripts/plotting/draw_local_quench_schematic.py
python3 scripts/plotting/plot_saved_diagnostics.py
```

The plotting commands write local files under `outputs/plots/`. The schematic is conceptual. The diagnostic plots visualize saved scientific inputs; the separate manuscript delivery retains the approved quantitative artwork. Generated figure files and manuscript binaries are not stored in this repository.

## Organization

| Directory | Contents |
|---|---|
| `data/` | Numerical enclosures, exact probability weights, source identities, plot inputs, and validation snapshots |
| `scripts/` | Standalone reducers, separate arithmetic checks, selected method code, and local plotting scripts |
| `results/` | Reviewed expected text/CSV/JSON outputs |
| `docs/` | Schemas, paper mapping, capabilities, and a neutral hashed file manifest |
| `outputs/` | Locally generated results and plots; excluded from repository publication |

The numerical collection is also supplied directly with the author-review materials, so recipients can run it without live repository access. Editable manuscript sources and the rebuildable typesetting package are separate deliverables.
