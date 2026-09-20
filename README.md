# Multiscale microbial network analysis

An English, reproducible reorganization of the original fungal-bacterial analysis.
The workflow imports the selected clustering result and performs directed network inference at two
scales, GLMY persistent path homology, curvature, and soil-trait importance.


## Two network scales

- **Coarse scale:** a node is a module; directed edges describe fitted effects
  between module abundance curves.
- **Fine scale:** a node is an OTU within a selected module; directed edges describe
  fitted effects among those OTUs. The procedure applies to every module. No
  particular module is used as the definition of fine-scale analysis.

The two scales use the same inference function. The default fine-scale setting is
`fine_modules = "all"`; change it to an integer vector to process any subset.
OTU identifiers are retained across treatments and analysis stages.

## Layout

```text
config.R                 Paths, filters, 40-module default, and fitting settings
run.R                    R workflow entry point
R/preprocess.R           Sample matching, GBM, CLR, scaling, power curves
R/clustering.R           Required clustering stage: load and validate k40.RData
R/network.R              Shared coarse/fine network inference and plots
R/importance.R           Elastic net, random forests, consensus, correlations
python/homology.py       External GLMY adapter and barcode export
python/plot_barcodes.py  H0-H2 barcode visualization
python/curvature.py      Signed-control Forman and projected LLY curvature
python/analyze_networks.py  Batch homology/curvature for both scales
data/                    Input workbooks (locally included; ignored by Git)
tools/                   Local GLMY executable (ignored by Git)
```

## Setup

Use R 4.5 or a compatible R release and Python 3.10 or later. Run commands from
this repository directory. Dependencies are installed explicitly, never at import.

```sh
Rscript scripts/install.R
python -m venv .venv
# Activate .venv using your shell's activation command, then:
python -m pip install -r requirements.txt
```

Place `fungi_otu.xlsx`, `bacteria_out.xlsx`, and `soil_properties.xlsx` in `data/`.
See [the data contract](data/README.md). The local delivery already contains copies
of the supplied workbooks and GLMY executable. The Git ignore rules let the code
be shared separately from these inputs. See [the GLMY note](tools/README.md).

## R workflow

```sh
Rscript run.R prepare
Rscript run.R clustering
Rscript run.R coarse
Rscript run.R fine
Rscript run.R importance
# Optional OTU importance within every configured module:
Rscript run.R fine-importance
```

`Rscript run.R all` runs preparation, module import, coarse and fine networks, and
coarse-scale importance. Fine-scale importance remains an explicit optional step.
The clustering stage in `R/clustering.R` is mandatory, not an optional component:
every downstream stage loads and validates the selected 40-module result before
analysis, including when a downstream command is run on its own. This stage reads
the existing result rather than estimating a new clustering.
An alternative configuration can be supplied as the second argument, for example
`Rscript run.R fine config_custom.R`. Use a separate `output_dir` for sensitivity
analyses. Re-run preparation when inputs or preprocessing settings change.

Module import verifies all 40 IDs and matches the 563 OTUs by identifier, never
by row position. The original module numbering is retained. The clustering file's
checksum and original BIC/likelihood fields are recorded for provenance, without
reinterpreting or recalculating them.

## Homology and curvature at either scale

Each network directory contains `edges.csv`, `nodes.csv`, `network.rds`, and
`decomposition.pdf`. The same Python commands accept a coarse or fine directory.

```sh
python python/curvature.py results/networks/coarse/CKN/edges.csv --output results/curvature/coarse/CKN
python python/homology.py results/networks/coarse/CKN/edges.csv --exe tools/GLMY.exe --edge-filter all --filtration signed --output results/homology/coarse/CKN/all
python python/plot_barcodes.py results/homology/coarse/CKN/all/barcodes.csv --output results/homology/coarse/CKN/all/barcodes.pdf
```

For fine networks, use `results/networks/fine/<treatment>/M<module>/edges.csv`.
Repeat homology with `--edge-filter positive` and `--edge-filter negative` and
separate output directories. `--shift` can fix a common shift across compared
networks; otherwise the smallest nonnegative shift is chosen separately. A
different shift changes edge entry times relative to vertex births at zero.
Inspect `metadata.json` before comparing barcodes across networks.

The GLMY adapter runs the executable in a new temporary directory, captures logs,
checks its exit code and JSON, and preserves isolated nodes. Its executable is
Windows-specific; other systems need a compatible build. The wrapper does not
substitute clique homology for directed path homology.

To process all exported networks in one batch (including all three edge-sign
filters and barcode plots), use:

```sh
python python/analyze_networks.py results/networks --output results --kind all --exe tools/GLMY.exe
```

Use `--kind curvature` for curvature alone or `--kind homology` for homology alone.
Batch homology of large fine networks may take substantial time; the timeout is
per network and filter and can be changed with `--timeout`.

## Soil indicators

Analyze CMF, NMF, and PMF from the supplied workbook. SMF is absent and is neither
imputed nor fabricated. Coarse predictors are module means on the retained real
samples; fine predictors are OTU values. Smoothed points are never used as soil
replicates. Consensus combines elastic-net stability and random-forest permutation
importance. These are exploratory predictive associations, not causal effects or
independent validation of network edges.

## Validation and scientific decisions

```sh
Rscript tests/test_core.R
python -m unittest discover -s tests -p "test_*.py"
```

Read [METHODS.md](docs/METHODS.md), [AUDIT.md](docs/AUDIT.md), and
[VALIDATION.md](docs/VALIDATION.md) before interpreting results. Numerical fixes
can change published/historical outputs; this is not a bitwise legacy replica.

## GitHub preparation

All maintained code, comments, configuration, and documentation are in English.
Raw workbook metadata and sample IDs are preserved. The source input directories
are unchanged. Commit this directory as its own repository; review the data and
binary ignore rules if you intend to distribute those files too. Add the actual
authors, study citation, and your chosen license before public release. No license
or authorship has been invented and no remote repository has been created.
