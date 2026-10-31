# Reproducing the interaction-tree analyses

Code for **Selective inference for pruned interaction trees with known unequal
variances in randomized trials**.

## Requirements

R (tested with version 4.6.1). For the simulation summaries and figures, install:

```r
install.packages(c("data.table", "ggplot2"))
```

The ACTG analysis and its figures use base R. Python and Ghostscript are not
required.

## Run

Open a terminal in this directory. Start with a small execution check:

```sh
Rscript reproduce.R quick
```

This runs a small simulation and the ACTG analysis, then creates tables and PDF
figures. The small simulation is an execution check, not the paper's experiment.

To reproduce the paper:

```sh
Rscript reproduce.R simulation --workers=4
Rscript reproduce.R actg
```

Use `--workers=1` on Windows. The default is one worker. Alternatively,
`Rscript reproduce.R all --workers=4` runs both analyses. The full simulation
is computationally intensive; it saves chunks so interrupted runs can resume.

Generated results are written to `results/`; tables and figures are written to
`results/` and `figures/`. Simulation outputs distinguish `main` from `quick`.
The full study uses 36,000 independent base trials, the paper's five effect
settings, sample sizes 100/200/400/800, variance ratios 0.25/1/4, maximum depths
1/2/3 and pruning penalties 0/4/8. Interval endpoints use the first 500 trial
indices as specified in the paper.

## Files

- `simulation/`: data generation, tree growth/pruning, selective inference,
  simulation execution, summaries and four figures.
- `actg/`: ACTG preparation, the five fixed-penalty analyses, two figures and
  the main results table. The training/holdout split is regenerated from seed 20261028.
- `data/`: the public ACTG175 data, variable dictionary and GPL-2 license.

The ACTG data come from **speff2trial 1.0.5**
([package source](https://cran.r-project.org/package=speff2trial)). The analysis
uses arms 1 and 2: 1,046 participants, with a 523/523 training/holdout split.
The original variable descriptions are in `data/ACTG175.Rd`; the data-package
license is in `data/LICENSE.GPL-2`.

The exact theory assumes known Gaussian covariance. The common-scale
substitution experiment and ACTG analysis with estimated arm variances have
the limited interpretations stated in the paper.
