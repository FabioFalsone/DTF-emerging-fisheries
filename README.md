# DTF emerging fisheries

This repository contains the reproducible materials associated with the manuscript:

> **A Decision-Tree Framework for the Sustainable Management of Emerging Fisheries and Fishing Innovations**  
> Falsone et al.

The Decision-Tree Framework (DTF) combines an automated landing-based screening pathway for eligible time series with an expert-supported pathway for shorter, non-consecutive, or otherwise data-limited series.

## Repository contents

```text
DTF-emerging-fisheries/
├── app/
│   ├── Shiny_App_DTF.R
│   └── example_data.csv
├── simulations/
│   └── DTF_simulation_reproducible_final.R
├── figures/
│   ├── Figure_1.png
│   ├── Figure_2.png
│   ├── Figure_3.png
│   ├── Figure_4.png
│   ├── Figure_5.png
│   ├── Figure_S1.png
│   ├── Figure_S2.png
│   ├── Figure_S3.png
│   ├── Figure_S4.png
│   └── Figure_S5.png
├── .gitignore
├── LICENSE
└── README.md
```

## Shiny application

The application implements the DTF questionnaire, the automated landing-based screening component, the expert-supported pathway, diagnostic outputs, decision-tree visualisation, and export functions.

### Required R packages

```r
install.packages(c(
  "shiny",
  "bslib",
  "visNetwork",
  "DT",
  "ggplot2",
  "gridExtra",
  "DiagrammeR",
  "DiagrammeRsvg",
  "rsvg"
))
```

### Run the application locally

From the repository root:

```r
source("app/Shiny_App_DTF.R")
```

The application also includes embedded example data. The same dataset is provided separately in `app/example_data.csv`.

A deployed version of the application is available at:

https://decisiontreetool.shinyapps.io/dtf-shiny/

## Simulation workflow

The simulation script reproduces the two-phase workflow described in Section 2.6 of the manuscript:

1. calibration using 20,000 simulated time series;
2. full-factorial performance evaluation using 100 stochastic replicates per combination.

The workflow evaluates alternative scenarios, time-series lengths, observation-noise levels, AR(1) autocorrelation, effect magnitudes, breakpoint positions, threshold configurations, and four composite-score weighting schemes.

### Required R packages

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr"))
```

The optional package `patchwork` is used only to assemble the three panels of main-text Figure 2 into a single image:

```r
install.packages("patchwork")
```

### Run the simulations

From the repository root:

```r
source("simulations/DTF_simulation.R")
```

The full workflow is computationally intensive. Output tables and figures are written to `tables/` and `figures/` in the active working directory.

## Input-data format

Landing data must be supplied as a CSV file containing at least:

```text
year,landing
2004,129
2005,237
```

Automated screening requires at least six unique and consecutive annual observations. Shorter or non-consecutive series are directed to the expert-supported pathway.

## Figures

The figures/ directory contains main-text Figures 1-5 and Supplementary Figures S1-S5.

## Reproducibility note

The code and figures in this repository accompany the manuscript and are intended to support transparency, independent evaluation, and future development of the DTF. The automated screening output is an early-warning diagnostic and is not a substitute for formal stock assessment.

## License

This repository is distributed under the MIT License. See `LICENSE`.
