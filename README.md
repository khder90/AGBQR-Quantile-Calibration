# AGBQR: Quantile-Specific Learning-Rate Calibration

Research code and analytical data supporting *Adaptive Gibbs Posteriors for Quantile Regression with Quantile-Specific Learning-Rate Calibration: Theory, Computation, and an R Implementation*.

## Contents

- `R/agbqr_core.R`: quantile loss, HAC covariance estimation, scalar learning-rate calibration, posterior sampling, interval construction, and simulation designs.
- `R/run_full_analysis.R`: Monte Carlo experiments, sensitivity analyses, Syrian inflation estimates, conditional expanding-window prediction, tables, and figures.
- `R/validate_mcmc.R`: four-chain numerical checks for 17 selected simulation and empirical settings.
- `R/smoke_test.R`: a short deterministic check of data integrity and core numerical operations.
- `data/Syrian_Inflation_1961_2025.xlsx`: annual analytical input in one worksheet named `Data`.
- `data/Syrian_Inflation_1961_2025.csv`: the same observations as a machine-readable CSV file.
- `validation/MCMC_Stability_FINAL_RESULTS.xlsx`: archived summary of the longer-chain diagnostic exercise supplied with the study materials.

## Input data

The workbook contains 65 annual observations, 1961–2025, and five variables:

| Column | Description |
| --- | --- |
| `YEAR` | Calendar year |
| `EXR` | Nominal exchange rate in the original input series |
| `TRDGDP` | Trade as a percentage of GDP |
| `INF` | Annual inflation rate, percent |
| `M2GDP` | Broad money as a percentage of GDP |

The original input workbook was `INFc.xlsx`, with headers `YEAR`, `EX`, `OPEN`, `INF`, and `M2\GDP`. Column names were standardized for the repository; the 65 rows and their numeric values were preserved. Data processing is performed by the analysis script: the nominal exchange rate is log-transformed, then the empirical response and regressors are standardized. Expanding-window comparisons use normalization constants estimated only from the training period.

The input workbook itself does not identify individual statistical releases or source vintages. The manuscript's data-source citations should therefore be checked against the underlying official publications independently of this repository.

## Setup

Install R and the required packages. From the repository root, run:

```r
install.packages(c(
  "readxl", "quantreg", "dplyr", "tidyr", "ggplot2",
  "openxlsx", "coda", "tseries", "posterior"
))
```

Run the quick check:

```r
source("R/smoke_test.R")
```

Recreate the numerical study:

```r
source("R/run_full_analysis.R")
```

`R/run_full_analysis.R` is configured for the published full study. It uses 300 Monte Carlo replications for each process–sample-size design, five quantiles, and four sample sizes, together with the stated MCMC settings. These computations are extensive. A short preliminary run can be performed by changing `RUN_MODE` to `"pilot"` in that script; pilot output must not be substituted for the full results.

The program writes tables and figures to the `results/` directory. Completed Monte Carlo blocks are retained so interrupted full runs can resume. The input and source files are not overwritten. Recorded seed settings, runtime settings, and an R session record are saved with the outputs.

For independently generated longer-chain diagnostics of the 17 selected settings, run:

```r
source("R/validate_mcmc.R")
```

The validation script writes chain files, diagnostic summaries, and empirical trace plots to `results/validation/`. Its seeds are fixed, but it generates its own Monte Carlo realizations; independent results need not be numerically identical to the archived summary. A numerical comparison should account for these differences rather than treating independent draws as an exact replay.

## Methods and interpretation

For each quantile, the analysis uses a check-loss Gibbs posterior and a prior-aware scalar learning-rate criterion. The calibration target is a coefficient-level HAC sandwich covariance matrix built from the same curvature and score covariance used for the Wald-HAC benchmark. MCMC is performed with a random-walk Metropolis kernel and burn-in-only proposal-scale adaptation.

The scalar criterion balances marginal variance discrepancies. It does not impose full covariance equality or nominal joint coverage. The time-ordered empirical test exercise is a *conditional* assessment using realized test-year regressors; it is not a real-time forecasting experiment.

## Software distinction

The published R package [AGBQR on CRAN](https://cran.r-project.org/package=AGBQR) and the expanded research scripts in this repository are distinct software resources. The scripts here implement the analysis described in the revised manuscript and do not represent a new CRAN release.

## Research contact

Khder Alakkari (corresponding author).
