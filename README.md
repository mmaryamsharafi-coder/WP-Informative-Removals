# WP-Informative-Removals
# Weibull-Poisson Distribution under Progressive Type-II Censoring with Informative Removals

This repository contains the data and R codes for the paper:

**Classical and Bayesian Inference for the Weibull-Poisson Distribution under Progressive Type-II Censoring with Informative Removals**  
by A. Falah Hasan and M. Sharafi (Corresponding author)

## 📁 Repository Structure
## Simulation Codes
- `R/02_simulation_n60_m15.R`: Simulation for n=60, m=15 with alpha=1.5, theta=1.0, lambda=3.0.
  ## Real Data Analysis
- `R/05_real_data_analysis.R`: Full analysis of bladder cancer data.
- Data: `data/bladder_cancer.csv`
## ⚙️ Requirements

- R version 4.0 or higher
- Required R packages:
  - `maxLik` (for MLE)
  - `MCMCpack` (for MCMC)
  - `coda` (for convergence diagnostics and HPD intervals)
  - `numDeriv` (for numerical derivatives)
  - `pracma` (for incomplete gamma functions)
  - `Rgof` (for goodness-of-fit tests)

You can install them with:
```r
install.packages(c("maxLik", "MCMCpack", "coda", "numDeriv", "pracma", "Rgof"))
