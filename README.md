
# CompetingRiskR

## Overview
`cmpp` is an R package designed to facilitate parametric modeling and inference for cumulative incidence functions in competing risks scenarios. The package implements methods discussed in seminal works by Jeong and Fine (2006, 2007), with efficient computation powered by Rcpp for high-performance applications.

## Key Features
- Direct parametric modeling of cumulative incidence functions using the Gompertz and Weibull distributions.
- Maximum likelihood estimation for parametric regression models.
- Support for goodness-of-fit tests and evaluation of proportional hazards or proportional odds assumptions.
- Fast computation using Rcpp integration.
- Application to clinical and epidemiological studies with competing risks.

## Background
Competing risks occur when multiple types of events prevent the observation of a particular event of interest. Traditional survival analysis often misrepresents such data, as it fails to account for competing risks. This package provides parametric methods for cumulative incidence functions (CIF), offering a more direct and interpretable analysis than traditional cause-specific hazard models.

## Implemented Methods
1. **Gompertz-Based Cumulative Incidence Models**:
   - Designed for scenarios with improper distributions (plateau in the cumulative incidence curve).
   - Applicable to long-term survival studies, such as breast cancer clinical trials.

2. **Regression Models**:
   - Proportional hazards and proportional odds models.
   - Flexible generalized odds rate models.
   - Covariate-adjusted CIF estimation.

3. **Likelihood-Based Inference**:
   - Full maximum likelihood estimation.
   - Delta method for confidence intervals and variance estimation.

## Example Use Case
The package is inspired by methods used to analyze data from clinical trials, such as the National Surgical Adjuvant Breast and Bowel Project (NSABP), where cumulative incidence functions are compared across treatment groups or adjusted for patient covariates like age and tumor size.

## Installation
To install the package:
```R
# Install from GitHub
devtools::install_github("stats9/cmpp")
```

## Usage
```R
library(cmpp)

# Load example data
data(dat)
names(dat)
feat <- dat[, -c(match(c('id', 'time', 'event', 'cause_burn', 'd1', 'd2', 'cause_hotObject3'), names(dat)))]  
timee <- dat[['time']]
d1 <- dat[['d1']]
d2 <- dat[['d2']]
feat2 <- feat |> data.matrix()

# Initialize the Cmpp model
Initialize(feat2, timee, d1, d2, 1e-10)

# Estimate parameters using LBFGS++
initial_params <- c(0.001, 0.001, 0.001, 0.001)
params <- estimate_parameters(initial_params)
print(params)

# Compute Hessian and other metrics
hessian <- compute_hessian(params$par)
print(hessian)

# Bootstrap variance estimation
results <- bootstrap_variance(feat2, timee, d1, d2, initial_params, n_bootstrap = 500)
print(results$variances)
print(results$bootstrap_estimates)

# Compute CIF results
cif_results <- CIF_res1(initial_params)
print(cif_results)

# Plot CIFs with confidence intervals
plot <- CIF_Figs(initial_params, timee)
print(plot)
```

## References
1. Jeong, J.-H., & Fine, J. (2006). Direct parametric inference for the cumulative incidence function. *Applied Statistics*, 55(2), 187-200.
2. Jeong, J.-H., & Fine, J. (2007). Parametric regression on cumulative incidence function. *Biostatistics*, 8(2), 184-196.
