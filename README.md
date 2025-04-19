# cmpp: Parametric Modeling for Competing Risks

![License](https://img.shields.io/badge/license-GPL%20(%3E%3D%202)-blue.svg)
![GitHub Actions Status](https://github.com/stats9/cmpp/actions/workflows/R-Check-Install.yaml/badge.svg)
![Platforms](https://img.shields.io/badge/tested%20on-Ubuntu%2C%20MacOS%2C%20Windows-blue)

## Overview
`cmpp` is an R package designed to facilitate parametric modeling and inference for cumulative incidence functions (CIFs) in competing risks scenarios. The package implements methods discussed in seminal works by Jeong and Fine (2006, 2007), with efficient computation powered by Rcpp for high-performance applications.

## Key Features
- **Parametric Modeling**:
  - Direct parametric modeling of cumulative incidence functions using Gompertz and Weibull distributions.
  - Support for Generalized Chance Models (GCM), Proportional Odds Models (POM), and Proportional Hazards Models (PHM).
- **Regression Models**:
  - Proportional hazards and proportional odds models.
  - Flexible generalized odds rate models.
  - Covariate-adjusted CIF estimation.
- **Likelihood-Based Inference**:
  - Full maximum likelihood estimation.
  - Delta method for confidence intervals and variance estimation.
- **Goodness-of-Fit Tests**:
  - Evaluation of proportional hazards or proportional odds assumptions.
- **Visualization**:
  - Plot cumulative incidence functions (CIFs) with confidence intervals.
  - Adjust CIFs based on covariate values.
- **High-Performance Computation**:
  - Fast computation using Rcpp and Eigen integration.

## Background
Competing risks occur when multiple types of events prevent the observation of a particular event of interest. Traditional survival analysis often misrepresents such data, as it fails to account for competing risks. This package provides parametric methods for cumulative incidence functions (CIFs), offering a more direct and interpretable analysis than traditional cause-specific hazard models.

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

4. **Visualization**:
   - Plot CIFs for competing risks with confidence intervals.
   - Adjust CIFs based on specific covariate values.

## Installation
To install the package:
```R
# Install from GitHub
devtools::install_github("your_github_username/cmpp")

# for windows 
install.packages(https://github.com/user-attachments/files/19817773/cmpp_0.0.1.zip, repos = NULL, type = "win-binary")

# for linux or mac
install.packages(https://github.com/user-attachments/files/19817772/cmpp_0.0.1.tar.gz, repos = NULL, type = "source")
```

## Usage

### Initialize the Cmpp Model
```R
library(cmpp)

# Example data
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
time <- rexp(100, rate = 1/10)

# Initialize the Cmpp model
Initialize(features, time, delta1, delta2, h = 1e-5)
```

### Estimate Parameters
```R
# Estimate parameters using the Generalized Chance Model (GCM)
initial_params <- rep(0.001, 2 * (ncol(features) + 3))
result <- estimate_parameters_GCM(initial_params)
print(result)

# Estimate parameters using the Proportional Odds Model (POM)
result_pom <- estimate_parameters_POM(initial_params)
print(result_pom)

# Estimate parameters using the Proportional Hazards Model (PHM)
result_phm <- estimate_parameters_PHM(initial_params)
print(result_phm)
```

### Compute Cumulative Incidence Functions (CIFs)
```R
# Compute CIFs for competing risks
cif_results <- CIF_res1(initial_params)
print(cif_results)

# Plot CIFs with confidence intervals
plot <- CIF_Figs(initial_params, time)
print(plot)
```

### Fine-Gray Model for Competing Risks
```R
# Fit a Fine-Gray model
result_fg <- FineGray_Model(
  featureNames = c("Feature1", "Feature2"),
  CovarNames = c("Covar1", "Covar2"),
  Failcode = 1,
  RiskNames = c("Event1", "Event2")
)
print(result_fg$Results)  # Summary of the Fine-Gray model
print(result_fg$Plot)     # Plot of the CIFs
```

### Adjusted CIFs Based on Covariates
```R
# Compute and plot adjusted CIFs
result_cif <- Cmpp_CIF(
  featureID = c(1, 2),
  featureValue = c(0.5, 1.2),
  RiskNames = c("Event1", "Event2"),
  TypeMethod = "GCM",
  predTime = seq(0, 10, by = 0.5)
)
print(result_cif$Plot$Plot_InputModel)  # Plot for the specified model
print(result_cif$CIF$CIFAdjusted)       # Adjusted CIF values
```

### Bootstrap Variance Estimation
```R
# Estimate variance of parameters using bootstrap
results <- bootstrap_variance(features, time, delta1, delta2, initial_params, n_bootstrap = 500)
print(results$variances)
print(results$bootstrap_estimates)
```

## Example Workflow
```R
library(cmpp)

# Load example data
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
time <- rexp(100, rate = 1/10)

# Initialize the Cmpp model
Initialize(features, time, delta1, delta2, h = 1e-5)

# Estimate parameters
initial_params <- rep(0.001, 2 * (ncol(features) + 3))
params <- estimate_parameters_GCM(initial_params)
print(params)

# Compute CIFs
cif_results <- CIF_res1(initial_params)
print(cif_results)

# Plot CIFs
plot <- CIF_Figs(initial_params, time)
print(plot)
```


## Example Of Data Inside cmpp package 

```R
data(dat)
names(dat)
feat <- dat[, -c(match(c('id', 'time', 'event', 'cause_burn', 'd1', 'd2', 'cause_hotObject3'), names(dat)))]  
timee <- dat[['time']]
d1 = dat[['d1']]
d2 = dat[['d2']]
feat2 = feat |> data.matrix()
Initialize(feat2, timee, d1, d2, 1e-10)

FineGray_Model()
(Res <- Cmpp_CIF())
Res$Plot$PlotNull_AllModels
Res$Plot$Plot_InputModel
Res$Plot$PlotAdjusted_AllModels

Res2 <- Cmpp_CIF(predTime = c(1, 2, 3), TypeMethod = "POM")
Res2$CIF
Res2$Plot$PlotAdjusted_AllModels
GetData()
Res2$Plot$PlotNull_AllModels
Res2$Plot$Plot_InputModel
```

## References
1. Jeong, J.-H., & Fine, J. (2006). Direct parametric inference for the cumulative incidence function. *Applied Statistics*, 55(2), 187-200.
2. Jeong, J.-H., & Fine, J. (2007). Parametric regression on cumulative incidence function. *Biostatistics*, 8(2), 184-196.