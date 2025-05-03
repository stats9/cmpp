# cmpp: Parametric Modeling for Competing Risks

![License](https://img.shields.io/badge/license-GPL%20(%3E%3D%202)-blue.svg)
![GitHub Actions Status](https://github.com/stats9/cmpp/actions/workflows/R-Check-Install.yaml/badge.svg)
![Platforms](https://img.shields.io/badge/tested%20on-Ubuntu%2C%20MacOS%2C%20Windows-blue)
[![CRAN Status](https://www.r-pkg.org/badges/version/cmpp)](https://CRAN.R-project.org/package=cmpp)


## Overview
`cmpp` is an R package designed to facilitate parametric modeling and inference for cumulative incidence functions (CIFs) in competing risks scenarios. The package implements methods discussed in seminal works by Jeong and Fine (2006, 2007), with efficient computation powered by Rcpp for high-performance applications.

## Key Features
- **Parametric Modeling**:
  - Direct parametric modeling of cumulative incidence functions using Gompertz and Weibull distributions.
  - Support for Generalized odds rates (GOR), Proportional Odds Models (POM), and Proportional Hazards Models (PHM).
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

# install from cran 
install.packages("cmpp")

# Install from GitHub
devtools::install_github("stats9/cmpp")

# for windows 
## after download cmpp zip file 
install.packages('cmpp_0.0.1.zip', repos = NULL, type = "win-binary")

# for linux or mac
## after download cmpp tar.gz file
install.packages('cmpp_0.0.1.tar.gz', repos = NULL, type = "source")
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
 library(cmpp)
 data("fertility_data")
 Nam <- names(fertility_data)
 fertility_data$Education
 datt <- make_Dummy(fertility_data, features = c("Education"))
 datt <- datt$New_Data 
 datt['Primary_Secondary'] <- datt$`Education:2`
 datt['Higher_Education'] <- datt$`Education:3`
 datt$`Education:2` <- datt$`Education:3` <- NULL
 datt2 <- make_Dummy(datt, features = 'Event')$New_Data
 d1 <- datt2$`Event:2`
 d2 <- datt2$`Event:3`
 feat <- datt2[c('age', 'Primary_Secondary', 'Higher_Education')] |> 
    data.matrix()
 timee <- datt2[['time']]
 Initialize(feat, timee, d1, d2, 1e-10)

# Estimate parameters using the Generalized odds rate (GOR)
initial_params <- rep(0.001, 2 * (ncol(features) + 3))
initial_params2 <- rep(0.001, 2 * (ncol(features) + 2))
result <- estimate_parameters_GOR(initial_params)
print(result)

# Estimate parameters using the Proportional Odds Model (POM)
result_pom <- estimate_parameters_POM(initial_params2)
print(result_pom)

# Estimate parameters using the Proportional Hazards Model (PHM)
result_phm <- estimate_parameters_PHM(initial_params2)
print(result_phm)

```

### Compute Cumulative Incidence Functions (CIFs)
```R
 library(cmpp)
 data("fertility_data")
 Nam <- names(fertility_data)
 fertility_data$Education
 datt <- make_Dummy(fertility_data, features = c("Education"))
 datt <- datt$New_Data 
 datt['Primary_Secondary'] <- datt$`Education:2`
 datt['Higher_Education'] <- datt$`Education:3`
 datt$`Education:2` <- datt$`Education:3` <- NULL
 datt2 <- make_Dummy(datt, features = 'Event')$New_Data
 d1 <- datt2$`Event:2`
 d2 <- datt2$`Event:3`
 feat <- datt2[c('age', 'Primary_Secondary', 'Higher_Education')] |> 
    data.matrix()
 timee <- datt2[['time']]
 Initialize(feat, timee, d1, d2, 1e-10)

# Compute CIFs for competing risks
cif_results <- CIF_res1(rep(0.001, 4))
print(cif_results)

# Plot CIFs with confidence intervals
Res <- CIF_Figs(rep(0.001, 4), timee)
print(Res)
```

### Fine-Gray Model for Competing Risks
```R
 library(cmpp)
 data("fertility_data")
 Nam <- names(fertility_data)
 fertility_data$Education
 datt <- make_Dummy(fertility_data, features = c("Education"))
 datt <- datt$New_Data 
 datt['Primary_Secondary'] <- datt$`Education:2`
 datt['Higher_Education'] <- datt$`Education:3`
 datt$`Education:2` <- datt$`Education:3` <- NULL
 datt2 <- make_Dummy(datt, features = 'Event')$New_Data
 d1 <- datt2$`Event:2`
 d2 <- datt2$`Event:3`
 feat <- datt2[c('age', 'Primary_Secondary', 'Higher_Education')] |> 
    data.matrix()
 timee <- datt2[['time']]
 Initialize(feat, timee, d1, d2, 1e-10)

# Fit a Fine-Gray model
result_fg <- FineGray_Model(
  CovarNames = c("Covar1", "Covar2", "Covar3"),
  Failcode = 1,
  RiskNames = c("Event1", "Event2")
)
print(result_fg$Results)  # Summary of the Fine-Gray model
print(result_fg$Plot)     # Plot of the CIFs
```

### Adjusted CIFs Based on Covariates
```R
 library(cmpp)
 data("fertility_data")
 Nam <- names(fertility_data)
 fertility_data$Education
 datt <- make_Dummy(fertility_data, features = c("Education"))
 datt <- datt$New_Data 
 datt['Primary_Secondary'] <- datt$`Education:2`
 datt['Higher_Education'] <- datt$`Education:3`
 datt$`Education:2` <- datt$`Education:3` <- NULL
 datt2 <- make_Dummy(datt, features = 'Event')$New_Data
 d1 <- datt2$`Event:2`
 d2 <- datt2$`Event:3`
 feat <- datt2[c('age', 'Primary_Secondary', 'Higher_Education')] |> 
    data.matrix()
 timee <- datt2[['time']]
 Initialize(feat, timee, d1, d2, 1e-10)
result_cif <- Cmpp_CIF(
  featureID = c(1, 2),
  featureValue = c(0.5, 1.2),
  RiskNames = c("Event1", "Event2"),
  TypeMethod = "GOR",
  predTime = seq(0, 10, by = 0.5)
)
print(result_cif$Plot$Plot_InputModel)  # Plot for the specified model
print(result_cif$CIF$CIFAdjusted)       # Adjusted CIF values
```

### Bootstrap Variance Estimation
```R
 library(cmpp)
 data("fertility_data")
 Nam <- names(fertility_data)
 fertility_data$Education
 datt <- make_Dummy(fertility_data, features = c("Education"))
 datt <- datt$New_Data 
 datt['Primary_Secondary'] <- datt$`Education:2`
 datt['Higher_Education'] <- datt$`Education:3`
 datt$`Education:2` <- datt$`Education:3` <- NULL
 datt2 <- make_Dummy(datt, features = 'Event')$New_Data
 d1 <- datt2$`Event:2`
 d2 <- datt2$`Event:3`
 feat <- datt2[c('age', 'Primary_Secondary', 'Higher_Education')] |> 
    data.matrix()
 timee <- datt2[['time']]
 Initialize(feat, timee, d1, d2, 1e-10)
datt <- GetData()

# Estimate variance of parameters using bootstrap
results <- bootstrap_variance(datt$features, datt$timee, 
    datt$delta1, datt$delta2, rep(0.001, 4), n_bootstrap = 500)
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
params <- estimate_parameters_GOR(initial_params)
print(params)

# Compute CIFs
cif_results <- CIF_res1(rep(0.001, 4))
print(cif_results)

# Plot CIFs
Res <- CIF_Figs(rep(0.01, 4), time)
print(Res)
```


## Example Of Data Inside cmpp package 

```R
 library(cmpp)
 data("fertility_data")
 Nam <- names(fertility_data)
 fertility_data$Education
 datt <- make_Dummy(fertility_data, features = c("Education"))
 datt <- datt$New_Data 
 datt['Primary_Secondary'] <- datt$`Education:2`
 datt['Higher_Education'] <- datt$`Education:3`
 datt$`Education:2` <- datt$`Education:3` <- NULL
 datt2 <- make_Dummy(datt, features = 'Event')$New_Data
 d1 <- datt2$`Event:2`
 d2 <- datt2$`Event:3`
 feat <- datt2[c('age', 'Primary_Secondary', 'Higher_Education')] |> 
    data.matrix()
 timee <- datt2[['time']]
 Initialize(feat, timee, d1, d2, 1e-10)
 data <- GetData()
 print(data$features)  # Feature matrix
 print(data$timee)      # Failure times
 print(data$delta1)    # Indicator for the first competing event
 print(data$delta2)    # Indicator for the second competing event
 result <- Cmpp_CIF(
   featureID = c(1, 2),
   featureValue = c(0.5, 1.2),
   RiskNames = c("Event1", "Event2"),
   TypeMethod = "GOR",
   predTime = seq(0, 10, by = 0.5)
 )
 print(result$Plot$Plot_InputModel)  # Plot for the specified model
 print(result$Plot$PlotAdjusted_AllModels)  # Adjusted CIFs for all models
 print(result$CIF$CIFAdjusted)  # Adjusted CIF values
```

## References
1. Jeong, J.-H., & Fine, J. (2006). Direct parametric inference for the cumulative incidence function. *Applied Statistics*, 55(2), 187-200.
2. Jeong, J.-H., & Fine, J. (2007). Parametric regression on cumulative incidence function. *Biostatistics*, 8(2), 184-196.