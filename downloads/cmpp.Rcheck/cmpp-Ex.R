pkgname <- "cmpp"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "cmpp-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('cmpp')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CIF_Figs")
### * CIF_Figs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CIF_Figs
### Title: Plot Cumulative Incidence Functions (CIF) with Confidence
###   Intervals
### Aliases: CIF_Figs

### ** Examples

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
initial_params <- c(0.001, 0.001, 0.001, 0.001)
result <- CIF_res1(initial_params)
print(result)
initial_params <- c(0.01, 0.01, 0.01, 0.01)
TimeFailure <- seq(0, 10, by = 0.1)
plot <- CIF_Figs(initial_params, TimeFailure)
print(plot)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CIF_Figs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CIF_res1")
### * CIF_res1

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CIF_res1
### Title: Compute Cumulative Incidence Function (CIF) Results
### Aliases: CIF_res1

### ** Examples

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
initial_params <- c(0.001, 0.001, 0.001, 0.001)
result <- CIF_res1(initial_params)
print(result)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CIF_res1", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Cleanup")
### * Cleanup

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Cleanup
### Title: Clean up memory by deleting the pointer to the Cmpp instance
### Aliases: Cleanup

### ** Examples

# Assuming you have previously initialized the Cmpp object with `Initialize()`
Cleanup()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Cleanup", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Cmpp_CIF")
### * Cmpp_CIF

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Cmpp_CIF
### Title: Compute and Plot Cumulative Incidence Functions (CIF) for
###   Competing Risks
### Aliases: Cmpp_CIF

### ** Examples

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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Cmpp_CIF", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_cdf_rcpp")
### * F_cdf_rcpp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_cdf_rcpp
### Title: Compute the CDF of the Parametric Generalized odds rate (GOR)
### Aliases: F_cdf_rcpp

### ** Examples

library(cmpp)
set.seed(321)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/2)
Initialize(features, x, delta1, delta2, h = 1e-3)
params <- rep(0.001, (ncol(features) + 3))
y <- 0.07
z <- features[1, ]
(cdf_value <- F_cdf_rcpp(params, z, y))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_cdf_rcpp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_cdf_rcpp2")
### * F_cdf_rcpp2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_cdf_rcpp2
### Title: Compute the CDF of the Parametric Proportional Odds Model (POM)
### Aliases: F_cdf_rcpp2

### ** Examples

library(cmpp)
set.seed(1984)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/2)
Initialize(features, x, delta1, delta2, h = 1e-5)
params <- rep(0.001, (ncol(features) + 2))
x <- 2
cdf_value <- F_cdf_rcpp2(params, features[1, ], x)
print(cdf_value)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_cdf_rcpp2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_cdf_rcpp3")
### * F_cdf_rcpp3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_cdf_rcpp3
### Title: Compute the CDF of the Parametric Proportional Hazards Model
###   (PHM)
### Aliases: F_cdf_rcpp3

### ** Examples

library(cmpp)
set.seed(1984)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/10)
Initialize(features, x, delta1, delta2, h = 1e-5)
params <- rep(0.001, (ncol(features) + 2))
x <- 5
cdf_value <- F_cdf_rcpp3(params, features[1, ], x)
print(cdf_value)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_cdf_rcpp3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("FineGray_Model")
### * FineGray_Model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: FineGray_Model
### Title: Fine-Gray Model for Competing Risks Data
### Aliases: FineGray_Model

### ** Examples

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
result <- FineGray_Model(
  CovarNames = c("Covar1", "Covar2", "Covar3"),
  Failcode = 1,
  RiskNames = c("Event1", "Event2")
)
print(result$Results)  # Summary of the Fine-Gray model
#print(result$Plot)     # Plot of the CIFs
print(result$CIF_Results)  # CIF data




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("FineGray_Model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("GetData")
### * GetData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: GetData
### Title: Retrieve Initialized Data from the Cmpp Model
### Aliases: GetData

### ** Examples

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("GetData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("GetDim")
### * GetDim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: GetDim
### Title: Get Dimensions of the Cmpp Object
### Aliases: GetDim

### ** Examples

# Initialize Cmpp object
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
# Get dimensions
dims <- GetDim()
dims$Nsamp    # Number of samples
dims$Nfeature # Number of features



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("GetDim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Initialize")
### * Initialize

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Initialize
### Title: Initialize Data for the Cmpp Model
### Aliases: Initialize

### ** Examples

library(cmpp)
features <- matrix(rnorm(100), ncol = 5)
x <- rnorm(20)
delta1 <- sample(0:1, 20, replace = TRUE)
delta2 <- 1 - delta1
Initialize(features, x, delta1, delta2, h = 1e-5)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Initialize", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LogLike1")
### * LogLike1

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LogLike1
### Title: Compute the Log-Likelihood for the Model
### Aliases: LogLike1

### ** Examples

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
param <- c(0.01, 0.01, 0.01, 0.01)
LogLike1(param)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LogLike1", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bootstrap_variance")
### * bootstrap_variance

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bootstrap_variance
### Title: Estimate Variance of Parameters Using Bootstrap Method
### Aliases: bootstrap_variance

### ** Examples

library(cmpp)
features <- matrix(rnorm(100), ncol = 5)
x <- rnorm(20)
delta1 <- sample(0:1, 20, replace = TRUE)
delta2 <- 1 - delta1
initial_params <- c(0.01, 0.01, 0.01, 0.01)
n_bootstrap <- 100
results <- bootstrap_variance(features, x, delta1, delta2, 
initial_params, n_bootstrap, optimMethod = "BFGS")
print(results$variances)
print(results$bootstrap_estimates)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bootstrap_variance", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cdf_gomp")
### * cdf_gomp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cdf_gomp
### Title: Compute the CDF of the Gompertz Distribution
### Aliases: cdf_gomp

### ** Examples

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
x <- c(1, 2, 3)
alpha <- 0.5
beta <- 0.1
lapply(x, cdf_gomp, alpha = alpha, beta = beta) |> unlist()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cdf_gomp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cmpp")
### * cmpp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cmpp
### Title: Direct Parametric Inference for the Cumulative Incidence
###   Function in Competing Risks
### Aliases: cmpp-package cmpp
### Keywords: cumulative incidence parametric regression risks survival

### ** Examples

## Example: Initialize the Cmpp model and compute CIFs
library(cmpp)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
timee <- rexp(100, rate = 1/10)
Initialize(features, timee, delta1, delta2, h = 1e-5)
# Initialize the Cmpp model
# Estimate parameters using the Generalized odds rate (GOR)
initial_params <- rep(0.001, 2 * (ncol(features) + 3))
result <- estimate_parameters_GOR(initial_params)
print(result)
# Compute CIFs for competing risks (without covariate effect | Not Regression model)
cif_results <- CIF_res1()
print(cif_results)
# Plot CIFs with confidence intervals
plot <- CIF_Figs(rep(0.01, 4), timee)
print(plot)
# Compute and plot adjusted CIFs
result_cif <- Cmpp_CIF(
featureID = c(1, 2),
featureValue = c(0.5, 1.2),
RiskNames = c("Event1", "Event2"),
TypeMethod = "GOR",
predTime = seq(0, 10, by = 0.5)
)
print(result_cif$Plot$Plot_InputModel) # Plot for the specified model
print(result_cif$CIF$CIFAdjusted) # Adjusted CIF values
# Fit a Fine-Gray model for competing risks
result_fg <- FineGray_Model(
CovarNames = c("Covar1", "Covar2", 'Covar3'),
Failcode = 1,
RiskNames = c("Event1", "Event2")
)
print(result_fg$Results)  # Summary of the Fine-Gray model
print(result_fg$Plot) # Plot of the CIFs

# Clean up memory
Cleanup()




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cmpp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_grad")
### * compute_grad

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_grad
### Title: Compute the Numerical Gradient of the Log-Likelihood
### Aliases: compute_grad

### ** Examples

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
param <- c(0.5, 0.1, 0.6, 0.2)
compute_grad(param)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_grad", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_hessian")
### * compute_hessian

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_hessian
### Title: Compute the Hessian Matrix of the Log-Likelihood
### Aliases: compute_hessian

### ** Examples

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
# Estimate model parameters using default initial values and the BFGS method
result <- estimate_parameters()
print(result)
param <- c(0.5, 0.1, 0.6, 0.2)
hessian <- compute_hessian(param)
print(hessian)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_hessian", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_log_f_gradient_rcpp")
### * compute_log_f_gradient_rcpp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_log_f_gradient_rcpp
### Title: Compute the Gradient of the Log-Likelihood Function Generalized
###   odds rate (GOR)
### Aliases: compute_log_f_gradient_rcpp

### ** Examples

library(cmpp)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
set.seed(1984)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/3)
Initialize(features, x, delta1, delta2, h = 1e-5)
params <- rep(0.001, 2 * (ncol(features) + 3))
gradient <- compute_log_f_gradient_rcpp(params)
print(gradient)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_log_f_gradient_rcpp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_log_f_gradient_rcpp2")
### * compute_log_f_gradient_rcpp2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_log_f_gradient_rcpp2
### Title: Compute the Gradient of the Log-Likelihood Function Proportional
###   Odds Model (POM)
### Aliases: compute_log_f_gradient_rcpp2

### ** Examples

library(cmpp)
set.seed(1984)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/5)
Initialize(features, x, delta1, delta2, h = 1e-5)
params <- rep(0.001, 2 * (ncol(features) + 2))
gradient <- compute_log_f_gradient_rcpp2(params)
print(gradient)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_log_f_gradient_rcpp2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_log_f_gradient_rcpp3")
### * compute_log_f_gradient_rcpp3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_log_f_gradient_rcpp3
### Title: Compute the Gradient of the Log-Likelihood Function Proportional
###   Hazards Model (PHM)
### Aliases: compute_log_f_gradient_rcpp3

### ** Examples

library(cmpp)
set.seed(1984)  
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/10)
Initialize(features, x, delta1, delta2, h = 1e-5)
params <- rep(0.001, 2 * (ncol(features) + 2))
gradient <- compute_log_f_gradient_rcpp3(params)
print(gradient)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_log_f_gradient_rcpp3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_log_f_hessian_rcpp")
### * compute_log_f_hessian_rcpp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_log_f_hessian_rcpp
### Title: Compute the Hessian Matrix of the Log-Likelihood Function
###   Generalized odds rate (GOR)
### Aliases: compute_log_f_hessian_rcpp

### ** Examples

library(cmpp)
set.seed(1984)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/7)
Initialize(features, x, delta1, delta2, h = 1e-4)
params <- rep(0.001, 2 * (ncol(features) + 3))
hessian <- compute_log_f_hessian_rcpp(params)
print(hessian)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_log_f_hessian_rcpp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("estimate_parameters")
### * estimate_parameters

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: estimate_parameters
### Title: Estimate Model Parameters Using Optimization
### Aliases: estimate_parameters

### ** Examples

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
# Estimate model parameters using default initial values and the BFGS method
result <- estimate_parameters()
print(result)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("estimate_parameters", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("estimate_parameters_GOR")
### * estimate_parameters_GOR

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: estimate_parameters_GOR
### Title: Estimate Parameters for the Generalized odds rate (GOR)
### Aliases: estimate_parameters_GOR

### ** Examples

library(cmpp)
# Example data
set.seed(371)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/4)
# Initialize the Cmpp model
Initialize(features, x, delta1, delta2, h = 1e-5)
# Define initial parameter values
initial_params <- rep(0.001, 2 * (ncol(features) + 3))
# Estimate parameters using the GOR
result <- estimate_parameters_GOR(initial_params)
print(result)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("estimate_parameters_GOR", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("estimate_parameters_PHM")
### * estimate_parameters_PHM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: estimate_parameters_PHM
### Title: Estimate Parameters for the Proportional Hazards Model (PHM)
### Aliases: estimate_parameters_PHM

### ** Examples

library(cmpp)
set.seed(1984)
# Example data
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/10)

# Initialize the Cmpp model
Initialize(features, x, delta1, delta2, h = 1e-5)

# Define initial parameter values
initial_params <- rep(0.001, 2 * (ncol(features) + 2))

# Estimate parameters using the PHM
result <- estimate_parameters_PHM(initial_params)
print(result)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("estimate_parameters_PHM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("estimate_parameters_POM")
### * estimate_parameters_POM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: estimate_parameters_POM
### Title: Estimate Parameters for the Proportional Odds Model (POM)
### Aliases: estimate_parameters_POM

### ** Examples

library(cmpp)
set.seed(1984)
# Example data
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/10)

# Initialize the Cmpp model
Initialize(features, x, delta1, delta2, h = 1e-5)

# Define initial parameter values
initial_params <- rep(0.001, 2 * (ncol(features) + 2))

# Estimate parameters using the POM
result <- estimate_parameters_POM(initial_params)
print(result)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("estimate_parameters_POM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("f_pdf_rcpp")
### * f_pdf_rcpp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: f_pdf_rcpp
### Title: Compute the PDF of the Parametric Generalized odds rate (GOR)
### Aliases: f_pdf_rcpp

### ** Examples

library(cmpp)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/10)
Initialize(features, x, delta1, delta2, h = 1e-5)
params <- rep(0.001, (ncol(features) + 3))
pdf_value <- f_pdf_rcpp(params, features[1, ], x[3])
print(pdf_value)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("f_pdf_rcpp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("f_pdf_rcpp2")
### * f_pdf_rcpp2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: f_pdf_rcpp2
### Title: Compute the PDF of the Parametric Proportional Odds Model (POM)
### Aliases: f_pdf_rcpp2

### ** Examples

library(cmpp)
set.seed(1984)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/9)
Initialize(features, x, delta1, delta2, h = 1e-4)
params <- rep(0.001, (ncol(features) + 2))
pdf_value <- f_pdf_rcpp2(params, features[1, ], x[3])
print(pdf_value)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("f_pdf_rcpp2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("f_pdf_rcpp3")
### * f_pdf_rcpp3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: f_pdf_rcpp3
### Title: Compute the PDF of the Parametric Proportional Hazards Model
###   (PHM)
### Aliases: f_pdf_rcpp3

### ** Examples

library(cmpp)
set.seed(21)
features <- matrix(rnorm(300, -1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1)
Initialize(features, x, delta1, delta2, h = 1e-5)
params <- rep(0.0001, (ncol(features) + 2))
pdf_value <- f_pdf_rcpp3(params, features[4, ], x[4])
print(pdf_value)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("f_pdf_rcpp3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fertility_data")
### * fertility_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fertility_data
### Title: Fertility History of Rural Women in Shiraz
### Aliases: fertility_data
### Keywords: datasets

### ** Examples

data(fertility_data)
head(fertility_data)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fertility_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("log_f_rcpp")
### * log_f_rcpp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: log_f_rcpp
### Title: Compute the Log-Likelihood Function Generalized odds rate (GOR)
### Aliases: log_f_rcpp

### ** Examples

library(cmpp)
set.seed(1984)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/4)
Initialize(features, x, delta1, delta2, h = 1e-5)
params <- rep(0.001, 2 * (ncol(features) + 3))
log_likelihood <- log_f_rcpp(params)
print(log_likelihood)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("log_f_rcpp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("log_f_rcpp2")
### * log_f_rcpp2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: log_f_rcpp2
### Title: Compute the Log-Likelihood Function Proportional Odds Model
###   (POM)
### Aliases: log_f_rcpp2

### ** Examples

library(cmpp)
set.seed(1984)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/8)
Initialize(features, x, delta1, delta2, h = 1e-5)
params <- rep(0.001, 2 * (ncol(features) + 2))
log_likelihood <- log_f_rcpp2(params)
print(log_likelihood)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("log_f_rcpp2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("log_f_rcpp3")
### * log_f_rcpp3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: log_f_rcpp3
### Title: Compute the Log-Likelihood Function Proportional Hazards Model
###   (PHM)
### Aliases: log_f_rcpp3

### ** Examples

library(cmpp)
set.seed(1984)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/10)
Initialize(features, x, delta1, delta2, h = 1e-5)
params <- rep(0.001, 2 * (ncol(features) + 2))
log_likelihood <- log_f_rcpp3(params)
print(log_likelihood)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("log_f_rcpp3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("makeMat")
### * makeMat

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: makeMat
### Title: Create a matrix of given size filled with a constant value
### Aliases: makeMat

### ** Examples

library(cmpp)
# Create a 3x3 matrix filled with 5
mat <- makeMat(3, 3, 5)
print(mat)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("makeMat", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("make_Dummy")
### * make_Dummy

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: make_Dummy
### Title: Create Dummy Variables
### Aliases: make_Dummy

### ** Examples

dat <- data.frame(sex = c('M', 'F', 'M'), cause_burn = c('A', 'B', 'A'))
result <- make_Dummy(Data = dat, features = c('sex', 'cause_burn'), reff = "first")
print(result$New_Data)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("make_Dummy", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pdf_gomp")
### * pdf_gomp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pdf_gomp
### Title: Compute the PDF of the Gompertz Distribution
### Aliases: pdf_gomp

### ** Examples

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
x <- c(1, 2, 3)
alpha <- 0.5
beta <- 0.1
lapply(x, pdf_gomp, alpha = alpha, beta = beta) |> unlist()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pdf_gomp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
