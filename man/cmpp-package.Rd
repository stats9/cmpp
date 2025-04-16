\name{cmpp-package}
\alias{cmpp-package}
\alias{cmpp}
\docType{package}
\title{
  Parametric Methods for Cumulative Incidence Functions in Competing Risks
}
\description{
  The `cmpp` package provides parametric modeling methods for analyzing cumulative incidence functions (CIFs) in the context of competing risks. It includes Gompertz-based models, regression techniques, and parametric approaches such as the Generalized Chance Model (GCM), Proportional Odds Model (POM), and Proportional Hazards Model (PHM). The package enables users to estimate and compare CIFs using maximum likelihood estimation, perform regression analysis, and visualize CIFs with confidence intervals. It also supports covariate adjustment and bootstrap variance estimation.
}
\details{
  The `cmpp` package offers functions for modeling cumulative incidence functions (CIFs) using the Gompertz distribution and generalized regression models. Key features of the package include:
  - Direct parametric modeling for cumulative incidence functions.
  - Maximum likelihood estimation of parameters.
  - Regression analysis with covariates, including treatment effects.
  - Functions for evaluating goodness-of-fit and performing model diagnostics.
  - Visualization of CIFs with confidence intervals.
  - Covariate-adjusted CIF estimation.
  - Bootstrap variance estimation for model parameters.

  Commonly used functions include:
  - `Initialize()`: Initializes the data for the Cmpp model.
  - `LogLike1()`: Computes the negative log-likelihood for the model.
  - `compute_grad()`: Computes the gradient of the log-likelihood.
  - `compute_hessian()`: Computes the Hessian matrix of the log-likelihood.
  - `estimate_parameters_GCM()`: Estimates parameters using the Generalized Chance Model (GCM).
  - `estimate_parameters_POM()`: Estimates parameters using the Proportional Odds Model (POM).
  - `estimate_parameters_PHM()`: Estimates parameters using the Proportional Hazards Model (PHM).
  - `CIF_res1()`: Computes CIF results for competing risks.
  - `CIF_Figs()`: Plots CIFs with confidence intervals.
  - `Cmpp_CIF()`: Computes and plots CIFs for competing risks using GCM, POM, and PHM.
  - `FineGray_Model()`: Fits a Fine-Gray model for competing risks data.
  - `bootstrap_variance()`: Estimates variance of parameters using the bootstrap method.
  - `GetData()`: Retrieves initialized data from the Cmpp model.
  - `Cleanup()`: Cleans up memory by deleting the Cmpp instance.

  This package is particularly useful in clinical trial analysis, especially in the presence of competing risks (e.g., death, recurrence).
}
\author{
  Habib Ezzat Abadi Pour (Stats9) \email{habibezati@outlook.com}
}
\references{
  Jeong, J.-H., & Fine, J. (2006). Direct parametric inference for the cumulative incidence function. *Applied Statistics*, 55(2), 187-200.  
  Jeong, J.-H., & Fine, J. (2007). Parametric regression on cumulative incidence function. *Biostatistics*, 8(2), 184-196.
}
\keyword{ package, survival analysis, competing risks, cumulative incidence, regression, Gompertz distribution, parametric modeling }
\seealso{
  \code{\link{Initialize}}, \code{\link{LogLike1}}, \code{\link{compute_grad}}, \code{\link{compute_hessian}}, \code{\link{estimate_parameters_GCM}}, \code{\link{estimate_parameters_POM}}, \code{\link{estimate_parameters_PHM}}, \code{\link{CIF_res1}}, \code{\link{CIF_Figs}}, \code{\link{Cmpp_CIF}}, \code{\link{FineGray_Model}}, \code{\link{bootstrap_variance}}, \code{\link{GetData}}, \code{\link{Cleanup}}
}
\examples{
  \dontrun{
    ## Example: Initialize the Cmpp model and compute CIFs
    library(cmpp)

    # Example data
    features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
    delta1 <- sample(c(0, 1), 100, replace = TRUE)
    delta2 <- 1 - delta1
    time <- rexp(100, rate = 1/10)

    # Initialize the Cmpp model
    Initialize(features, time, delta1, delta2, h = 1e-5)

    # Estimate parameters using the Generalized Chance Model (GCM)
    initial_params <- rep(0.001, 2 * (ncol(features) + 3))
    result <- estimate_parameters_GCM(initial_params)
    print(result)

    # Compute CIFs for competing risks
    cif_results <- CIF_res1(initial_params)
    print(cif_results)

    # Plot CIFs with confidence intervals
    plot <- CIF_Figs(initial_params, time)
    print(plot)

    # Compute and plot adjusted CIFs
    result_cif <- Cmpp_CIF(
      featureID = c(1, 2),
      featureValue = c(0.5, 1.2),
      RiskNames = c("Event1", "Event2"),
      TypeMethod = "GCM",
      predTime = seq(0, 10, by = 0.5)
    )
    print(result_cif$Plot$Plot_InputModel)  # Plot for the specified model
    print(result_cif$CIF$CIFAdjusted)      # Adjusted CIF values

    # Fit a Fine-Gray model for competing risks
    result_fg <- FineGray_Model(
      featureNames = c("Feature1", "Feature2"),
      CovarNames = c("Covar1", "Covar2"),
      Failcode = 1,
      RiskNames = c("Event1", "Event2")
    )
    print(result_fg$Results)  # Summary of the Fine-Gray model
    print(result_fg$Plot)     # Plot of the CIFs

    # Clean up memory
    Cleanup()
  }
}