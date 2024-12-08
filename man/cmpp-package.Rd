\name{cmpp-package}
\alias{cmpp-package}
\alias{cmpp}
\docType{package}
\title{
  Parametric Methods for Cumulative Incidence Functions in Competing Risks
}
\description{
  The `cmpp` package provides parametric modeling methods for analyzing cumulative incidence functions in the context of competing risks. It includes Gompertz-based models and regression techniques for survival analysis, enabling users to estimate and compare cumulative incidence functions using maximum likelihood estimation. The package also supports the calculation of variance, goodness-of-fit tests, and covariate adjustment.
}
\details{
  The `cmpp` package offers functions for modeling cumulative incidence functions (CIF) using the Gompertz distribution and generalized regression models. Key features of the package include:
  - Direct parametric modeling for cumulative incidence functions.
  - Maximum likelihood estimation of parameters.
  - Regression analysis with covariates, including treatment effects.
  - Functions for evaluating goodness-of-fit and performing model diagnostics.
  
  Commonly used functions include:
  - `gompertz_cif()`: Fits a Gompertz model for CIF.
  - `regression_cif()`: Performs regression analysis for CIF with covariates.
  - `predict_cif()`: Predicts cumulative incidence based on fitted models.
  - `check_fit()`: Evaluates the goodness-of-fit of the model.

  This package is useful in clinical trial analysis, especially in the presence of competing risks (e.g., death, recurrence).
}
\author{
  Habib Ezzat Abadi Pour (Stats9) \email{habibezati@outlook.com}
}
\references{
  Jeong, J.-H., & Fine, J. (2006). Direct parametric inference for the cumulative incidence function. *Applied Statistics*, 55(2), 187-200.  
  Jeong, J.-H., & Fine, J. (2007). Parametric regression on cumulative incidence function. *Biostatistics*, 8(2), 184-196.
}
\keyword{ package, survival analysis, competing risks, cumulative incidence, regression, Gompertz distribution }
\seealso{
  \code{\link{gompertz_cif}}, \code{\link{regression_cif}}, \code{\link{predict_cif}}, \code{\link{check_fit}}
}
\examples{
  \dontrun{
    ## Example: Fit a Gompertz model to a dataset and predict CIF
    data(dat)  # Load your dataset
    fit <- gompertz_cif(time = dat$time, event = dat$event)
    summary(fit)
    
    # Predict cumulative incidence for a specific time point
    predict_cif(fit, time = 5)
  }
}
