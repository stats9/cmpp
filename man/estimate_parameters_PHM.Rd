% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCode.r
\name{estimate_parameters_PHM}
\alias{estimate_parameters_PHM}
\title{Estimate Parameters for the Proportional Hazards Model (PHM)}
\usage{
estimate_parameters_PHM(initial_params, FeaturesNames = NULL, Method = "BFGS")
}
\arguments{
\item{initial_params}{A numeric vector of initial parameter values to start the optimization.
Default is \code{rep(0.001, 2 * (2 + ncol(features)))}, where \code{features} is the matrix of predictor variables.}

\item{FeaturesNames}{A character vector specifying the names of the features (covariates).
If \code{NULL}, default names (\code{beta1}, \code{beta2}, etc.) will be generated.}

\item{Method}{A character string specifying the optimization method to use.
The default is \code{"BFGS"}. See \code{?optim} for a list of available methods.}
}
\value{
A data frame containing:
\item{Parameter}{The parameter names, including \code{tau1}, \code{rho1}, \code{tau2}, \code{rho2}, and covariate coefficients (\code{beta1}, \code{beta2}, etc.).}
\item{Estimate}{The estimated parameter values.}
\item{S.E}{The standard errors of the estimated parameters.}
\item{PValue}{The p-values for the estimated parameters.}
}
\description{
This function estimates the parameters of the Proportional Hazards Model (PHM) using maximum likelihood estimation.
It computes the Hessian matrix, calculates standard errors, and derives p-values for the estimated parameters.
The function ensures that the diagonal elements of the covariance matrix are positive for valid variance estimates.
}
\details{
This function performs the following steps:
\itemize{
\item Estimates the model parameters using the \code{optim} function with the BFGS method.
\item Computes the gradient of the log-likelihood using the \code{compute_log_f_gradient_rcpp3} function.
\item Computes the Hessian matrix numerically using the \code{hessian} function from the \code{numDeriv} package.
\item Ensures that the diagonal elements of the covariance matrix are positive to avoid invalid variance estimates.
\item Calculates standard errors and p-values for the estimated parameters.
}

The Proportional Hazards Model (PHM) is a parametric model for cumulative incidence functions in competing risks analysis.
It uses Gompertz distributions to model the failure times for competing events.
}
\examples{
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

}
\seealso{
\link[stats:optim]{stats::optim},
\link{compute_log_f_gradient_rcpp3},
\link{log_f_rcpp3}.
}
