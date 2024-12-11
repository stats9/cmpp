#' Initialize Data for the Cmpp Model
#'
#' This function initializes the data used in the Cmpp model by passing the feature matrix, failure times, and the censoring indicators for competing risks and event indicators.
#' It sets up the necessary data structures for later computation.
#'
#' @param features A matrix of features (predictor variables) for each observation. Each row corresponds to an observation, and each column corresponds to a feature.
#' @param x A vector of failure times corresponding to the observations. The failure times can represent the time to an event or censoring.
#' @param delta1 A vector of Competing Risk indicator values. A binary vector where 1 indicates that the first event was observed and 0 means the event was censored or the second event occurred.
#' @param delta2 A vector of Event indicator values. A binary vector where 1 indicates that the second event was observed, and 0 means the event was censored or the first event occurred.
#' @param h A small step size used for numerical gradient computation. It is essential for calculating the gradient in optimization routines.
#'
#' @return This function does not return any value. It initializes the Cmpp model with the provided data, which is then used in subsequent computations.
#' 
#' @examples
#' # Example data initialization for Cmpp model
#' features <- matrix(rnorm(100), ncol = 5)  # 20 observations with 5 features
#' x <- rnorm(20)  # Failure times for 20 observations
#' delta1 <- sample(0:1, 20, replace = TRUE)  # Competing risk indicator
#' delta2 <- sample(0:1, 20, replace = TRUE)  # Event indicator
#' Initialize(features, x, delta1, delta2, h = 1e-5)
Initialize <- function(features, x, delta1, delta2, h) {
  .Call("Initialize", features, x, delta1, delta2, h)
}

#' Compute the CDF of the Gompertz Distribution
#'
#' This function calculates the cumulative distribution function (CDF) of the Gompertz distribution for a given vector of input values, shape, and scale parameters.
#'
#' @param x A numeric vector of input values (failure times, for example) at which to evaluate the CDF.
#' @param alpha A numeric value representing the shape parameter of the Gompertz distribution. It controls the growth rate.
#' @param beta A numeric value representing the scale parameter of the Gompertz distribution. It affects the scale of the distribution.
#'
#' @return A numeric vector of the CDF values evaluated at the given points `x` using the specified parameters `alpha` and `beta`.
#' 
#' @examples
#' # Example CDF calculation
#' x <- c(1, 2, 3)  # Failure times
#' alpha <- 0.5     # Shape parameter
#' beta <- 0.1      # Scale parameter
#' cdf_gomp(x, alpha, beta)
cdf_gomp <- function(x, alpha, beta) {
  .Call("cdf_gomp", x, alpha, beta)
}

#' Compute the PDF of the Gompertz Distribution
#'
#' This function calculates the probability density function (PDF) of the Gompertz distribution for a given vector of input values, shape, and scale parameters.
#'
#' @param x A numeric vector of input values (failure times, for example) at which to evaluate the PDF.
#' @param alpha A numeric value representing the shape parameter of the Gompertz distribution.
#' @param beta A numeric value representing the scale parameter of the Gompertz distribution.
#'
#' @return A numeric vector of the PDF values evaluated at the given points `x` using the specified parameters `alpha` and `beta`.
#' 
#' @examples
#' # Example PDF calculation
#' x <- c(1, 2, 3)  # Failure times
#' alpha <- 0.5     # Shape parameter
#' beta <- 0.1      # Scale parameter
#' pdf_gomp(x, alpha, beta)
pdf_gomp <- function(x, alpha, beta) {
  .Call("pdf_gomp", x, alpha, beta)
}

#' Compute Log-Likelihood for the Model
#'
#' This function computes the negative log-likelihood of the model, given a set of parameters. The log-likelihood is based on the Gompertz distributions and competing risks.
#'
#' @param Param A numeric vector of parameters (alpha, beta, etc.) to be optimized. The function assumes the first two parameters are related to the Gompertz distribution for the first event, and the next two are for the second event.
#'
#' @return The negative log-likelihood value for the given parameters.
#' 
#' @examples
#' # Example log-likelihood computation
#' Param <- c(0.5, 0.1, 0.6, 0.2)  # Example parameters
#' LogLike1(Param)
LogLike1 <- function(Param) {
  .Call("LogLike1", Param)
}

#' Compute the Gradient of the Log-Likelihood
#'
#' This function computes the gradient of the negative log-likelihood with respect to the model parameters using numerical differentiation.
#'
#' @param Param A numeric vector of parameters (alpha, beta, etc.) for which the gradient should be computed.
#'
#' @return A numeric vector representing the gradient of the log-likelihood function at the specified parameters.
#' 
#' @examples
#' # Example gradient computation
#' Param <- c(0.5, 0.1, 0.6, 0.2)  # Example parameters
#' GrLike(Param)
GrLike <- function(Param) {
  .Call("GrLike", Param)
}

#' Compute the Numerical Gradient of the Log-Likelihood
#'
#' This function computes the numerical gradient of the negative log-likelihood using finite differences.
#'
#' @param Param A numeric vector of parameters (alpha, beta, etc.) for which the numerical gradient should be computed.
#'
#' @return A numeric vector representing the gradient of the log-likelihood function at the specified parameters.
#' 
#' @examples
#' # Example numerical gradient computation
#' Param <- c(0.5, 0.1, 0.6, 0.2)  # Example parameters
#' compute_grad(Param)
compute_grad <- function(Param) {
  .Call("compute_grad", Param)
}

#' Estimate Parameters for Parametric CIF Using Optim
#'
#' This function uses the `optim` function to estimate the parameters of the model by minimizing the negative log-likelihood. 
#' Before using this function, ensure that the data is initialized with the `Initialize` function.
#'
#' @param initial_params A numeric vector of initial guesses for the model parameters. The default is a vector of 0.01 for each parameter.
#'
#' @return A list containing the optimization result, including the estimated parameters and other details.
#' 
#' @examples
#' # Example of estimating parameters using optim
#' features <- matrix(rnorm(100), ncol = 5)
#' x <- rnorm(20)
#' delta1 <- sample(0:1, 20, replace = TRUE)
#' delta2 <- sample(0:1, 20, replace = TRUE)
#' Initialize(features, x, delta1, delta2, h = 1e-5)
#' optim_result <- estimate_parameters(c(0.1, 0.2, 0.3, 0.4))
#' optim_result
estimate_parameters <- function(initial_params = rep(0.01, 4)) {  
  # Call optim to minimize the negative log-likelihood
  optim_result <- optim(
    par = initial_params,        # Initial parameters
    fn = LogLike1,               # Log-likelihood function to minimize
    gr = compute_grad,           # Gradient function
    method = "BFGS"              # Optimization method
  )
  return(optim_result)
}
