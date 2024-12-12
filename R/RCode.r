#' Initialize Data for the Cmpp Model
#'
#' This function initializes the data used in the Cmpp model by storing the feature matrix, failure times, 
#' and the competing risks indicators in the model environment. These are required for subsequent computations.
#'
#' @param features A numeric matrix of predictor variables. Each row corresponds to an observation.
#' @param x A numeric vector of failure times corresponding to observations.
#' @param delta1 A binary vector indicating the occurrence of the first competing event (1 for observed).
#' @param delta2 A binary vector indicating the occurrence of the second event (1 for observed).
#' @param h A numeric value specifying the step size for numerical gradient computations.
#'
#' @details This function does not return any value but sets up internal data structures required for model computation.
#' Ensure that `features`, `x`, `delta1`, and `delta2` have matching lengths or dimensions.
#'
#' @usage Initialize(features, x, delta1, delta2, h)
#'
#' @return This function returns `NULL`. The initialized data is stored in the package environment.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' features <- matrix(rnorm(100), ncol = 5)
#' x <- rnorm(20)
#' delta1 <- sample(0:1, 20, replace = TRUE)
#' delta2 <- sample(0:1, 20, replace = TRUE)
#' Initialize(features, x, delta1, delta2, h = 1e-5)
#' }

# Initialize <- function(features, x, delta1, delta2, h) {
#   .Call('cpp_Initialize', as.matrix(features), as.numeric(x), as.integer(delta1), as.integer(delta2), as.numeric(h))
# }

#' Compute the CDF of the Gompertz Distribution
#'
#' Calculates the cumulative distribution function (CDF) of the Gompertz distribution 
#' for given input values and parameters.
#'
#' @param x A numeric vector of non-negative input values (e.g., failure times).
#' @param alpha A positive numeric value representing the shape parameter.
#' @param beta A positive numeric value representing the scale parameter.
#'
#' @details The Gompertz distribution is commonly used in survival analysis and reliability studies.
#' Ensure that `alpha` and `beta` are positive for meaningful results.
#'
#' @usage cdf_gomp(x, alpha, beta)
#'
#' @return A numeric vector of the CDF values for each input in `x`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' x <- c(1, 2, 3)
#' alpha <- 0.5
#' beta <- 0.1
#' cdf_gomp(x, alpha, beta)
#' }

# cdf_gomp <- function(x, alpha, beta) {
#   .Call('cpp_cdf_gomp', as.numeric(x), as.numeric(alpha), as.numeric(beta))
# }

#' Compute the PDF of the Gompertz Distribution
#'
#' Calculates the probability density function (PDF) of the Gompertz distribution 
#' for given input values and parameters.
#'
#' @param x A numeric vector of non-negative input values (e.g., failure times).
#' @param alpha A positive numeric value representing the shape parameter.
#' @param beta A positive numeric value representing the scale parameter.
#'
#' @details The PDF provides the relative likelihood of a failure or event occurring at specific time points.
#' Ensure that `alpha` and `beta` are positive for meaningful results.
#'
#' @usage pdf_gomp(x, alpha, beta)
#'
#' @return A numeric vector of the PDF values for each input in `x`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' x <- c(1, 2, 3)
#' alpha <- 0.5
#' beta <- 0.1
#' pdf_gomp(x, alpha, beta)
#' }
# pdf_gomp <- function(x, alpha, beta) {
#   .Call('cpp_pdf_gomp', as.numeric(x), as.numeric(alpha), as.numeric(beta))
# }




#' Compute the Log-Likelihood for the Model
#'
#' Computes the negative log-likelihood of the Cmpp model given parameters and the initialized data. 
#' The log-likelihood considers Gompertz distributions for competing risks.
#'
#' @param Param A numeric vector of model parameters: [alpha1, beta1, alpha2, beta2], where 
#' the first two are for the first event and the next two are for the second event.
#'
#' @details This function requires the data to be initialized using `Initialize` before being called. 
#' The log-likelihood is based on survival probabilities derived from the Gompertz distributions.
#'
#' @usage LogLike1(Param)
#'
#' @return A single numeric value representing the negative log-likelihood.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Param <- c(0.01, 0.01, 0.01, 0.01)
#' LogLike1(Param)
#' }

# LogLike1 <- function(Param) {
#   .Call('cpp_LogLike1', as.numeric(Param))
# }


#' Compute the Numerical Gradient of the Log-Likelihood
#'
#' Calculates the gradient of the negative log-likelihood using finite differences. 
#' The function uses a small step size (`h`) defined during initialization.
#'
#' @param Param A numeric vector of parameters for which the gradient is calculated.
#'
#' @details This function approximates the gradient using central finite differences.
#' Ensure that `h` is appropriately set to avoid numerical instability.
#'
#' @usage compute_grad(Param)
#'
#' @return A numeric vector of the same length as `Param`, representing the gradient at the specified parameters.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Param <- c(0.5, 0.1, 0.6, 0.2)
#' compute_grad(Param)
#' }
# compute_grad <- function(Param) {
#   .Call('cpp_compute_grad', as.numeric(Param))
# }

#' Clean up memory by deleting the pointer to the Cmpp instance
#'
#' This function is used to clean up and delete the instance of the Cmpp class in
#' the C++ code. It ensures proper memory management and prevents memory leaks by
#' deleting the pointer to the `Cmpp` object when it is no longer needed. 
#' It is important to call this function after you are done with the `Cmpp` object
#' to ensure that no memory is leaked.
#'
#' @usage Cleanup()
#'
#' @details
#' The `cpp_Cleanup` function must be called after using the `Cmpp` object to clean up
#' the allocated memory in C++. Failure to call this function may result in memory 
#' leaks, as the memory allocated for the `Cmpp` object is not automatically freed.
#' 
#' @return NULL
#' @export
#' @examples
#' \dontrun{
#' # Assuming you have previously initialized the Cmpp object with cpp_Initialize()
#' Cleanup()
#' }

# Cleanup <- function() {
#   .Call("cpp_Cleanup")
# }



#' Create a matrix of given size filled with a constant value
#'
#' This function creates an `n x m` matrix of type `Eigen::MatrixXd`, where each
#' element is set to the specified constant value. This is useful for generating
#' matrices with uniform values for testing, initialization, or other purposes in
#' computational tasks where a matrix filled with a constant is needed.
#'
#' @usage makeMat(n, m, value)
#'
#' @details
#' The `cpp_makeMat` function generates a matrix with the given dimensions `n x m`
#' where all elements are initialized to the same constant value. It is useful
#' in scenarios where a specific value needs to be assigned to all elements of the
#' matrix, for example in machine learning algorithms, matrix manipulations, or tests.
#'
#' @param n An integer representing the number of rows in the matrix.
#' @param m An integer representing the number of columns in the matrix.
#' @param value A numeric value that will be used to fill the matrix.
#'
#' @return A numeric matrix of dimensions `n x m` filled with the specified value.
#' @export
#' @examples
#' \dontrun{
#' # Create a 3x3 matrix filled with 5
#' mat <- makeMat(3, 3, 5)
#' print(mat)
#' }

# makeMat <- function(n, m, value) {
#   .Call("cpp_makeMat", n, m, value)
# }



#' Get Dimensions of the Cmpp Object
#'
#' This function returns the number of samples and features stored in the Cmpp object.
#' It is primarily used to retrieve the dimensions of the data within the class.
#'
#' @return A list containing:
#'   \item{Nsamp}{Number of samples (rows in the feature matrix).}
#'   \item{Nfeature}{Number of features (columns in the feature matrix).}
#'
#' @export
#' @usage GetDim()
#' @details
#' The `GetDim` function allows users to access the internal dimensions of the
#' `Cmpp` class instance, such as the number of samples (`Nsamp`) and the number of features
#' (`Nfeature`). This is useful when working with large datasets, especially for
#' checking the size of the data without needing to manually access the underlying `Eigen::MatrixXd`
#' or `Eigen::VectorXd` objects directly.
#'
#' @examples
#' 
#' # Initialize Cmpp object
#' Initialize(features, x, delta1, delta2, h)
#' 
#' # Get dimensions
#' dims <- GetDim()
#' dims$Nsamp    # Number of samples
#' dims$Nfeature # Number of features
# GetDim <- function() {
#   # Call the C++ method via Rcpp
#   result <- .Call("cpp_GetDim")
#   return(result)
# }



#' Estimate Model Parameters Using Optimization
#'
#' This function estimates the parameters of a model by minimizing the negative
#' log-likelihood function using the specified optimization method. It utilizes
#' the `optim()` function in R, with the provided initial parameter values and
#' gradient computation. The optimization method can be specified, with "BFGS" being
#' the default.
#'
#' @usage estimate_parameters(initial_params = rep(0.01, 4), optimMethod = 'BFGS')
#'
#' @details
#' The `estimate_parameters` function performs parameter estimation by minimizing
#' the negative log-likelihood function using the chosen optimization method. The function
#' requires an initial guess for the parameters (a numeric vector) and will optimize the
#' log-likelihood function. The optimization also takes into account the gradient of the 
#' log-likelihood function, which is computed using the `compute_grad` function. The result 
#' of the optimization is an object of class `optim` containing the estimated parameters and
#' other details of the optimization process.
#'
#' The optimization method can be specified via the `optimMethod` argument. The default method
#' is "BFGS", but any method supported by R's `optim()` function (such as "Nelder-Mead", "CG", etc.)
#' can be used.
#'
#' @param initial_params A numeric vector of initial parameter values to start the optimization.
#'                        Default is a vector of four values, all set to 0.01.
#' @param optimMethod A character string specifying the optimization method to use. 
#'                    The default is `"BFGS"`. See `?optim` for a list of available methods.
#'
#' @return An `optim` object containing the optimization results, including the estimated
#'         parameters, value of the objective function at the optimum, and other optimization details.
#'
#' @seealso \link{optim} for more details on optimization methods and usage.
#'
#' @export
#' @examples
#' # Estimate model parameters using default initial values and the BFGS method
#' result <- estimate_parameters()
#' print(result)
#'
#' # Estimate model parameters with custom initial values and the Nelder-Mead method
#' result <- estimate_parameters(initial_params = c(0.1, 0.1, 0.1, 0.1), optimMethod = "Nelder-Mead")
#' print(result)
estimate_parameters <- function(initial_params = rep(0.01, 4), optimMethod = 'BFGS') { 
  optim_result <- optim(
    par = initial_params,        # Initial parameters
    fn = LogLike1,               # Log-likelihood function to minimize
    gr = compute_grad,           # Gradient function
    method = optimMethod         # Optimization method 
  )
  return(optim_result)
}
