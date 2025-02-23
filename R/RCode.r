#' @name Initialize
#' @title Initialize Data for the Cmpp Model
#'
#' @description This function initializes the data used in the Cmpp model by storing the feature matrix, failure times, 
#'    and the competing risks indicators in the model environment. These are required for subsequent computations.
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
NULL 

#' @name cdf_gomp
#' @title Compute the CDF of the Gompertz Distribution
#'
#' @description Calculates the cumulative distribution function (CDF) of the Gompertz distribution 
#'    for given input values and parameters.
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
NULL

#' @name pdf_gomp
#' @title Compute the PDF of the Gompertz Distribution
#'
#' @description Calculates the probability density function (PDF) of the Gompertz distribution 
#'    for given input values and parameters.
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
NULL 

#' @name LogLike1
#' @title Compute the Log-Likelihood for the Model
#'
#' @description Computes the negative log-likelihood of the Cmpp model given parameters and the initialized data. 
#'    The log-likelihood considers Gompertz distributions for competing risks.
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
NULL 

#' @name compute_grad
#' @title Compute the Numerical Gradient of the Log-Likelihood
#'
#' @description Calculates the gradient of the negative log-likelihood using finite differences. 
#'    The function uses a small step size (`h`) defined during initialization.
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
#' @examples
#' \dontrun{
#' Param <- c(0.5, 0.1, 0.6, 0.2)
#' compute_grad(Param)
#' }
NULL 

#' @name Cleanup
#' @title Clean up memory by deleting the pointer to the Cmpp instance
#'
#' @description This function is used to clean up and delete the instance of the Cmpp class in
#'    the C++ code. It ensures proper memory management and prevents memory leaks by
#'    deleting the pointer to the `Cmpp` object when it is no longer needed. 
#'    It is important to call this function after you are done with the `Cmpp` object
#'    to ensure that no memory is leaked.
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
NULL 

#' @name makeMat
#' @title Create a matrix of given size filled with a constant value
#'
#' @description This function creates an `n x m` matrix of type `Eigen::MatrixXd`, where each
#'    element is set to the specified constant value. This is useful for generating
#'    matrices with uniform values for testing, initialization, or other purposes in
#'    computational tasks where a matrix filled with a constant is needed.
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
NULL

#' @name GetDim
#' @title Get Dimensions of the Cmpp Object
#'
#' @description This function returns the number of samples and features stored in the Cmpp object.
#'    It is primarily used to retrieve the dimensions of the data within the class.
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
#' \dontrun{
#' # Initialize Cmpp object
#' Initialize(features, x, delta1, delta2, h)
#' # Get dimensions
#' dims <- GetDim()
#' dims$Nsamp    # Number of samples
#' dims$Nfeature # Number of features
#'}
NULL 

#' @name estimate_parameters
#' @title Estimate Model Parameters Using Optimization
#' @description  This function estimates the parameters of a model by minimizing the negative
#'    log-likelihood function using the specified optimization method. It utilizes
#'    the `optim()` function in R, with the provided initial parameter values and
#'    gradient computation. The optimization method can be specified, with "BFGS" being
#'    the default.
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
#' @examples
#'
#'\dontrun{
#' # Estimate model parameters using default initial values and the BFGS method
#' result <- estimate_parameters()
#' print(result)
#' }
#' 
#' @export
#'
estimate_parameters <- function(initial_params = rep(0.01, 4), optimMethod = 'BFGS') { 
  optim_result <- optim(
    par = initial_params,        # Initial parameters
    fn = LogLike1,               # Log-likelihood function to minimize
    gr = compute_grad,           # Gradient function
    method = optimMethod         # Optimization method 
  )
  return(optim_result)
}


#' @name compute_hessian
#' @title Compute the Hessian Matrix of the Log-Likelihood
#'
#' @description Calculates the Hessian matrix of the negative log-likelihood function using finite differences.
#'    This function is useful for understanding the curvature of the log-likelihood surface and for optimization purposes.
#'
#' @param Param A numeric vector of parameters for which the Hessian matrix is calculated.
#'
#' @details This function approximates the Hessian matrix using central finite differences.
#' Ensure that the step size `h` is appropriately set during initialization to avoid numerical instability.
#' The function requires the data to be initialized using `Initialize` before being called.
#'
#' @usage compute_hessian(Param)
#'
#' @return A numeric matrix representing the Hessian matrix at the specified parameters.
#'
#' @export
#' @examples
#' \dontrun{
#' Param <- c(0.5, 0.1, 0.6, 0.2)
#' hessian <- compute_hessian(Param)
#' print(hessian)
#' }
NULL

#' @name bootstrap_variance
#' @title Estimate Variance of Parameters Using Bootstrap Method
#'
#' @description This function estimates the variance of model parameters using the bootstrap method.
#'    It repeatedly samples the data with replacement, estimates the parameters for each sample,
#'    and computes the variance of the estimated parameters.
#'
#' @param features A numeric matrix of predictor variables. Each row corresponds to an observation.
#' @param x A numeric vector of failure times corresponding to observations.
#' @param delta1 A binary vector indicating the occurrence of the first competing event (1 for observed).
#' @param delta2 A binary vector indicating the occurrence of the second event (1 for observed).
#' @param initial_params A numeric vector of initial parameter values to start the optimization.
#' @param n_bootstrap An integer specifying the number of bootstrap samples.
#' @param optimMethod A character string specifying the optimization method to use. Default is `"BFGS"`.
#'
#' @details This function performs bootstrap sampling to estimate the variance of the model parameters.
#'    It requires the data to be initialized using `Initialize` before being called.
#'
#' @usage bootstrap_variance(features, x, delta1, delta2, initial_params, n_bootstrap, optimMethod)
#'
#' @return A list containing:
#'    - `variances`: A numeric vector representing the variance of the estimated parameters.
#'    - `bootstrap_estimates`: A matrix of parameter estimates for each bootstrap sample.
#'
#' @export
#' @examples
#' \dontrun{
#' features <- matrix(rnorm(100), ncol = 5)
#' x <- rnorm(20)
#' delta1 <- sample(0:1, 20, replace = TRUE)
#' delta2 <- sample(0:1, 20, replace = TRUE)
#' initial_params <- c(0.01, 0.01, 0.01, 0.01)
#' n_bootstrap <- 1000
#' results <- bootstrap_variance(features, x, delta1, delta2, initial_params, n_bootstrap)
#' print(results$variances)
#' print(results$bootstrap_estimates)
#' }
NULL

#' @name CIF_res1
#' @title Compute Cumulative Incidence Function (CIF) Results
#'
#' @description This function estimates the parameters of the model, computes the Hessian matrix, and calculates the variances and p-values for the parameters. It ensures that the diagonal elements of the covariance matrix are positive.
#'
#' @param initial_params A numeric vector of initial parameter values to start the optimization. Default is `rep(0.001, 4)`.
#'
#' @details This function performs the following steps:
#' \itemize{
#'   \item Estimates the model parameters using the `estimate_parameters` function.
#'   \item Computes the Hessian matrix using the `compute_hessian` function.
#'   \item Ensures that the diagonal elements of the covariance matrix are positive.
#'   \item Calculates the variances and p-values for the parameters.
#' }
#'
#' @return A data frame containing:
#' \item{Params}{The parameter names ("alpha1", "beta1", "alpha2", "beta2").}
#' \item{STD}{The standard deviations of the parameters.}
#'
#' @examples
#' \dontrun{
#' initial_params <- c(0.001, 0.001, 0.001, 0.001)
#' result <- CIF_res1(initial_params)
#' print(result)
#' }
#'
#' @export
#'
CIF_res1 <- function(initial_params = rep(0.001, 4)) {
  Params = estimate_parameters(initial_params = initial_params)$par
  hessian_mat <- compute_hessian(Params)
  Information_matrix <- - hessian_mat

  # for confident that Diagonal elements of covariance Matrix are positive 
  subMat1 <- Information_matrix[1:2, 1:2] |> solve()
  diag(subMat1) <- abs(diag(subMat1))
  subMat2 <- Information_matrix[3:4, 3:4] |> solve()
  diag(subMat2) <- abs(diag(subMat2))
  var_alpha1 <- subMat1[1, 1]
  var_beta1 <- subMat1[2, 2]
  var_alpha2 <- subMat2[1, 1]
  var_beta2 <- subMat2[2, 2]
  pval_alpha1 <- 2 * min(c(
  pnorm(Params[1], mean = 0, sd = sqrt(var_alpha1)), 
  pnorm(Params[1], mean = 0, sd = sqrt(var_alpha1), lower.tail = FALSE)
  ))

result <- data.frame(
Params = c("alpha1", "beta1", "alpha2", "beta2"), 
Estimation = Params, 
STD = c(sqrt(var_alpha1), sqrt(var_beta1), sqrt(var_alpha2), sqrt(var_beta2))
)

return (result)
}
NULL



#' @name CIF_Figs
#' @title Plot Cumulative Incidence Functions (CIF) with Confidence Intervals
#'
#' @description This function plots the cumulative incidence functions (CIF) for two competing risks based on the estimated parameters and their variances. It includes confidence intervals for the CIFs.
#'
#' @param initial_params A numeric vector of initial parameter values to start the optimization.
#' @param TimeFailure A numeric vector of failure times corresponding to observations.
#' @param OrderType A numeric vector indicating the order of the competing risks. Default is `c(2, 1)`.
#' @param RiskNames A character vector of names for the competing risks. Default is `NULL`.
#'
#' @details This function performs the following steps:
#' \itemize{
#'   \item Estimates the model parameters using the `estimate_parameters` function.
#'   \item Computes the Hessian matrix using the `compute_hessian` function.
#'   \item Ensures that the diagonal elements of the covariance matrix are positive.
#'   \item Computes the cumulative incidence functions (CIF) for two competing risks.
#'   \item Plots the CIFs along with their confidence intervals.
#' }
#'
#' @return A ggplot object showing the CIFs and their confidence intervals.
#'
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @export
#' @examples
#' \dontrun{
#' initial_params <- c(0.01, 0.01, 0.01, 0.01)
#' TimeFailure <- seq(0, 10, by = 0.1)
#' plot <- CIF_Figs(initial_params, TimeFailure)
#' print(plot)
#' }

CIF_Figs <- function(initial_params, TimeFailure, OrderType = c(2, 1), RiskNames = NULL) {
  Params = estimate_parameters(initial_params = initial_params)$par
  hessian_mat <- compute_hessian(Params)
  Information_matrix <- -hessian_mat
  
  # Ensure diagonal elements of covariance matrix are positive
  subMat1 <- solve(Information_matrix[1:2, 1:2])
  diag(subMat1) <- abs(diag(subMat1))
  subMat2 <- solve(Information_matrix[3:4, 3:4])
  diag(subMat2) <- abs(diag(subMat2))
  
  xtime <- TimeFailure
  diff_CIF <- function(x, a, b) {
    temp1 <- (b/a^2) * (1 + exp(a*x) * (a * x - 1)) * exp(b * (1 - exp(a * x)) / a)
    temp2 <- exp(b * (1-exp(a*x))/a) * (1/a) * (exp(a*x) - 1) * exp(b * (1 - exp(a * x))/a)
    return(c(temp1, temp2))
  }
  diff_CIF_1 <- function(x) diff_CIF(x, a = Params[1], b = Params[2])
  diff_CIF_2 <- function(x) diff_CIF(x, a = Params[3], b = Params[4])
  
  varHat_CIF1 <- lapply(xtime, FUN = function(x) {
    temp1 <- diff_CIF_1(x) 
    res <- temp1 %*% subMat1 %*% temp1 
    return(as.numeric(res))
  }) |> unlist() |> abs()
  
  varHat_CIF2 <- lapply(xtime, FUN = function(x) {
    temp1 <- diff_CIF_2(x) 
    res <- temp1 %*% subMat2 %*% temp1 
    return(as.numeric(res))
  }) |> unlist() |> abs()
  ## Define CIF 
  Fk <- function(x, a, b) {
    res <- 1 - exp(b * (1 - exp(a*x))/a)
    return(res)
  }
  
  CIF1 <- function(x) Fk(x = x, a = Params[1], b = Params[2])
  CIF2 <- function(x) Fk(x = x, a = Params[3], b = Params[4]) 
  
  Lower_Bonds_1 <- CIF1(xtime) - 1.96 * sqrt(varHat_CIF1)
  Lower_Bonds_1[Lower_Bonds_1 < 0] <- 0 
  
  Upper_Bonds_1 <- CIF1(xtime) + 1.96 * sqrt(varHat_CIF1)
  Upper_Bonds_1[Upper_Bonds_1 > 1] <- 1
  
  Lower_Bonds_2 <- CIF2(xtime) - 1.96 * sqrt(varHat_CIF2)
  Lower_Bonds_2[Lower_Bonds_2 < 0] <- 0 
  
  Upper_Bonds_2 <- CIF2(xtime) + 1.96 * sqrt(varHat_CIF2)
  Upper_Bonds_2[Upper_Bonds_2 > 1] <- 1
  
  cif1 <- CIF1(xtime)
  cif2 <- CIF2(xtime)
  cifdat <- data.frame( 
    Time = xtime,
    cif_1 = cif1, 
    LBond_1 = Lower_Bonds_1, 
    UBond_1 = Upper_Bonds_1, 
    cif_2 = cif2, 
    LBond_2 = Lower_Bonds_2, 
    UBond_2 = Upper_Bonds_2
  )
  
  dat11 <- cifdat |> 
              subset(select = c(Time, cif_1, cif_2))
  dat22 <- cifdat |> 
            subset(select = c(LBond_1, LBond_2))
  dat33 <- cifdat |> 
            subset(select = c(UBond_1, UBond_2))

  ldat1 <- dat11 |> 
            tidyr::pivot_longer(!Time,
            values_to = 'CIF', 
            names_to = "cif_fun") 

  if(is.null(RiskNames)) {
    Levelnames <- c("Event", "CompettingRisk")[OrderType]
  } else {
    Levelnames <- RiskNames 
  }

  ldat1 <- ldat1 |> 
            dplyr::mutate(groups = rep(Levelnames, nrow(ldat1)/2))

  ldat2 <- dat22 |> 
            tidyr::pivot_longer(cols = everything(), 
            values_to = 'Confidence_Lower', 
            names_to = "LBond") 

  ldat3 <- dat33 |> 
            tidyr::pivot_longer(cols = everything(), 
            values_to = 'Confidence_Upper', 
            names_to = "UBond") 

  new_dat <- cbind(ldat1, ldat2, ldat3)

  plot <- new_dat |> ggplot2::ggplot(ggplot2::aes(x = Time)) + 
    ggplot2::geom_line(ggplot2::aes(y = CIF, group = cif_fun, color = groups)) + 
    ggplot2::geom_line(ggplot2::aes(y = Confidence_Lower, group = LBond, color = groups), linetype = 2) + 
    ggplot2::geom_line(ggplot2::aes(y = Confidence_Upper, group = UBond, color = groups), linetype = 2) + 
    ggplot2::theme_bw() + 
    ggplot2::labs(title = "Cumulative Incidence Function with Confidence Interval")

  return(plot)
}

#' ################## add 2025-2-22 %%%%%%%%%%%%%%%%%%%%%%%%

#' Create Dummy Variables
#'
#' This function creates dummy variables for specified features in a dataset.
#'
#' @param Data A data frame containing the data.
#' @param features A character vector of feature names for which dummy variables are to be created.
#' @param reff A character string indicating the reference level. It can be either "first" or "last".
#' @return A list containing two elements: 
#' \item{New_Data}{A data frame with the original data and the newly created dummy variables.}
#' \item{Original_Data}{The original data frame.}
#' @examples
#' dat <- data.frame(sex = c('M', 'F', 'M'), cause_burn = c('A', 'B', 'A'))
#' result <- make_Dummy(Data = dat, features = c('sex', 'cause_burn'), reff = "first")
#' print(result$New_Data)
#' @export
make_Dummy <- function(Data = dat, features = c('sex', 'cause_burn'), reff = "first") {
    tempDat <- data.frame(Temp_column = rep(NA, nrow(Data)))
    for (j in features) {
        temp <- Data[[j]]
        temp2 <- as.factor(temp) |> levels()
        index <- ifelse(reff == "last", length(temp2), 1)
        for(h in temp2[-index]) {
            temp3 <- 1 * (as.factor(temp) == h)
            index2 <- which(as.factor(temp2) == h)
            name_temp <- paste(j, index2, sep = ":")
            tempDat[[name_temp]] <- temp3
        }
    }
    f_dat <- tempDat[, -1]
    index_name <- match(features, names(Data))
    dat2 <- Data |> 
        subset(select = -index_name)
    return_data <- cbind(dat2, f_dat) |> as.data.frame() 
    return(list(New_Data = return_data, Original_Data = Data))
}

#' @name f_pdf_rcpp
#' @title Compute the PDF of the Parametric Model
#' @description This function computes the probability density function (PDF) of the parametric model.
#' @param Params A numeric vector of parameters.
#' @param Z A numeric vector of covariates.
#' @param x A numeric value representing the time point.
#' @return A numeric value representing the PDF.
#' @export
#' @examples
#' \dontrun{
#' library(cmpp)
#' features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
#' delta1 <- sample(c(0, 1), 100)
#' delta2 <- 1 - delta1
#' x <- rexp(100, rate = 1/10)
#' Initialize(features, x, delta1, delta2, h = 1e-5)
#' params <- rep(0.001, 2 * (ncol(features) + 3))
#' pdf_value <- f_pdf_rcpp(params, Z[1, ], x[3])
#' print(pdf_value)
#' }
NULL

#' @name F_cdf_rcpp
#' @title Compute the CDF of the Parametric Model
#' @description This function computes the cumulative distribution function (CDF) of the parametric model.
#' @param Params A numeric vector of parameters.
#' @param Z A numeric vector of covariates.
#' @param x A numeric value representing the time point.
#' @return A numeric value representing the CDF.
#' @export
#' @examples
#' \dontrun{
#' library(cmpp)
#' features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
#' delta1 <- sample(c(0, 1), 100)
#' delta2 <- 1 - delta1
#' x <- rexp(100, rate = 1/10)
#' Initialize(features, x, delta1, delta2, h = 1e-5)
#' params <- rep(0.001, 2 * (ncol(features) + 3))
#' x <- 5
#' cdf_value <- F_cdf_rcpp(params, features[1, ], x)
#' print(cdf_value)
#' }
NULL

#' @name log_f_rcpp
#' @title Compute the Log-Likelihood Function
#' @description This function computes the log-likelihood function for the parametric model.
#' @param Params A numeric vector of parameters.
#' @return A numeric value representing the log-likelihood.
#' @export
#' @examples
#' \dontrun{
#' library(cmpp)
#' features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
#' delta1 <- sample(c(0, 1), 100)
#' delta2 <- 1 - delta1
#' x <- rexp(100, rate = 1/10)
#' Initialize(features, x, delta1, delta2, h = 1e-5)
#' params <- rep(0.001, 2 * (ncol(features) + 3))
#' log_likelihood <- log_f_rcpp(params)
#' print(log_likelihood)
#' }
NULL

#' @name compute_log_f_gradient_rcpp
#' @title Compute the Gradient of the Log-Likelihood Function
#' @description This function computes the gradient of the log-likelihood function for the parametric model.
#' @param Params A numeric vector of parameters.
#' @return A numeric vector representing the gradient of the log-likelihood.
#' @export
#' @examples
#' \dontrun{
#' features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
#' delta1 <- sample(c(0, 1), 100)
#' delta2 <- 1 - delta1
#' x <- rexp(100, rate = 1/10)
#' Initialize(features, x, delta1, delta2, h = 1e-5)
#' params <- rep(0.001, 2 * (ncol(features) + 3))
#' gradient <- compute_log_f_gradient_rcpp(params)
#' print(gradient)
#' }
NULL

#' @name compute_log_f_hessian_rcpp
#' @title Compute the Hessian Matrix of the Log-Likelihood Function
#' @description This function computes the Hessian matrix of the log-likelihood function for the parametric model.
#' @param Params A numeric vector of parameters.
#' @return A numeric matrix representing the Hessian matrix of the log-likelihood.
#' @export
#' @examples
#' \dontrun{
#' library(cmpp)
#' features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
#' delta1 <- sample(c(0, 1), 100)
#' delta2 <- 1 - delta1
#' x <- rexp(100, rate = 1/10)
#' Initialize(features, x, delta1, delta2, h = 1e-5)
#' params <- rep(0.001, 2 * (ncol(features) + 3))
#' hessian <- compute_log_f_hessian_rcpp(params)
#' print(hessian)
#' }
NULL


#' @name estimate_parameters2
#' @title Compute Cumulative Incidence Regression Results
#'
#' @description This function estimates the parameters of the model, computes the Hessian matrix, and calculates the variances and p-values for the parameters. It ensures that the diagonal elements of the covariance matrix are positive.
#'
#' @param initial_params A numeric vector of initial parameter values to start the optimization. Default is `rep(0.001, 2 * (3 + ncol(covars)))`.
#'
#' @details This function performs the following steps:
#' \itemize{
#'   \item Estimates the model parameters using the `optim` function and `log_f_rcpp` and `compute_log_f_gradient_rcpp`, `compute_log_f_hessian_rcpp`
#'   \item Computes the Hessian matrix using the `` function.
#'   \item Ensures that the diagonal elements of the covariance matrix are positive.
#'   \item Calculates the variances and p-values for the parameters.
#' }
#'
#' @return A data frame containing:
#' \item{Params}{The parameter names ('alpha1', "beta1", 'rho1', 'alpha2', 'beta2', 'rho2', 'beta11', ..., 'beta1k', 'beta21', ..., 'beta2k').}
#' \item{STD}{The standard deviations of the parameters.}
#' \item{Pvalue}{p-values.}
#'
#' @seealso \link{optim}, \link{compute_log_f_gradient_rcpp}, \link{log_f_rcpp}, \link{compute_log_f_hessian_rcpp}.
#'
#' @importFrom numDeriv hessian
#'
#' @examples
#' \dontrun{
#' library(cmpp)
#' features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
#' delta1 <- sample(c(0, 1), 100, replace = TRUE)
#' delta2 <- 1 - delta1
#' x <- rexp(100, rate = 1/10)
#' Initialize(features, x, delta1, delta2, h = 1e-5)
#' initial_params <- rep(0.001, 2 * (ncol(features) + 3))
#' result <- estimate_parameters2(initial_params)
#' print(result)
#' }
#'
#' @export
estimate_parameters2 <- function(initial_params) {
  # Optimize the log-likelihood function to estimate parameters
  tempk <- length(initial_params) - 6
  tempk <- tempk/2
  optim_result <- optim(
    par = initial_params,
    fn = log_f_rcpp,
    gr = compute_log_f_gradient_rcpp,
    method = "BFGS",
    control = list(fnscale = -1)
  )
  
  # Extract estimated parameters
  estimated_params <- optim_result$par
  
  # Compute the Hessian matrix at the estimated parameters
  # Define the objective function
  objective_function <- function(Params) {
    log_f_rcpp(Params)  # Call the Rcpp function that computes the log-likelihood
  }

  # Compute Hessian numerically
  compute_hessian_r <- function(Params) {
  hessian_matrix <- numDeriv :: hessian(func = objective_function, x = Params)
  return(hessian_matrix)
}

  hessian_matrix <- compute_hessian_r(initial_params)

  # Compute the standard deviations of the parameters
  std_devs <- diag(solve(-hessian_matrix)) |> abs() |> sqrt()
  
  # Compute the p-values for the parameters
  p_values <- 2 * (1 - pnorm(abs(estimated_params / std_devs)))
  
  # Create a data frame with parameter names, estimates, standard deviations, and p-values
  param_names <- c('alpha1', 'beta1', 'rho1', paste("beta1", 1:tempk, sep = ":"), 
  'alpha2', 'beta2', 'rho2', paste("beta2", 1:tempk, sep = ":"))
  result_df <- data.frame(
    Parameter = param_names,
    Estimate = estimated_params,
    StdDev = std_devs,
    PValue = p_values
  )
  
  return(result_df)
}
NULL
