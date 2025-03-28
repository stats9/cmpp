packageVersion('cmpp')
remove.packages('cmpp')
devtools :: document()
devtools :: install()
library(cmpp)
help(package = 'cmpp')


data(dat)
names(dat)
feat <- dat[, -c(match(c('id', 'time', 'event', 'cause_burn', 'd1', 'd2', 'cause_hotObject3'), names(dat)))]  
timee <- dat[['time']]
d1 = dat[['d1']]
d2 = dat[['d2']]
feat2 = feat |> data.matrix()
Initialize(feat2, timee, d1, d2, 1e-10)
optim(par = c(0.001, 0.001, 0.001, 0.001), 
    fn = LogLike1, 
    gr = compute_grad, 
    method = 'BFGS', hessian = TRUE)


(parr <- estimate_parameters(c(0.001, 0.001, 0.001, 0.001))$par)
tempp <- -compute_hessian(parr)
solve(-tempp)
tempp2 <- -numDeriv :: hessian(func = LogLike1, x = parr)

mat1 <- tempp[1:2, 1:2] 
mat2 <- tempp[3:4, 3:4]
mat1 |> solve()
mat2 |> solve()
mat1 |> eigen() |> _$values
mat2 |> eigen() |> _$values

t1 <- Sys.time()
ress <- bootstrap_variance(feat2, timee, d1, d2, 
    initial_params = rep(0.001, 4), n_bootstrap = 500, optimMethod = 'BFGS')
t2 <- Sys.time()
(elapsedTime <- difftime(t2, t1, units = 'secs'))

CIF_res1()
CIF_Figs(initial_params = rep(0.01, 4), timee) 

############ check install from github

devtools::install_github("stats9/cmpp")
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
feat2 |> dim()
# Estimate parameters using LBFGS++
initial_params <- c(0.001, 0.001, 0.001, 0.001)
params <- estimate_parameters(initial_params)
print(params)
colnames(feat2)
initial_params <- rep(0.001, 34)
(optim_result <- optim(
    par = initial_params,        # Initial parameters
    fn = log_f_rcpp,               # Log-likelihood function to minimize
    gr = compute_log_f_gradient_rcpp,           # Gradient function
    method = "BFGS",         # Optimization method 
    control = list(fnscale = -1)
))
log_f_rcpp(rep(0.001, 34))
estimate_parameters_GCM(initial_params)

initial_params <- rep(0.001, 32)
estimate_parameters_POM(initial_params)

log_f_rcpp(rep(0.001, 32))
estimate_parameters_PHM(initial_params)
library(cmpp)
help(package = "cmpp")# Load necessary package
if (!requireNamespace("numDeriv", quietly = TRUE)) {
  install.packages("numDeriv")
}
library(numDeriv)

# Define the objective function
objective_function <- function(Params) {
  log_f_rcpp(Params)  # Call the Rcpp function that computes the log-likelihood
}

# Compute Hessian numerically
compute_hessian_r <- function(Params) {
  hessian_matrix <- hessian(func = objective_function, x = Params)
  return(hessian_matrix)
}

# Example Usage
library(numDeriv)
hessian_matrix <- compute_hessian_r(initial_params)
print(hessian_matrix)
solve(-hessian_matrix)

# Compute Hessian and other metrics
hessian <- compute_hessian(params$par)
print(hessian)

# Bootstrap variance estimation
# results <- bootstrap_variance(feat2, timee, d1, d2, initial_params, n_bootstrap = 500, "BFGS")
print(results$variances)
print(results$bootstrap_estimates)

# Compute CIF results
cif_results <- CIF_res1(initial_params)
print(cif_results)

# Plot CIFs with confidence intervals
plot <- CIF_Figs(initial_params, timee)
print(plot)

dat2 <- readRDS('./Assisst_objects/PragnancyData.rds')
dat2 |> names()
dat2$edu2
feat <- dat2[, -c(match(c('id', 't', 'event',  'd1', 'd2'), names(dat2)))]  
feat |> as.data.frame() ->  feat
feat |> typeof()

feat$edu |> unlist() |> table()
timee <- dat2[['t']]
d1 <- dat[['d1']]
d2 <- dat[['d2']]
feat2 <- feat |> data.matrix()
feat2 |> _[, 1] |> table()
feat2 <- feat2[, -1]
Initialize(feat2, timee, d1, d2, 1e-10)

# Estimate parameters using LBFGS++
initial_params <- c(0.001, 0.001, 0.001, 0.001)
params <- estimate_parameters(initial_params)
print(params)
cif_results <- CIF_res1(initial_params)
print(cif_results)

tempPar <- c(params$par, 0.01, -0.01)
CIF_Figs(initial_params = rep(0.01, 4), timee) 
f_pdf_rcpp(params$par, feat2, 15)
F_cdf_rcpp(params$par, feat2, 1000)



F_cdf <- function(Params, Z, x) {
    Params = tempPar; Z = feat2[1, ]; x = 21
    alpha <- Params[1] 
    tau <- Params[2] 
    rho <- Params[3] 
    rho <- (rho < 0) * rho - (rho > 0) * rho
    Beta <- Params[4:length(Params)]
    tempval = crossprod(Z, Beta) |> as.numeric()
    temp1 <- 1 - ( 1 + alpha * 
        exp(tempval) * tau * (exp(rho * x) - 1)/rho)^(-1/alpha)
    return(temp1)
}

f_pdf <- function(Params = parr1, Z = z, x = xx) {
    Params = tempPar; Z = feat2[1, ]; x = 21
    alpha <- Params[1] 
    tau <- Params[2] 
    rho <- Params[3] 
    rho <- (rho < 0) * rho - (rho > 0) * rho
    Beta <- Params[4:length(Params)]
    tempval = crossprod(Z, Beta) |> as.numeric()
    (tau * (alpha * tau * ((exp(rho * x) - 1) * exp(tempval))/rho + 1)^(-1/alpha) * 
        (exp(tempval) * exp(rho * x))) / (((alpha * tau * (exp(rho * x) - 1) * 
            exp(tempval)) / rho ) + 1) -> result 
    names(result) <- NULL
    return(result)
}
F_cdf(tempPar, feat2, 21)
F_cdf_rcpp(tempPar, feat2[1, ], 21)
f_pdf_rcpp(tempPar, feat2[1, ], 21)
feat2 |> dim()


####### Add two other log-likelihood functions
