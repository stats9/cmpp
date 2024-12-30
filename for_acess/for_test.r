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

# Estimate parameters using LBFGS++
initial_params <- c(0.001, 0.001, 0.001, 0.001)
params <- estimate_parameters(initial_params)
print(params)

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
timee <- dat2[['t']]
d1 <- dat[['d1']]
d2 <- dat[['d2']]
feat2 <- feat |> data.matrix()
Initialize(feat2, timee, d1, d2, 1e-10)

# Estimate parameters using LBFGS++
initial_params <- c(0.001, 0.001, 0.001, 0.001)
params <- estimate_parameters(initial_params)
print(params)
cif_results <- CIF_res1(initial_params)
print(cif_results)

CIF_Figs(initial_params = rep(0.01, 4), timee) 
