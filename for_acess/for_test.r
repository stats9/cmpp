packageVersion('cmpp')
remove.packages('cmpp')
devtools :: document()
devtools :: install()
library(cmpp)
# help(package = 'cmpp')


data(dat)
names(dat)
feat <- dat[, -c(match(c('id', 'time', 'event', 'cause_burn', 'd1', 'd2', 'cause_hotObject3'), names(dat)))]  
timee <- dat[['time']]
d1 = dat[['d1']]
d2 = dat[['d2']]
feat2 = feat |> data.matrix()
Initialize(feat2, timee, d1, d2, 1e-10)
FineGray_Model()
Res <- Cmpp_CIF()
Res$Plot$PlotNull_AllModels
Res$Plot$Plot_InputModel
Res$Plot$PlotAdjusted_AllModels
Res2 <- Cmpp_CIF(, )
GetData()

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


############## add 2025-04-15 ############
Cmpp_CIF <- function(featureID = NULL, featreValue = NULL, RiskNames = NULL, 
    TypeMethod = "GOM", predTime = NULL) {
    
    if(!TypeMethod %in% c("GOM", "POM", "PHM")) {
      stop("TypeMethod must be one of `GOM` or `POM` or `PHM`.")    
    }
    if(is.null(featreValue)) {
       Features <- GetData()$features
      z <- colMeans(Features)
    } else {
  if(is.null(featureID)) {
  stop("To compute the CIF for specific values of certain features, include their indices in the featureID argument!")
  } else {
    if(length(featureID) != length(featreValue)) {
      stop("The length of featureID and featreValue must be the same!")
    } else {
      Features <- GetData()$features
      z <- colMeans(Features)
      for(i in 1:length(featureID)) {
        z[featureID[i]] <- featreValue[i]
      }
    }
  }
}
  if(is.null(RiskNames)) {
    RiskNames <- c("Risk1", "Risk2")
  } else {
    if(length(RiskNames) != 2) {
      stop("RiskNames must be a vector of length 2!")
    }
  }

  Par1 <- estimate_parameters_GCM(rep(0.01, 2*(3 + GetDim()$Nfeature)))
  Par2 <- estimate_parameters_POM(rep(0.01, 2*(2 + GetDim()$Nfeature)))
  Par3 <- estimate_parameters_PHM(rep(0.01, 2*(2 + GetDim()$Nfeature)))
  Par11 <- Par1[1:3, 2]
  Par12 <- Par1[(4 + GetDim()$Nfeature):(4 + GetDim()$Nfeature + 2), 2] 
  Par21 <- Par2[1:2, 2]
  Par22 <- Par2[(3 + GetDim()$Nfeature):(3 + GetDim()$Nfeature + 1), 2]
  Par31 <- Par3[1:2, 2]
  Par32 <- Par3[(3 + GetDim()$Nfeature):(3 + GetDim()$Nfeature + 1), 2]
  StoreTime <- predTime
  if(is.null(predTime)) {
    predTime <- GetData()$time  
  } 
  if(length(predTime) == 1) {
    tempTime <- GetData()$time 
    SD <- sd(tempTime)
    rangeTemp <- c(max(c(0, Time - 2*SD)), Time + 2*SD)
    timex <- seq(rangeTemp[1], rangeTemp[2], length.out = 100)
  } else {
    timex <- seq(min(preTime), max(preTime), length.out = 100)
  }
  timexNull <- seq(min(GetData()$time), max(GetData()$time), length.out = 100)
  zNull <- colMeans(GetData()$features)
  CIF11Null <- lapply(timexNull, \(timeVal) F_cdf_rcpp(Params = Par11, Z = zNull, x = timeVal)) |> unlist()
  CIF12Null <- lapply(timexNull, \(timeVal) F_cdf_rcpp(Params = Par12, Z = zNull, x = timeVal)) |> unlist()
  CIF21Null <- lapply(timexNull, \(timeVal) F_cdf_rcpp2(Params = Par21, Z = zNull, x = timeVal)) |> unlist()
  CIF22Null <- lapply(timexNull, \(timeVal) F_cdf_rcpp2(Params = Par22, Z = zNull, x = timeVal)) |> unlist()
  CIF31Null <- lapply(timexNull, \(timeVal) F_cdf_rcpp3(Params = Par31, Z = zNull, x = timeVal)) |> unlist()
  CIF32Null <- lapply(timexNull, \(timeVal) F_cdf_rcpp3(Params = Par32, Z = zNull, x = timeVal)) |> unlist()

  CIF11Fig <- lapply(timex, \(timeVal) F_cdf_rcpp(Params = Par11, Z = z, x = timeVal)) |> unlist()
  CIF12Fig <- lapply(timex, \(timeVal) F_cdf_rcpp(Params = Par12, Z = z, x = timeVal)) |> unlist()
  CIF21Fig <- lapply(timex, \(timeVal) F_cdf_rcpp2(Params = Par21, Z = z, x = timeVal)) |> unlist()
  CIF22Fig <- lapply(timex, \(timeVal) F_cdf_rcpp2(Params = Par22, Z = z, x = timeVal)) |> unlist()
  CIF31Fig <- lapply(timex, \(timeVal) F_cdf_rcpp3(Params = Par31, Z = z, x = timeVal)) |> unlist()
  CIF32Fig <- lapply(timex, \(timeVal) F_cdf_rcpp3(Params = Par32, Z = z, x = timeVal)) |> unlist()

  CIF11Val <- lapply(predTime, \(timeVal) F_cdf_rcpp(Params = Par11, Z = z, x = timeVal)) |> unlist()
  CIF12Val <- lapply(predTime, \(timeVal) F_cdf_rcpp(Params = Par12, Z = z, x = timeVal)) |> unlist()
  CIF21Val <- lapply(predTime, \(timeVal) F_cdf_rcpp2(Params = Par21, Z = z, x = timeVal)) |> unlist()
  CIF22Val <- lapply(predTime, \(timeVal) F_cdf_rcpp2(Params = Par22, Z = z, x = timeVal)) |> unlist()
  CIF31Val <- lapply(predTime, \(timeVal) F_cdf_rcpp3(Params = Par31, Z = z, x = timeVal)) |> unlist()
  CIF32Val <- lapply(predTime, \(timeVal) F_cdf_rcpp3(Params = Par32, Z = z, x = timeVal)) |> unlist()

  CIFnull <- data.frame(
    Model = rep(c("GOM", "POM", "PHM"), each = 2*length(timexNull)),
    CIF = c(CIF11Null, CIF12Null, CIF21Null, CIF22Null, CIF31Null, CIF32Null),
    Time = rep(timexNull, 6),
    Risk = rep(rep(RiskNames, each = length(timexNull)), 3)
  )

  CIFAdjustedFig <- data.frame(
    Model = rep(c("GOM", "POM", "PHM"), each = 2*length(timexNull)),
    CIFAdjusted = c(CIF11Fig, CIF12Fig, CIF21Fig, CIF22Fig, CIF31Fig, CIF32Fig), 
    Time = rep(timex, 6),
    Risk = rep(rep(RiskNames, each = length(timexNull)), 3)
    
  )

  CIFAdjustedVal <- data.frame(
    Model = rep(c("GOM", "POM", "PHM"), each = 2*length(predTime)),
    Time = rep(predTime, 6),
    CIFAdjusted = c(CIF11Val, CIF12Val, CIF21Val, CIF22Val, CIF31Val, CIF32Val),
    Risk = rep(rep(RiskNames, each = length(time)), 3) 
  )
  CIFnull <- CIFnull |> within(Model <- factor(Model, levels = c("GOM", "POM", "PHM")))
  CIFnull <- CIFnull |> within(Risk <- factor(Risk, levels = c(RiskNames[1], RiskNames[2])))
  CIFAdjustedFig <- CIFAdjustedFig |> within(Model <- factor(Model, levels = c("GOM", "POM", "PHM")))
  CIFAdjustedFig <- CIFAdjustedFig |> within(Risk <- factor(Risk, levels = c(RiskNames[1], RiskNames[2])))
  CIFAdjustedVal <- CIFAdjustedVal |> within(Model <- factor(Model, levels = c("GOM", "POM", "PHM")))
  CIFAdjustedVal <- CIFAdjustedVal |> within(Risk <- factor(Risk, levels = c(RiskNames[1], RiskNames[2])))

Plot_Adjusted <- CIFAdjustedFig |> 
    ggplot2::ggplot(ggplot2::aes(x = Time, y = CIFAdjusted, group = Risk, color = Risk)) + 
    ggplot2::geom_line(linewidth = 1) + 
    ggplot2::ylim(c(0, 1)) + 
    ggplot2::theme_bw() + 
    ggplot2::facet_wrap(~Model, scales = "fixed") + 
    ggplot2::labs(title = "Cumulative Incidence Function (CIF) for Competing Risks", 
    caption = "All Models | Adjusted by covariates")

Plot_NULL <- CIFnull |> 
    ggplot2::ggplot(ggplot2::aes(x = Time, y = CIF, group = Risk, color = Risk)) + 
    ggplot2::geom_line(linewidth = 1) + 
    ggplot2::ylim(c(0, 1)) + 
    ggplot2::theme_bw() + 
    ggplot2::facet_wrap(~Model, scales = "fixed") + 
    ggplot2::labs(title = "Cumulative Incidence Function (CIF) for Competing Risks", 
    caption = "All Models | Not Adjusted")

Plot_GOM <- CIFAdjustedFig |> 
    subset(subset = Model == "GOM") |> 
    ggplot2::ggplot(ggplot2::aes(x = Time, y = CIFAdjusted, group = Risk, color = Risk)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::theme_bw() + 
    ggplot2::labs(title = "Cumulative Incidence Function (CIF) for Competing Risks | GOM Model",
    caption = "Adjusted by covariates | GOM Model")

Plot_POM <- CIFAdjustedFig |> 
    subset(subset = Model == "POM") |> 
    ggplot2::ggplot(ggplot2::aes(x = Time, y = CIFAdjusted, group = Risk, color = Risk)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::theme_bw() + 
    ggplot2::labs(title = "Cumulative Incidence Function (CIF) for Competing Risks | POM Model",
    caption = "Adjusted by covariates | POM Model")

Plot_PHM <- CIFAdjustedFig |> 
    subset(subset = Model == "PHM") |> 
    ggplot2::ggplot(ggplot2::aes(x = Time, y = CIFAdjusted, group = Risk, color = Risk)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::theme_bw() + 
    ggplot2::labs(title = "Cumulative Incidence Function (CIF) for Competing Risks | PHM Model",
    caption = "Adjusted by covariates | PHM Model")

PlotType = switch(TypeMethod, 
  "GOM" = Plot_GOM, 
  "POM" = Plot_POM,
  "PHM" = Plot_PHM)

  OutPut <- list(Time = list(
  InputTime = StoreTime, 
  TimeForPlotAdjusted = timex, 
  TimeForPlotnull = timexNull
  ), 
  CIF = list(
  CIFNULL = CIFnull |> subset(subset = Model == TypeMethod) |> 
            subset(select = -c(Model)), 
  CIFAdjusted = CIFAdjustedVal |> 
      subset(subset = Model == TypeMethod) |> 
        subset(select = -c(Model))
  ),
  Plot = list(
    PlotNull_AllModels = Plot_NULL,
    PlotAdjusted_AllModels = Plot_Adjusted,
    Plot_InputModel = PlotType
  )
  )

  return(OutPut)
}


FineGray_Model()
cmpp_CIF()
