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
