library(cmpp)
data(dat)
dat
feat <- dat[, -c(match(c('id', 'time', 'event', 'cause_burn', 'd1', 'd2', 'cause_hotObject3'), names(dat)))]  
timee <- dat[['time']]
d1 = dat[['d1']]
d2 = dat[['d2']]
feat2 = feat |> data.matrix()
cpp_Initialize(feat2, timee, d1, d2, 1e-5)
optim(par = c(0.01, 0.01, 0.01, 0.01), 
    fn = cpp_LogLike1, 
    gr = cpp_compute_grad, 
    method = 'BFGS')
