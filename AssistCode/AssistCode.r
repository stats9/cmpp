packageVersion("cmpp")
remove.packages("cmpp")
devtools::clean_dll()
devtools::document()              # DESCRIPTION + NAMESPACE آپدیت شه
Rcpp::compileAttributes()         # RcppExports آپدیت شه
system("R CMD build .")           # ساخت tar.gz
system("R CMD INSTALL --build cmpp_0.0.1.tar.gz")  # ساخت zip

roxygen2::roxygenize(clean = TRUE)
devtools :: check()
devtools :: install()
devtools :: build(binary = TRUE)
devtools :: build()
?devtools::build
install.packages("./downloads/cmpp_0.0.1.zip", repos = NULL, type = "win-binary")
install.packages("./downloads/cmpp_0.0.1.tar.gz", repos = NULL, type = "source")
# Install from GitHub
devtools::install_github("stats9/cmpp")
library(cmpp)
data(dat)
names(dat)
feat <- dat[, -c(match(c('id', 'time', 'event', 'cause_burn', 'd1', 'd2', 'cause_hotObject3'), names(dat)))]  
timee <- dat[['time']]
d1 = dat[['d1']]
d2 = dat[['d2']]
feat2 = feat |> data.matrix()
Initialize(feat2, timee, d1, d2, 1e-10)

FineGray_Model()
(Res <- Cmpp_CIF())
Res$Plot$PlotNull_AllModels
Res$Plot$Plot_InputModel
Res$Plot$PlotAdjusted_AllModels

Res2 <- Cmpp_CIF(predTime = c(1, 2, 3), TypeMethod = "POM")
Res2$CIF
Res2$Plot$PlotAdjusted_AllModels
GetData()
Res2$Plot$PlotNull_AllModels
Res2$Plot$Plot_InputModel