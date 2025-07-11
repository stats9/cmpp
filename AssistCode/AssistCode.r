packageVersion("cmpp")
remove.packages("cmpp", lib = .libPaths()[1])
remove.packages("cmpp", lib = .libPaths()[2])
devtools::clean_dll()
devtools::document()              # DESCRIPTION + NAMESPACE آپدیت شه
roxygen2::roxygenize(clean = TRUE)
unlink("src/RcppExports.*")
Rcpp::compileAttributes()

system("R CMD build .")           # ساخت tar.gz
system("R CMD INSTALL --build cmpp_0.0.3.tar.gz")  # ساخت zip


devtools :: check()
devtools :: install()
devtools :: build(binary = TRUE)
devtools :: build()

install.packages("./downloads/cmpp_0.0.2.zip", repos = NULL, type = "win-binary")
install.packages("./downloads/cmpp_0.0.2.tar.gz", repos = NULL, type = "source")
# Install from GitHub
devtools::install_github("stats9/cmpp")
## create an example ############

library(cmpp)
data("fertility_data")
help(package = "cmpp")
(Nam <- names(fertility_data))
fertility_data$Education
datt <- make_Dummy(fertility_data, features = c("Education"))
datt <- datt$New_Data 
datt['Primary_Secondary'] <- datt$`Education:2`
datt['Higher_Education'] <- datt$`Education:3`
datt$`Education:2` <- datt$`Education:3` <- NULL
datt2 <- make_Dummy(datt, features = 'Event')$New_Data
d1 <- datt2$`Event:2`
d2 <- datt2$`Event:3`
feat <- datt2[c('age', 'Primary_Secondary', 'Higher_Education')] |> 
    data.matrix()
timee <- datt2[['time']]
Initialize(feat, timee, d1, d2, 1e-10)

FineGray_Model(, )
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

################### 2025-04-27 #####################

# library(haven)
# dat <- read_sav(file.choose())
# saveRDS(dat, file = "./data/Pregnancy.rds")
# dat <- readRDS("./data/Pregnancy.rds")
load("./data/fertility_data.rda")
ls()
dat <- fertility_data
dat |> names()
dat |> head()
dat$Education
# dat$edu |> table()
# edu <- dat$edu
# edu
# library(tidyverse)
# edu2 <- as.numeric(edu)
# edu2 |> table()
# edu2
# edu3 <- case_match(edu2,
#   1 ~ "1",
#   2 ~ "2",
#   3 ~ "2",
#   4 ~ "3",
#   .default = NA_character_
# ) 

# # attributes(edu3)
# edu4 <- edu3 |> factor(levels = 1:3) 
# levels(edu4)
# attributes(edu4) <- list(
#     label = "Education Levels",
#     levels = c("1", "2", "3"),
#     class = "factor", 
#     labels = c("Illiterate" = "1", "Primary/Secondary" = "2", "Higher Education" = "3")
# )
# dat |> head()
# dat2 <- within(dat, Education <- edu4)
# dat2$edu <- NULL
# dat2 |> head()
# event1 <- dat2$event
# event2 <- factor(event1, levels = 0:2)
# attributes(event2) <- list(
#     label = "Event",
#     levels = c("0", "1", "2"),
#     class = "factor", 
#     labels = c("Censored" = "0", "Live Birth" = "1", "Stillbirth" = "2")
# )
# dat3 <- within(dat2, Event <- event2)
# dat3$event <- NULL
# dat3 |> head()
# # save(dat3, file = "./data/Pregnancy.RData")
# Pregnancy <- dat3
# usethis::use_data(Pregnancy, overwrite = TRUE)
# Pregnancy |> head()
# fertility_data <- Pregnancy
# save(fertility_data, file = "./data/fertility_data.rda")

## Important Example
library(cmpp)
data("fertility_data")
Nam <- names(fertility_data)
fertility_data$Education
datt <- make_Dummy(fertility_data, features = c("Education"))
datt <- datt$New_Data 
datt['Primary_Secondary'] <- datt$`Education:2`
datt['Higher_Education'] <- datt$`Education:3`
datt$`Education:2` <- datt$`Education:3` <- NULL
datt2 <- make_Dummy(datt, features = 'Event')$New_Data
d1 <- datt2$`Event:2`
d2 <- datt2$`Event:3`
feat <- datt2[c('age', 'Primary_Secondary', 'Higher_Education')] |> 
   data.matrix()
timee <- datt2[['time']]
Initialize(feat, timee, d1, d2, 1e-10)
result <- Cmpp_CIF(
  featureID = c(1, 2),
  featureValue = c(0.5, 1.2),
  RiskNames = c("Event1", "Event2"),
  TypeMethod = "GOR",
  predTime = seq(0, 10, by = 0.5)
)
print(result$Plot$Plot_InputModel)  # Plot for the specified model
print(result$Plot$PlotAdjusted_AllModels)  # Adjusted CIFs for all models
print(result$CIF$CIFAdjusted)  # Adjusted CIF values
