% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCode.r
\name{GetDim}
\alias{GetDim}
\title{Get Dimensions of the Cmpp Object}
\usage{
GetDim()
}
\value{
A list containing:
\item{Nsamp}{Number of samples (rows in the feature matrix).}
\item{Nfeature}{Number of features (columns in the feature matrix).}
}
\description{
This function returns the number of samples and features stored in the Cmpp object.
It is primarily used to retrieve the dimensions of the data within the class.
}
\details{
The \code{GetDim} function allows users to access the internal dimensions of the
\code{Cmpp} class instance, such as the number of samples (\code{Nsamp}) and the number of features
(\code{Nfeature}). This is useful when working with large datasets, especially for
checking the size of the data without needing to manually access the underlying \code{Eigen::MatrixXd}
or \code{Eigen::VectorXd} objects directly.
}
\examples{
# Initialize Cmpp object
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
# Get dimensions
dims <- GetDim()
dims$Nsamp    # Number of samples
dims$Nfeature # Number of features
}
