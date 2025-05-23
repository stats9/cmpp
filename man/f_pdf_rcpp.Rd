% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCode.r
\name{f_pdf_rcpp}
\alias{f_pdf_rcpp}
\title{Compute the PDF of the Parametric Generalized odds rate (GOR)}
\arguments{
\item{Params}{A numeric vector of parameters.}

\item{Z}{A numeric vector of covariates.}

\item{x}{A numeric value representing the time point.}
}
\value{
A numeric value representing the PDF.
}
\description{
This function computes the probability density function (PDF) of the parametric model (GOR Approach).
}
\examples{
library(cmpp)
features <- matrix(rnorm(300, 1, 2), nrow = 100, ncol = 3)
delta1 <- sample(c(0, 1), 100, replace = TRUE)
delta2 <- 1 - delta1
x <- rexp(100, rate = 1/10)
Initialize(features, x, delta1, delta2, h = 1e-5)
params <- rep(0.001, (ncol(features) + 3))
pdf_value <- f_pdf_rcpp(params, features[1, ], x[3])
print(pdf_value)

}
