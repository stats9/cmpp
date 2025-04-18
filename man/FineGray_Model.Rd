% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCode.r
\name{FineGray_Model}
\alias{FineGray_Model}
\title{Fine-Gray Model for Competing Risks Data}
\usage{
FineGray_Model(
  featureNames = NULL,
  CovarNames = NULL,
  Failcode = 1,
  RiskNames = NULL
)
}
\arguments{
\item{featureNames}{A character vector of feature names for the covariates. If \code{NULL}, default names will be generated.}

\item{CovarNames}{A character vector of names for the covariates. If \code{NULL}, default names will be generated.}

\item{Failcode}{An integer specifying the event of interest (default is \code{1}).}

\item{RiskNames}{A character vector specifying the names of the competing risks. If \code{NULL}, default names ("Risk1" and "Risk2") will be used.}
}
\value{
A list containing:
\item{Results}{A summary of the Fine-Gray model fit.}
\item{Plot}{A ggplot object showing the cumulative incidence functions (CIFs) for the competing risks.}
\item{CIF_Results}{A data frame containing the CIFs for the competing risks, along with their corresponding time points.}
}
\description{
This function fits a Fine-Gray model for competing risks data using the \code{cmprsk} package.
It estimates the subdistribution hazard model parameters, computes cumulative incidence functions (CIFs),
and provides a summary of the results along with a plot of the CIFs.
}
\details{
This function retrieves the data initialized in the Cmpp model using the \code{GetData} function.
It uses the \code{crr} function from the \code{cmprsk} package to fit the Fine-Gray model for competing risks.
The function also computes cumulative incidence functions (CIFs) using the \code{cuminc} function and
generates a plot of the CIFs for the competing risks.
}
\examples{
\dontrun{
# Assuming the Cmpp model has been initialized
result <- FineGray_Model(
  featureNames = c("Feature1", "Feature2"),
  CovarNames = c("Covar1", "Covar2"),
  Failcode = 1,
  RiskNames = c("Event1", "Event2")
)
print(result$Results)  # Summary of the Fine-Gray model
print(result$Plot)     # Plot of the CIFs
print(result$CIF_Results)  # CIF data
}
}
