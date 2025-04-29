#' Fertility History of Rural Women in Shiraz
#'
#' This dataset includes fertility history information from a cross-sectional study
#' of 652 women aged 15–49 years in rural areas of Shiraz, Iran. 
#' It was used in the article "A parametric method for cumulative incidence modeling with a new four-parameter log-logistic distribution" 
#' to model the cumulative incidence of live births and competing risks (stillborn fetus or abortion).
#'
#' @format A data frame with 652 rows and the following variables:
#' \describe{
#'   \item{id}{Unique identifier for each case.}
#'   \item{time}{Time from marriage to event (live birth, competing event, or censoring).}
#'   \item{Event}{Event indicator: \code{0} = censored, \code{1} = live birth, \code{2} = stillborn fetus or abortion.}
#'   \item{age}{Age of the woman at the time of the event or censoring.}
#'   \item{Education}{Education level: \code{1} = Illiterate, \code{2} = Primary/Secondary, \code{3} = Higher Education.}
#' }
#'
#' 
#' @source <doi:10.1186/1742-4682-8-43>
#' https://doi.org/10.1186/1742-4682-8-43
#' @references
#' Shayan, Z., Ayatollahi, S. M. T., & Zare, N. (2011). 
#' "A parametric method for cumulative incidence modeling with a new four-parameter log-logistic distribution." 
#' Theoretical Biology and Medical Modelling, 8:43. <doi:10.1186/1742-4682-8-43>
#'
#'@note 
#' To view the article, follow this link: \url{https://tbiomed.biomedcentral.com/articles/10.1186/1742-4682-8-43}
#'
#' @examples
#' data(fertility_data)
#' head(fertility_data)
#' @keywords datasets
"fertility_data"
