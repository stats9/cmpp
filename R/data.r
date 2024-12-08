#' Example Dataset: dat
#'
#' This dataset contains survival data for a study involving burn injuries.
#'
#' @format A data frame with 968 rows and 21 variables:
#' \describe{
#'   \item{id}{Unique identifier for each individual.}
#'   \item{time}{Time to event or censoring.}
#'   \item{event}{Event indicator Variable:, 3 = event occurred (Discharge), 1 = Cardio Arrest (Competing Risk), 2 = Sepsis (Main Event).}
#'   \item{age}{Age of the individual.}
#'   \item{sex}{Sex of the individual (e.g., "Male", "Female").}
#'   \item{percent_burn}{Percentage of body burned.}
#'   \item{cause_burn}{Cause of burn (e.g., "Fire", "Hot Object").}
#'   \item{ca}{Calcium levels.}
#'   \item{k}{Potassium levels.}
#'   \item{mg}{Magnesium levels.}
#'   \item{na}{Sodium levels.}
#'   \item{pho}{Phosphorus levels.}
#'   \item{alb}{Albumin levels.}
#'   \item{bun}{Blood urea nitrogen levels.}
#'   \item{cr}{Creatinine levels.}
#'   \item{wbc}{White blood cell count.}
#'   \item{d1}{A Indicator variable to indicate the cause of death, which takes the value 1 if death is due to Sepsis and the value zero otherwise.}
#'   \item{d2}{A Indicator variable to indicate the cause of death, which takes the value 1 if death was due to a Cardio Arrest, and the value 0 otherwise.}
#'   \item{cause_fire1}{Indicator for burns caused by fire.}
#'   \item{cause_boiling2}{Indicator for burns caused by boiling liquids.}
#'   \item{cause_hotObject3}{Indicator for burns caused by hot objects.}
#' }
#' @source real-world dataset from a study on burn injuries (Shiraz Accident and Burn Hospital, Iran).
"dat"