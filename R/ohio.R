#' Ohio Children Wheeze Status
#'
#' The dataset is a subset of the six-city study, a longitudinal study of the health effects of air pollution.
#'
#' @docType data
#'
#' @usage data(ohio)
#'
#' @format The ohio data frame has 2148 rows and 4 columns
#' \itemize{
#' \item{resp}{an indicator of wheeze status (1=yes, 0=no).}
#' \item{id}{a numeric vector for subject id.}
#' \item{age}{a numeric vector of age, 0 is 9 years old}
#' \item{smoke}{an indicator of maternal smoking at the first year of the study}
#'}
#' @keywords datasets
#'
#'
#' @references Fitzmaurice, G.M. and Laird, N.M. (1993) A likelihood-based method for
#' analyzing longitudinal binary responses, Biometrika 80: 141–151.
#'
#' @examples
#' data(ohio)
#'
"ohio"
