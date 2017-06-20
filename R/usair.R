#' US Air Pollution
#'
#'
#' The data were collected to investigate the determinants of pollution by considering SO2
#' level as the dependent variable and the remaining variables being potential explanatory variables.
#'
#' @docType data
#'
#' @usage data(usair)
#'
#' @format A data frame with 7 columns and 41 rows. The columns are:
#' \itemize{
#' \item{SO2}{Sulphur dioxide content of air in micrograms per cubic meter.}
#' \item{negtemp}{Negative value of Average annual temperature in fahrenheit.}
#' \item{manuf}{Number of manufacturing enterprises employing 20 or more workers.}
#' \item{pop}{Population size (1970 census) in thousands.}
#' \item{wind}{Average annual wind speed in miles per hour.}
#' \item{precip}{Average annual precipitation in inches.}
#' \item{days}{Average number of days with precipitation per year}
#' }
#' @details Note that Neg.Temp represents the negative value of Average annual temperature.
#'
#' @keywords datasets
#'
#' @references Hand, D. J., Daly, F., Lunn, A. D., McConway, K. J. and Ostrowski, E. (1994),
#' A handbook of small data sets, Chapman and Hall, London.
#'
#' Brian S Everitt. An R and S-PLUS companion to multivariate analysis.
#' Springer Science & Business Media, 2006.
#'
#' @source Sokal and Rohlf (1981) Biometry: The Principles and Practices of Statistics in Biological Research
#'
#' @examples
#' data(usair)
"usair"
