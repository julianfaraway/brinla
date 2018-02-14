#' lowbwt
#'
#'
#' Low birth weight data that has been presented in Hosmer and Lemeshow (2004). The dataset contains information on 189 births to women seen in the obstetric clinic, where data were collected as part of a larger study at Baystate Medical Center in Springfield, Massachusetts.  
#'
#' @docType data
#'
#' @usage data(lowbwt)
#'
#' @format A data frame with 8 columns and 189 rows. The columns are:
#' \itemize{
#' \item{LOW: }{Indicator of low birth weight; 0 = "≥ 2500g"; 1 = "< 2500g"}
#' \item{AGE: }{Age of mother.}
#' \item{LWT: }{Weight of mother at last menstrual period.}
#' \item{RACE: }{Race of mother; 1 = white; 2 = black; 3 = other.}
#' \item{SMOKE: }{Smoking status during pregnancy; 0 = no; 1 = yes.}
#' \item{HT: }{History of hypertension; 0 = no; 1 = yes.}
#' \item{UI: }{Presence of uterine irritability; 0 = no; 1 = yes.}
#' \item{FTV: }{Number of physician visits during the first trimester.}
#' }
#' @details  The response variable is LOW indicating birth weight less than 2500 grams, which has been of concern to physicians for years.    
#'
#' @keywords datasets
#'
#' @references Hosmer, D. W. and S. Lemeshow (2004). Applied Logistic Regression. New York: John Wiley & Sons.
#'
#'
#' @examples
#' data(lowbwt)
"lowbwt"
