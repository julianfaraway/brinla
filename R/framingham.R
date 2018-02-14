#' framingham
#'
#'
#' The data set is from the Framingham Heart Study (Carroll et al., 2006).
#' 
#' @docType data
#'
#' @usage data(framingham)
#'
#' @format A data frame with 7 columns and 1615 rows. The columns are:
#' \itemize{
#' \item{AGE: }{Age at exam 2.}
#' \item{SBP21: }{First systolic blood pressure at exam 2.}
#' \item{SBP22: }{Second systolic blood pressure at exam 2.}
#' \item{SBP31: }{First systolic blood pressure at exam 3.}
#' \item{SBP32: }{Second systolic blood pressure at exam 3.}
#' \item{SMOKE: }{present smoking at exam 1; 0 = No; 1 = Yes.}
#' \item{FIRSTCHD: }{Indicator of first evidence of CHD occurring at exam 3 through 6; 0 = No; 1 = Yes.}
#' }
#' @details The Framingham study consists of a series of exams taken two years apart. There are 1615 males aged from 31 to 65 in this dataset, with the outcome indicating the occurrence of coronary heart disease (CHD) within an eight-year period following Exam 3. There were 128 total cases of CHD.      
#'
#' @keywords datasets
#'
#' @references Carroll, R. J., D. Ruppert, L. A. Stefanski, and C. Crainiceanu (2006). Measurement Error in Nonlinear Models: A Modern Perspective (2nd ed.). New York: Chapman & Hall/CRC Press.
#'
#'
#' @examples
#' data(framingham)
"framingham"
