#' tooth24
#'
#'
#' The data set is from the Signal Tandmobiel study, a prospective oral health study conducted in Belgium from 1996 to 2001.
#' 
#' @docType data
#'
#' @usage data(tooth24)
#'
#' @format A data frame with 5 columns and 4386 rows. The columns are:
#' \itemize{
#' \item{ID: }{Subject’s identification code.}
#' \item{LEFT: }{Lower limit of tooth emergence in years since age 5.}
#' \item{RIGHT: }{Upper limit of tooth emergence in years since age 5.}
#' \item{SEX: }{gender of the child; 0 = boy; 1 = girl.}
#' \item{DMF: }{Status of primary predecessor; 0 = sound; 1 = decayed, missing, or filled}
#' }
#' @details The study contains a cohort of randomly sampled schoolchildren who attended the first year of the primary school at the beginning of the study. The original dataset has been presented and studied by Bogaerts and Lesaffre (2004) and Gomez et al. (2009). Here we restrict the analysis to the age of the emergence time of the permanent upper left first premolars (tooth 24 in European dental notation). Since permanent teeth do not emerge before the age of 5, the origin time for all analyses is set at 5 years.      
#'
#' @keywords datasets
#'
#' @references Bogaerts, K. and E. Lesaffre (2004). A new, fast algorithm to find the regions of possible support for bivariate interval-censored data. Journal of Computational and Graphical Statistics 13(2), 330 – 340.
#' 
#' @references Gomez, G., M. L. Calle, R. Oller, and K. Langohr (2009). Tutorial on methods for interval-censored data and their implementation in R. Statistical Modelling 9(4), 259–297.
#'
#'
#' @examples
#' data(tooth24)
"tooth24"
