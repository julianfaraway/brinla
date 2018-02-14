#' larynx
#'
#'
#' The data set includes 90 males diagnosed with cancer of the larynx.
#' 
#' @docType data
#'
#' @usage data(larynx)
#'
#' @format A data frame with 5 columns and 90 rows. The columns are:
#' \itemize{
#' \item{stage: }{Stage of disease; 1 = stage1; 2 = stage2; 3 = stage3, 4 = stage4.}
#' \item{time: }{Time to death or on-study time (in months).}
#' \item{age: }{Age at diagnosis of larynx cancer (years).}
#' \item{diagyr: }{Year of diagnosis of larynx cancer (years).}
#' \item{delta: }{Death status; 0=alive, 1=dead.}
#' }
#' @details The dataset was first reported by Kardaun (1983), including 90 males diagnosed with caner of the larynx during the period 1970 – 1978 at a Dutch hospital. Patient survival times were recorded between first treatment and either death or the end of the study. Patients were classified into one of four stages using the American Joint Committee for Cancer Staging. Other variables also recorded include the patient’s age at the time of diagnosis, and the year of diagnosis.      
#'
#' @keywords datasets
#'
#' @references Kardaun, O. (1983). Statistical survival analysis of male larynx-cancer patients: A case study. Statistica Neerlandica 37(3), 103–125.
#' @references Klein, J. P. and M. L. Moeschberger (2005). Survival Analysis: Techniques for Censored and Truncated Data. New York: Springer.
#'
#'
#' @examples
#' data(larynx)
"larynx"
