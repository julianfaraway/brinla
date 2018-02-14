#' joint
#'
#'
#' A list contents both longitudinal and survival data for the AIDS clinical trial.
#' 
#' @docType data
#'
#' @usage data(joint)
#'
#' @format 
#' The longitudinal data set is a data frame with 7 columns and 2335 rows. The columns are:
#' \itemize{
#' \item{y: }{The square root of the CD4.}
#' \item{TIME: }{Time that CD4 counts were recorded.}
#' \item{TIMEDRUG: }{Time and drug interaction.}
#' \item{SEX: }{Sex; 1 = male; -1 = female.}
#' \item{PREVOI: }{Previous opportunistic infection (AIDS diagnosis) at study entry; 1 = yes; -1 = no}
#' \item{STRATUM: }{Failure or intolerance of zidovudine (AZT) therapy; 1 = failure; -1 = intolerance}
#' \item{ID: }{patient identification}
#' }
#' The survival data set is a data frame with 6 columns and 467 rows. The columns are:
#' \itemize{
#' \item{DRUG: }{receive either didanosine (ddI) or zalcitabine (ddC); 1 = ddI; 0 = ddC.}
#' \item{SEX: }{Sex; 1 = male; -1 = female.}
#' \item{PREVOI: }{Previous opportunistic infection (AIDS diagnosis) at study entry; 1 = yes; -1 = no.}
#' \item{STRATUM: }{Failure or intolerance of zidovudine (AZT) therapy; 1 = failure; -1 = intolerance.}
#' \item{SURVTIME: }{Time to death in months}
#' \item{CENSOR: }{Death status; 1 = death; 0 = censored}
#' }
#' 
#' @details In this AIDS study, both longitudinal and survival data were collected to compare the efficacy and safety of two antiretroviral drugs in treating patients who had failed or were intolerant of zidovudine (AZT) therapy. There were 467 HIV-infected patients who met entry con- ditions (either an AIDS diagnosis or two CD4 counts of 300 or fewer, and fulfilling specific criteria for AZT intolerance or failure). The patients were randomly assigned to receive either didanosine (ddI) or zalcitabine (ddC). CD4 counts were recorded at study entry, and again at the 2-, 6-, 12-, and 18-month visits. The times to death were also recorded.      
#'
#' @keywords datasets
#'
#' @references Guo, X. and B. P. Carlin (2004). Separate and joint modeling of longitudinal and event time data using standard computer packages. The American Statisti- cian 58(1), 16–24.
#' 
#'
#' @examples
#' data(joint)
#' str(joint$longitudinal)
#' str(joint$survival)
"joint"
