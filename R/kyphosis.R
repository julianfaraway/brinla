#' Children who have had Corrective Spinal Surgery
#'
#' Kyphosis is excessive outward curvature of the spine, causing hunching of the back.
#'
#' @docType data
#'
#' @usage data(kyphosis)
#'
#' @format A data frame with 4 columns and 81 rows. The columns are:
#' \itemize{
#' \item{Age}{in months.}
#' \item{NumVert}{The number of vertebrae involved.}
#' \item{StartVert}{the number of the first (topmost) vertebra operated on.}
#' \item{Kyphosis}{takes the value 0 if absent or 1 if present after surgery}}
#'
#' @keywords datasets
#'
#' @references John M. Chambers and Trevor J. Hastie eds. (1992)
#' Statistical Models in S, Wadsworth and Brooks/Cole, Pacific Grove, CA.
#'
#' @source Also found in the \code{rpart} package
#'
#' @examples
#' data(kyphosis)
"kyphosis"
