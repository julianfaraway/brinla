#' Meat spectrometry to determine fat content
#' 
#' A Tecator Infratec Food and Feed Analyzer working in the wavelength range
#' 850 - 1050 nm by the Near Infrared Transmission (NIT) principle was used to
#' collect data on samples of finely chopped pure meat. 215 samples were
#' measured. For each sample, the fat content was measured along with a 100
#' channel spectrum of absorbances. Since determining the fat content via
#' analytical chemistry is time consuming we would like to build a model to
#' predict the fat content of new samples using the 100 absorbances which can
#' be measured more easily.
#' 
#' 
#' @format Dataset contains the following variables \describe{
#' \item{list("V1-V100")}{ absorbances across a range of 100 wavelengths }
#' \item{list("fat")}{ fat content} }
#' @source H. H. Thodberg (1993) "Ace of Bayes: Application of Neural Networks
#' With Pruning", report no. 1132E, Maglegaardvej 2, DK-4000 Roskilde, Danmark
#' @keywords datasets
"meatspec"