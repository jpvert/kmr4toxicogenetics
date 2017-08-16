#' Kernel matrices for cell lines
#' 
#' A list of example kernel matrices of cell line features.
#' 
#' @format A list of 24 kernel matrices for cell lines of three categories of cell line features:
#' \describe{
#'   \item{Kcovariates}{4 kernel matrices of cell line covariates provided by the challenge, three of which correspond to a linear kernel of each attribute (sex, population, batch) plus one more combining all three attributes.}
#'   \item{KrnaseqRbf}{10 kernel matrices of normalized and NA-fixed RNA-seq counts of cell lines, corresponding to Gaussian RBF kernel with various bandwidth.}
#'   \item{KsnpRbf}{10 kernel matrices of genotypic SNP data of cell lines, corresponding to Gaussian RBF kernel with various bandwidth.}
#' }
#' 
#' @references 
#' Bernard, E., Jiao, Y., Scornet, E., Stoven, V., Walter, T., and Vert, J.-P. (2017). "Kernel multitask regression for toxicogenetics." \href{https://doi.org/10.1101/171298}{bioRxiv-171298}.
#' 
#' @source \url{https://doi.org/10.7303/syn1761567}
#' 

"kcell"