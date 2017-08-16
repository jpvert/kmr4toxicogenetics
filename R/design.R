#' The design matrix
#' 
#' The design matrix processed from the Dream 8 toxicogenetics challenge dataset.
#' 
#' @format A data frame with 884 cell lines in rows and 1,237 variables in columns of three categories:
#' \describe{
#'   \item{covariate}{16 variables binarized from 3 attributes (sex, population, batch) characterizing the nature of cell lines provided by the challenge.}
#'   \item{RNAgram}{337 variables corresponding to RNA Gram matrix (linear kernel) of normalized and NA-fixed RNA-seq counts of cell lines.}
#'   \item{SNPgram}{884 variables corresponding to genotypic SNP distance matrix (squared Euclidean distance) of SNP data of cell lines.}
#' }
#' 
#' @note RNA-seq data are available only for 337 cells from the challenge.
#' 
#' @source \url{https://doi.org/10.7303/syn1761567}
#' 

"design"