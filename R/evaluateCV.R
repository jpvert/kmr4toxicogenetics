#' Performance evaluation by cross-validation
#' 
#' Evaluate regression performance of a predictor by cross-validation.
#' 
#' @param mypredictor Character indicating a predictor function. Possible options are KMR (default), ElasticNet, Lasso and RF.
#' @param celllines Matrix of descriptors for \code{ncell} cell lines, of dimension \code{ncell x pcell}. Only used for ElasticNet, Lasso and RF.
#' @param celllinesKernel Kernel Gram matrix for \code{ncell} cell lines, of dimension \code{ncell x ncell}. Only used for KMR.
#' @param chemicals Matrix of descriptors for \code{nchem} chemicals, of dimension \code{nchem x pchem}. Only used for ElasticNet, Lasso and RF.
#' @param chemicalsKernel Kernel Gram matrix for \code{nchem} chemicals, of dimension \code{nchem x nchem}. Only used for KMR.
#' @param toxicity Matrix of toxicity values for \code{ncell} cell lines responding to \code{nchem} chemicals, of dimension \code{ncell x nchem}.
#' @param nfolds Number of folds for cross-validation. Default is 5.
#' @param nrepeats Number of times the k-fold cross-validation is performed. Default is 1.
#' @param seed A seed number for the random number generator (useful to have the same CV splits).
#' @param mc.cores Number of parallelable CPU cores to use.
#' @param ... Other arguments to pass to predictor function.
#' 
#' @return A list with matrices of cross-validation performance scores.
#' Each score matrix is of dimension \code{nexp x nchem} (per CV experiment, per chemical) where \code{nexp=nfolds*nrepeats}
#' and corresponds to one of the evaluation criteria:
#' \item{matrix.ci}{Concordance index.}
#' \item{matrix.rho}{Pearson correlation.}
#' 
#' @export
#' 
#' @seealso \code{\link{predictorKMR}}, \code{\link{predictorElasticNet}}, \code{\link{predictorLasso}}, \code{\link{predictorRF}}
#' 
#' @references 
#' Bernard, E., Jiao, Y., Scornet, E., Stoven, V., Walter, T., and Vert, J.-P. (2017). "Kernel multitask regression for toxicogenetics." \href{https://doi.org/10.1101/171298}{bioRxiv-171298}.
#' 
#' @importFrom parallel mclapply
#' @importFrom Hmisc rcorr.cens
#' @importFrom stats cor
#' 

evaluateCV <- function(mypredictor = c("predictorKMR", 
                                       "predictorElasticNet", 
                                       "predictorLasso", 
                                       "predictorRF"), 
                       celllines, 
                       celllinesKernel, 
                       chemicals, 
                       chemicalsKernel, 
                       toxicity, 
                       nfolds = 5, 
                       nrepeats = 10, 
                       seed = 47, 
                       mc.cores = 1, 
                       ...)
{
  mypredictor <- match.arg(mypredictor)
  
  # Set random number generator seed
  set.seed(seed)
  
  # Number of cell lines
  ncelllines <- if (!missing(celllines)) dim(celllines)[1] else dim(celllinesKernel)[1]
  
  # Number of chemicals
  nchemicals <- ncol(toxicity)
  
  # Make folds
  n <- ncelllines
  folds <- list()
  for (i in seq(nrepeats)) {
    folds <- c(folds, 
               split(sample(seq(n)), rep(1:nfolds, length = n)))
  }
  nexp <- length(folds)
  
  ## Main CV loop (parallelized):
  resCV <- parallel::mclapply(seq(nexp), function(iexp){
    # Set ytrain and ytest
    toxicityTrain <- toxicity[-folds[[iexp]],]
    toxicityTest <- toxicity[folds[[iexp]],]
    
    # Set xtrain and xtest
    if (mypredictor == "predictorKMR") {
      celllinesKernelTrain <- celllinesKernel[-folds[[iexp]], -folds[[iexp]], drop=F]
      celllinesKernelTest <- celllinesKernel[folds[[iexp]], -folds[[iexp]], drop=F]
      toxicityPred <- get(mypredictor, mode = "function")(celllinesKernelTrain, celllinesKernelTest, toxicityTrain, chemicalsKernel, ...)
    } else {
      celllinesTrain <- celllines[-folds[[iexp]], , drop=F]
      celllinesTest <- celllines[folds[[iexp]], , drop=F]
      toxicityPred <- get(mypredictor, mode = "function")(celllinesTrain, celllinesTest, toxicityTrain, ...)
    }
    
    # Evaluation criterion: c-index (per chemical)
    toxicityCI <- apply(rbind(toxicityPred, toxicityTest), 2, function(u){
      Hmisc::rcorr.cens(u[1:nrow(toxicityPred)], u[-(1:nrow(toxicityPred))], outx = FALSE)[[1]]
    })
    # Evaluation criterion: pearson correlation (per chemical)
    toxicityRho <- sapply(1:ncol(toxicityTest), function(i){
      stats::cor(toxicityTest[ ,i], toxicityPred[ ,i])
    })
    return(list(ci = toxicityCI, 
                rho = toxicityRho))
  }, mc.cores = mc.cores)
  
  matrix.ci <- do.call("rbind", lapply(resCV, function(rr) rr$ci))
  matrix.rho <- do.call("rbind", lapply(resCV, function(rr) rr$rho))
  colnames(matrix.ci) <- colnames(toxicity)
  colnames(matrix.rho) <- colnames(toxicity)
  
  return(list(matrix.ci = matrix.ci, 
              matrix.rho = matrix.rho))
}
