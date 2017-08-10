#' Wrapper function for kernel multitask regression
#' 
#' Wrapper function to perform kernel multitask regression with \code{cv.kmr} that trains a model on training set and then predicts on test set for multiple tasks.
#' 
#' @param patientsKernelTrain Precomputed kernel Gram matrix of \code{n} training patients, of dimension \code{n x n}.
#' @param patientsKernelTest Precomputed kernel Gram matrix of \code{m} test patients crossing \code{n} training patients, of dimension \code{m x n}.
#' @param response Matrix of observed toxicity values, of dimension \code{n x t}, for the \code{n} training patients responding to \code{t} drugs.
#' @param drugsKernel Kernel Gram matrix of the \code{t} drugs, of dimension \code{t x t}.
#' @param lambdas Sequence of lambdas that must be tested to fit a cross-validated KMR model. Default is exp(-15:25).
#' @param nfolds Number of folds for cross-validation. Default is 5.
#' @param nrepeats Number of times the k-fold cross-validation is performed. Default is 1.
#' 
#' @return A matrix of predicted toxicity values, of dimension \code{m x t}, for the \code{m} test patients responding to the \code{t} drugs.
#' 
#' @export
#' 
#' @note Multitask prediction is made, for which task relationships are encoded in \code{drugsKernel}.
#' 
#' @references 
#' Bernard, E., Jiao, Y., Scornet, E., Stoven, V., Walter, T., and Vert, J.-P. (2017). "Kernel multitask regression for toxicogenetics." \href{https://doi.org/10.1101/171298}{bioRxiv-171298}.
#' 
#' @seealso \code{\link[kmr]{cv.kmr}}
#' 
#' @importFrom kmr cv.kmr
#' 

predictorKMR <- function(patientsKernelTrain, 
                         patientsKernelTest, 
                         response, 
                         drugsKernel, 
                         lambdas = exp(-15:25), 
                         nfolds = 5, 
                         nrepeats = 1)
{
  patientsKernelTrain <- as.matrix(patientsKernelTrain)
  patientsKernelTest <- as.matrix(patientsKernelTest)
  
  # Train a mkr model
  cvobj <- kmr::cv.kmr(x = patientsKernelTrain, y = response, kx_type = "precomputed", 
                       kt_type = "precomputed", kt_option = list(kt = drugsKernel),
                       lambda = lambdas, nfolds = nfolds, nrepeats = nrepeats)
  # Make predictions
  pred <- predict(cvobj, patientsKernelTest)
  
  pred <- data.frame(pred)
  dimnames(pred) <- list(rownames(patientsKernelTest), colnames(response))
  return(pred)
}
