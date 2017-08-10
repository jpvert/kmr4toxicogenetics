#' Wrapper function for lasso regression
#' 
#' Wrapper function to perform lasso regression with \code{cv.glmnet} that trains a model on training set and then predicts on test set for multiple tasks.
#' 
#' @param patientsTrain Matrix of training descriptors, of dimension \code{n x p}, for \code{n} training patients with \code{p} descriptors.
#' @param patientsTest Matrix of test descriptors, of dimension \code{m x p}, for \code{m} test patients with the same set of \code{p} descriptors.
#' @param response Matrix of observed toxicity values, of dimension \code{n x t}, for the \code{n} training patients responding to \code{t} drugs.
#' 
#' @return A matrix of predicted toxicity values, of dimension \code{m x t}, for the \code{m} test patients responding to the \code{t} drugs.
#' 
#' @export
#' 
#' @note Prediction is made per task with no special treatment for multitask learning, nor are task features needed.
#' 
#' @seealso \code{\link[glmnet]{cv.glmnet}}, lasso implements a special case of elastic net \code{\link{predictorElasticNet}}
#' 

predictorLasso <- function(patientsTrain, 
                           patientsTest, 
                           response)
{
  predictorElasticNet(patientsTrain, patientsTest, response, alpha = 1)
}
