#' Wrapper function for elastic net regression
#' 
#' Wrapper function to perform elastic net regression with \code{cv.glmnet} that trains a model on training set and then predicts on test set for multiple tasks.
#' 
#' @param patientsTrain Matrix of training descriptors, of dimension \code{n x p}, for \code{n} training patients with \code{p} descriptors.
#' @param patientsTest Matrix of test descriptors, of dimension \code{m x p}, for \code{m} test patients with the same set of \code{p} descriptors.
#' @param response Matrix of observed toxicity values, of dimension \code{n x t}, for the \code{n} training patients responding to \code{t} drugs.
#' @param alpha The elasticnet mixing parameter. \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty. Default is 0.5. All other arguments are taken by default implementation of \code{randomForest}.
#' 
#' @return A matrix of predicted toxicity values, of dimension \code{m x t}, for the \code{m} test patients responding to the \code{t} drugs.
#' 
#' @export
#' 
#' @note Prediction is made per task with no special treatment for multitask learning, nor are task features needed.
#' 
#' @seealso \code{\link[glmnet]{cv.glmnet}}
#' 
#' @importFrom glmnet cv.glmnet
#' 

predictorElasticNet <- function(patientsTrain, 
                                patientsTest, 
                                response, 
                                alpha = 0.5)
{
  npatientsTest <- dim(patientsTest)[1]
  npatientsTrain <- dim(patientsTrain)[1]
  nfeatcell <- dim(patientsTrain)[2]
  nchemicals <- dim(response)[2]
  
  patientsTrain <- as.matrix(patientsTrain)
  patientsTest <- as.matrix(patientsTest)
  
  pred <- matrix(data = 0, nrow = npatientsTest, ncol = nchemicals)
  
  # Treat chemicals one by one
  for (i in seq(nchemicals)) {
    # Train a lasso model
    cvob1 <- glmnet::cv.glmnet(patientsTrain, response[ ,i] , alpha = alpha)
    # Make predictions
    pred[ ,i] <- predict(cvob1, patientsTest, s = "lambda.min")
  }
  
  pred <- data.frame(pred)
  dimnames(pred) <- list(rownames(patientsTest), colnames(response))
  return(pred)
}
