#' Wrapper function for random forest regression
#' 
#' Wrapper function to perform random forest regression with \code{randomForest} that trains a model on training set and then predicts on test set for multiple tasks.
#' 
#' @param patientsTrain Matrix of training descriptors, of dimension \code{n x p}, for \code{n} training patients with \code{p} descriptors.
#' @param patientsTest Matrix of test descriptors, of dimension \code{m x p}, for \code{m} test patients with the same set of \code{p} descriptors.
#' @param response Matrix of observed toxicity values, of dimension \code{n x t}, for the \code{n} training patients responding to \code{t} drugs.
#' @param ntree Number of trees to grow a random forest. Default is 500. All other arguments are taken by default implementation of \code{randomForest}.
#' 
#' @return A matrix of predicted toxicity values, of dimension \code{m x t}, for the \code{m} test patients responding to the \code{t} drugs.
#' 
#' @export
#' 
#' @note Prediction is made per task with no special treatment for multitask learning, nor are task features needed.
#' 
#' @seealso \code{\link[randomForest]{randomForest}}
#' 
#' @importFrom randomForest randomForest
#' 

predictorRF <- function(patientsTrain, 
                        patientsTest, 
                        response, 
                        ntree = 500)
{
  npatientsTest <- dim(patientsTest)[1]
  npatientsTrain <- dim(patientsTrain)[1]
  nfeatcell <- dim(patientsTrain)[2]
  nchemicals <- dim(response)[2]
  
  pred <- matrix(data = 0, nrow = npatientsTest, ncol = nchemicals)
  
  # Treat chemicals one by one
  for (i in seq(nchemicals)) {
    # Train a random forest model
    model_rf <- randomForest::randomForest(patientsTrain, response[ ,i], ntree = ntree)
    # Make predictions
    pred[ ,i] <- predict(model_rf, patientsTest)
  }
  
  pred <- data.frame(pred)
  dimnames(pred) <- list(rownames(patientsTest), colnames(response))
  return(pred)
}
