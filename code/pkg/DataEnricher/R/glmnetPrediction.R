######################################################################
# Data Sharing for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee T. 
# Email: asiae002@umn.edu
######################################################################

myElasticNetRegression <- function(data, givenLambda){
  response <- as.matrix(data[,"y", drop=FALSE])
  # predictors <- as.matrix(data[, -ncol(data)])
  predictors <- as.matrix(data[ , !(names(data) %in% c("y", "colors"))])
  
  best.predictors <- predictors
  print("Starting Elastic Net Regression ...")
  if(missing(givenLambda)){
    print("There is no Lambda! Let's cross-validate!")
    tmpCv1005 <- cv.glmnet(best.predictors, response, nfolds=10, alpha=.5)
    print("Done with the cross-validation of elastic net!")
  }
  else{
    print("There is a Lambda! Let's use it!")
    tmpCv1005 <- glmnet(best.predictors, response, alpha=.5, lambda = givenLambda)
    print("Done with the elastic net!")
  }
  return (tmpCv1005)
}
