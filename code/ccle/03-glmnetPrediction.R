######################################################################
# Data Sharing for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee T. 
# Email: asiae002@umn.edu
######################################################################
source("00-paths.R")
library(glmnet)
set.seed(123)
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


dataDir <- file.path(paths$clean, "xy") 

allFiles <- list.files(dataDir)
perDrugBestPerformanceMean <- c()
perDrugBestPerformanceSds <- c()
for(file in allFiles){
  # if(file == "README.md")
  #   next
  print(paste("Loading data from file ",  file ,"..."))
  load(file.path(dataDir,file))
  
  data <- XyActArea 
  data <- data[,!(names(data) %in% "cancerId")] #removing group id, not necessary for the elastic net. 
  
  print("Data loaded!")
  
  cv1005 <- myElasticNetRegression(data)
  savedLambda.min <- cv1005$lambda.min	
  # Save the best performance.
  perDrugBestPerformanceMean <- append(perDrugBestPerformanceMean, cv1005$cvm[cv1005$lambda == cv1005$lambda.min])
  perDrugBestPerformanceSds <- append(perDrugBestPerformanceSds, cv1005$cvsd[cv1005$lambda == cv1005$lambda.min])
    # par(mar=c(4,4,4,4))
    # plot(log(cv1005$lambda),cv1005$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv1005$name)
    # print("Done plotting the results!")
  #   
  ############################
  # This is for important feature selection. Disabled for now.
  ############################

#   print("Starting the bootstrap!")
#   runningSumOfImpIndex <- matrix(0L, nrow = ncol(data), ncol = 1)
#   numBootStrap <- 50
#   for(i in 1:numBootStrap){
#     print(paste("Resample number ", i))
#     bootStrapSampleIndex <- sample(1:nrow(data), nrow(data), replace=TRUE)
#     bootStrapData <- data[bootStrapSampleIndex,]
#     bootStrapCv1005 <- myElasticNetRegression(bootStrapData, savedLambda.min)
#     runningSumOfImpIndex <- runningSumOfImpIndex + as.integer(abs(as.matrix(coef(bootStrapCv1005))) >0.00000001)
#   }
  
  #readline(print("Stop!"))
 
  #avgFreqOfImpIndex <- runningSumOfImpIndex / numBootStrap
  #impFeatureIndex <- (avgFreqOfImpIndex >= .5)
  
  #valueOfCoef <- as.matrix(coef(cv1005, s="lambda.min"))
  #valueOfImpCoef <- valueOfCoef[impFeatureIndex]
  #nameOfImpCoef <- rownames(valueOfCoef)[impFeatureIndex]
  #sortedValueOfImpCoef <- sort(valueOfImpCoef, decreasing = TRUE, index.return=TRUE)
  #print(nameOfImpCoef[sortedValueOfImpCoef$ix])
  ############################
  
  # readline(print("Should I go to the next drug?!"))
}
save(perDrugBestPerformanceMean, file=file.path(paths$scratch, "perDrugBestPerformance-Mean-ElasticNet.RData"))
save(perDrugBestPerformanceSds, file=file.path(paths$scratch, "perDrugBestPerformance-Sd-ElasticNet.RData"))
