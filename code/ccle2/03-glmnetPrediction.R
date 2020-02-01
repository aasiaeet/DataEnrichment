######################################################################
# Data Enrichment for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee 
# Email: asiae002@umn.edu
######################################################################
# Remove these two lines at the end. 
# focusedCancerTypes <- c("BREAST", "OVARY", "SKIN")
# # focusedCancerTypes <- c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "LUNG")
# corrThresh <- .3
###########################
# if(!file.exists(file.path(paths$scratch, paste("cvMeans_", paste(focusedCancerTypes, collapse = "_") ,".RData")))){
source("00-paths.R")
library(glmnet)
set.seed(123)
dataDir <- file.path(paths$clean, "xy") 
allFiles <- list.files(dataDir)
nfolds <- 5
numLambda <- 10 

#Prepare result storage.
drugList <- c()
for(file in allFiles){
  ## Recovering drug name
  drugName <- strsplit(file, "_")[[1]][2]
  drugName <- strsplit(drugName, "[.]")[[1]][1]
  drugList <- append(drugList, drugName)
}
cvMeans <- matrix(NA, nrow = length(drugList), ncol = numLambda)
rownames(cvMeans) <- drugList
cvSds <- matrix(NA, nrow = length(drugList), ncol = numLambda)
rownames(cvSds) <- drugList
cvLambdas <- matrix(NA, nrow = length(drugList), ncol = numLambda)
rownames(cvLambdas) <- drugList
cvMinLambdas <- matrix(NA, nrow = length(drugList), ncol = 1)
rownames(cvMinLambdas) <- drugList


for(file in allFiles){
  load(file.path(dataDir,file))
  data <- XyActArea 
  # Selecting cancer types
  data <- data[(data$cancerType %in% focusedCancerTypes), ] 
  response <- data[,"y", drop=TRUE]
  ## Removing group id, not necessary for the elastic net. 
  predictors <- data[, !(names(data) %in% c("y", "colors", "cancerType"))]
  ## Removing features with variance zero
  # predictors <- predictors[, !(sapply(predictors, var) <= 0.05)]
  ## Recovering drug name
  drugName <- strsplit(file, "_")[[1]][2]
  drugName <- strsplit(drugName, "[.]")[[1]][1]
  print(paste("Data for", drugName, "loaded!"))
  
  # Cross-validation
  nTrain <- length(response)
  set.seed(123)
  foldsId <- sample(rep(1:nfolds, length.out = nTrain))
  # Determining lambda.
  lambdaMax <- .5 #max(sapply(predictors, function(x,y) abs(sum(x * y)), response)) / nTrain
  lambdaMin <- .001 * lambdaMax
  myLambda <- exp(1)^seq(log(lambdaMin), log(lambdaMax), length.out = 10)
  
  cvResults <- matrix(NA, nrow = nfolds, ncol = numLambda)
  lowestErrs <- c()
  for (k in 1:nfolds) {
    print(paste("Fold", k, "started."))
    # ptm <- proc.time()
    testId <- which(foldsId == k)
    trainX <- predictors[-testId, ]
    trainY <- response[-testId]
    testX <- predictors[testId, ]
    testY <- response[testId]
    
    # Using train data to trim features: Don't use test, bc you'll leak info from test to train. 
    ## Removing features with variance zero
    zeroVarFeatureIndex <- (sapply(trainX, var) <= .02^2)
    trainX <- trainX[, !zeroVarFeatureIndex]
    testX <- testX[, !zeroVarFeatureIndex]
    ## Removing less relevant features 
    correlations <- sapply(trainX, cor, y=trainY)
    # correlations[sapply(correlations, is.na)] <- 0
    # print(paste("Here is number of NA in corr:", sum(sapply(correlations, is.na))))
    mask <- (abs(correlations) >= corrThresh) 
    bestTrainX <- trainX[,mask]
    bestTestX <- testX[,mask]
    print(dim(bestTrainX))
    
    
    # glmnetFit <- glmnet(x = as.matrix(bestTrainX), y = trainY, alpha=2/3, lambda = myLambda)
    glmnetFit <- glmnet(x = as.matrix(bestTrainX), y = trainY, alpha=1, lambda = myLambda)
    yHat <- predict(glmnetFit, newx=as.matrix(bestTestX),s=myLambda) # make predictions
    err <- sapply(as.data.frame(yHat), function(x,y) mean((x - y)^2), testY)
    cvResults[k, ] <- err
    # print(paste("Done with the elastic net on fold", k))
  }
  cvMeans[drugName,] <-  colMeans(cvResults)
  cvSds[drugName,] <- apply(cvResults, 2, sd)
  cvMinLambdas[drugName] <- myLambda[which.min(cvMeans[drugName,])]
  cvLambdas[drugName,] <- myLambda
  print(min(cvMeans[drugName,]))
  
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
save(cvMeans, file=file.path(paths$scratch, paste("new_cvMeans_", paste(focusedCancerTypes, corrThresh, collapse = "_") ,".RData")))
save(cvSds, file=file.path(paths$scratch, paste("new_cvSds_", paste(focusedCancerTypes, corrThresh, collapse = "_") ,".RData")))
save(cvLambdas, file=file.path(paths$scratch, paste("new_cvLambdas_", paste(focusedCancerTypes, corrThresh, collapse = "_") ,".RData")))
save(cvMinLambdas, file=file.path(paths$scratch, paste("new_cvMinLambdas_", paste(focusedCancerTypes, corrThresh, collapse = "_") ,".RData")))



# }


