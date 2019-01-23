######################################################################
# Data Enrichment for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee 
# Email: asiae002@umn.edu
######################################################################
# Remove these two lines at the end. 
# focusedCancerTypes <- c("BREAST", "OVARY", "SKIN")
focusedCancerTypes <- c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "LUNG")
corrThresh <- .3
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
    mask <- (abs(correlations) >= corrThresh) 
    bestTrainX <- trainX[,mask]
    bestTestX <- testX[,mask]
    print(dim(bestTrainX))
    
    
    
    cvResults[k, ] <- err
  }
  cvMeans[drugName,] <-  colMeans(cvResults)
  cvSds[drugName,] <- apply(cvResults, 2, sd)
  cvMinLambdas[drugName] <- myLambda[which.min(cvMeans[drugName,])]
  cvLambdas[drugName,] <- myLambda
  print(min(cvMeans[drugName,]))
  
}
save(cvMeans, file=file.path(paths$scratch, paste("new_cvMeans_", paste(focusedCancerTypes, collapse = "_") ,".RData")))
save(cvSds, file=file.path(paths$scratch, paste("new_cvSds_", paste(focusedCancerTypes, collapse = "_") ,".RData")))
save(cvLambdas, file=file.path(paths$scratch, paste("new_cvLambdas_", paste(focusedCancerTypes, collapse = "_") ,".RData")))
save(cvMinLambdas, file=file.path(paths$scratch, paste("new_cvMinLambdas_", paste(focusedCancerTypes, collapse = "_") ,".RData")))

minIndex <- apply(cvMeans, 1, which.min)
minMean <- sapply(1:24, function(i) cvMeans[i, minIndex[i]])
minSds <- sapply(1:24, function(i) cvSds[i, minIndex[i]])



######################################################################
# Data Enrichment for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee
# Email: asiae002@umn.edu
######################################################################
# if(!file.exists(file=file.path(paths$scratch, paste("cvDeMeans_", paste(focusedCancerTypes, collapse = "_") ,".RData")))){
set.seed(123)
dataDir <- file.path(paths$clean, "xy") 
allFiles <- list.files(dataDir)
nfolds <- 5
nGroup <- length(focusedCancerTypes)
numBaseTaus <- 5

#Prepare result storage.
drugList <- c()
for(file in allFiles){
  ## Recovering drug name
  drugName <- strsplit(file, "_")[[1]][2]
  drugName <- strsplit(drugName, "[.]")[[1]][1]
  drugList <- append(drugList, drugName)
}
cvDeMeans <- matrix(NA, nrow = length(drugList), ncol = numBaseTaus ^ (nGroup + 1))
rownames(cvDeMeans) <- drugList
cvDeSds <- matrix(NA, nrow = length(drugList), ncol = numBaseTaus ^ (nGroup + 1))
rownames(cvDeSds) <- drugList
cvDeTaus <- array(rep(NA, length(drugList) * numBaseTaus^(nGroup + 1) * (nGroup + 1)),
                  dim = c(numBaseTaus^(nGroup + 1), (nGroup + 1), length(drugList)),
                  dimnames = list(c(),c(),drugList)
)
cvMinLambdas <- matrix(NA, nrow = length(drugList), ncol = (nGroup + 1))
rownames(cvMinLambdas) <- drugList

counter <- 0
for(file in allFiles){
  load(file.path(dataDir,file))
  counter <- counter + 1
  # if(counter == 5)
  #   break
  data <- XyActArea 
  # Selecting cancer types
  data <- data[(data$cancerType %in% focusedCancerTypes), ] 
  response <- data[,"y", drop=TRUE]
  ## Separating group id, necessary for the elastic net. 
  groupsId <- as.numeric(as.factor(data$cancerType))
  predictors <- data[, !(names(data) %in% c("y", "colors", "cancerType"))]
  ## Recovering drug name
  drugName <- strsplit(file, "_")[[1]][2]
  drugName <- strsplit(drugName, "[.]")[[1]][1]
  print(paste("Data for", drugName, "loaded!"))
  
  # Cross-validation
  nTrain <- length(response)
  set.seed(123)
  foldsId <- sample(rep(1:nfolds, length.out = nTrain))
  ## Determining taus. 
  # Determining lambda.
  lambdaMax <- .5 #max(sapply(predictors, function(x,y) abs(sum(x * y)), response)) / nTrain
  lambdaMin <- .001 * lambdaMax
  myLambda <- exp(1)^seq(log(lambdaMin), log(lambdaMax), length.out = 10)
  myGroupLambda <- exp(1)^seq(log(lambdaMin), log(.1), length.out = 5)
  
  cvResults <- matrix(NA, nrow = nfolds, ncol = dim(gridTaus)[1])
  lowestErrs <- c()
  for (k in 1:nfolds) {
    glmnetFits <- list()
    print(paste("Fold", k, "started."))
    testId <- which(foldsId == k)
    trainX <- predictors[-testId, ]
    trainG <- groupsId[-testId]
    trainY <- response[-testId]
    testX <- predictors[testId, ]
    testG <- groupsId[testId]
    testY <- response[testId]
    
    # Using train data to trim features: Don't use test, bc you'll leak info from test to train. 
    ## Removing features with variance zero
    zeroVarFeatureIndex <- (sapply(trainX, var) <= .02^2)
    trainX <- trainX[, !zeroVarFeatureIndex]
    testX <- testX[, !zeroVarFeatureIndex]
    ## Removing less relevant features 
    correlations <- sapply(trainX, cor, y=trainY)
    # correlations[sapply(correlations, is.na)] <- 0
    mask <- (abs(correlations) >= corrThresh) 
    bestTrainX <- trainX[,mask]
    bestTestX <- testX[,mask]
    print(dim(bestTrainX))
    
    glmnetFit0 <- glmnet(x = as.matrix(bestTrainX), y = trainY, alpha=1, lambda = myLambda)
    yHatTrain <- predict(glmnetFit, newx=as.matrix(bestTrainX),s=myLambda) # make predictions
    trResidual <- sapply(as.data.frame(yHatTrain), function(x,y) (y - x), trainY)
    for(g in 1:nGroup){
      groupIndex <- trainG == g
      glmnetFits[[g]] <- list()
      for(i in 1:dim(trResidual)[2]){
        glmnetFits[[g]][[i]] <- glmnet(x = as.matrix(bestTrainX[groupIndex,]), y = trResidual[groupIndex, i], alpha=1, lambda = myGroupLambda)
      }
    }
    dealerFit <- dealer(x = as.matrix(bestTrainX), y = trainY, g = as.factor(trainG), gamma = 1, grdTaus = gridTaus, normalize = TRUE)
    yHat <- predictDeal(dealerFit, newx=as.matrix(bestTestX),newg = as.factor(testG), gamma = 1, s=gridTaus) 
    err <- sapply(as.data.frame(yHat), function(x,y) mean((x - y)^2), testY)
    cvResults[k, ] <- err
  }
  cvDeMeans[drugName,] <-  colMeans(cvResults)
  print(min(cvDeMeans[drugName,]))
  cvDeSds[drugName,] <- apply(cvResults, 2, sd)
  cvMinLambdas[drugName,] <- gridTaus[which.min(cvDeMeans[drugName,]), ]
  # cvMinLambdas <- gridTaus[which.min(cvDeMeans[drugName,]),]
  dealerFit <- dealer(x = as.matrix(bestTrainX), y = trainY, g = as.factor(trainG), gamma = 1, grdTaus = t(cvMinLambdas[drugName,]), normalize = TRUE)
  cvDeTaus[,,drugName] <- gridTaus
  
}
save(cvDeMeans, file=file.path(paths$scratch, paste("newer_cvDeMeans_", paste(focusedCancerTypes, collapse = "_") ,".RData")))
save(cvDeSds, file=file.path(paths$scratch, paste("newer_cvDeSds_", paste(focusedCancerTypes, collapse = "_") ,".RData")))
save(cvDeTaus, file=file.path(paths$scratch, paste("newer_cvDeTaus", paste(focusedCancerTypes, collapse = "_") ,".RData")))

save(cvMinLambdas, file=file.path(paths$scratch, paste("newer_cvMinLambdas", paste(focusedCancerTypes, collapse = "_") ,".RData")))



# }


