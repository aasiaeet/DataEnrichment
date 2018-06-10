######################################################################
# Data Sharing for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee T. 
# Email: asiae002@umn.edu
######################################################################
library("Matrix")
source("projectOntoElasticNet.R")
# source("projectOntoL1.R")
# library(glmnet)
source("00-paths.R")
# source("dataSharing.R")

source("projectOntoElasticNet.R")
# source("projectOntoL1.R")
# library(glmnet)
nStepsGD <- 10
stoppingCriteria <- .01 #1 / (20 * nrow(x))#.0001


norm2 <- function(x) sqrt(sum(x^2))
norm1 <- function(x) sum(abs(x))
normInf <- function(x) max(abs(x))


predict.dataSharing <- function(fittedModel, x, groups, taus, addIntercept = TRUE){
  if(addIntercept == TRUE){
    x <- cbind(x, rep(1, nrow(x))) #intercept
  }
  learnedBeta <- matrix(fittedModel, ncol=3)
  predictedY <- rep(0, nrow(x))
  for(g in 1:2){
    groupIndex <- groups == g
    predictedY[groupIndex] <-  x[groupIndex,] %*% (learnedBeta[,g] + learnedBeta[,3])
  }
  return(predictedY)
}

dataSharing <- function(x, y, groups, gamma, taus){
  # print(paste("Testing taus:", taus))
  eta <- 1/(10 * nrow(x))
  x <- cbind(x, rep(1, nrow(x))) #intercept
  BetaPrv <- matrix(0, ncol(x), 3)
  marginalImp <- rep(0, length(y))
  BetaNex <- BetaPrv
  meanSqErrOld <- 0
  for(i in 1:nStepsGD){
    # if(i %% 25 == 0){
    #   print(paste("This is iteration", i, "of SPGD"))
    # }
    eta <- (1 / (1+ i)) * eta
    gradient <- 0
    for(g in 1:2){
      groupIndex <- groups == g
      # (2/nrow(x)) * 
      scaledGradient <- sqrt(eta * (1/sum(groupIndex))) * t(x[groupIndex,]) %*% (y[groupIndex] - x[groupIndex,] %*% (BetaPrv[,g] + BetaPrv[,3]))
      # marginalImp[groupIndex] <- y[groupIndex] - x[groupIndex,] %*% (BetaPrv[,g] + BetaPrv[,3])
      BetaNex[,g] <-  BetaPrv[,g] + scaledGradient
      gradient <- gradient + scaledGradient
      # BetaNex[,g] <- projectOntoL1(BetaNex[,g], as.double(taus[g]))
      BetaNex[,g] <- projectOntoElasticNet(BetaNex[,g], gamma, as.double(taus[g]))
    }
    BetaNex[,3] <-  BetaPrv[,3] + gradient #2/(nrow(x)) * t(x) %*% marginalImp
    # BetaNex[,3] <- projectOntoL1(BetaNex[,3], as.double(taus[3]))
    BetaNex[,3] <- projectOntoElasticNet(BetaNex[,3], gamma, as.double(taus[3]))
    ### Stopping criteria. 
    predictedY <- predict.dataSharing(BetaNex, x, groups, taus, addIntercept = FALSE)
    meanSqErrNew <- mean((y - predictedY)^2)
    # print(paste("Mean squared error is: ", meanSqErrNew))
    if(abs(meanSqErrNew - meanSqErrOld) < stoppingCriteria){
      # print(paste("breaking bad! at step", i, "diff in meanSqErr is", abs(meanSqErrNew - meanSqErrOld)))
      break 
    }
    # if(i == nStepsGD){
    #   print("hey")
    #   print(paste("After", nStepsGD, "steps of GD, I reached ", normInf(BetaPrv - BetaNex), "consecutive change in beta_{ij}."))
    # }
    meanSqErrOld <- meanSqErrNew
    BetaPrv <- BetaNex
  }
  # print(paste("Training MSE was", meanSqErrOld))
  return(BetaNex)
}



# Cross-validation for data sharing.
cv.dataSharing <- function(predictors, response, groupsId, nfolds=10, gamma=2, tausGrid){
  stopifnot(nrow(predictors) == length(response))
  nTrain <- length(response)
  foldsId <- sample(rep(1:nfolds, length.out = nTrain))
  
  cvResults <- matrix(NA, nrow = nfolds, ncol = nrow(tausGrid))
  for (k in 1:nfolds) {
    # print(paste("Fold", k, "started."))
    # ptm <- proc.time()
    testId <- which(foldsId == k)
    
    trainX <- predictors[-testId, ]
    trainY <- response[-testId]
    trainG <- groupsId[-testId]
    testX <- predictors[testId, ]
    testY <- response[testId]
    testG <- groupsId[testId]
    
    fittedModels <- apply(t(tausGrid), 2, function(myTaus) dataSharing(trainX, trainY, trainG, gamma, myTaus))
    
    preditions <- apply(t(1:nrow(tausGrid)), 2, function(i) predict.dataSharing(fittedModels[,i], testX, testG, tausGrid[i,]) ) 
    
    cvResults[k, ] <- apply(preditions, 2, function(yHat) mean((testY - yHat)^2))
    print(paste("Fold", k, "best taus test MSE is", min(cvResults[k,])))
    # readline("pause")
  }
  cv <- list()
  cv$cvm <- colMeans(cvResults)
  cv$cvsd <- apply(cvResults, 2, sd)
  cv$taus.min <- tausGrid[which.min(cv$cvm),]
  #Fit all data and get predictor for each set of taus.
  #TODO: return for any desired taus
  # cv$dataSharing.fit <- as(apply(t(tausGrid), 2, function(myTaus) dataSharing(predictors, response, groupsId, gamma, myTaus)), "dgCMatrix")
  # cv$dataSharing.fit <- as(apply(cv$taus.min, 1, function(myTaus) dataSharing(predictors, response, groupsId, gamma, myTaus)), "dgCMatrix")
  cv$dataSharing.fit <- dataSharing(predictors, response, groupsId, gamma, cv$taus.min)
  return(cv)
}


coef <- function(cv){
  #TODO: return for any desired taus
  # return (cv$dataSharing.fit[,which.min(cv$cvm)])
  return (cv$dataSharing.fit)
}

myDataSharingSparse <- function(data, tausGrid, givenTau){
  #TODO: make the groupId first/one to last column. Or enforce  "y" for the response col name. 
  groupsId <- data$cancerId
  data <- data[,!(names(data) %in% "cancerId")] #removing group id after saving it. 
  response <- as.matrix(data[,"y", drop=FALSE])
  best.predictors <- as.matrix(data[, !(names(data) %in% c("y", "colors"))])
  
  # best.predictors <- predictors
  
  print("Starting SPGD for Data Sharing Regression ...")
  if(missing(givenTau)){
    print("There is no Tau! Let's cross-validate!")
    tmpCv1005 <- cv.dataSharing(best.predictors, response, groupsId, nfolds=10, gamma=10, tausGrid)
    print("Done with the cross-validation of data sharing!")
  }
  else{
    print("There is a Tau! Let's use it!")
    tmpCv1005 <- dataSharing(best.predictors, response, groupsId, gamma=10, tau = givenTau)
    print("Done with the data sharing!")
  }
  return (tmpCv1005)
}



# TODO: Save them as an object!!!
load(file.path(paths$scratch, "listOfRefinedCV-DS-Sparse.RData"))
tausMinAll <- matrix(NA, 24, 3)
for(i in 1:96){
    if(i %% 4 == 3){
      tausMinAll[i / 4 + 1, ] <- listOfRefinedCV[[i]]
    }
}

dataDir <- file.path(paths$clean, "xy") 

allFiles <- list.files(dataDir)
set.seed(123)
counter <- 5 # 5 is for the Erlotinib, 3 AZD530
for(file in allFiles){
  # counter <- counter + 1
  print(paste("Loading data from file ",  file ,"..."))
  load(paste(dataDir,file, sep="/"))
  data <- XyActArea 
  print(paste("Data loaded!"))
  
  
  print("Starting the bootstrap!")
  runningSumOfImpIndex <- matrix(0L, nrow = ncol(data)-2, ncol = 3) #-2 for y, colors, and cancer/groupId which is in data
  runningsumOfCoeffsEffect <- matrix(0L, nrow = ncol(data)-2, ncol = 3)
  numBootStrap <- 10
  for(i in 1:numBootStrap){
    print(paste("Resample number ", i))
    bootStrapSampleIndex <- sample(1:nrow(data), nrow(data), replace=TRUE)
    bootStrapData <- data[bootStrapSampleIndex,]
    bootStrapCv1005 <- myDataSharingSparse(bootStrapData, givenTau = tausMinAll[counter,])
    # bootStrapCv1005 <- myDataSharingSparse(bootStrapData, givenTau = c(1.9,2.2,3.2))
    temp <- 1 * (abs(bootStrapCv1005[-nrow(bootStrapCv1005),]) > 0.0005) #remove the intercept
    runningSumOfImpIndex <- runningSumOfImpIndex + temp
    
    # bootStrapCv1005 <- bootStrapCv1005[-nrow(bootStrapCv1005),]
    # ids <- bootStrapData$cancerId
    # bootStrapData <- bootStrapData[, !(names(bootStrapData) %in% c("y", "color", "cancerId"))]
    # meansOfFeaturesForLung <- colMeans(bootStrapData[ids == 1,])
    # meansOfFeaturesForBlood <- colMeans(bootStrapData[ids == 2,])
    # bootStrapCv1005[,1] <- (bootStrapCv1005[,1]+bootStrapCv1005[,3]) * meansOfFeaturesForLung
    # bootStrapCv1005[,2] <- (bootStrapCv1005[,2]+bootStrapCv1005[,3]) * meansOfFeaturesForBlood
    # bootStrapCv1005[,3] <- rep(1, length(bootStrapCv1005[,3])) #bootStrapCv1005[,3] * colMeans(bootStrapData)
    # runningsumOfCoeffsEffect <- runningsumOfCoeffsEffect + 1 * (bootStrapCv1005 < 0)
  }
  
  # avgOfCoeffsEffect <- runningsumOfCoeffsEffect / numBootStrap
  avgFreqOfImpIndex <- runningSumOfImpIndex / numBootStrap
  # impFeatureIndex <- (avgFreqOfImpIndex >= .8 & avgOfCoeffsEffect > 0.8) #avgOfCoeffsEffect > 0 select those features that increase y on average.
  # print(sum(impFeatureIndex) / 3)
  impFeatureIndex <- (avgFreqOfImpIndex >= .8)
  print(sum(impFeatureIndex) / 3)
  
  
  # cv1005 <- myDataSharingSparse(data, givenTau = tausMinAll[counter,])
  # valueOfCoef <- as.matrix(cv1005)
  # valueOfCoef <- valueOfCoef[-nrow(valueOfCoef), ] #remove the intercept
  # valueOfImpCoef <- valueOfCoef[impFeatureIndex]
  fileConnGroup1<-file(file.path(paths$scratch, paste(file, "group1.txt")))
  fileConnGroup2<-file(file.path(paths$scratch, paste(file, "group2.txt")))
  fileConnShared<-file(file.path(paths$scratch, paste(file, "shared.txt")))
  data <- data[, !(names(data) %in% c("y", "colors", "cancerId"))]
  
  for(i in 1:3){
    nameOfImpCoef <- names(data)[impFeatureIndex[,i]]
    nameOfImpCoef <- sub("-cna", "", nameOfImpCoef)
    nameOfImpCoef <- sub("-exp", "", nameOfImpCoef)
    nameOfImpCoef <- sub("-hybridMut", "", nameOfImpCoef)
    nameOfImpCoef <- sub("-oncoMut", "", nameOfImpCoef)
    
    if (i == 1){
      writeLines(nameOfImpCoef, fileConnGroup1)
    }
    if (i == 2){
      writeLines(nameOfImpCoef, fileConnGroup2)
    }
    if (i == 3){
      writeLines(nameOfImpCoef, fileConnShared)
    }
  }
  
  # sortedValueOfImpCoef <- sort(valueOfImpCoef, decreasing = TRUE, index.return=TRUE)
  # print(nameOfImpCoef[sortedValueOfImpCoef$ix])
  
  
  close(fileConnShared)
  close(fileConnGroup1)
  close(fileConnGroup2)
  
  ############################
  # readline(print("Should I go to the next drug?!"))
}


# currentMean <- 0
# currentSd <- 0
# minIndex <- 0
# sds <- c()
# means <- c()
# sds <- c()
# for(i in 1:96){
#   if(i %% 4 == 1){
#     currentMean <- listOfRefinedCV[[i]]
#     minIndex <- which.min(currentMean)
#   }
#   else if(i %% 4 == 2){
#     currentSd <- listOfRefinedCV[[i]]
#     means <- append(means, currentMean[minIndex])
#     sds <- append(sds, currentSd[minIndex])
#   }
# }

# 
# library(ggplot2)
# gList <- list()
# apply(t(1:3), 2,
#       function(i){
#         toPlot <- as.data.frame(matrix(NA,14,1))
#         toPlot$means <- apply(expRes[[1]], 2, mean)
#         toPlot$sds <- apply(expRes[[i]], 2, sd)
#         index <- i
#         if(i==3){
#           toPlot$n <- numSamples
#           index <- 0
#         }
#         else
#           toPlot$n <- ((3 - (i%%2)) / 5) * numSamples
#         gList[[i]] <- ggplot(data=toPlot,aes(x = n,y=means)) +
#           theme_bw() +
#           geom_errorbar(data=toPlot,aes(ymin=means-sds,ymax=means+sds), width=80, size=.5, color="black")+
#           geom_line(aes(y=means), color="red", size=.8) +
#           theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), text = element_text(size=20))+
#           xlab(substitute(n[nn], list(nn=index))) +
#           ylab(substitute("||"~delta[nn]~"||=||"~hat(beta)[nn] - beta[nn]^{"*"}~"||", list(nn=index)))
#       }
# )
# 



# nameDescMapping <- as.data.frame(expressions$Description)
# row.names(nameDescMapping) <- row.names(expressions)
# featureNames <- names(data)[1:(ncol(data)-2)]
# featureNames[impFeatureIndex[,1]]
# nameDescMapping[featureNames[impFeatureIndex[,1]],]
