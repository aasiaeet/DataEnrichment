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
