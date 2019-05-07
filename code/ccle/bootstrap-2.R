######################################################################
# Data Sharing for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee T. 
# Email: asiae002@umn.edu
######################################################################
source("00-paths.R")
norm2 <- function(x) sqrt(sum(x^2))
norm1 <- function(x) sum(abs(x))
normInf <- function(x) max(abs(x))

source("projectOntoElasticNet.R")
source("projectOntoL1.R")
###########################
ssinglePredict <- function(fittedBetas, x, gs){#, addIntercept){
  # if(addIntercept == TRUE){
  groups <- as.numeric(gs)
  x <- cbind(rep(1, nrow(x)), x) #intercept is the first element
  # }
  predictedY <- rep(0, nrow(x))
  nGroup <- dim(fittedBetas)[2] - 1
  for(g in 1:nGroup){
    groupIndex <- groups == g
    predictedY[groupIndex] <-  x[groupIndex,] %*% (fittedBetas[,g] + fittedBetas[,nGroup + 1])
  }
  return(predictedY)
}

smultiDeal <- function(.Object){
  # There is no auto-complete for objects! Simplify:
  x <- .Object@x
  y <- .Object@y
  gs <- .Object@g
  gamma <- .Object@gamma
  grdTaus <- .Object@grdTaus
  
  nGroup <- dim(grdTaus)[2] - 1
  groups <- as.numeric(gs)
  # Scaling both to avoid intercept. 
  p <- dim(x)[2]
  avgX <- matrix(NA, nGroup, p)
  avgY <- rep(NA, nGroup)
  ng <- rep(0, nGroup)
  n <- length(y)
  # Saving averages for future intercept computation:
  for(g in 1:nGroup){
    groupIndex <- groups == g
    if(!any(groupIndex)){
      print(paste("There should be at least one elements in group", levels(groups)[g]))
      return(NULL)
    }
    ng[g] <- sum(groupIndex)
    avgX[g,] <- colMeans(x[groupIndex,])
    avgY[g] <- mean(y[groupIndex])
    x[groupIndex,] <- x[groupIndex,] - avgX[g,] #scale(x[groupIndex,]) : scale will generate NaN
    y[groupIndex] <- y[groupIndex] - avgY[g]
  }
  
  BetaNex <- matrix(0, p, nGroup + 1) #warm start beta. 
  
  # etas <- etas * 10
  nStepsGD <- 50
  stoppingCriteria <- .001 #1 / (20 * nrow(x))#.0001
  # For each tuning parameter:
  for(i in 1:dim(grdTaus)[1]){
    # .Object@betaList[[i]] <- deal(.Object@x, .Object@y, .Object@g, .Object@gamma, .Object@grdTaus[i,])
    # marginalImp <- rep(0, n)
    ## Walk nStepGD number of gradient descent steps:
    etas <- sapply(1:(nGroup + 1), function(g) ifelse(g != nGroup + 1, 1/(sqrt(n * ng[g])) ,1/n))
    stableErrCounter <- 0
    meanSqErrNew <- Inf 
    BetaNex <- matrix(0, p, nGroup + 1) #no warm start beta.
    for(j in 1:nStepsGD){
      BetaPrv <- BetaNex # Warm start from previous set of parameters. 
      meanSqErrOld <- meanSqErrNew 
      accumMult <- c(rep(0, n))
      ### Update per group parameters:
      for(g in 1:nGroup){
        groupIndex <- groups == g
        if(!any(groupIndex)){
          print(paste("There should be at least one elements in group", levels(groups)[g]))
          return(NULL)
        }
        recycle <- as.numeric(x[groupIndex,] %*% BetaPrv[,g])
        accumMult[groupIndex] <- recycle# to recycle the computation
        BetaNex[,g] <-  BetaPrv[,g] + etas[g] * t(x[groupIndex,]) %*% (y[groupIndex] - recycle - x[groupIndex,] %*% BetaPrv[,nGroup + 1])
        BetaNex[,g] <- projectOntoL1(BetaNex[,g], as.double(grdTaus[i,g]))
        # BetaNex[,g] <- projectOntoElasticNet(BetaNex[,g], gamma, as.double(grdTaus[i,g]))
      }
      ### Update the shared parameter:
      BetaNex[,nGroup + 1] <-  BetaPrv[,nGroup + 1] + etas[nGroup + 1] * t(x) %*% (y - x %*% BetaPrv[,nGroup + 1] - accumMult)
      # BetaNex[,nGroup + 1] <- projectOntoElasticNet(BetaNex[,nGroup + 1], gamma, as.double(grdTaus[i,nGroup + 1]))
      BetaNex[,nGroup + 1] <- projectOntoL1(BetaNex[,nGroup + 1], as.double(grdTaus[i,nGroup + 1]))
      ### Stopping criteria.
      # predictedY <- singlePredict(rbind(compIntercept(BetaNex, avgY, avgX, ng, n), BetaNex), x, groups)
      predictedY <- singlePredict(rbind(0, BetaNex), x, groups)
      
      meanSqErrNew <- mean((y - predictedY)^2)
      # print(paste("New mean square error:", meanSqErrNew))
      
      if(meanSqErrNew > meanSqErrOld){
        stableErrCounter <- stableErrCounter + 1
        if(stableErrCounter == 5){
          # print(paste("Enough learning! at step", j, "no marginal improvement")) # of meanSqErr is", abs(meanSqErrNew - meanSqErrOld)))
          break
        }
        BetaNex <- BetaPrv #do not update the parameters, keep the best one
        meanSqErrNew <- meanSqErrOld #and err will stay the same.
        etas <- etas / 1.5 #choose smaller step sizes over time. 
        # print(paste("Mean squared error is: ", meanSqErrNew))
      }else{stableErrCounter <- 0}
      # if(i == nStepsGD){
      #   print("hey")
      #   print(paste("After", nStepsGD, "steps of GD, I reached ", normInf(BetaPrv - BetaNex), "consecutive change in beta_{ij}."))
      # }
      
      
    }
    ## Computing intercepts
    # .Object@betaList[[i]] <- rbind(compIntercept(BetaNex, avgY, avgX, ng, n), BetaNex)
    .Object@betaList[[i]] <- rbind(0, BetaNex)
  }
  return(.Object@betaList)
}

sdealer <- setClass(
  # Set the name for the class
  "sdealer",
  
  # Define the slots
  slots = c(
    x = "matrix",
    y = "numeric",
    g = "factor",
    gamma = "numeric",
    grdTaus = "matrix",
    betaList = "list",
    meanX = "numeric",
    sdX = "numeric",
    meanY = "numeric",
    sdY = "numeric"
  )
)

setGeneric(name="spredictDeal",
           def=function(.Object, newx, newg, gamma, s)
           {
             standardGeneric("spredictDeal")
           }
)

setMethod(f="spredictDeal",
          signature="sdealer",
          definition=function(.Object, newx, newg, gamma, s)
          {
            # TODO: Assert size of s matches @grdTaus
            # Currently we are not learning gamma by CV so ignore it for now.  
            predictedRes <- matrix(NA, nrow(newx), dim(s)[1])
            ## Normalize by trainX information
            newx <- sweep(newx, 2, FUN = "-", .Object@meanX) #subtract mean train features
            newx <- sweep(newx, 2, FUN = "/", .Object@sdX) #divide by sd of train features
            for(i in 1:dim(s)[1]){
              # Finding the closest recorded tau to the queried one
              tausIndex <- which.min(apply(sweep(.Object@grdTaus, 2, FUN = "-", s[i,]), 1, norm2))
              predictedRes[,i] <- ssinglePredict(.Object@betaList[[tausIndex]], newx, newg)#, addIntercept = TRUE)
            }
            # Add back mean and sd of the response from trainY
            return(predictedRes * .Object@sdY + .Object@meanY)
          }
)

setMethod(f = "initialize",
          signature = "sdealer",
          definition = function(.Object, x, y, g, gamma, grdTaus, normalize = TRUE)
          {
            if(!normalize){
              # TODO: here :) 
              # Assuming data is normalized. TODO: check this! still experimental. 
              .Object@x <- x
              .Object@y <- y  
              .Object@meanX <- rep(0, dim(x)[2])
              .Object@sdX <- rep(1, dim(x)[2])
              .Object@meanY <- 0
              .Object@sdY <- 1
            }else{
              # Keeping mean and sd of the test data to compute z-score for test data in prediction time
              # Note that even though we are normalizing, we will need intercepts. 
              .Object@meanY <- mean(y)
              .Object@sdY <- sd(y)
              .Object@y <- as.numeric(scale(y))
              
              xSd <- apply(x, 2, sd)
              names(xSd) <- NULL
              constFeatureIndex <- xSd == 0
              variableFeatureIndex <- xSd != 0
              actualP <- dim(x)[2]
              
              x <- x[,variableFeatureIndex]
              p <- dim(x)[2]
              .Object@sdX <- xSd <- apply(x, 2, sd)
              .Object@meanX <- colMeans(x)
              .Object@x <- scale(x)
            }
            .Object@gamma <- gamma 
            .Object@grdTaus <- grdTaus
            .Object@g <- g
            
            # .Object@nGroup <- length(levels(g))
            tempBetaList <- multiDeal(.Object)
            tempBeta <- matrix(0,actualP,dim(grdTaus)[2])
            .Object@betaList <- rep(list(matrix(NA,actualP + 1,dim(grdTaus)[2])), length(tempBetaList)) #+1 for the intercept
            for(i in 1:length(tempBetaList)){
              tempBeta[variableFeatureIndex,] <- tempBetaList[[i]][2:(p+1),]
              .Object@betaList[[i]] <- rbind(tempBetaList[[i]][1,], tempBeta) #add intercept
            }
            return(.Object)
          }
)

###########################
singlePredict <- function(fittedBetas, x, gs){#, addIntercept){
  # if(addIntercept == TRUE){
  groups <- as.numeric(gs)
  x <- cbind(rep(1, nrow(x)), x) #intercept is the first element
  # }
  predictedY <- rep(0, nrow(x))
  nGroup <- dim(fittedBetas)[2] - 1
  for(g in 1:nGroup){
    groupIndex <- groups == g
    predictedY[groupIndex] <-  x[groupIndex,] %*% (fittedBetas[,g] + fittedBetas[,nGroup + 1])
  }
  return(predictedY)
}

multiDeal <- function(.Object){
  # There is no auto-complete for objects! Simplify:
  x <- .Object@x
  y <- .Object@y
  gs <- .Object@g
  gamma <- .Object@gamma
  grdTaus <- .Object@grdTaus
  
  nGroup <- dim(grdTaus)[2] - 1
  groups <- as.numeric(gs)
  # Scaling both to avoid intercept. 
  p <- dim(x)[2]
  avgX <- matrix(NA, nGroup, p)
  avgY <- rep(NA, nGroup)
  ng <- rep(0, nGroup)
  n <- length(y)
  # Saving averages for future intercept computation:
  for(g in 1:nGroup){
    groupIndex <- groups == g
    if(!any(groupIndex)){
      print(paste("There should be at least one elements in group", levels(groups)[g]))
      return(NULL)
    }
    ng[g] <- sum(groupIndex)
    avgX[g,] <- colMeans(x[groupIndex,])
    avgY[g] <- mean(y[groupIndex])
    x[groupIndex,] <- x[groupIndex,] - avgX[g,] #scale(x[groupIndex,]) : scale will generate NaN
    y[groupIndex] <- y[groupIndex] - avgY[g]
  }
  
  BetaNex <- matrix(0, p, nGroup + 1) #warm start beta. 
  
  # etas <- etas * 10
  nStepsGD <- 50
  stoppingCriteria <- .001 #1 / (20 * nrow(x))#.0001
  # For each tuning parameter:
  for(i in 1:dim(grdTaus)[1]){
    # .Object@betaList[[i]] <- deal(.Object@x, .Object@y, .Object@g, .Object@gamma, .Object@grdTaus[i,])
    # marginalImp <- rep(0, n)
    ## Walk nStepGD number of gradient descent steps:
    etas <- sapply(1:(nGroup + 1), function(g) ifelse(g != nGroup + 1, 1/(sqrt(n * ng[g])) ,1/n))
    stableErrCounter <- 0
    meanSqErrNew <- Inf 
    BetaNex <- matrix(0, p, nGroup + 1) #no warm start beta.
    for(j in 1:nStepsGD){
      BetaPrv <- BetaNex # Warm start from previous set of parameters. 
      meanSqErrOld <- meanSqErrNew 
      accumMult <- c(rep(0, n))
      ### Update per group parameters:
      for(g in 1:nGroup){
        groupIndex <- groups == g
        if(!any(groupIndex)){
          print(paste("There should be at least one elements in group", levels(groups)[g]))
          return(NULL)
        }
        recycle <- as.numeric(x[groupIndex,] %*% BetaPrv[,g])
        accumMult[groupIndex] <- recycle# to recycle the computation
        BetaNex[,g] <-  BetaPrv[,g] + etas[g] * t(x[groupIndex,]) %*% (y[groupIndex] - recycle - x[groupIndex,] %*% BetaPrv[,nGroup + 1])
        BetaNex[,g] <- projectOntoL1(BetaNex[,g], as.double(grdTaus[i,g]))
        # BetaNex[,g] <- projectOntoElasticNet(BetaNex[,g], gamma, as.double(grdTaus[i,g]))
      }
      ### Update the shared parameter:
      BetaNex[,nGroup + 1] <-  BetaPrv[,nGroup + 1] + etas[nGroup + 1] * t(x) %*% (y - x %*% BetaPrv[,nGroup + 1] - accumMult)
      # BetaNex[,nGroup + 1] <- projectOntoElasticNet(BetaNex[,nGroup + 1], gamma, as.double(grdTaus[i,nGroup + 1]))
      BetaNex[,nGroup + 1] <- projectOntoL1(BetaNex[,nGroup + 1], as.double(grdTaus[i,nGroup + 1]))
      ### Stopping criteria.
      predictedY <- singlePredict(rbind(compIntercept(BetaNex, avgY, avgX, ng, n), BetaNex), x, groups)
      meanSqErrNew <- mean((y - predictedY)^2)
      # print(paste("New mean square error:", meanSqErrNew))
      
      if(meanSqErrNew > meanSqErrOld){
        stableErrCounter <- stableErrCounter + 1
        if(stableErrCounter == 5){
          # print(paste("Enough learning! at step", j, "no marginal improvement")) # of meanSqErr is", abs(meanSqErrNew - meanSqErrOld)))
          break
        }
        BetaNex <- BetaPrv #do not update the parameters, keep the best one
        meanSqErrNew <- meanSqErrOld #and err will stay the same.
        etas <- etas / 1.5 #choose smaller step sizes over time. 
        # print(paste("Mean squared error is: ", meanSqErrNew))
      }else{stableErrCounter <- 0}
    }
    ## Computing intercepts
    .Object@betaList[[i]] <- rbind(compIntercept(BetaNex, avgY, avgX, ng, n), BetaNex)
  }
  return(.Object@betaList)
}


compIntercept <- function(Beta, avgY, avgX, ng, n){
  nGroup <- dim(Beta)[2] - 1
  intercept <- rep(NA, nGroup + 1)
  for(g in 1:nGroup){
    intercept[g] <- avgY[g] - avgX[g, ] %*% (Beta[,g] + Beta[,nGroup + 1])
  }
  intercept[nGroup + 1] <- as.numeric(intercept[1:nGroup] %*% ng) / n
  return(intercept)
}



dealer <- setClass(
  # Set the name for the class
  "dealer",
  
  # Define the slots
  slots = c(
    x = "matrix",
    y = "numeric",
    g = "factor",
    gamma = "numeric",
    grdTaus = "matrix",
    betaList = "list",
    meanX = "numeric",
    sdX = "numeric",
    meanY = "numeric",
    sdY = "numeric"
  )
)

setGeneric(name="predictDeal",
           def=function(.Object, newx, newg, gamma, s)
           {
             standardGeneric("predictDeal")
           }
)

setMethod(f="predictDeal",
          signature="dealer",
          definition=function(.Object, newx, newg, gamma, s)
          {
            # TODO: Assert size of s matches @grdTaus
            # Currently we are not learning gamma by CV so ignore it for now.  
            predictedRes <- matrix(NA, nrow(newx), dim(s)[1])
            ## Normalize by trainX information
            newx <- sweep(newx, 2, FUN = "-", .Object@meanX) #subtract mean train features
            newx <- sweep(newx, 2, FUN = "/", .Object@sdX) #divide by sd of train features
            for(i in 1:dim(s)[1]){
              # Finding the closest recorded tau to the queried one
              tausIndex <- which.min(apply(sweep(.Object@grdTaus, 2, FUN = "-", s[i,]), 1, norm2))
              predictedRes[,i] <- singlePredict(.Object@betaList[[tausIndex]], newx, newg)#, addIntercept = TRUE)
            }
            # Add back mean and sd of the response from trainY
            return(predictedRes * .Object@sdY + .Object@meanY)
          }
)

setMethod(f = "initialize",
          signature = "dealer",
          definition = function(.Object, x, y, g, gamma, grdTaus, normalize = TRUE)
          {
            if(!normalize){
              # TODO: here :) 
              # Assuming data is normalized. TODO: check this! still experimental. 
              .Object@x <- x
              .Object@y <- y  
              .Object@meanX <- rep(0, dim(x)[2])
              .Object@sdX <- rep(1, dim(x)[2])
              .Object@meanY <- 0
              .Object@sdY <- 1
            }else{
              # Keeping mean and sd of the test data to compute z-score for test data in prediction time
              # Note that even though we are normalizing, we will need intercepts. 
              .Object@meanY <- mean(y)
              .Object@sdY <- sd(y)
              .Object@y <- as.numeric(scale(y))
              
              xSd <- apply(x, 2, sd)
              names(xSd) <- NULL
              constFeatureIndex <- xSd == 0
              variableFeatureIndex <- xSd != 0
              actualP <- dim(x)[2]
              
              x <- x[,variableFeatureIndex]
              p <- dim(x)[2]
              .Object@sdX <- xSd <- apply(x, 2, sd)
              .Object@meanX <- colMeans(x)
              .Object@x <- scale(x)
            }
            .Object@gamma <- gamma 
            .Object@grdTaus <- grdTaus
            .Object@g <- g
            
            # .Object@nGroup <- length(levels(g))
            tempBetaList <- multiDeal(.Object)
            tempBeta <- matrix(0,actualP,dim(grdTaus)[2])
            .Object@betaList <- rep(list(matrix(NA,actualP + 1,dim(grdTaus)[2])), length(tempBetaList)) #+1 for the intercept
            for(i in 1:length(tempBetaList)){
              tempBeta[variableFeatureIndex,] <- tempBetaList[[i]][2:(p+1),]
              .Object@betaList[[i]] <- rbind(tempBetaList[[i]][1,], tempBeta) #add intercept
            }
            return(.Object)
          }
)
###########################
library("glmnet")
source("00-paths.R")
source("projectOntoElasticNet.R")
source("projectOntoL1.R")
# source("04-dataEnrichmentPrediction.R")


corrThresh <- 0.2
focusedCancerTypes <- c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "LUNG")
# focusedCancerTypes <- c("BREAST", "OVARY", "SKIN")


load(file=file.path(paths$scratch, paste("newer_cvMinLambdas", paste(focusedCancerTypes, collapse = "_") ,".RData")))
bestTaus <- cvMinLambdas
# load(file=file.path(paths$scratch, paste("new_cvMinLambdas_", paste(focusedCancerTypes, collapse = "_") ,".RData")))
# bestLambdas <- cvMinLambdas
# bestLambdas <- bestLambdas[25:48]#why the first 24 is garbage :)
load(file=file.path(paths$scratch, paste("newester_cvMinLambdas", paste(focusedCancerTypes, collapse = "_") ,".RData")))
bestSharedTaus <- cvMinLambdas



# minDeMeanIndex <- apply(cvDeMeans, 1, which.min)
# bestTaus <- sapply(1:24, function(i) cvDeTaus[minDeMeanIndex[i],,i])
nGroup <- dim(bestTaus)[2] - 1

set.seed(123)
dataDir <- file.path(paths$clean, "xy") 
allFiles <- list.files(dataDir)

errFull <- matrix(NA, nrow = 24, ncol = 100)
errShared <- matrix(NA, nrow = 24, ncol = 100)

drugCnt <- 0
drugs <- c()
for(file in allFiles){
  load(file.path(dataDir,file))
  drugCnt <- drugCnt + 1
  data <- XyActArea 
  # Selecting cancer types
  data <- data[(data$cancerType %in% focusedCancerTypes), ] 
  response <- data[,"y", drop=TRUE]
  ## Separating group id, necessary for the elastic net. 
  groupsId <- as.factor(data$cancerType)
  predictors <- data[, !(names(data) %in% c("y", "colors", "cancerType"))]
  ## Recovering drug name
  drugName <- strsplit(file, "_")[[1]][2]
  drugName <- strsplit(drugName, "[.]")[[1]][1]
  print(paste("Data for", drugName, "loaded!"))
  drugs[drugCnt] <- drugName
  
  # Using data to trim features.
  ## Removing features with variance zero
  zeroVarFeatureIndex <- (sapply(predictors, var) <= .02^2)
  predictors <- predictors[, !zeroVarFeatureIndex]
  ## Removing less relevant features 
  correlations <- sapply(predictors, cor, y=response)
  # correlations[sapply(correlations, is.na)] <- 0
  # print(paste("Here is number of NA in corr:", sum(sapply(correlations, is.na))))
  mask <- (abs(correlations) >= corrThresh) 
  bestPredictors <- predictors[,mask]
  p <- ncol(bestPredictors)
  n <- nrow(bestPredictors)
  print(dim(bestPredictors))
  print("Starting the bootstrap!")
  runningSumOfImpIndex <- matrix(0L, nrow = p, ncol = nGroup + 1) #-2 for y, colors, and cancer/groupId which is in data
  # runningsumOfCoeffsEffect <- matrix(0L, nrow = p, ncol = nGroup)
  numBootStrap <- 100
  flag <- FALSE
  for(i in 1:numBootStrap){
    print(paste("Resample number ", i))
    bsIndex <- sample(1:n, n, replace=TRUE)
    
    # Check if there is zero sample in any group of the resampled data
    for(g in 1:nGroup){
      groupIndex <- as.numeric(groupsId[bsIndex]) == g
      if(!any(groupIndex)){
        flag <- TRUE
        break
      }
    }
    if(flag){
      i <- i - 1
      next 
    }
    
    
    fullTaus <- t(as.matrix(c(.5,.25,1)))
    # fullTaus <- t(as.matrix(bestTaus[drugName,]))
    dealerFit <- dealer(x = as.matrix(bestPredictors[bsIndex,]), y = response[bsIndex], g = groupsId[bsIndex], gamma = 1, grdTaus = fullTaus, normalize = TRUE)
    sharedTaus  <- t(as.matrix(c(0,0,.1)))
    # sharedTaus <- t(as.matrix(bestSharedTaus[drugName,]))
    dealerFitShared <- sdealer(x = as.matrix(bestPredictors[bsIndex,]), y = response[bsIndex], g = groupsId[bsIndex], gamma = 1, grdTaus = sharedTaus, normalize = TRUE)
    
    # t-test for errors
    yHatFull <- predictDeal(dealerFit, newx=as.matrix(bestPredictors[bsIndex,]),newg = groupsId[bsIndex], gamma = 1, s=fullTaus)
    yHatShared <- spredictDeal(dealerFitShared, newx=as.matrix(bestPredictors[bsIndex,]),newg = groupsId[bsIndex], gamma = 1, s=sharedTaus) 
    # Here yHats are just vectors
    errFull[drugCnt, i] <- sapply(as.data.frame(yHatFull), function(x,y) mean((x - y)^2), response[bsIndex])
    errShared[drugCnt, i] <- sapply(as.data.frame(yHatShared), function(x,y) mean((x - y)^2), response[bsIndex])
    
    # p-val support
    temp <- 1 * (abs(dealerFit@betaList[[1]][-1,]) > 0) #remove the intercept and very small coeffs.
    runningSumOfImpIndex <- runningSumOfImpIndex + temp
  }
  avgFreqOfImpIndex <- runningSumOfImpIndex / numBootStrap
  # impFeatureIndex <- (avgFreqOfImpIndex >= .8 & avgOfCoeffsEffect > 0.8) #avgOfCoeffsEffect > 0 select those features that increase y on average.
  # print(sum(impFeatureIndex) / 3)
  impFeatureIndex <- (avgFreqOfImpIndex >= .5)
  print(sum(impFeatureIndex) / (nGroup + 1))
  
  for(i in 1:(nGroup+1)){
    fileConnGroup <- file(file.path(paths$scratch, paste(file, "group-small-2-manual",i,".txt")))
    nameOfImpCoef <- names(bestPredictors)[impFeatureIndex[,i]]
    nameOfImpCoef <- sub("-cna", "", nameOfImpCoef)
    nameOfImpCoef <- sub("-exp", "", nameOfImpCoef)
    nameOfImpCoef <- sub("-hybridMut", "", nameOfImpCoef)
    nameOfImpCoef <- sub("-oncoMut", "", nameOfImpCoef)
    writeLines(nameOfImpCoef, fileConnGroup)
    close(fileConnGroup)
  }
}


save(errFull, file=file.path(paths$scratch, paste("errDE", paste(focusedCancerTypes, collapse = "_") ,".RData")))
save(errShared, file=file.path(paths$scratch, paste("errDEShared", paste(focusedCancerTypes, collapse = "_") ,".RData")))

pdf("errorBoxPlots.pdf")
cnt <- 0 
par(mfrow=c(3,2),mai = c(.2, 0.2, 0.2, 0.2), oma = c(1, 1, 2, 2))
for(i in 1:24){
  myData <- data.frame(method = rep(c("DE", "LASSO"), each=100), err = c(errFull[i,], errShared[i,]))  
  boxplot(err~method, data=myData, col=c("red","blue"), xlab = "Method", ylab="MSE")
  testRes <- t.test(errFull[i,], errShared[i,], paired=TRUE, alternative = "less")
  if(testRes$p.value < 0.05){
    cnt <- cnt+1
  }
  title(paste(drugs[i], ":Mean diff:", round(testRes$estimate,2), ",p-val.:", formatC(testRes$p.value, format = "e", digits = 2)))
}
dev.off()

