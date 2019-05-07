######################################################################
# Data Enrichment for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee
# Email: asiae002@umn.edu
######################################################################
# Remove these two lines at the end. 
# focusedCancerTypes <- c("BREAST", "OVARY", "SKIN")
# focusedCancerTypes <- c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "LUNG")
# corrThresh <- .3
###########################
source("00-paths.R")
# source("dataEnricher.R")
norm2 <- function(x) sqrt(sum(x^2))
norm1 <- function(x) sum(abs(x))
normInf <- function(x) max(abs(x))

source("projectOntoElasticNet.R")
source("projectOntoL1.R")

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
      # if(i == nStepsGD){
      #   print("hey")
      #   print(paste("After", nStepsGD, "steps of GD, I reached ", normInf(BetaPrv - BetaNex), "consecutive change in beta_{ij}."))
      # }
      
      
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
  
  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  # validity=function(.Object)
  # {
  #   if((dim(.Object@x)[1] != length(.Object@y)) || (length(.Object@y) != length(.Object@g))) {
  #     return("Dimensions of x, y, or g are not matching.")
  #   }
  #   if(length(unique(.Object@g)) + 1 != dim(.Object@grdTaus)[2]){
  #     return("Number of taus is not equal to number of groups.")
  #   }
  #   return(TRUE)
  # }
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
              .Object@meanX <- colMeans(x)
              .Object@sdX <- apply(x, 2, sd)
              .Object@meanY <- mean(y)
              .Object@sdY <- sd(y)
              .Object@x <- scale(x)
              .Object@y <- as.numeric(scale(y))
            }
            .Object@g <- g
            # .Object@nGroup <- length(levels(g))
            .Object@gamma <- gamma 
            .Object@grdTaus <- grdTaus
            .Object@betaList <- multiDeal(.Object)
            # validObject(.Object)
            return(.Object)
          }
)
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

# cvDeMeans <- matrix(NA, nrow = length(drugList), ncol = numBaseTaus)
# rownames(cvDeMeans) <- drugList
# cvDeSds <- matrix(NA, nrow = length(drugList), ncol = numBaseTaus)
# rownames(cvDeSds) <- drugList
# cvDeTaus <- array(rep(NA, length(drugList) * numBaseTaus^(nGroup + 1) * (nGroup + 1)), 
#                   dim = c(numBaseTaus, (nGroup + 1), length(drugList)),
#                   dimnames = list(c(),c(),drugList)
# )
# cvMinLambdas <- matrix(NA, nrow = length(drugList), ncol = (nGroup + 1))
# rownames(cvMinLambdas) <- drugList


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
  groupsId <- data$cancerType
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
  tauMax <- .4#(max(sapply(predictors, function(x,y) abs(sum(x * y)), response)) / nTrain) / 10
  tauMin <- .01 #* (tauMax) 
  # myTaus <- exp(1)^seq(log(tauMin), log(tauMax), length.out = numBaseTaus)
  myTaus <- seq(tauMin, tauMax, length.out = numBaseTaus)
  gridTaus <- as.matrix(expand.grid(replicate(nGroup + 1, myTaus, simplify = F)))
  # gridTaus <- t(sapply(1:numBaseTaus, function(x) c(myTaus[x],0,0)))
  
  cvResults <- matrix(NA, nrow = nfolds, ncol = dim(gridTaus)[1])
  lowestErrs <- c()
  for (k in 1:nfolds) {
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
    # print(paste("Here is number of NA in corr:", sum(sapply(correlations, is.na))))
    mask <- (abs(correlations) >= corrThresh) 
    bestTrainX <- trainX[,mask]
    bestTestX <- testX[,mask]
    print(dim(bestTrainX))
    
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
save(cvDeMeans, file=file.path(paths$scratch, paste("newer_cvDeMeans_", paste(focusedCancerTypes, corrThresh, collapse = "_") ,".RData")))
save(cvDeSds, file=file.path(paths$scratch, paste("newer_cvDeSds_", paste(focusedCancerTypes, corrThresh, collapse = "_") ,".RData")))
save(cvDeTaus, file=file.path(paths$scratch, paste("newer_cvDeTaus", paste(focusedCancerTypes, corrThresh, collapse = "_") ,".RData")))

save(cvMinLambdas, file=file.path(paths$scratch, paste("newer_cvMinLambdas", paste(focusedCancerTypes,corrThresh, collapse = "_") ,".RData")))






# perpareTaus <- function(myTaus, nGroups = 2){
#   # One for the shared component and two for each of Lung and Blood cancers. 
#   gridTaus <- expand.grid(replicate(nGroups + 1, myTaus, simplify = F))
#   return(gridTaus)
#   # # print("We are looking for a good range of parameters ...")
#   # cv102 <- myDataSharingSparse(data, gridTaus)
#   # 
#   # bestTaus <- gridTaus[sort(cv102$cvm, index.return=TRUE)$ix[1:10],]
#   # thirdQuantile <- apply(bestTaus, 2, quantile)[4,]
#   # firstQuantile <- apply(bestTaus, 2, quantile)[2,]
#   # tausMax <- apply(bestTaus, 2, max)
#   # tausMin <- apply(bestTaus, 2, min)
#   # 
#   # for(i in 1:3){
#   #   if(max(taus) == thirdQuantile[i]){  # We are on the margine of the search space, let's expand it. 
#   #     tausMax[i] = 10
#   #   } 
#   #   if(min(taus) == firstQuantile[i]){
#   #     tausMin[i] = 0.001
#   #   }
#   # }
#   # # print(paste("We found the following range", tausMin, tausMax))
#   # refinedTaus <- apply(t(1:3), 2, function(i) (seq(from=tausMin[i], to=tausMax[i], length.out = 50)))
#   
# }
# }