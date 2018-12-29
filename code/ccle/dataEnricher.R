######################################################################
# Data Enrichment for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee
# Email: asiae002@umn.edu
######################################################################
norm2 <- function(x) sqrt(sum(x^2))
norm1 <- function(x) sum(abs(x))
normInf <- function(x) max(abs(x))

source("projectOntoElasticNet.R")

singlePredict <- function(fittedBetas, x, groups, gamma, taus){#, addIntercept){
  # if(addIntercept == TRUE){
  x <- cbind(x, rep(1, nrow(x))) #intercept
  # }
  predictedY <- rep(0, nrow(x))
  nGroup <- length(taus) - 1
  for(g in 1:nGroup){
    groupIndex <- groups == g
    predictedY[groupIndex] <-  x[groupIndex,] %*% (fittedBetas[,g] + fittedBetas[,nGroup + 1])
  }
  return(predictedY)
}

multiDeal <- function(.Object){
  for(i in 1:dim(.Object@grdTaus)[1]){
    .Object@betaList[[i]] <- deal(.Object@x, .Object@y, .Object@g, .Object@gamma, .Object@grdTaus[i,])
  }
  return(.Object@betaList)
}


deal <- function(x, y, gs, gamma, taus){#, normalized = FALSE){
  nGroup <- length(taus) - 1
  groups <- as.numeric(gs)
  # Scaling both to avoid intercept. 
  p <- dim(x)[2]
  avgX <- matrix(NA, nGroup, p)
  avgY <- rep(NA, nGroup)
  ng <- rep(0, nGroup)
  n <- length(y)
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
  
  
  etas <- sapply(1:(nGroup + 1), function(g) ifelse(g != nGroup + 1, 1/(sqrt(n * ng[g])) ,1/n))
  # etas <- etas * 10
  nStepsGD <- 5 
  stoppingCriteria <- .001 #1 / (20 * nrow(x))#.0001
  
  # x <- cbind(x, rep(1 / nrow(x), nrow(x))) #no need for intercept computation
  BetaPrv <- matrix(0, p, nGroup + 1) #last column is dedicated for \beta_0
  # marginalImp <- rep(0, n)
  BetaNex <- BetaPrv
  meanSqErrOld <- 0
  for(i in 1:nStepsGD){
    # if(i %% 5 == 0){
    #   print(paste("This is iteration", i, "of SPGD"))
    # }
    accumMult <- c(rep(0, n))
    for(g in 1:nGroup){
      groupIndex <- groups == g
      if(!any(groupIndex)){
        print(paste("There should be at least one elements in group", levels(groups)[g]))
        return(NULL)
      }
      recycle <- as.numeric(x[groupIndex,] %*% BetaPrv[,g])
      accumMult[groupIndex] <- recycle# to recycle the computation
      # marginalImp[groupIndex] <- y[groupIndex] - x[groupIndex,] %*% (BetaPrv[,g] + BetaPrv[,nGroup + 1])
      BetaNex[,g] <-  BetaPrv[,g] + etas[g] * t(x[groupIndex,]) %*% (y[groupIndex] - recycle - x[groupIndex,] %*% BetaPrv[,nGroup + 1])
      # BetaNex[,g] <- projectOntoL1(BetaNex[,g], as.double(taus[g]))
      BetaNex[,g] <- projectOntoElasticNet(BetaNex[,g], gamma, as.double(taus[g]))
    }
    # BetaNex[,nGroup + 1] <-  BetaPrv[,nGroup + 1] + gradient #2/(nrow(x)) * t(x) %*% marginalImp
    BetaNex[,nGroup + 1] <-  BetaPrv[,nGroup + 1] + etas[nGroup + 1] * t(x) %*% (y - x %*% BetaPrv[,nGroup + 1] - accumMult)
    # BetaNex[,nGroup + 1] <- projectOntoL1(BetaNex[,nGroup + 1], as.double(taus[nGroup + 1]))
    BetaNex[,nGroup + 1] <- projectOntoElasticNet(BetaNex[,nGroup + 1], gamma, as.double(taus[nGroup + 1]))
    ### Stopping criteria.
    
    # predictedY <- singlePredict(BetaNex, x, groups, gamma = 2, taus)#, addIntercept = FALSE)
    # meanSqErrNew <- mean((y - predictedY)^2)
    # 
    
    # if(abs(meanSqErrNew - meanSqErrOld) < stoppingCriteria)){
    #   # print(paste("Mean squared error is: ", meanSqErrNew))
    #   print(paste("Enough learning! at step", i, "marginal improvement of meanSqErr is", abs(meanSqErrNew - meanSqErrOld)))
    #   break
    # }
    # # if(i == nStepsGD){
    # #   print("hey")
    # #   print(paste("After", nStepsGD, "steps of GD, I reached ", normInf(BetaPrv - BetaNex), "consecutive change in beta_{ij}."))
    # # }
    # meanSqErrOld <- meanSqErrNew
    # 
    
    BetaPrv <- BetaNex
  }
  # print(paste("Training MSE was", meanSqErrOld))
  # Computing intercepts
  intercept <- rep(NA, nGroup + 1)
  for(g in 1:nGroup){
    intercept[g] <- avgY[g] - avgX[g, ] %*% (BetaNex[,g] + BetaNex[,nGroup])
  }
  intercept[nGroup + 1] <- as.numeric(intercept[1:nGroup] %*% ng) / n
  BetaNex <- rbind(intercept, BetaNex)
  
  return(BetaNex)
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
    betaList = "list"
  ),
  
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
            predictedRes <- matrix(NA, nrow(newx), dim(s)[1])
            for(i in 1:dim(s)[1]){
              # Finding the closest recorded tau to the queried one
              tausIndex <- which.min(apply(sweep(.Object@grdTaus, 2, FUN = "-", s[i,]), 1, norm2))
              predictedRes[,i] <- singlePredict(.Object@betaList[[tausIndex]], newx, newg, gamma, .Object@grdTaus[tausIndex,])#, addIntercept = TRUE)
            }
            return(predictedRes)
          }
)

setMethod(f = "initialize",
          signature = "dealer",
          definition = function(.Object, x, y, g, gamma, grdTaus)
          {
            .Object@x <- x
            .Object@y <- y
            .Object@g <- g
            .Object@gamma <- gamma
            .Object@grdTaus <- grdTaus
            .Object@betaList <- multiDeal(.Object)
            # validObject(.Object)
            return(.Object)
          }
)


# print("start")
# dealerFit <- dealer(x = as.matrix(bestTrainX), y = trainY, g = as.factor(trainG), gamma = 2, grdTaus = gridTaus)
# print("stop")
# set.seed(123)
# xx <- matrix(1:20, 4, 5)
# yy <- c(0,9,8,7)
# gg <- c(1,2,1,2)
# gTaus <- matrix(runif(30, 0, 1), 10, 3)
# de <- dealer(x = xx, y = yy, g = as.factor(gg), gamma = 2, grdTaus = gTaus)



# dealerFit <- dealer(x = as.matrix(bestTrainX), y = trainY, g = as.factor(trainG), gamma = 2, grdTaus = gridTaus) #, normalized = FALSE

