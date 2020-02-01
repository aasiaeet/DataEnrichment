#############################
# Synthetic Experiement for High Dimensional Data Sharing
# Author: asiae002@umn.edu
#############################



#############################
# Other sources.
#############################
source("projectOntoL1.R")
library(MASS)

#############################
# Helper functions. 
#############################
norm2 <- function(x) sqrt(sum(x^2))
norm1 <- function(x) sum(abs(x))
normInf <- function(x) max(abs(x))

#############################
# Parameter setup. 
#############################
set.seed(123)

#############################
# Optimization parameters.
nStepsGD <- 100
stoppingCriteria <- .001 

#############################
# Statistical parameters.
numSamples <- c(seq(from=100, to=1000, by=100))
# expRes <- lapply(1:3, function(x) matrix(NA, 5, length(numSamples)))
p <- 1000 #the dimension
beta0 <- c(rep(1, 20), rep(0, p - 20))
beta1 <- c(rep(0,p))
beta1[51:60] <- 2
beta2 <- c(rep(0,p))
beta2[96:100] <- -2
tau1 <- norm1(beta1)
tau2 <- norm1(beta2)
tau0 <- norm1(beta0)

#############################
# Beta is a p-by-(G+1) matrix where column g corresponds to the gth group.
# The last column is for the common component. 
# The last row is for the intercept. 
Beta <- cbind(beta1, beta2, beta0)
Beta <- rbind(Beta, rep(0, 3))


#############################
# Data Sharing code.
#############################

#############################
# Predict from the fitted model for the new data point x. 
predict.dataSharing <- function(fittedModel, x, groups, taus){
  x <- cbind(x, rep(1, nrow(x))) #intercept
  learnedBeta <- matrix(fittedModel, ncol=3)
  predictedY <- rep(0, nrow(x))
  for(g in 1:2){
    groupIndex <- groups == g
    predictedY[groupIndex] <-  x[groupIndex,] %*% (learnedBeta[,g] + learnedBeta[,3])
  }
  return(predictedY)
}

#############################
# Learn parameters of DS model using the input data. 
dataSharing <- function(x, y, groups, gamma, taus){
  print(paste("Testing taus:", taus))
  commonStepSize <- 1/nrow(x)
  x <- cbind(x, rep(1, nrow(x))) #intercept
  BetaNex <- BetaPrv <- matrix(0, ncol(x), 3)
  marginalImp <- rep(0, length(y))
  for(i in 1:nStepsGD){
    if(i %% 50 == 0){
      print(paste("This is iteration", i, "of PGD"))
    }
    for(g in 1:2){
      groupIndex <- groups == g
      groupStepSize <- 1/sqrt(sum(groupIndex) * nrow(x))  
      marginalImp[groupIndex] <- y[groupIndex] - x[groupIndex,] %*% (BetaPrv[,g] + BetaPrv[,3])
      BetaNex[,g] <-  BetaPrv[,g] + groupStepSize * t(x[groupIndex,]) %*% marginalImp[groupIndex]
      BetaNex[,g] <- projectOntoL1(BetaNex[,g], taus[g])
    }
    BetaNex[,3] <-  BetaPrv[,3] + commonStepSize * t(x) %*% marginalImp
    BetaNex[,3] <- projectOntoL1(BetaNex[,3], taus[3])
    ### Stopping criteria. 
    if(normInf(BetaPrv - BetaNex) < stoppingCriteria){
      print(paste("breaking bad! at step", i))
      break 
    }
    BetaPrv <- BetaNex
  }
  return(BetaNex)
}

#############################
# Cross validation for the DS model
# Here we only seek best taus and ignore gammas: assume gammas are zero and we are only project onto l1 ball. 
cv.dataSharing <- function(predictors, response, groupsId, nfolds=10, gamma=2){
  stopifnot(nrow(predictors) == length(response))
  nTrain <- length(response)
  foldsId <- sample(rep(1:nfolds, length.out = nTrain))
  
  #tau is the set of parameters that we want to test. 
  taus <- seq(from=5, to=25, by=5)
  # One for the shared component and the other two for groups. 
  tausGrid <- expand.grid(taus,taus,taus)
  
  cvResults <- matrix(NA, nrow = nfolds, ncol = nrow(tausGrid))
  for (k in 1:nfolds) {
    print(paste("Fold", k, "started."))
    testId <- which(foldsId == k)
    trainX <- predictors[-testId, ]
    trainY <- response[-testId]
    trainG <- groupsId[-testId]
    testX <- predictors[testId, ]
    testY <- response[testId]
    testG <- groupsId[testId]
    
    fittedModels <- apply(t(tausGrid), 2, function(myTaus) dataSharing(trainX, trainY, trainG, gamma, myTaus))
    # Each col of the fitted model is the unraveled representation of the Beta matrix with length p * (G+1)
    preditions <- apply(t(1:nrow(tausGrid)), 2, function(i) predict.dataSharing(fittedModels[,i], testX, testG, tausGrid[i,]) ) 
    
    cvResults[k, ] <- apply(preditions, 2, function(yHat) mean((testY - yHat)^2))
  }
  cv <- list()
  cv$cvm <- colMeans(cvResults)
  cv$cvsd <- apply(cvResults, 2, sd)
  cv$cvm.min <- min(cv$cvm)
  cv$taus.min <- tausGrid[which.min(cv$cvm),]
  cv$taus.grid <- tausGrid 
  return(cv)
}

########################
# Synthetic experiment.
########################
rounds <- 30
numSamples <- c(seq(from=100, to=1000, by=100))
expRes <- lapply(1:3, function(x) matrix(NA, rounds, length(numSamples)))

for(e in 1:rounds){
  j <- 0
  for(n in numSamples){
    j <- j + 1
    X <- mvrnorm(n, rep(0, p), 0.3 * diag(p))
    epsilon <- rnorm(n, 0, .1)
    n1 <- 2 * (n/5)
    n2 <- n - n1
    y <- rep(0, n)
    y[1:n1] <- X[1:n1,] %*% (beta0 + beta1) + epsilon[1:n1]
    y[(n1+1):n] <- X[(n1+1):n,] %*% (beta0 + beta2) + epsilon[(n1+1):n]
    g <- rep(0, n)
    g[1:n1] <-1
    g[(n1+1):n] <- 2

    # Here we know the actual taus, so no need to cross-validation. 
    inferredBeta <- dataSharing(X, y, g, gamma=0, c(tau1, tau2, tau0))
    print(paste("n=", n, apply(abs(inferredBeta - Beta),  2, function(x) sqrt(sum(x^2)) )))
    currentResults <- apply(abs(inferredBeta - Beta),  2, function(x) sqrt(sum(x^2)) )
    for(i in 1:3)
      expRes[[i]][e,j] <- currentResults[i]
    # readline("pause")
  }
}
save(expRes, file = "../../../scratch/syntheticRes.RData")

library(ggplot2)
gList <- list()
apply(t(1:3), 2,
      function(i){
        toPlot <- as.data.frame(matrix(NA,length(numSamples),1))
        toPlot$means <- apply(expRes[[i]], 2, mean)
        toPlot$sds <- apply(expRes[[i]], 2, sd)
        index <- i
        if(i==3){
          toPlot$n <- numSamples
          index <- 0
        }
        else
          toPlot$n <- ((3 - (i%%2)) / 5) * numSamples
        gList[[i]] <- ggplot(data=toPlot,aes(x = n,y=means)) +
          theme_bw() +
          geom_errorbar(data=toPlot,aes(ymin=means-sds,ymax=means+sds), width=60, size=.5, color="gray")+
          geom_line(aes(y=means), color="red", size=1) +
          theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), text = element_text(size=20))+
          xlab(substitute(n[nn], list(nn=index))) +
          ylab(substitute("||"~delta[nn]~"||=||"~hat(beta)[nn] - beta[nn]^{"*"}~"||", list(nn=index)))
      }
)


#############################
# Main part of the synthetic experiment. 
#############################
# for(n in numSamples){
#   j <- j + 1
#   X <- mvrnorm(n, rep(0, p), 0.3 * diag(p))
#   epsilon <- rnorm(n, 0, .1)
#   n1 <- 2 * (n/5)
#   n2 <- n - n1
#   y <- rep(0, n)
#   y[1:n1] <- X[1:n1,] %*% (beta0 + beta1) + epsilon[1:n1]
#   y[(n1+1):n] <- X[(n1+1):n,] %*% (beta0 + beta2) + epsilon[(n1+1):n]
#   g <- rep(0, n)
#   g[1:n1] <-1 
#   g[(n1+1):n] <- 2
#   learned.cv <-cv.dataSharing(X, y, g, nfolds=10, gamma=0)
#   
# #   inferredBeta <- dataSharing(X, y, g, gamma=0.01, c(tau1, tau2, tau0))
# #   print(paste("n=", n, apply(abs(inferredBeta - Beta),  2, function(x) sqrt(sum(x^2)) )))
# #   currentResults <- apply(abs(inferredBeta - Beta),  2, function(x) sqrt(sum(x^2)) )
#   # readline("pause")
# }
# save(expRes, file = "syntheticRes-CV.RData")




