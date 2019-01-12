######################################################################
# Data Sharing for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee T. 
# Email: asiae002@umn.edu
######################################################################
library("glmnet")
source("00-paths.R")
source("projectOntoElasticNet.R")
source("projectOntoL1.R")
source("04-dataEnrichmentPrediction.R")


corrThresh <- 0.25
focusedCancerTypes <- c("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "LUNG")
load(file=file.path(paths$scratch, paste("cvDeMeans_", paste(focusedCancerTypes, collapse = "_") ,".RData")))
load(file=file.path(paths$scratch, paste("cvDeSds_", paste(focusedCancerTypes, collapse = "_") ,".RData")))
load(file=file.path(paths$scratch, paste("cvDeTaus", paste(focusedCancerTypes, collapse = "_") ,".RData")))

minDeMeanIndex <- apply(cvDeMeans, 1, which.min)
bestTaus <- sapply(1:24, function(i) cvDeTaus[minDeMeanIndex[i],,i])
nGroup <- dim(bestTaus)[1]

set.seed(123)
dataDir <- file.path(paths$clean, "xy") 
allFiles <- list.files(dataDir)

drugCnt <- 0
for(file in allFiles){
  load(file.path(dataDir,file))
  drugCnt <- drugCnt + 1
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
  runningSumOfImpIndex <- matrix(0L, nrow = p, ncol = nGroup) #-2 for y, colors, and cancer/groupId which is in data
  # runningsumOfCoeffsEffect <- matrix(0L, nrow = p, ncol = nGroup)
  numBootStrap <- 100
  for(i in 1:numBootStrap){
    print(paste("Resample number ", i))
    bsIndex <- sample(1:n, n, replace=TRUE)
    dealerFit <- dealer(x = as.matrix(bestPredictors[bsIndex,]), y = response[bsIndex], g = as.factor(groupsId[bsIndex]), gamma = 1, grdTaus = t(as.matrix(bestTaus[,drugCnt])), normalize = TRUE)
    temp <- 1 * (abs(dealerFit@betaList[[1]][-1,]) > 0.0005) #remove the intercept and very small coeffs.
    runningSumOfImpIndex <- runningSumOfImpIndex + temp
  }
  avgFreqOfImpIndex <- runningSumOfImpIndex / numBootStrap
  # impFeatureIndex <- (avgFreqOfImpIndex >= .8 & avgOfCoeffsEffect > 0.8) #avgOfCoeffsEffect > 0 select those features that increase y on average.
  # print(sum(impFeatureIndex) / 3)
  impFeatureIndex <- (avgFreqOfImpIndex >= .8)
  print(sum(impFeatureIndex) / nGroup)
  
  for(i in 1:nGroup){
    fileConnGroup <- file(file.path(paths$scratch, paste(file, "group",i,".txt")))
    nameOfImpCoef <- names(data)[impFeatureIndex[,i]]
    nameOfImpCoef <- sub("-cna", "", nameOfImpCoef)
    nameOfImpCoef <- sub("-exp", "", nameOfImpCoef)
    nameOfImpCoef <- sub("-hybridMut", "", nameOfImpCoef)
    nameOfImpCoef <- sub("-oncoMut", "", nameOfImpCoef)
    writeLines(nameOfImpCoef, fileConnGroup)
    close(fileConnGroup)
  }
  
  # sortedValueOfImpCoef <- sort(valueOfImpCoef, decreasing = TRUE, index.return=TRUE)
  # print(nameOfImpCoef[sortedValueOfImpCoef$ix])
  
  
  
  
  
}

# 
# dataDir <- file.path(paths$clean, "xy") 
# 
# allFiles <- list.files(dataDir)
# set.seed(123)
# counter <- 5 # 5 is for the Erlotinib, 3 AZD530
# for(file in allFiles){
#   # counter <- counter + 1
#   print(paste("Loading data from file ",  file ,"..."))
#   load(paste(dataDir,file, sep="/"))
#   data <- XyActArea 
#   print(paste("Data loaded!"))
#   
#   
#   print("Starting the bootstrap!")
#   runningSumOfImpIndex <- matrix(0L, nrow = ncol(data)-2, ncol = 3) #-2 for y, colors, and cancer/groupId which is in data
#   runningsumOfCoeffsEffect <- matrix(0L, nrow = ncol(data)-2, ncol = 3)
#   numBootStrap <- 10
#   for(i in 1:numBootStrap){
#     print(paste("Resample number ", i))
#     bootStrapSampleIndex <- sample(1:nrow(data), nrow(data), replace=TRUE)
#     bootStrapData <- data[bootStrapSampleIndex,]
#     bootStrapCv1005 <- myDataSharingSparse(bootStrapData, givenTau = tausMinAll[counter,])
#     # bootStrapCv1005 <- myDataSharingSparse(bootStrapData, givenTau = c(1.9,2.2,3.2))
#     temp <- 1 * (abs(bootStrapCv1005[-nrow(bootStrapCv1005),]) > 0.0005) #remove the intercept
#     runningSumOfImpIndex <- runningSumOfImpIndex + temp
#     
#     # bootStrapCv1005 <- bootStrapCv1005[-nrow(bootStrapCv1005),]
#     # ids <- bootStrapData$cancerId
#     # bootStrapData <- bootStrapData[, !(names(bootStrapData) %in% c("y", "color", "cancerId"))]
#     # meansOfFeaturesForLung <- colMeans(bootStrapData[ids == 1,])
#     # meansOfFeaturesForBlood <- colMeans(bootStrapData[ids == 2,])
#     # bootStrapCv1005[,1] <- (bootStrapCv1005[,1]+bootStrapCv1005[,3]) * meansOfFeaturesForLung
#     # bootStrapCv1005[,2] <- (bootStrapCv1005[,2]+bootStrapCv1005[,3]) * meansOfFeaturesForBlood
#     # bootStrapCv1005[,3] <- rep(1, length(bootStrapCv1005[,3])) #bootStrapCv1005[,3] * colMeans(bootStrapData)
#     # runningsumOfCoeffsEffect <- runningsumOfCoeffsEffect + 1 * (bootStrapCv1005 < 0)
#   }
#   
#   # avgOfCoeffsEffect <- runningsumOfCoeffsEffect / numBootStrap
#   avgFreqOfImpIndex <- runningSumOfImpIndex / numBootStrap
#   # impFeatureIndex <- (avgFreqOfImpIndex >= .8 & avgOfCoeffsEffect > 0.8) #avgOfCoeffsEffect > 0 select those features that increase y on average.
#   # print(sum(impFeatureIndex) / 3)
#   impFeatureIndex <- (avgFreqOfImpIndex >= .8)
#   print(sum(impFeatureIndex) / 3)
#   
#   
#   # cv1005 <- myDataSharingSparse(data, givenTau = tausMinAll[counter,])
#   # valueOfCoef <- as.matrix(cv1005)
#   # valueOfCoef <- valueOfCoef[-nrow(valueOfCoef), ] #remove the intercept
#   # valueOfImpCoef <- valueOfCoef[impFeatureIndex]
#   fileConnGroup1<-file(file.path(paths$scratch, paste(file, "group1.txt")))
#   fileConnGroup2<-file(file.path(paths$scratch, paste(file, "group2.txt")))
#   fileConnShared<-file(file.path(paths$scratch, paste(file, "shared.txt")))
#   data <- data[, !(names(data) %in% c("y", "colors", "cancerId"))]
#   
#   for(i in 1:3){
#     nameOfImpCoef <- names(data)[impFeatureIndex[,i]]
#     nameOfImpCoef <- sub("-cna", "", nameOfImpCoef)
#     nameOfImpCoef <- sub("-exp", "", nameOfImpCoef)
#     nameOfImpCoef <- sub("-hybridMut", "", nameOfImpCoef)
#     nameOfImpCoef <- sub("-oncoMut", "", nameOfImpCoef)
#     
#     if (i == 1){
#       writeLines(nameOfImpCoef, fileConnGroup1)
#     }
#     if (i == 2){
#       writeLines(nameOfImpCoef, fileConnGroup2)
#     }
#     if (i == 3){
#       writeLines(nameOfImpCoef, fileConnShared)
#     }
#   }
#   
#   # sortedValueOfImpCoef <- sort(valueOfImpCoef, decreasing = TRUE, index.return=TRUE)
#   # print(nameOfImpCoef[sortedValueOfImpCoef$ix])
#   
#   
#   close(fileConnShared)
#   close(fileConnGroup1)
#   close(fileConnGroup2)
#   
#   ############################
#   # readline(print("Should I go to the next drug?!"))
# }
# 
# 
# # currentMean <- 0
# # currentSd <- 0
# # minIndex <- 0
# # sds <- c()
# # means <- c()
# # sds <- c()
# # for(i in 1:96){
# #   if(i %% 4 == 1){
# #     currentMean <- listOfRefinedCV[[i]]
# #     minIndex <- which.min(currentMean)
# #   }
# #   else if(i %% 4 == 2){
# #     currentSd <- listOfRefinedCV[[i]]
# #     means <- append(means, currentMean[minIndex])
# #     sds <- append(sds, currentSd[minIndex])
# #   }
# # }
# 
# # 
# # library(ggplot2)
# # gList <- list()
# # apply(t(1:3), 2,
# #       function(i){
# #         toPlot <- as.data.frame(matrix(NA,14,1))
# #         toPlot$means <- apply(expRes[[1]], 2, mean)
# #         toPlot$sds <- apply(expRes[[i]], 2, sd)
# #         index <- i
# #         if(i==3){
# #           toPlot$n <- numSamples
# #           index <- 0
# #         }
# #         else
# #           toPlot$n <- ((3 - (i%%2)) / 5) * numSamples
# #         gList[[i]] <- ggplot(data=toPlot,aes(x = n,y=means)) +
# #           theme_bw() +
# #           geom_errorbar(data=toPlot,aes(ymin=means-sds,ymax=means+sds), width=80, size=.5, color="black")+
# #           geom_line(aes(y=means), color="red", size=.8) +
# #           theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), text = element_text(size=20))+
# #           xlab(substitute(n[nn], list(nn=index))) +
# #           ylab(substitute("||"~delta[nn]~"||=||"~hat(beta)[nn] - beta[nn]^{"*"}~"||", list(nn=index)))
# #       }
# # )
# # 
# 
# 
# 
# # nameDescMapping <- as.data.frame(expressions$Description)
# # row.names(nameDescMapping) <- row.names(expressions)
# # featureNames <- names(data)[1:(ncol(data)-2)]
# # featureNames[impFeatureIndex[,1]]
# # nameDescMapping[featureNames[impFeatureIndex[,1]],]
