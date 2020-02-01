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


resRghl <- c()
for(file in allFiles){
  load(file.path(dataDir,file))
  data <- XyActArea 
  # Selecting cancer types
  data <- data[(data$cancerType %in% focusedCancerTypes), ] 
  response <- data[,"y", drop=TRUE]
  rg <- range(response)
  resRghl[file] <- max(rg) - min(rg)
  drugName <- strsplit(file, "_")[[1]][2]
  drugName <- strsplit(drugName, "[.]")[[1]][1]
  print(paste("Data for", drugName, "loaded!"))
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