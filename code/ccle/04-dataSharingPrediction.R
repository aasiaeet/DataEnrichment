######################################################################
# Data Sharing for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee T. 
# Email: asiae002@umn.edu
######################################################################
source("00-paths.R")
source("dataSharing.R")

source("projectOntoElasticNet.R")
# source("projectOntoL1.R")

dataDir <- file.path(paths$clean, "xy") 

allFiles <- list.files(dataDir)
perDrugBestPerformanceDSSparse <- c()
listOfRefinedCV <- c()
set.seed(123)
counter <- 0
for(file in allFiles){
  if(file == "README.md")
    next

  print(paste("Loading data from file ",  file ,"..."))
  load(paste(dataDir,file, sep="/"))
  data <- XyActArea 
  print(paste("Data loaded!"))
  
  #tau is the set of parameters that we want to test. 
  taus <- seq(from=1, to=5, by=1)
  # One for the shared component and two for each of Lung and Blood cancers. 
  tausGrid <- expand.grid(taus,taus, taus)
  
  # print("We are looking for a good range of parameters ...")
  cv102 <- myDataSharingSparse(data, tausGrid)

  bestTaus <- tausGrid[sort(cv102$cvm, index.return=TRUE)$ix[1:10],]
  thirdQuantile <- apply(bestTaus, 2, quantile)[4,]
  firstQuantile <- apply(bestTaus, 2, quantile)[2,]
  tausMax <- apply(bestTaus, 2, max)
  tausMin <- apply(bestTaus, 2, min)
  
  for(i in 1:3){
    if(max(taus) == thirdQuantile[i]){  # We are on the margine of the search space, let's expand it. 
      tausMax[i] = 10
    } 
    if(min(taus) == firstQuantile[i]){
      tausMin[i] = 0.001
    }
  }
  # print(paste("We found the following range", tausMin, tausMax))
  refinedTaus <- apply(t(1:3), 2, function(i) (seq(from=tausMin[i], to=tausMax[i], length.out = 50)))
  
  # print("We are looking for the best parameters in the range ...")
  cv102Refined <- myDataSharingSparse(data, refinedTaus)

  listOfRefinedCV <- append(listOfRefinedCV, cv102Refined)
  # savedTau.min <- cv102Refined$taus.min
  # Save the best performance.
  perDrugBestPerformanceDSSparse <- append(perDrugBestPerformanceDSSparse, min(cv102Refined$cvm))
  ############################
  # This is for feature extraction. Disabled for now.
  ############################
  #   par(mar=c(4,4,4,4))
  #   plot(log(cv1005$lambda),cv1005$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv1005$name)
  #   print("Done plotting the results!")
  #   
  
  # readline(print("Stop!"))

  # # Checking beta_0, beta_1, beta_2 non-zero component (not that informative)
  # bestBeta <- matrix(coef(cv102Refined), ncol = 3)
  # print(paste("first", bestBeta[which(bestBeta[,1]>0),1], collapse = ""))
  # print(paste("second", bestBeta[which(bestBeta[,2]>0),2], collapse = ""))
  # print(paste("third", bestBeta[which(bestBeta[,3]>0),3], collapse = ""))
  
  # readline()
  print(perDrugBestPerformanceDSSparse)
  # readline(print("Should I go to the next drug?!"))
}
save(perDrugBestPerformanceDSSparse, file=paste(file.path(paths$scratch, "perDrugBestPerformance-DS-Sparse.RData")))

save(listOfRefinedCV, file=paste(file.path(paths$scratch,  "listOfRefinedCV-DS-Sparse.RData", sep="")))


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
