######################################################################
# Data Sharing for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee T. 
# Email: asiae002@umn.edu
######################################################################
dataDir <- "../../../data/ccle/clean/processedRaw2RData/"
outputDir <- "../../../data/ccle/clean/xy/"
cutOffCnt <- 50

###################################
# Preprocessing
# Run this once and store the objects for future loading.
###################################
# expressions <- read.csv(paste(dataDir,"expression.csv", sep=""), row.names = 1, check.names = FALSE, skip = 2)
# expressions[["Description"]] <- NULL
# expressions <- scale(t(expressions))
# doseResponse <- read.csv(paste(dataDir,"doseResponse.csv", sep=""), check.names =  FALSE)
# oncoMutDf <- read.csv(paste(dataDir,"oncomapmut.csv", sep=""), check.names = FALSE)
# oncoMutList = levels(oncoMutDf[["Hugo_Symbol"]])
# hybridMutDf <- read.csv(paste(dataDir,"hybridmut.csv", sep=""), check.names = FALSE)
# hybridMutList = levels(hybridMutDf[["Hugo_Symbol"]])
# save(expressions, file=paste(dataDir,"expressions.RData", sep=""))
# save(doseResponse, file=paste(dataDir,"doseResponse.RData", sep=""))
# save(oncoMutDf, file=paste(dataDir,"oncoMutDf.RData", sep="")) #Definition of mutations captured by oncoMute
# save(oncoMutList, file=paste(dataDir,"oncoMutList.RData", sep=""))
# save(hybridMutDf, file=paste(dataDir,"hybridMutDf.RData", sep=""))
# save(hybridMutList, file=paste(dataDir,"hybridMutList.RData", sep=""))

###################################
# Loading processed data
###################################
load(paste(dataDir,"expressions.RData", sep=""))
load(paste(dataDir,"doseResponse.RData", sep=""))
load(paste(dataDir,"oncoMutDf.RData", sep=""))
load(paste(dataDir,"oncoMutList.RData", sep=""))
load(paste(dataDir,"hybridMutDf.RData", sep=""))
load(paste(dataDir,"hybridMutList.RData", sep=""))


# Run 24 experiments, one per drug. 
drugs <- doseResponse[["Compound"]]
for(drug in levels(drugs)){
  #Conceptually cell line names are not factor, so we keep them as character. 
  ccleNames <- as.character(doseResponse[drugs == drug, "﻿CCLE Cell Line Name"])
  ccleNames <- intersect(row.names(expressions), ccleNames)
  ###################################
  #Extract x: expression part
  ###################################
  XExpression <- expressions[ccleNames,]
  
  ###################################
  #Extract x: oncomap/hybrid mutation part
  ###################################
  XOncoMut <- data.frame(matrix(0, nrow=length(ccleNames), ncol=length(oncoMutList)), row.names = ccleNames, check.rows = FALSE)
  XHybridMut <- data.frame(matrix(0, nrow=length(ccleNames), ncol=length(hybridMutList)), row.names = ccleNames, check.rows = FALSE)  
  
  #Column names can not be set in the data.frame constructor. 
  names(XOncoMut) <- oncoMutList 
  names(XHybridMut) <- hybridMutList 
  
  
  for(name in ccleNames){
    #There mayebe multiple mutation per gene so we need unique. 
    #Keep indices as character, factor indecies are converted to number.  
    XOncoMut[name, as.character(unique(oncoMutDf[oncoMutDf[["Tumor_Sample_Barcode"]] == name, "Hugo_Symbol"]))] <- 1
    XHybridMut[name, as.character(unique(hybridMutDf[hybridMutDf[["Tumor_Sample_Barcode"]] == name, "Hugo_Symbol"]))] <- 1
  }
  #remove features (mutations) that didn't happen in the whole data set.
  XOncoMut <- XOncoMut[, colSums(XOncoMut) != 0]  
  XHybridMut <- XHybridMut[, colSums(XHybridMut) != 0]
  
  ###################################
  #Bind x parts to make glmnet input
  ###################################
  X <- merge(as.data.frame(XOncoMut), as.data.frame(XHybridMut), by='row.names', all=TRUE, sort=TRUE, suffixes=c("-oncoMut","-hybridMut"))
  #Clean up the mess of merge for the next merge!
  row.names(X) <- X[["Row.names"]]
  X[["Row.names"]] <- NULL
  
  X <- merge(as.data.frame(X), as.data.frame(XExpression), by='row.names', all=TRUE, sort=TRUE, suffixes=c("", "-exp"))
  row.names(X) <- X[["Row.names"]]
  X[["Row.names"]] <- NULL
  X <- X[complete.cases(X), ]
  X <- as.data.frame(scale(X))
  completeCasesCL <- row.names(X)
  ###################################
  #Extract y
  ###################################
  completeCasesDrugIndex <- (doseResponse$Compound == drug) & (doseResponse$`﻿CCLE Cell Line Name` %in% completeCasesCL)
  yActArea <- doseResponse[completeCasesDrugIndex , "ActArea"]
  names(yActArea) <- doseResponse[completeCasesDrugIndex , "﻿CCLE Cell Line Name"]
  
  print("Removing less relevant features ...")
  correlations <- sapply(X, cor, y=yActArea) 
  mask <- (abs(correlations) >= .1)
  best.X <- X[,mask]  
  print("Features are ready!")
  
  
  XyActArea <- merge(best.X, yActArea, by='row.names', all=TRUE, sort=TRUE)
  row.names(XyActArea) <- XyActArea[["Row.names"]]
  XyActArea[["Row.names"]] <- NULL
  
  ##################################
  # Remove cancers with less than 10 cell lines. 
  ##################################
  # Putting similar cancers consecutively.
  rownames(XyActArea) <- lapply(rownames(XyActArea), function(x) paste(rev(strsplit(x, NULL)[[1]]), collapse = ''))
  XyActArea <- XyActArea[ order(row.names(XyActArea)), ]
  rownames(XyActArea) <- lapply(rownames(XyActArea), function(x) paste(rev(strsplit(x, NULL)[[1]]), collapse = ''))
  XyActArea$cancerId <- rep(0, nrow(XyActArea)) #ids start from 1.
  
  
  # saveXyAct <- XyActArea
  
  cellLineNames <- rownames(XyActArea)
  separatorPlaces <- gregexpr(pattern = '_', cellLineNames)
  
  cancerIdCnt <- 1 #running counter of group id for cancers.
  startRange <- 0
  
  removeIndex <- c()
  oldCellLineName <- ""
  # Counting similar cancers. 
  for(i in 1:nrow(XyActArea)){
    #print(paste("Working on ", cellLineNames[i]))
    #First occurance of the _ till the end is the cancer name. 
    newCellLineName <- substr(cellLineNames[i], start=separatorPlaces[[i]][1]+1, stop=nchar(cellLineNames[i]))
    if (newCellLineName == "Lung"){
      XyActArea$cancerId[i] <- 1
    }
    else if (newCellLineName == "OVARY"){
      XyActArea$cancerId[i] <- 2
    }
    else{
      XyActArea$cancerId[i] <- 0
    }
    ### Commented for now. This is for the general case. 
#     if (newCellLineName != oldCellLineName){
#       #print(paste("Done reading ", oldCellLineName, "from", startRange, "to", i-1))
#       if(i > 1){
#         cnt <- i - startRange 
#         if(cnt < cutOffCnt){
#           removeIndex <- c(removeIndex, startRange:(i-1))
#         }
#         else{
#           XyActArea$cancerId[startRange:(i-1)] <- cancerIdCnt
#           cancerIdCnt <- cancerIdCnt + 1
#         }
#       }
#       oldCellLineName <- newCellLineName
#       startRange <- i
#     }
  }
  #Currently we only pick lung and blood
  XyActArea <- XyActArea[!(XyActArea$cancerId == 0), ]
#   XyActArea <- XyActArea[-removeIndex, ]
  
  save(XyActArea, file=paste(outputDir, "XyActArea-", drug, ".RData", sep=""))
  print(paste("Drug", drug, "is processed."))
  # readline("pause")
}



