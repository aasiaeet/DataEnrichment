######################################################################
# Data Sharing for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee T. 
# Email: asiae002@umn.edu
######################################################################
library("ggplot2")
source("00-paths.R")

###################################
# Pre-processing the data
# Comment it if already done once. 
###################################
# source("01-preprocess.R")

# cutOffCnt <- 50



###################################
# Loading processed data
###################################
load(file.path(paths$clean,"expressions.RData"))
load(file.path(paths$clean,"doseResponse.RData"))
load(file.path(paths$clean,"oncoMut.RData"))
load(file.path(paths$clean,"hybridMut.RData"))
load(file.path(paths$clean,"copyNumbers.RData"))

# Run 24 experiments, one per drug. 
drugs <- doseResponse[["Compound"]]
for(drug in levels(drugs)){
  if(as.character(drug) != "AZD0530"){
    next 
  }
  # Conceptually cell line names are not factor, so we keep them as character.
  # Here we keep the cell lines for which we have all of the informations. 
  # Filtering out all of the others early on makes the code more efficient. 
  ccleNames <- as.character(doseResponse[drugs == drug, "CellLineName"])
  ccleNames <- intersect(row.names(expressions), ccleNames)
  ccleNames <- intersect(row.names(copyNumbers), ccleNames)
  
  ###################################
  #Extract x: expression and cpNumber parts
  ###################################
  XExpression <- expressions[ccleNames,]
  XCopyNumber <- copyNumbers[ccleNames,]
  
  ###################################
  #Extract x: oncomap/hybrid mutation part
  ###################################
  
  # Initialization
  oncoMutList <- levels(oncoMut[["Hugo_Symbol"]])
  hybridMutList <- levels(hybridMut[["Hugo_Symbol"]])
  XOncoMut <- data.frame(matrix(0, nrow=length(ccleNames), ncol=length(oncoMutList)), row.names = ccleNames, check.rows = FALSE)
  XHybridMut <- data.frame(matrix(0, nrow=length(ccleNames), ncol=length(hybridMutList)), row.names = ccleNames, check.rows = FALSE)   #Column names can not be set in the data.frame constructor. 
  names(XOncoMut) <- oncoMutList 
  names(XHybridMut) <- hybridMutList 
  
  # Recording the mutations
  for(name in ccleNames){
    #There mayebe multiple mutation per gene so we need unique. 
    #Keep indices as character because factor indecies are converted to number.  
    XOncoMut[name, as.character(unique(oncoMut[oncoMut[["Tumor_Sample_Barcode"]] == name, "Hugo_Symbol"]))] <- 1
    XHybridMut[name, as.character(unique(hybridMut[hybridMut[["Tumor_Sample_Barcode"]] == name, "Hugo_Symbol"]))] <- 1
  }
  #remove features (mutations) that didn't happen in the whole data set.
  XOncoMut <- XOncoMut[, colSums(XOncoMut) != 0]  
  XHybridMut <- XHybridMut[, colSums(XHybridMut) != 0]
  
  ###################################
  #Bind x parts to make the input
  ###################################
  X <- merge(as.data.frame(XOncoMut), as.data.frame(XHybridMut), by='row.names', all=TRUE, sort=TRUE, suffixes=c("-oncoMut","-hybridMut"))
  #Clean up the mess of merge for the next merge!
  row.names(X) <- X[["Row.names"]]
  X[["Row.names"]] <- NULL
  
  X <- merge(as.data.frame(X), as.data.frame(XExpression), by='row.names', all=TRUE, sort=TRUE, suffixes=c("", "-exp"))
  row.names(X) <- X[["Row.names"]]
  X[["Row.names"]] <- NULL
  
  X <- merge(as.data.frame(X), as.data.frame(XCopyNumber), by='row.names', all=TRUE, sort=TRUE, suffixes=c("", "-cna"))
  row.names(X) <- X[["Row.names"]]
  X[["Row.names"]] <- NULL

  X <- X[complete.cases(X), ]
  completeCasesCL <- row.names(X)
  ###################################
  #Extract y
  ###################################
  completeCasesDrugIndex <- (doseResponse$Compound == drug) & (doseResponse$CellLineName %in% completeCasesCL)
  yActArea <- doseResponse[completeCasesDrugIndex , "ActArea"]
  names(yActArea) <- doseResponse[completeCasesDrugIndex , "CellLineName"]
  
  print("Removing less relevant features ...")
  correlations <- sapply(X, cor, y=yActArea) 
  mask <- (abs(correlations) >= .1)
  best.X <- X[,mask]  
  print("Features are ready!")
  
  
  XyActArea <- merge(best.X, yActArea, by='row.names', all=TRUE, sort=TRUE)
  row.names(XyActArea) <- XyActArea[["Row.names"]]
  XyActArea[["Row.names"]] <- NULL
  XyActArea <- as.data.frame(scale(XyActArea))
  ##################################
  # Remove less frequent cancers
  # Keeping only Ovarian, Lung, and Blood
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
    if(newCellLineName == "LUNG"){
      XyActArea$cancerId[i] <- 1
    }
    else if (newCellLineName == "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"){# (newCellLineName == "OVARY"){
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
  
  # XyActArea$colors[XyActArea$cancerId == 1] <- "Lung"
  # XyActArea$colors[XyActArea$cancerId == 2] <- "Blood"
  # # pdf(file=file.path(paths$scratch, paste(drug, "Response-Lung-Blood.pdf", sep="-")))
  # # ggplot(XyActArea, aes(y, fill = colors)) + geom_histogram(binwidth = 0.15, alpha = 0.5, aes(y = ..density..), position = 'identity') +
  # #   ggtitle(drug)
  # ggplot(XyActArea, aes(y, fill = colors)) + geom_density(alpha = 0.5, aes(y = ..density..)) +
  #   ggtitle("Saracatinib") +
  #   xlab("Normalized Activity Area") +
  #   theme_bw() +
  #   ylab("Density of Responded Patients") +
  #   theme(legend.position = c(.85, .8)) 
  # # dev.off()
  # ggsave(file.path(paths$scratch, paste(drug, "Response-Lung-Blood.pdf", sep="-")))
  # 
  
  save(XyActArea, file=file.path(paths$clean,"xy", paste("XyActArea-", drug, ".RData")))
  
  print(paste("Drug", drug, "is processed."))
  # readline("pause")
}


