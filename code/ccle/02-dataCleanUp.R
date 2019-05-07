######################################################################
# Data Enrichment for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee
# Email: asiae002@umn.edu
######################################################################
library("ggplot2")

###################################
# Loading processed data
###################################
loadIfNotExists <- function(dataName){
  if(!exists(dataName)) {
    dataName <- paste(dataName, ".RData", sep = "")
    load(file.path(paths$clean,dataName), parent.env(environment()))
  }  
}
sapply(c("expressions","doseResponse","oncoMut","hybridMut","copyNumbers"), loadIfNotExists)
dir.create(file.path(paths$clean, "xy"), showWarnings = FALSE)


# Run 24 times, one per drug. 
drugs <- doseResponse[["Compound"]]
for(drug in levels(drugs)){
  if(!file.exists(file.path(paths$clean,"xy", paste("XyActArea-", drug, ".RData")))){
    # Filtering celllines with missing features. 
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
    
    XyActArea <- merge(X, yActArea, by='row.names', all=TRUE, sort=TRUE)
    row.names(XyActArea) <- XyActArea[["Row.names"]]
    XyActArea[["Row.names"]] <- NULL
    # XyActArea <- as.data.frame(scale(XyActArea))
    XyActArea <- as.data.frame(XyActArea)
    
    ##################################
    # Setting cancer types. 
    ##################################
    cellLineNames <- rownames(XyActArea)
    separatorPlaces <- gregexpr(pattern = '_', cellLineNames)
    
    # Filling out the cancer id column for similar cancers. 
    for(i in 1:nrow(XyActArea)){
      #First occurance of the _ till the end is the cancer name. 
      XyActArea$cancerType[i] <- substr(cellLineNames[i], start=separatorPlaces[[i]][1]+1, stop=nchar(cellLineNames[i]))
    }     
    
    
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
    
    save(XyActArea, file=file.path(paths$clean,"xy", paste("XyActArea_", drug, ".RData", sep="")))
    
    print(paste("Drug", drug, "is processed."))
    # readline("pause")
  }
}


