######################################################################
# Data Sharing for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee T. 
# Email: asiae002@umn.edu
# Description: Preprocessing
# Run this once and store the objects for future loading.
######################################################################
source("00-paths.R")

######################################################################
# Pre-process gene expression data. 
######################################################################
expressions <- read.csv(file.path(paths$raw,"expression.csv", sep=""), row.names = 1, check.names = FALSE, skip = 2)
IndexOfNoNameGenes <- which(duplicated(expressions[["Description"]])) #These are 90 genes 87 of them are without name, so we set them as AmirGenei
noNameGenes <- expressions[["Description"]][IndexOfNoNameGenes]
expressions[["Description"]] <- lapply(expressions[["Description"]], as.character)
expressions[["Description"]][IndexOfNoNameGenes] <- paste('AmirGene', 1:length(noNameGenes), noNameGenes)
row.names(expressions) <- expressions[["Description"]]
expressions[["Description"]] <- NULL
expressions <- scale(t(expressions))

doseResponse <- read.csv(file.path(paths$raw,"doseResponse.csv", sep=""), check.names =  FALSE)
doseResponse <- doseResponse[,c("Compound", "CellLineName","ActArea")]


oncoMut <- read.csv(file.path(paths$raw,"oncomapmut.csv", sep=""), check.names =  FALSE)
oncoMut <- oncoMut[,c("Hugo_Symbol", "Tumor_Sample_Barcode")]

hybridMut <- read.csv(file.path(paths$raw,"hybridmut.csv", sep=""), check.names =  FALSE)
hybridMut<- hybridMut[,c("Hugo_Symbol", "Tumor_Sample_Barcode")]

# oncoMutDf <- read.csv(paste(dataDir,"oncomapmut.csv", se=p""), check.names = FALSE)
# oncoMutList <- levels(oncoMutDf[["Hugo_Symbol"])]
# hybridMutDf <- read.csv(paste(dataDir,"hybridmut.csv", sep=""), check.names = FALSE)
# hybridMutList <- levels(hybridMutDf[["Hugo_Symbol"]])
copyNumbers <- read.csv(file.path(paths$raw,"copyNumbers.csv", sep=""), check.names =  FALSE)
IndexOfNoNameGenes <- which(duplicated(copyNumbers[["SYMBOL"]])) #These are 2 genes "1-Mar" and "2-Mar" that are problematic, so we set them as AmirGenei - 1-Mar, etc. 
noNameGenes <- copyNumbers[["SYMBOL"]][IndexOfNoNameGenes]
copyNumbers[["SYMBOL"]] <- lapply(copyNumbers[["SYMBOL"]], as.character)
copyNumbers[["SYMBOL"]][IndexOfNoNameGenes] <- paste('AmirGene', 1:length(noNameGenes), noNameGenes)
row.names(copyNumbers) <- copyNumbers[["SYMBOL"]]
copyNumbers[["SYMBOL"]] <- NULL
copyNumbers[["CHR"]] <- NULL
copyNumbers[["CHRLOC"]] <- NULL
copyNumbers[["CHRLOCEN"]] <- NULL
copyNumbers <- scale(t(copyNumbers))



save(expressions, file=file.path(paths$clean,"expressions.RData"))
save(doseResponse, file=file.path(paths$clean,"doseResponse.RData"))
save(oncoMut, file=file.path(paths$clean,"oncoMut.RData")) #Definition of mutations captured by oncoMute
# save(oncoMutList, file=paste(dataDir,"oncoMutList.RData", sep=""))
save(hybridMut, file=file.path(paths$clean,"hybridMut.RData"))
# save(hybridMutList, file=paste(dataDir,"hybridMutList.RData", sep=""))
save(copyNumbers, file=file.path(paths$clean,"copyNumbers.RData"))





