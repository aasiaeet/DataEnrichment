######################################################################
# Data Enrichment for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee
# Email: asiae002@umn.edu
######################################################################

source("00-paths.R")

######################################################################
# Pre-process cell line data. 
######################################################################


if (!file.exists(file.path(paths$clean,"expressions.RData"))) {
  expressions <- read.csv(file.path(paths$raw,"expression.csv", sep=""), row.names = 1, check.names = FALSE, skip = 2)
  IndexOfNoNameGenes <- which(duplicated(expressions[["Description"]])) #These are 90 genes 87 of them are without name, so we set them as GGenei
  noNameGenes <- expressions[["Description"]][IndexOfNoNameGenes]
  expressions[["Description"]] <- lapply(expressions[["Description"]], as.character)
  expressions[["Description"]][IndexOfNoNameGenes] <- paste('GGene', 1:length(noNameGenes), noNameGenes)
  row.names(expressions) <- expressions[["Description"]]
  expressions[["Description"]] <- NULL
  # expressions <- scale(t(expressions))
  expressions <- t(expressions)
  save(expressions, file=file.path(paths$clean,"expressions.RData"))
}

if (!file.exists(file.path(paths$clean,"doseResponse.RData"))) {
  doseResponse <- read.csv(file.path(paths$raw,"doseResponse.csv", sep=""), check.names =  FALSE)
  doseResponse <- doseResponse[,c("Compound", "CellLineName","ActArea")]
  save(doseResponse, file=file.path(paths$clean,"doseResponse.RData"))
}

if (!file.exists(file.path(paths$clean,"oncoMut.RData"))) {
  oncoMut <- read.csv(file.path(paths$raw,"oncomapmut.csv", sep=""), check.names =  FALSE)
  oncoMut <- oncoMut[,c("Hugo_Symbol", "Tumor_Sample_Barcode")]
  save(oncoMut, file=file.path(paths$clean,"oncoMut.RData")) #Definition of mutations captured by oncoMute
}  

if (!file.exists(file.path(paths$clean,"hybridMut.RData"))) {
  hybridMut <- read.csv(file.path(paths$raw,"hybridmut.csv", sep=""), check.names =  FALSE)
  hybridMut<- hybridMut[,c("Hugo_Symbol", "Tumor_Sample_Barcode")]
  save(hybridMut, file=file.path(paths$clean,"hybridMut.RData"))
}

if (!file.exists(file.path(paths$clean,"copyNumbers.RData"))) {
  copyNumbers <- read.csv(file.path(paths$raw,"copyNumbers.csv", sep=""), check.names =  FALSE)
  IndexOfNoNameGenes <- which(duplicated(copyNumbers[["SYMBOL"]])) #These are 2 genes "1-Mar" and "2-Mar" that are problematic, so we set them as Genei - 1-Mar, etc. 
  noNameGenes <- copyNumbers[["SYMBOL"]][IndexOfNoNameGenes]
  copyNumbers[["SYMBOL"]] <- lapply(copyNumbers[["SYMBOL"]], as.character)
  copyNumbers[["SYMBOL"]][IndexOfNoNameGenes] <- paste('Gene', 1:length(noNameGenes), noNameGenes)
  row.names(copyNumbers) <- copyNumbers[["SYMBOL"]]
  copyNumbers[["SYMBOL"]] <- NULL
  copyNumbers[["CHR"]] <- NULL
  copyNumbers[["CHRLOC"]] <- NULL
  copyNumbers[["CHRLOCEN"]] <- NULL
  # copyNumbers <- scale(t(copyNumbers))
  copyNumbers <- t(copyNumbers) #don't scale until the end
  save(copyNumbers, file=file.path(paths$clean,"copyNumbers.RData"))
}





