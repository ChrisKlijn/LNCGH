# mariekeNL_diffSegAnalysis.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Analysis of the segmented delta profiles
# -------------------------------------------------------------------

# WD - medoid
setwd('/home/klijn/data/smallproj/mariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/codeChris/generalFunctionsR/chris_DNAcopy_utils.R')

load('diffDataMarieke.Rda')
load('marieke_diffResults.Rda')

# Functions

getSegCount <- function (segFrame, originalFrame) {
  
  segCount <- as.data.frame(table(segFrame$ID, segFrame$genotype), stringsAsFactors=F)
  segCount <- segCount[-which(segCount$Freq == 0),]
  colnames(segCount) <- c('TumorID','Genotype', 'NumSeg')
  
  # Add missing tumor IDs
  
  AllID <- unique(originalFrame$ID)
  
  missingID <- which(!(AllID %in% segFrame$ID))
  
  for (m in missingID) {
      rowInsert <- data.frame(TumorID=AllID[m], 
        Genotype=originalFrame$genotype[originalFrame$ID==AllID[m]][1], 
        NumSeg=0)
      segCount <- rbind(segCount, rowInsert)    
  }
  
  
  return(segCount)

}

sampleInfoTumor <- sampleInfo[sampleInfo$Type == 'Tumor',]
output <- KCsegDiff$output
output <- output[-which(output$num.mark < 10),] 
output <- output[-which(abs(output$seg.mean) < .2),]



tempTable <- table(output$ID)
tumorSegStats <- data.frame(ID=names(tempTable), nSeg=as.numeric(tempTable))
tumorSegStats$Class <- sampleInfoTumor$BRCAness[match(tumorSegStats$ID, sampleInfoTumor$NR)]

