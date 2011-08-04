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
  
  segCount <- as.data.frame(table(segFrame$ID, segFrame$subtype), stringsAsFactors=F)
  segCount <- segCount[-which(segCount$Freq == 0),]
  colnames(segCount) <- c('TumorID','subtype', 'NumSeg')
  
  # Add missing tumor IDs
  
  AllID <- unique(originalFrame$ID)
  
  missingID <- which(!(AllID %in% segFrame$ID))
  
  for (m in missingID) {
      rowInsert <- data.frame(TumorID=AllID[m], 
        subtype=originalFrame$subtype[originalFrame$ID==AllID[m]][1], 
        NumSeg=0)
      segCount <- rbind(segCount, rowInsert)    
  }
  
  
  return(segCount)

}

outputAll <- extractSeg(segResult=KCsegDiff, minMark=10, cutoff=NULL)
outputAll$subtype <- sampleInfoTumor$Mol_subtype[match(as.numeric(outputAll$ID), sampleInfoTumor$NR)]
  
outputFilter <- extractSeg(segResult=KCsegDiff, minMark=10, cutoff=.2, higher='both')
outputFilter$subtype <- sampleInfoTumor$Mol_subtype[match(as.numeric(outputFilter$ID), sampleInfoTumor$NR)]

segCount <- getSegCount(outputFilter, outputAll)

tempTable <- table(output$ID)
tumorSegStats <- data.frame(ID=names(tempTable), nSeg=as.numeric(tempTable))
tumorSegStats$Class <- sampleInfoTumor$BRCAness[match(tumorSegStats$ID, sampleInfoTumor$NR)]

