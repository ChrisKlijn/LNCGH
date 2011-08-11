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

sampleInfoTumor <- sampleInfo[sampleInfo$Type == 'Tumor',]

outputAll <- extractSeg(segResult=KCsegDiff, minMark=10, cutoff=NULL)
outputAll$subtype <- sampleInfoTumor$Mol_subtype[match(as.numeric(outputAll$ID), sampleInfoTumor$NR)]
  
outputFilter <- extractSeg(segResult=KCsegDiff, minMark=10, cutoff=.2, higher='both')
outputFilter$subtype <- sampleInfoTumor$Mol_subtype[match(as.numeric(outputFilter$ID), sampleInfoTumor$NR)]

segCount <- getSegCount(outputFilter, outputAll)

postscript(file='Figures/deltaSeg.eps', width=4, height=6, paper='special', horizontal=F)
  qplot(data=segCount, x=subtype, y=NumSeg, position=position_jitter(w=0.2, h=0), 
    color=subtype, size=I(2)) + opts(title='Number of segments in delta profile')
dev.off()

# Output all segments in delta profiles





