# mariekeNL_diffSegAnalysis.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Analysis of the segmented delta profiles - 
#              quantile normalized version
# -------------------------------------------------------------------

# WD - medoid
setwd('/home/klijn/data/smallproj/mariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/codeChris/generalFunctionsR/chris_DNAcopy_utils.R')

library(ggplot2)
library(DNAcopy)

load('marieke_diffResults_quant_UD2.Rda')

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

# Load clinical data
sampleInfo <- read.delim('Clin_data_tumorLN.txt', stringsAsFactors=F)
sampleInfo <- sampleInfo[order(sampleInfo$File_name),]
sampleInfo <- sampleInfo[-which(sampleInfo$NR %in% c(592, 456)),]
sampleInfoTumor <- sampleInfo[sampleInfo$Type == 'Tumor',]


outputAll <- extractSeg(segResult=KCsegDiff, minMark=10, cutoff=NULL)

outputAll$subtype <- sampleInfoTumor$Mol_subtype[
  match(outputAll$ID, sampleInfoTumor$NR)]
  
outputFilter <- extractSeg(segResult=KCsegDiff, minMark=10, cutoff=.2, higher='both')
outputFilter$subtype <- sampleInfoTumor$Mol_subtype[
  match(outputFilter$ID, sampleInfoTumor$NR)]

segCount <- getSegCount(outputFilter, outputAll)

postscript(file='Figures/deltaSeg_quant_UD2.eps', width=4, height=6, paper='special', horizontal=F)
  qplot(data=segCount, x=subtype, y=NumSeg, position=position_jitter(w=0.2, h=0), 
    color=subtype, size=I(2)) + opts(title='Number of segments (UD) in delta profile')
dev.off()

# Output all segments in delta profiles

write.table(x=outputFilter, file="deltaSegTable_UD2.txt", col.names=TRUE,  quote=FALSE, row.names=FALSE, sep='\t')

a <- subset(KCsegDiff, samplelist=c('355', '322', '811'))

pdf(file='Figures/segmentDiff_UD2.pdf', width=20, height=20)
plot(KCsegDiff, plot.type='chrombysample', cex.axis=10)
dev.off()


