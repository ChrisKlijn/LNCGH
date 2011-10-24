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
source('~/codeChris/smallProjects/MariekeLN/mariekeNL_overlapSeg_functions.R')

library(ggplot2)
library(DNAcopy)
library(GenomicRanges)

load('marieke_diffResults_quant_UD2.Rda')
load('DNACopy/mariekeData.Rda')

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

outputAll <- outputAll[order(outputAll$chrom, outputAll$loc.start),]
outputAll$subtype <- sampleInfoTumor$Mol_subtype[
  match(outputAll$ID, sampleInfoTumor$NR)]
  
outputFilter <- extractSeg(segResult=KCsegDiff, minMark=10, cutoff=.2, higher='both')
outputFilter <- outputFilter[order(outputFilter$chrom, outputFilter$loc.start),]
outputFilter$subtype <- sampleInfoTumor$Mol_subtype[
  match(outputFilter$ID, sampleInfoTumor$NR)]

segCount <- getSegCount(outputFilter, outputAll)

postscript(file='Figures/deltaSeg_quant_UD2.eps', width=4, height=6, paper='special', horizontal=F)
  qplot(data=segCount, x=subtype, y=NumSeg, position=position_jitter(w=0.2, h=0), 
    color=subtype, size=I(2)) + opts(title='Number of segments (UD) in delta profile')
dev.off()

# Output all segments in delta profiles

write.table(x=outputFilter, file="deltaSegTable_UD2.txt", 
  col.names=TRUE, quote=FALSE, row.names=FALSE, sep='\t')

pdf(file='Figures/segmentDiff_UD2.pdf', width=20, height=20)
plot(KCsegDiff, plot.type='chrombysample', cex.axis=10)
dev.off()

# Calculate overlaps

filterSegs <- with(outputFilter, GRanges(
  seqnames = Rle(paste('chr', chrom)),
  ranges=IRanges(start=loc.start,end=loc.end), 
  seg.mean=seg.mean, 
  num.mark=num.mark,
  ID=ID))

overlMat <- as.matrix(findOverlaps(filterSegs))
overlMat <- overlMat[(overlMat[,1] - overlMat[,2]) != 0,]
filterSegsOverlap <- filterSegs[unique(overlMat[,1])]
outputFilterOverlap <- outputFilter[unique(overlMat[,1]),]

for (chr in unique(outputFilterOverlap$chrom)) {

  chrsubset <- subset(outputFilterOverlap, chrom==chr)
  uniqID <- unique(chrsubset$ID)
  IDcolors <- rich.colors(n=length(uniqID))
  names(IDcolors) <- unique(chrsubset$ID)

  png(file=paste('Figures/overlSeg', chr, '.png', sep=''), 
    width=1000, height=700)

  layout(rbind(matrix(c(
    seq(1,length(uniqID)),
    seq(1,length(uniqID))+length(uniqID)), nrow=2), 
    rep(length(uniqID)*2+1, times=length(uniqID))))

  for (n in uniqID) {
    
    IDsubs <- subset(chrsubset, ID==n)
    plotSegsID(IDsubs, KC, Type='Tumor', IDcolors)
    plotSegsID(IDsubs, KC, Type='LN', IDcolors)

  }

  plotSegsPerChrom(outputFilterOverlap, chr, KCdiff)

  dev.off()

}

testOverlap <- function (segM) {
    if (segM > .2) {
      return('gain')
    }
    else {
      if (segM < -.2) {
        return('loss')
      }
      else {
        return('none')
      }
    }
  }

for (i in 1:nrow(outputFilterOverlap)) {

  tumFile <- sampleInfo$File_name[sampleInfo$NR == 
    outputFilterOverlap$ID[i] & sampleInfo$Type == 'Tumor']
  LNFile <- sampleInfo$File_name[sampleInfo$NR == 
  outputFilterOverlap$ID[i] & sampleInfo$Type == 'LN']
  segProbes <- KC$maploc > outputFilterOverlap$loc.start[i] &
    KC$maploc < outputFilterOverlap$loc.end[i] &
    KC$chrom == outputFilterOverlap$chrom[i]
  outputFilterOverlap$tumState[i] <- 
    testOverlap(mean(KC[segProbes, tumFile]))
  outputFilterOverlap$LNState[i] <- 
    testOverlap(mean(KC[segProbes, LNFile]))

}

outputFilterOverlap$combineState <- paste(outputFilterOverlap$tumState,
  outputFilterOverlap$LNState, sep='->')

write.table(x=outputFilterOverlap, file="deltaSegOvelapTable_UD2.txt", 
  col.names=TRUE, quote=FALSE, row.names=FALSE, sep='\t')

uniqueChromIDCombs <- !duplicated(paste(outputFilterOverlap$ID, outputFilterOverlap$chrom, sep='-'))

overlapChangeTable <- table(outputFilterOverlap$combineState[
  uniqueChromIDCombs], outputFilterOverlap$chrom[uniqueChromIDCombs])

write.table(x=overlapChangeTable, file="overlapChangeTable.txt", 
  col.names=TRUE,  quote=FALSE, row.names=T, sep='\t')

