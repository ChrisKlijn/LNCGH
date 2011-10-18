# mariekeNL_overlapSeg.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Instead of using the delta profiles, find non-
#              overlapping segments in between the tumor and lymph
#              node metas
# -------------------------------------------------------------------

# WD - medoid
setwd('/home/klijn/data/smallproj/mariekeLN/')

library(DNAcopy)
library(GenomicRanges)
library(multicore)
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/codeChris/generalFunctionsR/chris_delta_functions.R')
source('~/codeChris/smallProjects/MariekeLN/mariekeNL_overlapSeg_functions.R')

# Load KC data
load('~/data/smallproj/mariekeLN/DNACopy/mariekeData.Rda')
# Load segmentation results
load('marieke_segResult.Rda')

# Body Code

# Load clinical data
sampleInfo <- read.delim('Clin_data_tumorLN.txt', stringsAsFactors=F)
sampleInfo <- sampleInfo[order(sampleInfo$File_name),]

# Remove samples 592 and 456

KC <- KC[,c('chrom', 'maploc', sampleInfo$File_name[which(!sampleInfo$NR %in% c(592, 456))])]
KCseg <- subset(KCseg, samplelist=sampleInfo$File_name[which(!sampleInfo$NR %in% c(592, 456))])
sampleInfo <- sampleInfo[-which(sampleInfo$NR %in% c(592, 456)),]

# Check if the samples are ordered the same in the dataframe, the sampleinfo and the segdata
all.equal(colnames(KC)[3:ncol(KC)], sampleInfo$File_name)
all.equal(names(KCseg$data)[3:ncol(KC)], sampleInfo$File_name)
all.equal( unique(KCseg$output$ID), sampleInfo$File_name)

# Testing

probeRanges <- GRanges(seqnames = Rle(paste('chr', KC$chrom)),
  ranges=IRanges(start=KC$maploc, width=60))

tumSamp <- sampleInfo$File_name[sampleInfo$NR == 44 & 
  sampleInfo$Type == 'Tumor']
LNSamp <- sampleInfo$File_name[sampleInfo$NR == 44 & 
  sampleInfo$Type == 'LN']

pairRanges <- getPairRanges(tumSamp=tumSamp, LNSamp=LNSamp, KCseg=KCseg)
pairRanges.remove <- lapply(pairRanges, removeSegments, 
  thres=.2, minProbes=10)
pairRanges.spec <- removeOverlappingSegments(pairRanges.remove)

patientList <- as.list(unique(sampleInfo$NR))
names(patientList) <- paste('P', unique(sampleInfo$NR), sep='')

segList <- mclapply(patientList, nonOverlapWrapper, KCseg, sampleInfo)






par(mfrow=c(2,4))
for (i in 1:length(tumSpecificRanges)) {
  plotProbes <- matchMatrix(findOverlaps(subject=tumSpecificRanges[i],
    query=probeRanges))[,'query']
  plot(KC$maploc[plotProbes], KC[plotProbes,tumSamp], ylim=c(-2,2), 
    pch=19, cex=.5, col='blue')
  points(KC$maploc[plotProbes], KC[plotProbes,LNSamp], ylim=c(-2,2), 
    pch=19, cex=.5, col='orange')
}

png(file="Figures/testOverlapSeg.png", width=2000, height=1600)
par(mfrow=c(4,5))
for (i in 1:length(LNSpecificRanges)) {
  plotProbes <- matchMatrix(findOverlaps(subject=LNSpecificRanges[i],
    query=probeRanges))[,'query']
  plot(KC$maploc[plotProbes], KC[plotProbes,tumSamp], ylim=c(-2,2), 
    pch=19, cex=.6, col=rgb(0, 0, 1, .4), type='n')
  # points(KC$maploc[plotProbes], KC[plotProbes,LNSamp], ylim=c(-2,2), 
  #   pch=19, cex=.6, col=rgb(1, .65, 0, .4))
  lines(KC$maploc[plotProbes], runmed(x=KC[plotProbes,tumSamp], 
    k=min(50, length(plotProbes))), 
    col=rgb(0, 0, 1, 1), lwd=3)
  lines(KC$maploc[plotProbes], runmed(x=KC[plotProbes,LNSamp], 
     k=min(50, length(plotProbes))), 
    col=rgb(1, .65, 0, 1), lwd=3)
}
dev.off()

