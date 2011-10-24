# mariekeNL_pairNorm_quantNorm.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Do quantile norm to create delta profiles and
#              segment those.
# -------------------------------------------------------------------

# WD - medoid
setwd('/home/klijn/data/smallproj/mariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/codeChris/generalFunctionsR/chris_delta_functions.R')

# Load KC data
load('~/data/smallproj/mariekeLN/DNACopy/mariekeData.Rda')

# Load segmentation results
load('marieke_segResult.Rda')

library(KCsmart)
library(preprocessCore)
library(multicore)
library(DNAcopy)
data(hsMirrorLocs)

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

uniqPatients <- unique(sampleInfo$NR)
sampleCombs <- matrix('', nrow=length(uniqPatients), ncol=2)

for (i in 1:length(uniqPatients)) {
  
  sampleCombs[i, 1] <- sampleInfo$File_name[
    sampleInfo$NR == uniqPatients[i] &
    sampleInfo$Type =='LN']
  sampleCombs[i, 2] <- sampleInfo$File_name[
    sampleInfo$NR == uniqPatients[i] &
    sampleInfo$Type =='Tumor']

}

KCdiff <- deltaQuant(KC, sampleCombs)

CNA.KCdiff<- CNA(as.matrix(KCdiff[,3:ncol(KCdiff)]), KCdiff$chrom, KCdiff$maploc, data.type=c("logratio"), 
  sampleid=colnames(KCdiff)[3:ncol(KCdiff)])
CNA.KCdiff.smoothed <- smooth.CNA(CNA.KCdiff)

#KCsegDiff <- segment(CNA.KCdiff.smoothed, verbose=1, 
#  undo.splits='sdundo')
KCsegDiff <- segment(CNA.KCdiff.smoothed, verbose=1, 
  undo.splits='sdundo',undo.SD=2)

names(KCsegDiff$data) <- gsub('X', '', names(KCsegDiff$data))
KCsegDiff$output$ID <- gsub('X', '', KCsegDiff$output$ID)

uniqPat <- unique(sampleInfo$NR)

for (n in uniqPat) {

  tumSamp <- sampleInfo$File_name[
    sampleInfo$NR == n &
    sampleInfo$Type =='Tumor']

  KCsegDiff$output$ID[grep(tumSamp, KCsegDiff$output$ID)] <- 
    as.character(n)
  colnames(KCdiff)[grep(tumSamp, colnames(KCdiff))] <-
    as.character(n)
  colnames(KCdiff)[grep(tumSamp, colnames(KCdiff))] <-
    as.character(n)
  colnames(KCsegDiff$data)[grep(tumSamp, colnames(KCsegDiff$data))] <-
    as.character(n)
  
}

save(file='marieke_diffResults_quant_UD2.Rda', 
  list=c('KCdiff','KCsegDiff'))

# Do quantile norm plots

# Reorder for nice visualization of Tumor and LN pairs
sampleInfo <- sampleInfo[order(sampleInfo$NR, sampleInfo$Type),]
KC <- KC[,c('chrom', 'maploc', sampleInfo$File_name)]

dataMatrix <- as.matrix(KC[,3:ncol(KC)])
dataMatrix <- normalize.quantiles(dataMatrix)
KCnorm <- KC
KCnorm[,3:ncol(KCnorm)] <- dataMatrix

png(file="Figures/preNormBoxplot.png", width=1000, height=800)
boxplot(KC[,3:ncol(KC)], las=2, main='boxplots beform qnorm', 
  names=paste(sampleInfo$NR, sampleInfo$Type))
dev.off()

png(file="Figures/postNormBoxplot.png", width=1000, height=800)
boxplot(KCnorm[,3:ncol(KCnorm)], las=2, main='boxplots after qnorm', 
  names=paste(sampleInfo$NR, sampleInfo$Type))
dev.off()

library(multicore)
mclapply(as.list(seq(1, nrow(sampleInfo), 2)), plotTwoSamples, sampleInfoSet=sampleInfo, KCset=KC)

plotTwoSamples <- function(i, sampleInfoSet, KCset) {
  fileName <- paste("Figures/dotplot_", sampleInfoSet$NR[i], ".png", sep='')
  png(file=fileName, width=1200, height=800)    
  par(mfrow=c(2,1))
  plotRawCghDotPlot(KCdataSet=KCset, mirrorLocs=hsMirrorLocs, samples=i,
   doFilter=T, plotTitle=paste(sampleInfoSet$Type[i], sampleInfoSet$NR[i]),
   setylim=c(-1,1), setcex=1)
  plotRawCghDotPlot(KCdataSet=KCset, mirrorLocs=hsMirrorLocs, samples=i+1,
   doFilter=T, plotTitle=paste(sampleInfoSet$Type[i+1], sampleInfoSet$NR[i+1]),
   setylim=c(-1,1), setcex=1)
  dev.off()
}

