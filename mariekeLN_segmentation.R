# mariekeLN_segmentation.R
# -------------------------------------------------------------------
# Copyright 2010 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Segmentation of the tumors to find differentially
#              changed region in the aCGH data.
# -------------------------------------------------------------------

# WD
setwd('/data/people/klijn/marieke')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')

# Load data
load('~/marieke/DNACopy/mariekeData.Rda')

# Load clinical data
sampleInfo <- read.delim('Clin_data_tumorLN.txt', stringsAsFactors=F)
sampleInfo <- sampleInfo[order(sampleInfo$File_name),]

# Check if the samples are ordered the same in the dataframe and the sampleinfo
all.equal(colnames(KC)[3:ncol(KC)], sampleInfo$File_name)

# Remove sample 592

KC <- KC[,-(which(sampleInfo$NR == 592) + 2)]
sampleInfo <- sampleInfo[-which(sampleInfo$NR == 592),]

# Check if the samples are ordered the same in the dataframe and the sampleinfo
all.equal(colnames(KC)[3:ncol(KC)], sampleInfo$File_name)

library(DNAcopy)

CNA.KC <- CNA(as.matrix(KC[,3:ncol(KC)]), KC$chrom, KC$maploc, data.type=c("logratio"), 
  sampleid=colnames(KC)[3:ncol(KC)])
CNA.KC.smoothed <- smooth.CNA(CNA.KC)
KCseg <- segment(CNA.KC.smoothed, verbose=1, undo.splits='sdundo')

KCseg_filter <- KCseg$output[KCseg$output$num.mark > 9,]
KCseg_filter$ID <- gsub('X', '', KCseg_filter$ID)

segResultFrame <- data.frame(ID=unique(KCseg_filter$ID), 
  segcount5=0,segcount2=0, segcountAll=0,
  stringsAsFactors=F)


for (r in 1:nrow(segResultFrame)) {
    segResultFrame$segcount5[r] <- 
      sum(abs(KCseg_filter$seg.mean[KCseg_filter$ID == segResultFrame$ID[r]]) > .5)
    segResultFrame$segcount2[r] <- 
      sum(abs(KCseg_filter$seg.mean[KCseg_filter$ID == segResultFrame$ID[r]]) > .2)
    segResultFrame$segcountAll[r] <- 
      sum(KCseg_filter$ID == segResultFrame$ID[r])
}

sampleInfo <- merge(x=sampleInfo, y=segResultFrame, by.x='File_name', by.y='ID')

# Clean up Sampleinfo
sampleInfo <- sampleInfo[order(sampleInfo$NR, sampleInfo$Type),]
sampleInfo$Mol_subtype[seq(1, nrow(a), 2)] <- sampleInfo$Mol_subtype[seq(2, nrow(a), 2)]
sampleInfo$Treatment[seq(1, nrow(a), 2)] <- sampleInfo$Treatment[seq(2, nrow(a), 2)]
sampleInfo$LN_10[seq(1, nrow(a), 2)] <- sampleInfo$LN_10[seq(2, nrow(a), 2)]
sampleInfo$BRgrade[seq(1, nrow(a), 2)] <- sampleInfo$BRgrade[seq(2, nrow(a), 2)]
sampleInfo$recurrencedeath[seq(1, nrow(a), 2)] <- sampleInfo$recurrencedeath[seq(2, nrow(a), 2)]
sampleInfo$osstat[seq(1, nrow(a), 2)] <- sampleInfo$osstat[seq(2, nrow(a), 2)]
sampleInfo$B1_classifier[seq(1, nrow(a), 2)] <- sampleInfo$B1_classifier[seq(2, nrow(a), 2)]
sampleInfo$B2_classifier[seq(1, nrow(a), 2)] <- sampleInfo$B2_classifier[seq(2, nrow(a), 2)]
sampleInfo$BRCAness[seq(1, nrow(a), 2)] <- sampleInfo$BRCAness[seq(2, nrow(a), 2)]
sampleInfo$ER_new[seq(1, nrow(a), 2)] <- sampleInfo$ER_new[seq(2, nrow(a), 2)]
sampleInfo$PR_new[seq(1, nrow(a), 2)] <- sampleInfo$PR_new[seq(2, nrow(a), 2)]
sampleInfo$Basal_Nielsen_IHC[seq(1, nrow(a), 2)] <- sampleInfo$Basal_Nielsen_IHC[seq(2, nrow(a), 2)]

# Plotting

library(ggplot2)
postscript(file='Figures/boxplot_seg2_brcaness.eps', paper='special', horizontal=F, width=5, height=4)
qplot(data=sampleInfo, x=BRCAness, y=segcount2, geom='boxplot', fill=BRCAness) + theme_bw()
dev.off()

