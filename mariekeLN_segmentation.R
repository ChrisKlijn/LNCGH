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

# Overlapping segments per tumor-lymphnode pairs

# Run local, to use correct version of R and Bioconductor

setwd('/home/cklijn/work/Matlab/Data/NKI Data/Marieke analysis/LN')
load('marieke_segResult.Rda')
library(GenomicRanges)

# Do analysis only on segments > .2 or < -.2
KCseg_filter_2 <- KCseg_filter[abs(KCseg_filter$seg.mean) > .2,]
uniqNR <- unique(sampleInfo$NR)
resultCount <- vector(mode='numeric', length=length(uniqNR))

for (i in 1:length(uniqNR)) {
  tempSegTum <- KCseg_filter_2[KCseg_filter_2$ID %in% sampleInfo$File_name[(sampleInfo$NR == uniqNR[i]) &
    (sampleInfo$Type == 'Tumor')],]
  tempSegLN <- KCseg_filter_2[KCseg_filter_2$ID %in% sampleInfo$File_name[(sampleInfo$NR == uniqNR[i]) &
    (sampleInfo$Type == 'LN')],]
  resultCount[i] <- nrow(tempSegLN) - nrow(tempSegTum)
}

countResultFrame <- as.data.frame(resultCount)
countResultFrame$NR <- sampleInfo$NR[seq(1,nrow(sampleInfo),2)] 
countResultFrame$BRCAness <- sampleInfo$BRCAness[seq(1,nrow(sampleInfo),2)]
countResultFrame$Mol_subtype <- sampleInfo$Mol_subtype[seq(1,nrow(sampleInfo),2)]

## SANDBOX AREA

tempIRTum <- IRanges(start=tempSegTum$loc.start, end=tempSegTum$loc.end)
tempIRLN <- IRanges(start=tempSegLN$loc.start, end=tempSegLN$loc.end)

tempGRangesTum <- GRanges(seqnames=tempSegTum$chrom, ranges=tempIRTum, strand=rep('+', nrow(tempSegTum)),
  values=tempSegTum$ID, score=tempSegTum$seg.mean)

tempGRangesLN <- GRanges(seqnames=tempSegLN$chrom, ranges=tempIRLN, strand=rep('+', nrow(tempSegLN)),
  values=tempSegLN$ID, score=tempSegLN$seg.mean)

## SANDBOX AREA

# Calculate differences between tumor and ln

for (i in 1:length(uniqNR)) {
  tempSegTum <- KCseg_filter_2[KCseg_filter_2$ID %in% sampleInfo$File_name[(sampleInfo$NR == uniqNR[i]) &
    (sampleInfo$Type == 'Tumor')],]
  tempSegLN <- KCseg_filter_2[KCseg_filter_2$ID %in% sampleInfo$File_name[(sampleInfo$NR == uniqNR[i]) &
    (sampleInfo$Type == 'LN')],]
}

# -------------------------------------------------------------------
# Plotting

library(ggplot2)
postscript(file='Figures/boxplot_seg2_brcaness.eps', paper='special', horizontal=F, width=5, height=4)
qplot(data=sampleInfo, x=BRCAness, y=segcount2, geom='boxplot', fill=BRCAness) + theme_bw() +
  opts(title='Segment counts - BRCAness')
dev.off()

postscript(file='Figures/boxplot_seg2_molsubtype.eps', paper='special', horizontal=F, width=5, height=4)
qplot(data=sampleInfo, x=Mol_subtype, y=segcount2, geom='boxplot', fill=Mol_subtype) + theme_bw() + 
  opts(title='Segment counts - Molecular Subtype')
dev.off()

# Segment differences

postscript(file='Figures/scatter_segmentdiff_2.eps', paper='special', horizontal=F, width=5, height=4)
qplot(data=countResultFrame, x=BRCAness, y=resultCount, position=position_jitter(w=0.3, h=0), 
  size=I(3), col=BRCAness) + opts(title='Verschil in segmenten > < .2 - T en LN') + theme_bw()
dev.off()

postscript(file='Figures/scatter_segmentdiff_2_molsub.eps', paper='special', horizontal=F, width=5, height=4)
qplot(data=countResultFrame, x=Mol_subtype, y=resultCount, position=position_jitter(w=0.3, h=0), 
  size=I(3), col=Mol_subtype) + opts(title='Verschil in segmenten > < .2 - T en LN') + theme_bw()
dev.off()

