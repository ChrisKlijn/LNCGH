# mariekeNL_compKC.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke Lymph node vs primary
# Description: Comparative KC smart to do correlation afterwards.
# -------------------------------------------------------------------

# WD - medoid
setwd('/home/klijn/data/smallproj/mariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')

# Load data
load('~/data/smallproj/mariekeLN/DNACopy/mariekeData.Rda')

# Load clinical data
sampleInfo <- read.delim('Clin_data_tumorLN.txt', stringsAsFactors=F)
sampleInfo <- sampleInfo[order(sampleInfo$File_name),]

#------------------------------------------------------------------
# KCsmart comparative
#
# Exclude sample 592 and 456
# Run compKC on all samples, only the Triple Negatives and only
# the luminals
#------------------------------------------------------------------

library(KCsmart)
data(hsMirrorLocs)

# Exclude sample 592 and 456 

KC <- KC[,c('chrom', 'maploc', with(sampleInfo, File_name[!NR %in% c(592, 456)]))]
sampleInfo <- subset(sampleInfo, !(NR %in% c(592, 456)))

# Check if the samples are ordered the same in the dataframe and the sampleinfo
all.equal(colnames(KC)[3:ncol(KC)], sampleInfo$File_name)

# Separate the data

KCtum <- KC[,c(1,2,which(sampleInfo$Type == 'Tumor')+2)]
KCtumTN <-  KC[,c(1,2,which(sampleInfo$Type == 'Tumor' & sampleInfo$Mol_subtype == 'TN')+2)]
KCtumLum <-  KC[,c(1,2,which(sampleInfo$Type == 'Tumor' & sampleInfo$Mol_subtype == 'luminal')+2)]
KCLN <- KC[,c(1,2,which(sampleInfo$Type == 'LN')+2)]
KCLNTN <- KC[,c(1,2,which(sampleInfo$Type == 'LN' & 
  (sampleInfo$NR %in% sampleInfo$NR[sampleInfo$Mol_subtype == 'TN']))+2)]
KCLNLum <- KC[,c(1,2,which(sampleInfo$Type == 'LN' & 
  (sampleInfo$NR %in% sampleInfo$NR[sampleInfo$Mol_subtype == 'luminal']))+2)]

# Run comparative analysis

KCcollTLN <- calcSpmCollection(data=KCtum, mirrorLocs=hsMirrorLocs, data2=KCLN)
KCcollTLNcomp <- compareSpmCollection(KCcollTLN, nperms=1000)
sigReg <- getSigRegionsCompKC(KCcollTLNcomp, fdr=.05)

KCcollTLN_TN <- calcSpmCollection(data=KCtumTN, mirrorLocs=hsMirrorLocs, data2=KCLNTN)
KCcollTLN_TNcomp <- compareSpmCollection(KCcollTLN_TN, nperms=1000)
sigReg_TN <- getSigRegionsCompKC(KCcollTLN_TNcomp, fdr=.05)

KCcollTLN_Lum <- calcSpmCollection(data=KCtumLum, mirrorLocs=hsMirrorLocs, data2=KCLNLum)
KCcollTLN_Lumcomp <- compareSpmCollection(KCcollTLN_Lum, nperms=1000)
sigReg_Lum <- getSigRegionsCompKC(KCcollTLN_Lumcomp, fdr=.05)

# Plot the results

postscript(file='Figures/compKC_all.eps', width=15, height=5, paper='special', horizontal=F)
plot(KCcollTLNcomp, sigRegions=sigReg, col1=colors()[122], col2=colors()[148])
abline(h=0)
dev.off()

postscript(file='Figures/compKC_TN.eps', width=15, height=5, paper='special', horizontal=F)
plot(KCcollTLN_TNcomp, sigRegions=sigReg_TN, col1=colors()[122], col2=colors()[148])
abline(h=0)
dev.off()

postscript(file='Figures/compKC_Lum.eps', width=15, height=5, paper='special', horizontal=F)
plot(KCcollTLN_Lumcomp, sigRegions=sigReg_Lum, col1=colors()[122], col2=colors()[148])
abline(h=0)
dev.off()
