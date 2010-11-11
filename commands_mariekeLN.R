##----------------------------------------------------------------------
## Commands Marieke Lymph node to primary tumor comparison
## Data: 720k Nimblegen aCGH
## 
## Christiaan Klijn
##----------------------------------------------------------------------

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

#------------------------------------------------------------------
# Exploratory Plots
#------------------------------------------------------------------

library(KCsmart)
data(hsMirrorLocs)

pdf(file='dotFigures/T_LN_comp_t44.pdf', width=12, height=10)
par(mfrow=c(2,1))
plotRawCghDotPlot(KCdataSet=KC, mirrorLocs=hsMirrorLocs, samples=1, doFilter=T, plotTitle='Tumor 44')
plotRawCghDotPlot(KCdataSet=KC, mirrorLocs=hsMirrorLocs, samples=2, doFilter=T, plotTitle='Lymph Node 44')
dev.off()

png(file='dotFigures/T_LN_comp_t44.png', width=1200, height=800)
par(mfrow=c(2,1))
plotRawCghDotPlot(KCdataSet=KC, mirrorLocs=hsMirrorLocs, samples=1, doFilter=T, plotTitle='Tumor 44')
plotRawCghDotPlot(KCdataSet=KC, mirrorLocs=hsMirrorLocs, samples=2, doFilter=T, plotTitle='Lymph Node 44')
dev.off()

#------------------------------------------------------------------
# KCsmart comparative
#
# All samples, Tumor vs. Lymph Node
#------------------------------------------------------------------

library(KCsmart)
data(hsMirrorLocs)

# Separate the data

KCtum <- KC[,c(1,2,which(sampleInfo$Type == 'Tumor')+2)]
KCtumTN <-  KC[,c(1,2,which(sampleInfo$Type == 'Tumor' & sampleInfo$Mol_subtype == 'TN')+2)]
KCtumLum <-  KC[,c(1,2,which(sampleInfo$Type == 'Tumor' & sampleInfo$Mol_subtype == 'luminal')+2)]
KCLN <- KC[,c(1,2,which(sampleInfo$Type == 'LN')+2)]
KCLNTN <- KC[,c(1,2,which(sampleInfo$Type == 'LN' & 
  (sampleInfo$NR %in% sampleInfo$NR[sampleInfo$Mol_subtype == 'TN']))+2)]
KCLNLum <- KC[,c(1,2,which(sampleInfo$Type == 'LN' & 
  (sampleInfo$NR %in% sampleInfo$NR[sampleInfo$Mol_subtype == 'luminal']))+2)]


KCcollTLN <- calcSpmCollection(data=KCtum, mirrorLocs=hsMirrorLocs, data2=KCLN)
KCcollTLNcomp <- compareSpmCollection(KCcollTLN, nperms=1000, fdr=0.1)
sigReg <- getSigRegionsCompKC(KCcollTLNcomp, fdr=.01)

KCcollTLN_TN <- calcSpmCollection(data=KCtumTN, mirrorLocs=hsMirrorLocs, data2=KCLNTN)
KCcollTLN_TNcomp <- compareSpmCollection(KCcollTLN_TN, nperms=1000, fdr=0.1)
sigReg_TN <- getSigRegionsCompKC(KCcollTLN_TNcomp, fdr=.01)

KCcollTLN_Lum <- calcSpmCollection(data=KCtumLum, mirrorLocs=hsMirrorLocs, data2=KCLNLum)
KCcollTLN_Lumcomp <- compareSpmCollection(KCcollTLN_Lum, nperms=1000, fdr=0.1)
sigReg_TN <- getSigRegionsCompKC(KCcollTLN_Lumcomp, fdr=.01)

#------------------------------------------------------------------
# Per Sample delta profile
#------------------------------------------------------------------

library(preprocessCore)

# Quantile normalization
dataMatrix <- as.matrix(KC[,3:ncol(KC)])
dataMatrix <- normalize.quantiles(dataMatrix)
KCnorm <- KC
KCnorm[,3:ncol(KCnorm)] <- dataMatrix

# Instantiate KC frame for the delta profile
diffKC <- KCnorm[,c(1,2)]

# Loop over the tumor numbers - making a delta profile per Tumor-LN pair
uniqTumNum <- unique(sampleInfo$NR)
for (t in 1:length(uniqTumNum)) {
  diffKC$temp <- KCnorm[,sampleInfo$File_name[
    (sampleInfo$NR == uniqTumNum[t] & sampleInfo$Type == 'LN')]] - KCnorm[,sampleInfo$File_name[
    (sampleInfo$NR == uniqTumNum[t] & sampleInfo$Type == 'Tumor')]]
  colnames(diffKC) <- gsub('temp', paste('Sample', uniqTumNum[t], sep='-'), colnames(diffKC))
}

plotRawCghDotPlot(KCdataSet=diffKC, mirrorLocs=hsMirrorLocs, samples=1, doFilter=T, plotTitle='Tumor 44')


