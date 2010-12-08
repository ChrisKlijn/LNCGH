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
KCcollTLNcomp <- compareSpmCollection(KCcollTLN, nperms=1000)
sigReg <- getSigRegionsCompKC(KCcollTLNcomp, fdr=.05)

KCcollTLN_TN <- calcSpmCollection(data=KCtumTN, mirrorLocs=hsMirrorLocs, data2=KCLNTN)
KCcollTLN_TNcomp <- compareSpmCollection(KCcollTLN_TN, nperms=1000)
sigReg_TN <- getSigRegionsCompKC(KCcollTLN_TNcomp, fdr=.05)

KCcollTLN_Lum <- calcSpmCollection(data=KCtumLum, mirrorLocs=hsMirrorLocs, data2=KCLNLum)
KCcollTLN_Lumcomp <- compareSpmCollection(KCcollTLN_Lum, nperms=1000)
sigReg_Lum <- getSigRegionsCompKC(KCcollTLN_Lumcomp, fdr=.05)

pdf(file='Figures/compKC_all.pdf', width=15, height=5)
plot(KCcollTLNcomp, sigRegions=sigReg, col1=colors()[122], col2=colors()[148])
abline(h=0)
dev.off()

pdf(file='Figures/compKC_TN.pdf', width=15, height=5)
plot(KCcollTLN_TNcomp, sigRegions=sigReg_TN, col1=colors()[122], col2=colors()[148])
abline(h=0)
dev.off()

pdf(file='Figures/compKC_LUM.pdf', width=15, height=5)
plot(KCcollTLN_Lumcomp, sigRegions=sigReg_Lum, col1=colors()[122], col2=colors()[148])
abline(h=0)
dev.off()

#------------------------------------------------------------------
# KCsmart comparative
#
# All samples, excluding sample Nr. 592 (due to mixed histology)
#------------------------------------------------------------------

KCtum_excl592 <- KCtum[,-which(colnames(KCtum) %in% sampleInfo$File_name[sampleInfo$NR == 592])]
KCLN_excl592 <- KCLN[,-which(colnames(KCLN) %in% sampleInfo$File_name[sampleInfo$NR == 592])]
KCtumLum_excl592 <- KCTumLum[,-which(colnames(KCTumLum) %in% sampleInfo$File_name[sampleInfo$NR == 592])]
KCLNLum_excl592 <- KCLNLum[,-which(colnames(KCLNLum) %in% sampleInfo$File_name[sampleInfo$NR == 592])]


KCcollTLN_excl592 <- calcSpmCollection(data=KCtum_excl592, mirrorLocs=hsMirrorLocs, data2=KCLN_excl592)
KCcollTLNcomp_excl592 <- compareSpmCollection(KCcollTLN_excl592, nperms=1000)
sigReg_all_excl592 <- getSigRegionsCompKC(KCcollTLNcomp_excl592, fdr=.05)


KCcollTLN_Lum_excl592 <- calcSpmCollection(data=KCtumLum_excl592, mirrorLocs=hsMirrorLocs, 
  data2=KCLNLum_excl592)
KCcollTLN_Lum_comp_excl592 <- compareSpmCollection(KCcollTLN_Lum_excl592, nperms=1000)
sigReg_Lum_excl592 <- getSigRegionsCompKC(KCcollTLN_Lum_comp_excl592, fdr=.05)


#------------------------------------------------------------------
# Sandbox
#------------------------------------------------------------------

plot(KCtum[order(KCtum$chrom, KCtum$maploc, decreasing=F),1], pch='.')

par(mfrow=c(2,1))
plotRawCghDotPlot(KCdataSet=KC, mirrorLocs=hsMirrorLocs, 
    samples=which(sampleInfo$NR == 104 & sampleInfo$Type == 'Tumor'), 
    doFilter=T, plotTitle=paste('Tumor - Sample', uniqTumNum[t]),
    chromosomes=11)
plotRawCghDotPlot(KCdataSet=KC, mirrorLocs=hsMirrorLocs, 
    samples=which(sampleInfo$NR == 245 & sampleInfo$Type == 'Tumor'), 
    doFilter=T, plotTitle=paste('Tumor - Sample', uniqTumNum[t]),
    chromosomes=11)

#------------------------------------------------------------------
# Correlation Plots
#------------------------------------------------------------------

library(gplots)

sampNames <- c(colnames(KCtum)[3:ncol(KCtum)], colnames(KCLN)[3:ncol(KCLN)])
labs <- sampleInfo$NR[match(sampNames, sampleInfo$File_name)]
TLNlabs <- as.factor(sampleInfo$Type[match(sampNames, sampleInfo$File_name)]))
sampCorMat <- cor(KCcollTLNcomp@spmCollection@data, use='na.or.complete')

pdf(file='Figures/corrMat_T_LN.pdf', width=7, height=7)
heatmap(sampCorMat, scale='none', labRow=labs, labCol=labs, 
  col=colorpanel(265, low='blue', high='yellow'), 
  ColSideColors=colors()[c(122, 148)][as.numeric(TLNlabs)],
  main='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

#------------------------------------------------------------------
# Per sample difference
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

# KCsmart of the difference

diffSpm <- calcSpm(data=diffKC, mirrorLocs=hsMirrorLocs)
pdf(file='Figures/diffKCsmart.pdf', width=15, height=12)
plot(diffSpm)
dev.off()

#------------------------------------------------------------------
# Plot Per Pair Three plots
#------------------------------------------------------------------

# Originals

uniqTumNum <- unique(sampleInfo$NR)

for (t in 1:length(uniqTumNum)) {
  png(file=paste('dotFigures/oriSample_', uniqTumNum[t], '.png', sep=''), width=1280, height=800)
  par(mfrow=c(2,1))
  plotRawCghDotPlot(KCdataSet=KC, mirrorLocs=hsMirrorLocs, 
    samples=which(sampleInfo$NR == uniqTumNum[t] & sampleInfo$Type == 'Tumor'), 
    doFilter=T, plotTitle=paste('Tumor - Sample', uniqTumNum[t]))
  plotRawCghDotPlot(KCdataSet=KC, mirrorLocs=hsMirrorLocs, 
    samples=which(sampleInfo$NR == uniqTumNum[t] & sampleInfo$Type == 'LN'), 
    doFilter=T, plotTitle=paste('Lymph Node - Sample', uniqTumNum[t]))
  dev.off()
}

# Quantile Normalized

uniqTumNum <- unique(sampleInfo$NR)

for (t in 1:length(uniqTumNum)) {
  png(file=paste('dotFigures/sample_', uniqTumNum[t], '.png', sep=''), width=1280, height=1024)
  par(mfrow=c(3,1))
  plotRawCghDotPlot(KCdataSet=KCnorm, mirrorLocs=hsMirrorLocs, 
    samples=which(sampleInfo$NR == uniqTumNum[t] & sampleInfo$Type == 'Tumor'), 
    doFilter=T, plotTitle=paste('Tumor - Sample', uniqTumNum[t]))
  plotRawCghDotPlot(KCdataSet=KCnorm, mirrorLocs=hsMirrorLocs, 
    samples=which(sampleInfo$NR == uniqTumNum[t] & sampleInfo$Type == 'LN'), 
    doFilter=T, plotTitle=paste('Lymph Node - Sample', uniqTumNum[t]))
  plotRawCghDotPlot(KCdataSet=diffKC, mirrorLocs=hsMirrorLocs, 
    samples=grep(paste(uniqTumNum[t], '$', sep=''), colnames(diffKC)) - 2, 
    doFilter=T, plotTitle=paste('Difference - Sample', uniqTumNum[t]))
  dev.off()
}





