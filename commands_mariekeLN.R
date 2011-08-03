##----------------------------------------------------------------------
## Commands Marieke Lymph node to primary tumor comparison
## Data: 720k Nimblegen aCGH
## 
## Christiaan Klijn
##----------------------------------------------------------------------

# WD - medoid
setwd('/home/klijn/data/smallproj/mariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')

# Load data
load('~data/smallproj/mariekeLN/DNACopy/mariekeData.Rda')

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
# Requested tumor plots
#------------------------------------------------------------------

plotSamp <- which(sampleInfo$NR == 44 & sampleInfo$Type == 'Tumor')
postscript(file='dotFigures/t44.eps', width=12, height=5, paper='special', horizontal=F)
plotRawCghDotPlot(KCdataSet=KC, mirrorLocs=hsMirrorLocs, samples=plotSamp, doFilter=T, plotTitle='Tumor 44')
dev.off()

plotSamp <- which(sampleInfo$NR == 565 & sampleInfo$Type == 'Tumor')
postscript(file='dotFigures/t565.eps', width=12, height=5, paper='special', horizontal=F)
plotRawCghDotPlot(KCdataSet=KC, mirrorLocs=hsMirrorLocs, samples=plotSamp, doFilter=T, plotTitle='Tumor 565')
dev.off()

#------------------------------------------------------------------
# KCsmart comparative
#
# All samples, Tumor vs. Lymph Node
#------------------------------------------------------------------

library(KCsmart)
data(hsMirrorLocs)

# Remove sample 592

KC <- KC[,-(which(sampleInfo$NR == 592) + 2)]
sampleInfo <- sampleInfo[-which(sampleInfo$NR == 592),]

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

load('allSampl_TNcomp.Rda')
#pdf(file='Figures/compKC_TN.pdf', width=15, height=5)
postscript(file='Figures/compKC_TN.eps', width=15, height=5, paper='special', horizontal=F)
plot(KCcollTLN_TNcomp, sigRegions=sigReg_TN, col1=colors()[122], col2=colors()[148])
abline(h=0)
dev.off()

load('allSampl_Lumcomp.Rda')
#pdf(file='Figures/compKC_LUM.pdf', width=15, height=5)
postscript(file='Figures/compKC_Lum.eps', width=15, height=5, paper='special', horizontal=F)
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

# pdf(file='Figures/corrMat_T_LN.pdf', width=7, height=7)
postscript(file='Figures/corrMat_T_LN.eps', paper='special', horizontal=F, width=7, height=7)
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

# ----------------------------------------------------------------------------___
# SANDBOX and old code
# ----------------------------------------------------------------------------___


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

shortSampleInfo <- sampleInfo[,c('File_name', 'Type', 'NR')]
combinedSeg <- merge(x=KCseg_filter_2, y=shortSampleInfo, by.x='ID', 
  by.y='File_name')


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



# -------------------------------------------------------------------

# Quantile normalization
dataMatrix <- as.matrix(KC[,3:ncol(KC)])
dataMatrix <- normalize.quantiles(dataMatrix)
KCnorm <- KC
KCnorm[,3:ncol(KCnorm)] <- dataMatrix

# -------------------------------------------------------------------
# Muck code, linear regression
# -------------------------------------------------------------------

smallKC <- KC[,c(1,2,3,4)]
smallSeg <- subset(KCseg, samplelist=paste('X', colnames(KC)[c(3,4)], sep=''))
smallFreq <- glFrequency(smallSeg, .2)

ind <- smallFreq$gain == 0 & smallFreq$loss == 0
g0 <- smallFreq$gain == 0

par(mfrow=c(2,2))
plot(smallKC[!ind, 2], smallKC[!ind,3] + 2*smallKC[!ind, 1], pch='.')
abline(h=2*unique(smallKC$chrom))
plot(smallKC[ind, 2], smallKC[ind,3] + 2*smallKC[ind, 1], pch='.')
abline(h=2*unique(smallKC$chrom))
plot(smallKC[!ind, 2], smallKC[!ind,4] + 2*smallKC[!ind, 1], pch='.')
abline(h=2*unique(smallKC$chrom))
plot(smallKC[ind, 2], smallKC[ind,4] + 2*smallKC[ind, 1], pch='.')
abline(h=2*unique(smallKC$chrom))

par(mfrow=c(1,2))
plot(smallKC[smallFreq$gain == 1, 2], smallKC[smallFreq$gain == 1,3] + 3*smallKC[smallFreq$gain == 1, 1], 
  pch='.', col='red')
points(smallKC[smallFreq$loss == -1, 2], smallKC[smallFreq$loss == -1,3] + 3*smallKC[smallFreq$loss == -1, 1], 
  pch='.', col='green')

abline(h=2*unique(smallKC$chrom))
plot(smallKC[ind, 2], smallKC[ind,3] + 2*smallKC[ind, 1], pch='.')
abline(h=2*unique(smallKC$chrom))

ind2 <- smallFreq$gain == 1 | smallFreq$loss == -1
plot(smallKC[ind2,3], smallKC[ind2,4], pch='.')

par(mfrow=c(1,2))
plot(smallKC[!ind2, 2], smallKC[!ind2,3] + 2*smallKC[!ind2, 1], pch='.')
abline(h=2*unique(smallKC$chrom))
plot(smallKC[ind2, 2], smallKC[ind2,3] + 2*smallKC[ind2, 1], pch='.')
abline(h=2*unique(smallKC$chrom))

par(mfrow=c(1,2))

fit <- lm(smallKC[ind2,4] ~ smallKC[ind2,3])
plot(smallKC[ind2,3], smallKC[ind2,4], pch='.', col='gray')
abline(a=coefficients(fit)[[1]], b=coefficients(fit)[[2]], col='gray')
points((smallKC[ind2,3]+ coefficients(fit)[[1]]) * coefficients(fit)[[2]], smallKC[ind2,4] , pch='.', 
  col='darkgray')
abline(a=0, b=1, col='black')

fitrob <- lmrob(smallKC[ind2,4] ~ smallKC[ind2,3])
plot(smallKC[ind2,3], smallKC[ind2,4], pch='.', col='gray')
abline(a=coefficients(fitrob)[[1]], b=coefficients(fitrob)[[2]], col='gray')
points((smallKC[ind2,3] + coefficients(fitrob)[[1]]) * coefficients(fitrob)[[2]], smallKC[ind2,4], pch='.', 
  col='darkgray')
abline(a=0, b=1, col='black')

diffKC <- smallKC[,c(1,2)]
diffKCnorm <- smallKC[,c(1,2)]
diffKCquant <- smallKC[,c(1,2)]
diffKCnorm$diff <- (smallKC[,3] + coefficients(fitrob)[[1]]) * coefficients(fitrob)[[2]] - smallKC[,4]
diffKC$diff <- smallKC[,3] - smallKC[,4]
diffKCquant$diff <- KCnorm[,3] - KCnorm[,4]

par(mfrow=c(2,1))

plotRawCghDotPlot(KCdataSet=smallKC, mirrorLocs=hsMirrorLocs, 
  samples=1, doFilter=T, plotTitle='LN')
plotRawCghDotPlot(KCdataSet=smallKC, mirrorLocs=hsMirrorLocs, 
  samples=2, doFilter=T, plotTitle='Tumor')

par(mfrow=c(3,1))  
  
plotRawCghDotPlot(KCdataSet=diffKC, mirrorLocs=hsMirrorLocs, 
  samples=1, doFilter=T, plotTitle='Diff non Norm')
plotRawCghDotPlot(KCdataSet=diffKCnorm, mirrorLocs=hsMirrorLocs, 
    samples=1, doFilter=T, plotTitle='Diff Norm')
plotRawCghDotPlot(KCdataSet=diffKCquant, mirrorLocs=hsMirrorLocs, 
    samples=1, doFilter=T, plotTitle='Diff quant Norm')

combMat <- matrix(1:nrow(sampleInfo), ncol=2)
combList <- vector(mode='list', length=nrow(sampleInfo)/2)
names(combList) <- sampleInfo$NR[duplicated(sampleInfo$NR)]

for (cl in 1:length(combList)) {
  tempName <- as.numeric(names(combList))[cl]
  combList[[cl]] <- with(sampleInfo, 
    c(which(Type=='LN' & NR==tempName), which(Type=='Tumor' & NR==tempName)))
}

normRobustLinear <- function(KC, KCseg combMat) {
  
  require(multicore)
  
  diffList <- mclapply(combList, doLinRob, KC, KCseg)

}

doLinRob <- function(comb, KC, KCseg) {
  
  require(robustbase)
  
  smallKC <- KC[,c(1,2, comb)]
  smallSeg <- subset(KCseg, samplelist=colnames(smallKC)[c(3,4)])
  smallFreq <- glFrequency(smallSeg, .2)
  
  ind <- smallFreq$gain == 1 | smallFreq$loss == -1
  
  fitrob <- lmrob(smallKC[ind,4] ~ smallKC[ind,3])
  
  diffNorm <- (smallKC[,3] + coefficients(fitrob)[[1]]) * coefficients(fitrob)[[2]] - smallKC[,4]
  

}




