# mariekeNL_pairNormTests.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Tests on how to pairwise normalise the tumor LN pairs
# -------------------------------------------------------------------

# WD - medoid
setwd('/home/klijn/data/smallproj/mariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')

# Load KC data
load('~/data/smallproj/mariekeLN/DNACopy/mariekeData.Rda')

# Load segmentation results
load('marieke_segResult.Rda')

library(KCsmart)
library(preprocessCore)
data(hsMirrorLocs)

# Load clinical data
sampleInfo <- read.delim('Clin_data_tumorLN.txt', stringsAsFactors=F)
sampleInfo <- sampleInfo[order(sampleInfo$File_name),]

# Check if the samples are ordered the same in the dataframe and the sampleinfo
all.equal(colnames(KC)[3:ncol(KC)], sampleInfo$File_name)

# -------------------------------------------------------------------

# Quantile normalization
dataMatrix <- as.matrix(KC[,3:ncol(KC)])
dataMatrix <- normalize.quantiles(dataMatrix)
KCnorm <- KC
KCnorm[,3:ncol(KCnorm)] <- dataMatrix

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
  smallSeg <- subset(KCseg, samplelist=paste('X', colnames(smallKC)[c(3,4)], sep=''))
  smallFreq <- glFrequency(smallSeg, .2)
  
  ind <- smallFreq$gain == 1 | smallFreq$loss == -1
  
  fitrob <- lmrob(smallKC[ind,4] ~ smallKC[ind,3])
  
  diffNorm <- (smallKC[,3] + coefficients(fitrob)[[1]]) * coefficients(fitrob)[[2]] - smallKC[,4]
  

}






