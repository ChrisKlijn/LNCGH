# mariekeNL_correlation.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke Lymph node vs primary
# Description: Data to do correlation plots between tumors
# -------------------------------------------------------------------

# WD - medoid
setwd('/home/klijn/data/smallproj/mariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')

# Load data
load('~/data/smallproj/mariekeLN/DNACopy/mariekeData.Rda')

# Load clinical data
sampleInfo <- read.delim('Clin_data_tumorLN.txt', stringsAsFactors=F)
sampleInfo <- sampleInfo[order(sampleInfo$File_name),]

# Check if the samples are ordered the same in the dataframe and the sampleinfo
all.equal(colnames(KC)[3:ncol(KC)], sampleInfo$File_name)

#------------------------------------------------------------------
# Correlation Plots
#------------------------------------------------------------------

library(gplots)

# Load comparative analysis
load('marieke_compKC.Rda')

sampNames <- colnames(KC[3:ncol(KC)])
labs <- sampleInfo$NR[match(sampNames, sampleInfo$File_name)]
TLNlabs <- as.factor(sampleInfo$Type[match(sampNames, sampleInfo$File_name)])
sampCorMat <- cor(KCcollTLN@data, use='na.or.complete')

# pdf(file='Figures/corrMat_T_LN.pdf', width=7, height=7)
postscript(file='Figures/corrMat_T_LN.eps', paper='special', horizontal=F, width=7, height=7)
heatmap(sampCorMat, scale='none', labRow=labs, labCol=labs, 
  col=colorpanel(265, low='blue', high='yellow'), 
  ColSideColors=colors()[c(122, 148)][as.numeric(TLNlabs)],
  main='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

#------------------------------------------------------------------
# Dotplots of non-concordant samples
#------------------------------------------------------------------

library(KCsmart)
data(hsMirrorLocs)

# Patient samples 592 and 456 show disconcordance - plot them in a figure

KC592 <- KC[,c(1,2, which(sampleInfo$NR == 592)+2)]

png(file='Figures/Patient592.png', width=1600, height=800*(ncol(KC592)-2))
  par(mfrow=c(ncol(KC592)-2, 1))
  for (p in 1:(ncol(KC592)-2)) {
    sampleNr <- which(sampleInfo$File_name == colnames(KC592)[p+2])
    plotRawCghDotPlot(KCdataSet=KC592, samples=p,
      mirrorLocs=hsMirrorLocs, doFilter=T, 
      plotTitle=paste(sampleInfo[sampleNr, c('File_name','NR','Type') ], collapse=' '))
  }
dev.off()


KC456 <- KC[,c(1,2, which(sampleInfo$NR == 456)+2)]
png(file='Figures/Patient456.png', width=1600, height=800*(ncol(KC456)-2))
  par(mfrow=c(ncol(KC456)-2, 1))
  for (p in 1:(ncol(KC456)-2)) {
    sampleNr <- which(sampleInfo$File_name == colnames(KC456)[p+2])
    plotRawCghDotPlot(KCdataSet=KC456, samples=p,
      mirrorLocs=hsMirrorLocs, doFilter=T, 
      plotTitle=paste(sampleInfo[sampleNr, c('File_name','NR','Type') ], collapse=' '))
  }
dev.off()



