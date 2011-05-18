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


