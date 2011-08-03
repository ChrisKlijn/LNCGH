# mariekeLN_segmentation.R
# -------------------------------------------------------------------
# Copyright 2010 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Segmentation of the tumors to find differentially
#              changed region in the aCGH data.
# -------------------------------------------------------------------

# WD
setwd('~/data/smallproj/mariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')

# Libraries
library(DNAcopy)

# Load data
load('~/data/smallproj/mariekeLN/DNACopy/mariekeData.Rda')

# Load clinical data
sampleInfo <- read.delim('Clin_data_tumorLN.txt', stringsAsFactors=F)
sampleInfo <- sampleInfo[order(sampleInfo$File_name),]

# Check if the samples are ordered the same in the dataframe and the sampleinfo
all.equal(colnames(KC)[3:ncol(KC)], sampleInfo$File_name)

CNA.KC <- CNA(as.matrix(KC[,3:ncol(KC)]), KC$chrom, KC$maploc, data.type=c("logratio"), 
  sampleid=colnames(KC)[3:ncol(KC)])
CNA.KC.smoothed <- smooth.CNA(CNA.KC)
KCseg <- segment(CNA.KC.smoothed, verbose=1, undo.splits='sdundo')

names(KCseg$data) <- gsub('X', '', names(KCseg$data))
KCseg$output$ID <- gsub('X', '', KCseg$output$ID)

save(file='marieke_segResult.Rda', list=c('KCseg')) 
