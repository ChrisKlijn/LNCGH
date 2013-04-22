# mariekeLN_segmentation.R
# -------------------------------------------------------------------
# Copyright 2010 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Segmentation of the tumors to find differentially
#              changed region in the aCGH data.
# -------------------------------------------------------------------

# WD
projectDir <- '~/Projects/NKIProjects/LNCGH/'
setwd(projectDir)
source('~/OldNKI/code Backup/codeChris/generalFunctionsR/chris_cghdata_analysis.R')

# Libraries
library(DNAcopy)

# Load data
load(file.path(projectDir, 'newDataImport.RData'))

# Check if the samples are ordered the same in the dataframe and the sampleinfo
all.equal(colnames(rawData)[3:ncol(rawData)], sampleInfo$File_name)

CNA.rawData <- CNA(as.matrix(rawData[,3:ncol(rawData)]), rawData$chrom, rawData$maploc, data.type=c("logratio"), 
  sampleid=colnames(rawData)[3:ncol(rawData)])
CNA.rawData.smoothed <- smooth.CNA(CNA.rawData)
rawDataseg <- segment(CNA.rawData.smoothed, verbose=1, undo.splits='sdundo')

names(rawDataseg$data) <- gsub('X', '', names(rawDataseg$data))
rawDataseg$output$ID <- gsub('X', '', rawDataseg$output$ID)

save(file='marieke_segResult.Rda', list=c('rawDataseg')) 
