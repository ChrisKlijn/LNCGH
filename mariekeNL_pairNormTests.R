# mariekeNL_pairNormTests.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Tests on how to pairwise normalise the tumor LN pairs
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

# -------------------------------------------------------------------

a <- KC[,c(1,2,3,4)]
zeroInd <- (abs(a[,3]) < .2) & (abs(a[,4]) < .2)
plot(a[!zeroInd,3], a[!zeroInd,4], pch='.')
