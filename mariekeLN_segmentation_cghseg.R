# mariekeLN_segmentation_cghseg.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Segmentation using the cghseg package
# -------------------------------------------------------------------


# WD
setwd('~/data/smallproj/mariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghsegmentation_functions.R')

# Libraries
library(cghseg)
library(multicore)
library(GenomicRanges)

# Load data
load('~/data/smallproj/mariekeLN/DNACopy/mariekeData.Rda')

# Load clinical data
sampleInfo <- read.delim('Clin_data_tumorLN.txt', stringsAsFactors=F)
sampleInfo <- sampleInfo[order(sampleInfo$File_name),]

# Check if the samples are ordered the same in the dataframe and the sampleinfo
all.equal(colnames(KC)[3:ncol(KC)], sampleInfo$File_name)

chromList <- unique(KC$chrom)
names(chromList) <- paste('chr', unique(KC$chrom), sep='')

segmentList <- mclapply(chromList, doCghSeg, allKC=KC, chrom=KC$chrom)

save(file='marieke_segResult_cghseg.Rda', list=c('segmentList'))
