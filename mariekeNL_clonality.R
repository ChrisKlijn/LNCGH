# mariekeNL_correlation.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke Lymph node vs primary
# Description: apply the clonality package to the data of Marieke
# -------------------------------------------------------------------

library(Clonality)

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

# Remove X and Y
KC_nXY <- KC[-which(KC$chrom %in% c(23,24)),]

# Put data in the CNA format
dataCNA <- CNA(genomdat=KC_nXY[,3:ncol(KC_nXY)], chrom=KC_nXY$chrom, 
  maploc=KC_nXY$maploc, data.type='logratio', sampleid=sampleInfo$File_name)

# Average the data so that 10,000 probes remain (per documentation)
# 720k array, 678k after XY removal, so average 68 probes
dataAve <- ave.adj.probes(dataCNA, 68)

# Set maploc to kbp instead of bp
dataAve$maploc <- floor(dataAve$maploc/1000)

dataAve$chrom <- splitChromosomes(dataAve$chrom, dataAve$maploc)

# Run the clonality analysis

clonResult <- clonality.analysis(dataAve, sampleInfo$NR, nmad = 1.25, 
  reference=T, allpairs=F)

  
