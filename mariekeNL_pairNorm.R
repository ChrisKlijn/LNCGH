# mariekeNL_pairNorm.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Do robust regression to create delta profiles and
#              segment those.
# -------------------------------------------------------------------

# WD - medoid
setwd('/home/klijn/data/smallproj/mariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/codeChris/generalFunctionsR/chris_delta_functions.R')

# Load KC data
load('~/data/smallproj/mariekeLN/DNACopy/mariekeData.Rda')

# Load segmentation results
load('marieke_segResult.Rda')

library(KCsmart)
library(preprocessCore)
library(multicore)
data(hsMirrorLocs)

# Body Code

# Load clinical data
sampleInfo <- read.delim('Clin_data_tumorLN.txt', stringsAsFactors=F)
sampleInfo <- sampleInfo[order(sampleInfo$File_name),]

# Remove samples 592 and 456

KC <- KC[,c('chrom', 'maploc', sampleInfo$File_name[which(!sampleInfo$NR %in% c(592, 456))])]
KCseg <- subset(KCseg, samplelist=sampleInfo$File_name[which(!sampleInfo$NR %in% c(592, 456))])
sampleInfo <- sampleInfo[-which(sampleInfo$NR %in% c(592, 456)),]

# Check if the samples are ordered the same in the dataframe, the sampleinfo and the segdata
all.equal(colnames(KC)[3:ncol(KC)], sampleInfo$File_name)
all.equal(names(KCseg$data)[3:ncol(KC)], sampleInfo$File_name)
all.equal( unique(KCseg$output$ID), sampleInfo$File_name)

# Make the combination numbers, always in the order LN - Tumor

combList <- vector(mode='list', length=length(unique(sampleInfo$NR)))
names(combList) <- as.character(unique(sampleInfo$NR))


for (nr in unique(sampleInfo$NR)) {
  tempComb <- c(0,0)
  tempComb[1] <- sampleInfo$File_name[sampleInfo$NR == nr & sampleInfo$Type == 'LN']
  tempComb[2] <- sampleInfo$File_name[sampleInfo$NR == nr & sampleInfo$Type == 'Tumor']
  
  combList[[as.character(nr)]] <- tempComb
}

diffList <- mclapply(combList, deltaLinear, KC, KCseg)
KCdiff <- cbind(KC[,c('chrom', 'maploc')], as.data.frame(diffList, optional=T))


CNA.KCdiff<- CNA(as.matrix(KCdiff[,3:ncol(KCdiff)]), KCdiff$chrom, KCdiff$maploc, data.type=c("logratio"), 
  sampleid=colnames(KCdiff)[3:ncol(KCdiff)])
CNA.KCdiff.smoothed <- smooth.CNA(CNA.KCdiff)
KCsegDiff <- segment(CNA.KCdiff.smoothed, verbose=1, undo.splits='sdundo')

names(KCsegDiff$data) <- gsub('X', '', names(KCsegDiff$data))
KCsegDiff$output$ID <- gsub('X', '', KCsegDiff$output$ID)

save(file='marieke_diffResults.Rda', list=c('KCdiff','KCsegDiff')






