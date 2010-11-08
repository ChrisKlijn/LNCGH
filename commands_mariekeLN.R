##----------------------------------------------------------------------
## Commands Marieke Lymph node to primary tumor comparison
## Data: 720k Nimblegen aCGH
## 
## Christiaan Klijn
##----------------------------------------------------------------------

# WD
setwd('/data/people/klijn/marieke')

# Load data
load('~/marieke/DNACopy/mariekeData.Rda')

# Load clinical data
sampleInfo <- read.delim('Clin_data_tumorLN.txt', stringsAsFactors=F)
sampleInfo <- sampleInfo[order(sampleInfo$File_name),]

# Check if the samples are ordered the same in the dataframe and the sampleinfo
all.equal(colnames(KC)[3:ncol(KC)], sampleInfo$File_name)

#------------------------------------------------------------------
# Exploratory Plots
#------------------------------------------------------------------




#------------------------------------------------------------------
# KCsmart comparative
#
# All samples, Tumor vs. Lymph Node
#------------------------------------------------------------------

library(KCsmart)
data(hsMirrorLocs)

KCcollTLN <- calcSpmCollection(data=KC, mirrorLocs=hsMirrorLocs, cl=sampleInfo$File_name)

