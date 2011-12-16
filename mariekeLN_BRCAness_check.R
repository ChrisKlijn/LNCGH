# mariekeLN_BRCAness_check.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@gmail.com>
# Project: Marieke LN
# Description: Check if BRCAness status is changed in Tumor vs LN
# -------------------------------------------------------------------

# WD - medoid
setwd('/home/klijn/data/smallproj/mariekeLN/philipBRCAclassification/MariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/codeChris/generalFunctionsR/chris_DNAcopy_utils.R')

library(ggplot2)

# load clindata with brcaness info
cdBRCA <- read.delim("ClinPCS151211.txt", header = TRUE,  
  sep = "\t",  stringsAsFactors = FALSE)

table(cdBRCA$NR, cdBRCA$BRCAnessPCS)


