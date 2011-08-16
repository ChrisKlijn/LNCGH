# mariekeNL_plotDiff.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Plot rainbow plots with the difference profiles
# -------------------------------------------------------------------

# WD - medoid
setwd('/home/klijn/data/smallproj/mariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/codeChris/generalFunctionsR/chris_DNAcopy_utils.R')

library(ggplot2)
library(multicore)
library(KCsmart)

load('marieke_diffResults.Rda')
load('DNACopy/mariekeData.Rda')

plotProfilesDiff <- function (NR, KC, KCdiff, chr=NA, setcex=1) {
  
  require(KCsmart)
  data(hsMirrorLocs)
  
  tumName <- sampleInfo$File_name[sampleInfo$NR == as.numeric(NR) & sampleInfo$Type == 'Tumor']
  LNName <- sampleInfo$File_name[sampleInfo$NR == as.numeric(NR) & sampleInfo$Type == 'LN']
    
  KCTLN <- KC[,c('chrom', 'maploc', tumName, LNName)]
  
  KCTLN <- cbind(KCTLN, KCdiff[,NR])
  
  if (is.na(chr)) {
    fileName <- paste('Figures/', NR, 'diffgroup.png', sep='')
  }
  else {
   fileName <- paste('Figures/', NR, 'diffgroup', chr,'.png', sep='')
  }
  
  png(file=fileName, width=800, height=1000)
  par(mfrow=c(3,1))
    plotRawCghDotPlot(KCdataSet=KCTLN, mirrorLocs=hsMirrorLocs, samples=1, doFilter=T, 
      plotTitle=paste('Tumor', NR), setylim=c(-1,1), chromosomes=chr, setcex=setcex)
    plotRawCghDotPlot(KCdataSet=KCTLN, mirrorLocs=hsMirrorLocs, samples=2, doFilter=T, 
      plotTitle=paste('LN', NR), setylim=c(-1,1),, chromosomes=chr, setcex=setcex)
    plotRawCghDotPlot(KCdataSet=KCTLN, mirrorLocs=hsMirrorLocs, samples=3, doFilter=T, 
      plotTitle=paste('Diff', NR), setylim=c(-1,1),, chromosomes=chr, setcex=setcex)
  dev.off()
  
}


# Load clinical data
sampleInfo <- read.delim('Clin_data_tumorLN.txt', stringsAsFactors=F)
sampleInfo <- sampleInfo[order(sampleInfo$File_name),]

# Remove samples 592 and 456
KC <- KC[,c('chrom', 'maploc', sampleInfo$File_name[which(!sampleInfo$NR %in% c(592, 456))])]
sampleInfo <- sampleInfo[-which(sampleInfo$NR %in% c(592, 456)),]

NRlist <- as.list(colnames(KCdiff)[3:ncol(KCdiff)])

mclapply(NRlist, plotProfilesDiff, KC, KCdiff)
mclapply(NRlist, plotProfilesDiff, KC, KCdiff, chr=22, setcex=3)
mclapply(NRlist, plotProfilesDiff, KC, KCdiff, chr=23, setcex=3)
mclapply(NRlist, plotProfilesDiff, KC, KCdiff, chr=12, setcex=3)

