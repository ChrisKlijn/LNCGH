# mariekeNL_plotBiphasicSamples.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Get high-res plots of the biphasic tumors
# -------------------------------------------------------------------

# WD - medoid
setwd('/home/klijn/data/smallproj/mariekeLN/')
source('~/codeChris/generalFunctionsR/chris_cghdata_analysis.R')
source('~/codeChris/generalFunctionsR/chris_DNAcopy_utils.R')

library(ggplot2)
library(multicore)
library(KCsmart)

data(hsMirrorLocs)

load('marieke_diffResults.Rda')
load('DNACopy/mariekeData.Rda')

samplesToPlot <- c('455138A02', '455952A01', '455952A02')
names(samplesToPlot) <- c('Meta','Expansive','Diffuse')

plotKC <- KC[,c('chrom','maploc', samplesToPlot)]

for (i in 1:3) {
  tiff(file=paste('Figures/', names(samplesToPlot)[i], '.tiff', 
    sep=''), width=3000, height=1000, pointsize=30)
  plotRawCghDotPlot(KCdataSet=plotKC, mirrorLocs=hsMirrorLocs, samples=i,
    doFilter=T, plotTitle=names(samplesToPlot)[i], setylim=c(-1,1),
    setcex=3, filterSize=20)
  dev.off()
}


