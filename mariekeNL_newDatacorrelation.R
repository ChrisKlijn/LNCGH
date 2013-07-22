# mariekeNL_correlation.R
# -------------------------------------------------------------------
# Copyright 2013 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke Lymph node vs primary
# Description: Data to do correlation plots between tumors
# -------------------------------------------------------------------

projectDir <- '~/Projects/NKIProjects/LNCGH/newData/'
setwd(projectDir)

load(file=file.path(projectDir, '../newDataImport.RData'))
sampleInfoTumor <- sampleInfo[sampleInfo$Type %in% "Tumor",]

# Exclude sample 592 and 456 as they failed qc

KC    <- segData[,!grepl('592|456', colnames(segData))]
noAnnCols <- grep('Tumor|LN|chrom|maploc', colnames(KC), invert=TRUE)
colnames(KC)[noAnnCols] <- paste0(colnames(KC)[noAnnCols], 'Tumor')

KCLN  <- KC[,grepl('maploc|chrom|LN', colnames(KC))]
KCTum <- KC[,!grepl('LN', colnames(KC))]
colnames(KCLN)  <-  gsub('LN', '', colnames(KCLN))
colnames(KCTum) <-  gsub('Tumor', '', colnames(KCTum))

all.equal(colnames(KCLN),colnames(KCTum))
KCcollTLN <- calcSpmCollection(data=KCTum, mirrorLocs=hsMirrorLocs, data2=KCLN)

#------------------------------------------------------------------
# Correlation Plots
#------------------------------------------------------------------

library(gplots)

sampNames <- c(paste0(colnames(KCTum[3:ncol(KCTum)]), 'Tumor'),
               paste0(colnames(KCLN[3:ncol(KCLN)])  , 'LN'))

TLNlabs <- as.factor(gsub('[0-9]*', '', sampNames))
sampCorMat <- cor(KCcollTLN@data, use='na.or.complete')

str(KCcollTLN)

# pdf(file='Figures/corrMat_T_LN.pdf', width=7, height=7)
pdf(file='Figures/corrMat_T_LN_newData.pdf', width=10, height=10)
heatmap(sampCorMat, scale='none', labRow=sampNames, labCol=sampNames, 
  col=colorpanel(265, low='blue', high='yellow'), 
  ColSideColors=colors()[c(122, 148)][as.numeric(TLNlabs)],
  main='Clustered Correlation Matrix - per-sample KCsmart curve')
dev.off()

#------------------------------------------------------------------
# Dotplots of non-concordant samples
#------------------------------------------------------------------

library(KCsmart)
data(hsMirrorLocs)

# Patient samples 592 and 456 show disconcordance - plot them in a figure

KC592 <- KC[,c(1,2, which(sampleInfo$NR == 592)+2)]

png(file='Figures/Patient592.png', width=1600, height=800*(ncol(KC592)-2))
  par(mfrow=c(ncol(KC592)-2, 1))
  for (p in 1:(ncol(KC592)-2)) {
    sampleNr <- which(sampleInfo$File_name == colnames(KC592)[p+2])
    plotRawCghDotPlot(KCdataSet=KC592, samples=p,
      mirrorLocs=hsMirrorLocs, doFilter=T, 
      plotTitle=paste(sampleInfo[sampleNr, c('File_name','NR','Type') ], collapse=' '))
  }
dev.off()


KC456 <- KC[,c(1,2, which(sampleInfo$NR == 456)+2)]
png(file='Figures/Patient456.png', width=1600, height=800*(ncol(KC456)-2))
  par(mfrow=c(ncol(KC456)-2, 1))
  for (p in 1:(ncol(KC456)-2)) {
    sampleNr <- which(sampleInfo$File_name == colnames(KC456)[p+2])
    plotRawCghDotPlot(KCdataSet=KC456, samples=p,
      mirrorLocs=hsMirrorLocs, doFilter=T, 
      plotTitle=paste(sampleInfo[sampleNr, c('File_name','NR','Type') ], collapse=' '))
  }
dev.off()



