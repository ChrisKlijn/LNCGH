# Code for revision manuscript
# Load data
# Christiaan Klijn

projectDir <- '~/Projects/NKIProjects/LNCGH/newData/'
setwd(projectDir)

dir()

segData    <- read.delim('130405LN-T_SegmentedRatios.txt', 
                      stringsAsFactors=FALSE)
rawData    <- read.delim('130405LN-T_LogRatios.txt', 
                      stringsAsFactors=FALSE)
sampleInfo <- read.delim('130405LN-T_Clinical.txt',
                         stringsAsFactors=FALSE)
colnames(segData) <- gsub('X', '', colnames(segData))
colnames(rawData) <- gsub('X', '', colnames(rawData))

save(file='../newDataImport.RData', list=c('segData',
                                           'rawData',
                                           'sampleInfo'))
