setwd('~/Projects/NKIProjects/LNCGH/DNAcopyData/')

files <- dir(pattern='txt')

tempData <- read.delim(files[1], stringsAsFactors=FALSE)
probes <- tempData$PROBE_ID

dataList <- vector(mode='list', length=length(files))
names(dataList) <- unlist(lapply(strsplit(
	files, split='_'), function (x) {return(x[1])}))

for (i in 1:length(files)) {
	dataList[[i]] <- read.delim(files[i], stringsAsFactors=FALSE)[,4]
}

dataFr <- as.data.frame(do.call('cbind', dataList))
dataFr <- cbind(probes, dataFr)

write.table(x=dataFr, file='../GEOprobedata.txt', sep='\t', 
	quote=FALSE, row.names=FALSE)

