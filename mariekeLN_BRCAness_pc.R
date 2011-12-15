## Philip Schouten December 13th 2011
## map 720k nimblegen profiles to 135k nimblegen probes. 135k are a subset of the 720k. 
## create nimblescan output files of the mapped profiles for classification with the nimblegen
## array CGH classifier.

## adjust to folder
setwd('D:/MariekeLN')

# load KCsmart for plotting function
require(KCsmart)
data(hsMirrorLocs)

## unavg DNAcopy file of 135k profile
tmp135 <- read.table('ng135kprofile_unavgDNAcopy.txt',
  sep="\t", header=T, stringsAsFactors=F)

# create identifier
map135 <- paste(tmp135$CHROMOSOME, tmp135$CHR_POSITION, sep='-')

# clinical data to match array to sample number
clin <- read.table('Clin_data_tumorLN-131211.txt', sep="\t", header=T, stringsAsFactors=F)
# for plotting etc
clin$id <- paste(clin$NR, clin$Type, sep="")

# folder with DNAcopy files
setwd('DNACopy')

# get all txt files
files <- dir(pattern='.txt' )
# get short file name to  match clin file
files <- cbind(files, short= gsub('(^[0-9NGA]{4,9}).*', '\\1', files))
# bind name of the sample
lookup <- cbind(files, clin$id[match(files[,'short'], clin$File_name)] ) 
# create folder for mapped profiles and for plots of the mapped profiles
dir.create('mapped')
dir.create('plots')

# loop over all txt files
for (i in 1:nrow(lookup)) {
  # read 720k profile
  tmp <- read.table(lookup[i,1], sep="\t", header=T, stringsAsFactors=F)
  
  # create map of the loaded 720k profiles
  map720 <- paste(tmp$CHROMOSOME, tmp$CHR_POSITION, sep='-')

  # check whether all probes are there should be ...
  print(sum(map720 %in% map135)==length(map135))

  print(all.equal(map720[map720%in%map135], map135))

  print(all.equal(map720[match(map135, map720)], map135))
  
  # the input files do not resemble the nimblescan output files that have to be loaded in the classifier.
  # the added columns do not influence classifier
  mapped <- with(tmp[match(map135, map720),  ], data.frame(CHROMOSOME=CHROMOSOME, CHR_POSITION=CHR_POSITION, RAW_DATAPOINTS=1,
      POSITION_COUNT=1,   RATIO_CORRECTED=RATIO_CORRECTED, WINDOW_SIZE=1, stringsAsFactors=F))

  print(unique(mapped$CHROMOSOME))
  # write the 135k profile extracted from the 720k profile in unavgDNAcopy format
  write.table(mapped, file=paste('mapped/', lookup[i,2], '.txt', sep=''), sep="\t", quote=F, row.names=F)

  # create dataframe for plotting the profiles
  plotDF<-  data.frame(chrom=tmp$CHROMOSOME, maploc=tmp$CHR_POSITION, tmp$RATIO_CORRECTED, stringsAsFactors=F)
  plotDF$chrom <- gsub('chr', '', plotDF$chrom)
  plotDF$chrom <- gsub('X', 23, plotDF$chrom)
  plotDF$chrom <- gsub('Y', 24, plotDF$chrom)
  plotDF$chrom <- as.numeric(plotDF$chrom)
  colnames(plotDF)[3] <- lookup[i,3]
  # and order on chrom and maploc  
  plotDF <- plotDF[order(plotDF$chrom, plotDF$maploc),]

  # plot profiles
  png (filename=paste('plots/', lookup[i,2], '.png', sep=""), height=768, width=1024)
  rawCghDotPlot(plotDF, mirrorLocs=hsMirrorLocs,
      samples=1)
      dev.off()
  
      }

# Perform classification of nimblegen profiles. For this I've loaded the mapped files into
# the classifier tool. There are 3 models currently which are exactly the same as the BRCA1
# 191 (published by Simon) and 371 (published by Marieke) classifers and the BRCA2 classifier.
# Profiles are translated to BAC CGH profiles by averaging the nimblegen probes per BAC clone. 
# The translated profiles are in the folder classification/bacProfiles. In the 
# classification/nimbProfiles is the segmentation of the 135k profiles. 

# 191 classifier:  21318   (cutoff: ?)
# 371 classifier:  1300    (cuoff: 0.63)
# B2 classifier:   10402    (cutoff:0.5)

# 3 quality scores:
# NoiseVariance: variance of ratios around segments
# DynamicRange: at the moment not useful as there is no outlier detection in place. Outliers heavily
# influence this measure
# Signal2NoiseRatio:  segments / noiseVar

# no decision has been made on a good score. Flatliners have low variance and low signal2noise. If noiseVar
# increases in a flatliner it indicates bad profile quality.
# Signal2Noise ration is high in genomic instable tumors, however this also increase noisevariance. Most
# likely the criterium witll end up looking like :  NoiseVar < x OR (signal2Noise > x AND  noiseVar < x )
# for now just plotted both the 135 profiles with segmentation and the mapped BAC profiles.

setwd('../')

# sort folders by filename:
bacFiles <- dir('classification/bacProfiles/', full.names=T)
ratio135Files <- dir('DNAcopy/mapped', pattern='.txt', full.names=T)
seg135Files <- dir('classification/nimbProfiles/', full.names=T)


# BAC nki platform to override genomic position column in translated bac profiles
plf <- read.table('platformnki.txt', sep="\t", header=T, stringsAsFactors=F)

# get BAC ratios and segments
bacratios <- matrix(ncol=length(bacFiles), nrow=3277)
bacsegments <-  matrix(ncol=length(bacFiles), nrow=3277)
colnames(bacratios) <-  gsub('.*/([0-9]{6}A[0-9]{2}|[NGA0-9]{5}).csv$', '\\1', bacFiles)
colnames(bacsegments) <- gsub('.*/([0-9]{6}A[0-9]{2}|[NGA0-9]{5}).csv$', '\\1', bacFiles)

for (i in 1:length(bacFiles)) {
  tmp <- read.csv(bacFiles[i], stringsAsFactors=F, skip=9)
  #watch out columns are 1 off !   
  bacratios[,i] <- tmp[,4]
  bacsegments[,i] <- tmp[,6]
  }

ratios135 <- matrix(ncol=length(ratio135Files), nrow=134937)
segments135 <- matrix(ncol=length(seg135Files), nrow=134937)
colnames(ratios135) <- gsub('.*/([0-9]{6}A[0-9]{2}|[NGA0-9]{5}).txt$', '\\1', ratio135Files)
colnames(segments135) <- gsub('.*/([0-9]{6}A[0-9]{2}|[NGA0-9]{5}).txt$', '\\1', seg135Files)

all.equal(colnames(ratios135), colnames(segments135))
all.equal(colnames(bacratios), colnames(ratios135))
all.equal(colnames(bacsegments), colnames(bacratios))
bacratios <- bacratios[, match(colnames(bacratios), colnames(ratios135))]

for (i in 1:length(seg135Files)) {
  # again columns are off due to read delim
  tmp <- read.delim(seg135Files[i], row.names=NULL, stringsAsFactors=F)
  
  segments135[,i] <- tmp[,3]
  }
  
for (i in 1:length(ratio135Files)) {
    # again columns are off due to read delim
  tmp <- read.delim(ratio135Files[i], header=T, stringsAsFactors=F)
  
  ratios135[,i] <- tmp$RATIO_CORRECTED
  }

rat135 <- data.frame(chrom=tmp135$CHROMOSOME, maploc=tmp135$CHR_POSITION, ratios135, stringsAsFactors=F)
seg135 <- data.frame(chrom=tmp135$CHROMOSOME, maploc=tmp135$CHR_POSITION, segments135, stringsAsFactors=F)
rat3277 <- data.frame(chrom=plf$chrom, maploc=plf$maploc, bacratios, stringsAsFactors=F)
seg3277 <- data.frame(chrom=plf$chrom, maploc=plf$maploc, bacsegments, stringsAsFactors=F)

rat135$chrom <- fixChrom(rat135$chrom)
seg135$chrom <- fixChrom(seg135$chrom)

dir.create('translationPlots')
  
for (i in 3:(ncol(rat135))) {
  png(paste('translationPlots/', colnames(rat135)[i],'.png', sep=""), height=768, width=1024)
  par(mfrow=c(2,1))
  rawCghDotPlot(rat135, samples=(i-2), mirrorLocs=hsMirrorLocs, doFilter=T, filterSize=21)
  addCghDotPlot(seg135, samples=(i-2), mirrorLocs=hsMirrorLocs, setcex=2, rainbow=F)
  
  rawCghDotPlot(rat3277, samples=(i-2), mirrorLocs=hsMirrorLocs)
  addCghDotPlot(seg3277, samples=(i-2), mirrorLocs=hsMirrorLocs, setcex=2, rainbow=F) 
  
  dev.off()
  }

# plots look good.

# read the output of the classification (converted txt to csv because read.delim doesn't get the number
# of columns right and shifts the table):
classification <- read.csv('classification/clsfr_results.csv', header=T, stringsAsFactors=F)

classification$array <- gsub('.*\\\\([0-9]{6}A[0-9]{2}|[NGA0-9]{5}).txt$', '\\1', classification$FileName)
clin <- cbind(clin, classification[match(clin$File_name, classification$array), c(5,6,7)])

# either B1 like (371 classifier) or B2 like.
clin$BRCAnessPCS <- clin$B1.B1Train1300. > 0.63 |  clin$B2.B2Train10402. >0.5


with(clin[clin$BRCAness!="",], table(BRCAness, BRCAnessPCS))

write.table(clin, 'ClinPCS151211.txt', quote=F, row.names=F, sep="\t")


#################################################################################
#################################################################################
###         FUNCTIONS:
#################################################################################


# function to draw CGH plot. Adjusted from Chris Klijn
rawCghDotPlot <- function (KCdataSet, mirrorLocs, samples=1, doFilter=F, filterSize=10, 
chromosomes=NA, setcex=1,setpch='.', plotTitle=NA, rainbow=T, color='red', useRange=c(-3,5), yLab='log2') {


  library('KCsmart')
  library('calibrate')
  library('fields')
  data(hsMirrorLocs)

# First remove missing values
	if (any(is.na(KCdataSet[, samples+2]))) {
		KCdataSet <- KCdataSet[-which(is.na(KCdataSet[,samples+2])),]
	}
## This function plot the raw ratios on a genomic axis

## Standard color for every chromosome
  if (rainbow==T) {
    chromCols <- rainbow(length(unique(KCdataSet$chrom)))
  }
  
  if (rainbow==F) {
    chromCols <- rep(color,24)
  }

## If chromosomes are selected, use those. Otherwise use all chromosomes
## in the dataset
	if (any(is.na(chromosomes))) {
		uniqChroms <- unique(KCdataSet$chrom)
	}
	else {
		uniqChroms <- chromosomes
	}
	
	uniqChroms <- sort(uniqChroms)
	mirrorLocs <- mirrorLocs[uniqChroms]
	attr(mirrorLocs, 'chromNames') <- as.character(uniqChroms)
	
## Get the end coordinates of the chromosomes, cumulative
	chromEnds <- cumsum(unlist(lapply(mirrorLocs, max)))
## Construct a conlinear vector
	KCdataSet <- KCdataSet[KCdataSet$chrom %in% uniqChroms,]
	KCdataSet <- KCdataSet[order(KCdataSet$chrom, KCdataSet$maploc),]
	maplocLin <- KCdataSet$maploc
	

	for (i in 1:length(uniqChroms)) {
		maplocLin[KCdataSet$chrom == uniqChroms[i]] <- 
			maplocLin[KCdataSet$chrom == uniqChroms[i]] + c(0,chromEnds)[i]
	}

## Use a running mean filter if selected, with supplied (or standard)
## filter size
  if (useRange == T) {
    if (doFilter) {
		  filteredData <-runmed(KCdataSet[,samples+2], filterSize)
		  dataRange <- range(filteredData[!is.na(filteredData)])
    }
    else {
		  dataRange <- range(KCdataSet[,samples+2])
      }

  }

  if (! useRange == T) {
    dataRange <- useRange
  }

  if (is.na(plotTitle)) {
    plotTitle=colnames(KCdataSet)[samples+2]
  }
  

    
  plot(c(0, max(maplocLin)), dataRange,type='n', xlab='Genomic Position (bp)', ylab=yLab, main=plotTitle)


	
	for (i in 1:length(uniqChroms)) {
		if (doFilter) {
			points(maplocLin[KCdataSet$chrom == uniqChroms[i]], runmed(KCdataSet[KCdataSet$chrom == uniqChroms[i],samples+2], 
				filterSize), col='black', pch=setpch, cex=setcex)
		}
		else {
			points(maplocLin[KCdataSet$chrom == uniqChroms[i]], KCdataSet[KCdataSet$chrom == uniqChroms[i],samples+2], 
				col='black', pch=setpch, cex=setcex)
		}
	}
	abline(v=chromEnds, col=colors()[615])
	abline(h=0, col='black')
	textX <- chromEnds - (unlist(lapply(mirrorLocs, max)))/2
	text(textX, max(dataRange)-0.05*max(dataRange), labels=attr(mirrorLocs, 'chromNames'))

}

##########################################################################################
##########################################################################################
# function to add to a CGH plot. Adapted from Chris Klijn.

addCghDotPlot <- function (KCdataSet, mirrorLocs, samples=1, doFilter=F, filterSize=10, 
chromosomes=NA, setcex=1, plotTitle=NA, rainbow=T, color='red', ylims=c(-3,5)) {

  library('KCsmart')
  library('calibrate')
  library('fields')
  data(hsMirrorLocs)




# First remove missing values
	if (any(is.na(KCdataSet[, samples+2]))) {
		KCdataSet <- KCdataSet[-which(is.na(KCdataSet[,samples+2])),]
	}
## This function plot the raw ratios on a genomic axis

## Standard color for every chromosome
	if (rainbow==T) {
    chromCols <- rainbow(length(unique(KCdataSet$chrom)))
  }
  
  if (rainbow==F) {
    chromCols <- rep(color,24)
  }
  
## If chromosomes are selected, use those. Otherwise use all chromosomes
## in the dataset
	if (any(is.na(chromosomes))) {
		uniqChroms <- unique(KCdataSet$chrom)
	}
	else {
		uniqChroms <- chromosomes
	}
	
	uniqChroms <- sort(uniqChroms)
	mirrorLocs <- mirrorLocs[uniqChroms]
	attr(mirrorLocs, 'chromNames') <- as.character(uniqChroms)
	
## Get the end coordinates of the chromosomes, cumulative
	chromEnds <- cumsum(unlist(lapply(mirrorLocs, max)))
## Construct a conlinear vector
	KCdataSet <- KCdataSet[KCdataSet$chrom %in% uniqChroms,]
	KCdataSet <- KCdataSet[order(KCdataSet$chrom, KCdataSet$maploc),]
	maplocLin <- KCdataSet$maploc
	

	for (i in 1:length(uniqChroms)) {
		maplocLin[KCdataSet$chrom == uniqChroms[i]] <- 
			maplocLin[KCdataSet$chrom == uniqChroms[i]] + c(0,chromEnds)[i]
	}

## Use a running mean filter if selected, with supplied (or standard)
## filter size
	if (doFilter) {
		filteredData <- runmed(KCdataSet[,samples+2], filterSize)
		dataRange <- range(filteredData[!is.na(filteredData)])
	}
	else {
		dataRange <- range(KCdataSet[,samples+2])
	}
	if (is.na(plotTitle)) {
    plotTitle=colnames(KCdataSet)[samples+2]
  }
#	plot(c(0, max(maplocLin)), dataRange,type='n', xlab='Genomic Position (bp)', ylab='log2', main=plotTitle)


	
	for (i in 1:length(uniqChroms)) {
		if (doFilter) {
			points(maplocLin[KCdataSet$chrom == uniqChroms[i]], runmed(KCdataSet[KCdataSet$chrom == uniqChroms[i],samples+2], 
				filterSize), col=chromCols[uniqChroms[i]], pch='.', cex=setcex)
		}
		else {
			points(maplocLin[KCdataSet$chrom == uniqChroms[i]], KCdataSet[KCdataSet$chrom == uniqChroms[i],samples+2], 
				col=chromCols[uniqChroms[i]], pch='.', cex=setcex)
		}
	}
	abline(v=chromEnds, col=colors()[615])
	abline(h=0, col='black')
	#textX <- chromEnds - (unlist(lapply(mirrorLocs, max)))/2
	#text(textX, max(dataRange)-0.05*max(dataRange), labels=attr(mirrorLocs, 'chromNames'))

}

################################################################################################################
################################################################################################################


# from cgh.r
fixChrom <- function(chrom) {

         chrom <- gsub('chr', '', chrom)
         chrom <- gsub('X', 23, chrom)
         chrom <- gsub('Y', 24, chrom)
         chrom <- as.numeric(chrom)
         return(chrom)
}

#################################################################################################################