##----------------------------------------------------------------------
## Commands Marieke Lymph node to primary tumor comparison
## 
## Data loading
## 
## Christiaan Klijn
##----------------------------------------------------------------------

# Data was preprocessed by bash script in the same folder
# This was done to extract only the relevant columns from the DNA copy files

# Working directory - work on Centroid

setwd('/data/people/klijn/marieke/DNACopy')

# Load files 
fileNamesData <- dir(pattern='reduce.txt')

for (i in 1:length(fileNamesData)) {
	cat(fileNamesData[i], '\n')
	tempFrame <- read.delim(fileNamesData[i], stringsAsFactors=F)
	if (i == 1) {
		resultFrame <- tempFrame
		colnames(resultFrame)[4] <- fileNamesData[i]
	}
	resultFrame[,i+3] <- tempFrame[,4]
	colnames(resultFrame)[i+3] <- fileNamesData[i]
}

# Remove all info after the first underscore from the colum names
colnames(resultFrame) <- sub('_.+$', '', colnames(resultFrame))

# Make a KC data frame of the data
KC <- resultFrame[,-1]
colnames(KC) <- c('chrom', 'maploc', colnames(KC)[3:ncol(KC)])
KC$chrom <- gsub('chr', '', KC$chrom)
KC$chrom <- gsub('X', '23', KC$chrom)
KC$chrom <- gsub('Y', '24', KC$chrom)
KC$chrom <- as.numeric(KC$chrom)

save(file='mariekeData.Rda', list=c('KC'))
