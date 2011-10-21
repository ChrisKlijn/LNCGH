# mariekeNL_overlapSeg_functions.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: Marieke LN
# Description: Functions to support finding non-overlapping copy
#              number changes
# -------------------------------------------------------------------

getPairRanges <- function(tumSamp, LNSamp, KCseg) {
  # Obtain two ranges objects from segmented output

  require(GenomicRanges)

  tumOutp <- KCseg$output[(KCseg$output$ID %in% tumSamp),]
  LNOutp <- KCseg$output[(KCseg$output$ID %in% LNSamp),]

  tumRanges <- GRanges(seqnames = Rle(paste('chr', tumOutp$chrom)),
    ranges=IRanges(start=tumOutp$loc.start, end=tumOutp$loc.end), 
    seg.mean=tumOutp$seg.mean, num.mark=tumOutp$num.mark)
  LNRanges <- GRanges(seqnames = Rle(paste('chr', LNOutp$chrom)),
    ranges=IRanges(start=LNOutp$loc.start, end=LNOutp$loc.end), 
    seg.mean=LNOutp$seg.mean, num.mark=LNOutp$num.mark)

  returnList <- list(tum=tumRanges, LN=LNRanges)

}

removeSegments <- function(rangeObj, thres=.2, minProbes=10) {
  
  metaData <-  elementMetadata(rangeObj)
  keepVect <- metaData$num.mark > 10 & abs(metaData$seg.mean) > .2
  return(subset(rangeObj, keepVect))

}

removeOverlappingSegments <- function(pairRangeList) {
  
  attach(pairRangeList)

  overlapMat <- matchMatrix(findOverlaps(subject=tum, query=LN))
     
  tumSpec <- tum[seq(1,length(tum))[-(overlapMat[,'subject'])]]
  LNSpec  <- LN[seq(1,length(LN))[-(overlapMat[,'query'])]]
  
  detach(pairRangeList)

  return(list(tum=tumSpec, LN=LNSpec))
}

nonOverlapWrapper <- function(patNR, KCseg, sampleInfo, 
  thres=.2, minprobes=10) {
  # Wrapper function to find pairwise non-overlapping segments

  tumSamp <- sampleInfo$File_name[sampleInfo$NR == patNR & 
    sampleInfo$Type == 'Tumor']
  LNSamp <- sampleInfo$File_name[sampleInfo$NR == patNR & 
    sampleInfo$Type == 'LN']
  pairRanges <- getPairRanges(tumSamp=tumSamp, LNSamp=LNSamp, KCseg=KCseg)
  pairRanges.remove <- lapply(pairRanges, removeSegments, 
    thres=.2, minProbes=10)
  pairRanges.spec <- removeOverlappingSegments(pairRanges.remove)

  return(list(nonO=pairRanges.spec, all=pairRanges))
}

plotNonOverlap <- function(patientNR, segList, KC, 
  sampleInfo, probeRanges) {
  
  tumSamp <- sampleInfo$File_name[sampleInfo$NR == patientNR & 
    sampleInfo$Type == 'Tumor']
  LNSamp <- sampleInfo$File_name[sampleInfo$NR == patientNR & 
    sampleInfo$Type == 'LN']

  pairList <- segList[[paste('P', patientNR, sep='')]]$nonO
  allRanges <- Reduce(pairList, f='c')
  fileName <- paste('Figures/', patientNR, 'segs', '.pdf', sep='')

  pdf(file=fileName, width=5, height=5)
  for (i in 1:length(allRanges)) {
    plotProbes <- matchMatrix(findOverlaps(subject=allRanges[i],
      query=probeRanges))[,'query']
    rangeString <- paste(KC$chrom[plotProbes[1]], 
      KC$maploc[plotProbes[1]], 
      KC$maploc[plotProbes[length(plotProbes)]], sep='-')
    plot(KC$maploc[plotProbes], KC[plotProbes,tumSamp], type='n',
      main=rangeString)
    
    # Plot segments
    lines(KC$maploc[plotProbes], runmed(x=KC[plotProbes,tumSamp], 
      k=min(11, length(plotProbes))), 
      col='blue', lwd=1)
    abline(h=mean(KC[plotProbes,tumSamp]), col=colors()[107])
    lines(KC$maploc[plotProbes], runmed(x=KC[plotProbes,LNSamp], 
      k=min(11, length(plotProbes))), 
      col='orange', lwd=1)
    abline(h=mean(KC[plotProbes,LNSamp]), col=colors()[571])
  }
  dev.off()

}

plotSegsPerChrom <- function (outputFrame, chromNum, KC) {
  
  require(gplots)

  chrsubset <- subset(outputFrame, chrom==chromNum)
#  IDcolors <- brewer.pal(length(unique(chrsubset$ID)),"Dark2")
  IDcolors <- rich.colors(n=length(unique(chrsubset$ID)))
  names(IDcolors) <- unique(chrsubset$ID)

  # initialize plot

  with(chrsubset, plot(0, 0, type='n', 
    xlim=c(min(loc.start), max(loc.end)),
    ylim=c((min(seg.mean, 0)-abs(min(seg.mean))*.2), 
      (max(seg.mean, 0) + abs(max(seg.mean)) * .2)),
    main=paste('Chromosome', chrom[1]), adj=1))
  abline(h=0, col='black')

  for (i in 1:nrow(chrsubset)) {
    
    with(chrsubset, lines(
      x=c(loc.start[i], loc.end[i]),
      y=c(seg.mean[i], seg.mean[i]),
      col=IDcolors[as.character(ID[i])],
      lwd=5))
    probeInd <- with(chrsubset, 
      loc.start[i] < KC$maploc & loc.end[i] > KC$maploc & 
      chrom[i]==KC$chrom)
    with(chrsubset, points(
      x=KC$maploc[probeInd], y=KC[probeInd, as.character(ID[i])],
      col=IDcolors[as.character(ID[i])], pch='.'))

  }

  #Grab the plotting region dimensions
  rng <- par("usr")
  #Call your legend with plot = FALSE to get its dimensions
  lg <- legend(rng[1],rng[2], names(IDcolors), lwd=5, col = IDcolors,
    plot = FALSE, horiz=T)
  #Once you have the dimensions in lg, use them to adjust
  # the legend position
  #Note the use of xpd = NA to allow plotting outside plotting region             
  legend(rng[1],rng[4] + lg$rect$h, names(IDcolors), lwd=5, 
    col = IDcolors,plot = TRUE, xpd = NA, horiz=T)


}

plotSegsID <- function (outputFrame, KC, Type='Tumor', IDcolors) {
  
  thisID <- outputFrame$ID[1]

  hybID <- sampleInfo$File_name[sampleInfo$NR == as.numeric(thisID) &
    sampleInfo$Type == Type]
 
  with(chrsubset, plot(0, 0, type='n', 
    xlim=c(min(loc.start), max(loc.end)),
    ylim=c(-1, 1),
    main=paste(Type, thisID), adj=0))
  abline(h=0, col='black')
  
  for (i in 1:nrow(outputFrame)) {
    probeInd <- with(outputFrame, 
      loc.start[i] < KC$maploc & loc.end[i] > KC$maploc & 
      chrom[i]==KC$chrom)    
    lines(x=c(outputFrame$loc.start[i], outputFrame$loc.end[i]),
      y=c(mean(KC[probeInd, hybID]), mean(KC[probeInd, hybID])),
      col=IDcolors[thisID],
      lwd=5)
    points(x=KC$maploc[probeInd], y=KC[probeInd, hybID],
      col=IDcolors[thisID], pch='.')

  }

}

