# MIT License
#
# Copyright 2017 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# plotSingleOrganism, plotPairOrganism, and their dependencies.

plotSingleOrganism<-function (outPlot, alignmentQualityFile, meanQualityAllFile, meanQualityAlignedFile, exonIntronFile,
                              exonIntronPerCellFile, startTagTrimFile, polyATagTrimFile, cellBCCountsFile, 
                              readsPerCellBarcodeFile, estimatedNumCells, estimatedNumBeads, selectedCellsFile, 
                              digitalExpressionSummaryFile, basePctMatrixCellFile, basePctMatrixMolecularFile,
                              alignmentQualityByCellFile, beadSynthesisErrorDetailFile,
                              molecularBarcodeDistributionByGeneFile, point.cex=1) {
	options(scipen=10000)
	pdf(outPlot, compress=T)

    plotAlignmentMetrics(
            outPlot=NULL,
            alignmentQualityFile=alignmentQualityFile,
            meanQualityAllFile=meanQualityAllFile,
            meanQualityAlignedFile=meanQualityAlignedFile,
            exonIntronFile=exonIntronFile,
            exonIntronPerCellFile=exonIntronPerCellFile,
            startTagTrimFile=startTagTrimFile,
            polyATagTrimFile=polyATagTrimFile,
            cellBCCountsFile=cellBCCountsFile,
            estimatedNumCells=estimatedNumCells,
            estimatedNumBeads=estimatedNumBeads,
            digitalExpressionSummaryFile=digitalExpressionSummaryFile,
            basePctMatrixCellFile=basePctMatrixCellFile,
            basePctMatrixMolecularFile=basePctMatrixMolecularFile,
            alignmentQualityByCellFile=alignmentQualityByCellFile,
            beadSynthesisErrorDetailFile=beadSynthesisErrorDetailFile,
            selectedCellsFile=selectedCellsFile,
            point.cex=point.cex)

	transcriptDownsamplingFile=makeTranscriptDownsamplingPath(outPlot)
	transcriptQuantileFile=makeTranscriptQuantilePath(outPlot)

	downsampleTranscriptsAndQuantiles(
		outDownsamplingFile=transcriptDownsamplingFile,
		outQuantileFile=transcriptQuantileFile,
		molecularBarcodeDistributionByGeneFile=molecularBarcodeDistributionByGeneFile,
		estimatedNumCells=estimatedNumCells)

    plotStandardAnalysisSingleOrganism(
        outPlot=NULL,
        digitalExpressionSummaryFile=digitalExpressionSummaryFile,
        molecularBarcodeDistributionByGeneFile=molecularBarcodeDistributionByGeneFile,
		numCells=estimatedNumCells,
		transcriptDownsamplingFile=transcriptDownsamplingFile,
		transcriptQuantileFile=transcriptQuantileFile,
        point.cex=point.cex)

	dev.off()

}

plotAlignmentMetrics<-function(outPlot, alignmentQualityFile, meanQualityAllFile, meanQualityAlignedFile,
    exonIntronFile, exonIntronPerCellFile=NULL, startTagTrimFile, polyATagTrimFile, cellBCCountsFile,
    estimatedNumCells, estimatedNumBeads, digitalExpressionSummaryFile=NULL, basePctMatrixCellFile, basePctMatrixMolecularFile,
    alignmentQualityByCellFile, beadSynthesisErrorDetailFile, organism1=NULL, organism2=NULL,
    readsPerCellBCOrganismFile1=NULL, readsPerCellBCOrganismFile2=NULL,
    digitalExpressionSummaryFile1=NULL, digitalExpressionSummaryFile2=NULL,
    exonIntronPerCellFile1=NULL, exonIntronPerCellFile2=NULL,
    selectedCellsFile=NULL, cellTypesFile=NULL, point.cex=1) {

    prevScipen = options(scipen=10000)[[1]]
    if (!is.null(outPlot)) {
        pdf(outPlot, compress=T)
    }
    readsPerCellBarcodeFile = cellBCCountsFile

    totalNumReads=plotAlignmentQuality(alignmentQualityFile)[["Total"]]
    plotMeanQualityPerCycle(meanQualityAllFile, meanQualityAlignedFile)

    plotExonIntronFraction(exonIntronFile, "all reads")

    plotStartTagTrimming(startTagTrimFile, totalNumReads)
    plotPolyATagTrimming(polyATagTrimFile, totalNumReads)

	plotNumCellBarcodes(cellBCCountsFile, readsPerCellBarcodeFile, xlimit=c(1, estimatedNumCells*2), selectedCellsFile)
	plotNumCellBarcodes(cellBCCountsFile, readsPerCellBarcodeFile, xlimit=c(1, estimatedNumBeads*2))

    if (!is.null(digitalExpressionSummaryFile)) {
        plotGenesPerCellBarcodeTwoRanges(digitalExpressionSummaryFile, xlimits=c(estimatedNumCells, estimatedNumBeads), point.cex=point.cex)
    }

    if (!is.null(digitalExpressionSummaryFile1)) {
        plotGenesPerCellBarcodeTwoRanges(digitalExpressionSummaryFile1, xlimits=c(estimatedNumCells, estimatedNumBeads), point.cex=point.cex, organism=organism1)
    }

    if (!is.null(digitalExpressionSummaryFile2)) {
        plotGenesPerCellBarcodeTwoRanges(digitalExpressionSummaryFile2, xlimits=c(estimatedNumCells, estimatedNumBeads), point.cex=point.cex, organism=organism2)
    }

    plotNumReadsPerCellBarcode(readsPerCellBarcodeFile, estimatedNumCells, point.cex=point.cex)

    plotBarcodeDegeneratePlots(beadSynthesisErrorDetailFile)
    plotBarcodeBaseDist(basePctMatrixCellFile, bases=c("A", "C", "G", "T"), strTitle="Cell Barcodes")
    plotBarcodeBaseDist(basePctMatrixMolecularFile, bases=c("A", "C", "G", "T"), strTitle="Molecular Barcodes")
    plotPerCellAlignment(alignmentQualityByCellFile, readsPerCellBarcodeFile, estimatedNumCells= estimatedNumCells, titleStr="", cellTypesFile=NULL, organism=NULL)

    if (!(is.null(organism1) | is.null(organism2) | is.null(readsPerCellBCOrganismFile1) | is.null(readsPerCellBCOrganismFile2) |
            is.null(exonIntronPerCellFile1) | is.null(exonIntronPerCellFile2))) {
        df=getNumReadsPerCellBarcodeByOrganismPair(readsPerCellBCOrganismFile1, readsPerCellBCOrganismFile2, organism1, organism2)
        plotNumReadsPerCellBarcodeByOrganismPair(df, estimatedNumCells, organism1, organism2, plotRemainder=F, reportDir=NULL, point.cex=point.cex)
        plotNumReadsPerCellBarcodeByOrganismPair(df, estimatedNumBeads, organism1, organism2, plotRemainder=F, reportDir=NULL, point.cex=point.cex)
        rm (df)
        # Since the exonIntronPerCellFile is only selected cells, plot all of them.
        plotPerCellAnnos(exonIntronPerCellFile1, readsPerCellBarcodeFile, titleStr="", cellTypesFile= cellTypesFile, organism=organism1)
        plotPerCellAnnos(exonIntronPerCellFile2, readsPerCellBarcodeFile, titleStr="", cellTypesFile= cellTypesFile, organism=organism2)
    } else {
        # Since the exonIntronPerCellFile is only selected cells, plot all of them.
        plotPerCellAnnos(exonIntronPerCellFile, readsPerCellBarcodeFile, titleStr="", cellTypesFile=NULL, organism=NULL)
    }
    if (!is.null(outPlot)) {
        dev.off()
    }
    options(scipen=prevScipen)
}


plotPairOrganism<-function (outPlot, alignmentQualityFile, meanQualityAllFile, meanQualityAlignedFile, exonIntronFile,
                            exonIntronPerCellFile, startTagTrimFile, polyATagTrimFile, cellBCCountsFile, 
                            readsPerCellBarcodeFile, estimatedNumCells, estimatedNumBeads, selectedCellsFile, 
                            digitalExpressionSummaryFile, basePctMatrixCellFile, basePctMatrixMolecularFile,
                            alignmentQualityByCellFile, beadSynthesisErrorDetailFile,
                            digitalExpressionSummaryFile1, digitalExpressionSummaryFile2,
                            readsPerCellBCOrganismFile1, readsPerCellBCOrganismFile2, cellTypesFile,
                            organism1, organism2,
                            molecularBarcodeDistributionByGeneFile1, molecularBarcodeDistributionByGeneFile2, point.cex=1) {
	options(scipen=10000)
	pdf(outPlot, compress=T)

    plotAlignmentMetrics(
        outPlot=NULL,
        alignmentQualityFile=alignmentQualityFile,
        meanQualityAllFile= meanQualityAllFile,
        meanQualityAlignedFile=meanQualityAlignedFile,
        exonIntronFile=exonIntronFile,
        exonIntronPerCellFile=exonIntronPerCellFile,
        startTagTrimFile=startTagTrimFile,
        polyATagTrimFile=polyATagTrimFile,
        cellBCCountsFile=cellBCCountsFile,
        estimatedNumCells=estimatedNumCells,
        estimatedNumBeads=estimatedNumBeads,
        digitalExpressionSummaryFile=digitalExpressionSummaryFile,
        basePctMatrixCellFile=basePctMatrixCellFile,
        basePctMatrixMolecularFile=basePctMatrixMolecularFile,
        alignmentQualityByCellFile=alignmentQualityByCellFile,
        beadSynthesisErrorDetailFile=beadSynthesisErrorDetailFile,
        organism1=organism1,
        organism2=organism2,
        readsPerCellBCOrganismFile1=readsPerCellBCOrganismFile1,
        readsPerCellBCOrganismFile2=readsPerCellBCOrganismFile2,
        selectedCellsFile=selectedCellsFile,
        cellTypesFile=cellTypesFile,
        point.cex=point.cex)

    plotStandardAnalysisPairOrganism(
        outPlot=NULL,
		digitalExpressionSummaryFile1=digitalExpressionSummaryFile1,
		digitalExpressionSummaryFile2=digitalExpressionSummaryFile2,
        molecularBarcodeDistributionByGeneFile1=molecularBarcodeDistributionByGeneFile1,
        molecularBarcodeDistributionByGeneFile2=molecularBarcodeDistributionByGeneFile2,
        organism1=organism1,
        organism2=organism2,
        cellTypesFile=cellTypesFile,
        estimatedNumCells=estimatedNumCells,
        estimatedNumBeads=estimatedNumBeads,
        point.cex=point.cex)

	dev.off()
}


plotAlignmentQuality<-function (alignmentQualityFile) {
	
	summaryStats=read.table(alignmentQualityFile, header=T, stringsAsFactors=F, sep="\t", fill=T, nrows=1)
	readDist=read.table(alignmentQualityFile, header=T, stringsAsFactors=F, sep="\t", fill=T, skip=4)
	
	d=c(summaryStats$totalReads, summaryStats$mappedReads, summaryStats$hqMappedReads, summaryStats$hqMappedReadsNoPCRDupes)
	names (d)=c("Total","Mapped", "HQ", "HQ No Dupes")
	dd=d/1e6
	z=barplot(dd, col="light blue", ylim=c(0,max(dd)*1.1), cex.names=0.65, main="Alignment Quality", ylab="# Reads [millions]")
	text(x=z, y=(as.numeric(dd)*1.05), labels=prettyNum(d,big.mark=",",scientific=FALSE))
	text(x=z, y=(as.numeric(dd)*0.5), labels=paste(round(d/d[1]*100,1),"%", sep=""))
	
	plot(readDist$read.quality, readDist$num.reads/sum(as.numeric(readDist$num.reads))*100, type='l', col="blue", lwd=3, xlab="map quality", ylab="pct of reads", main="Map quality distribution")
	return (d)
	
}


plotMeanQualityPerCycle<-function (meanQualityAllFile, meanQualityAlignedFile) {
	a=read.table(meanQualityAllFile, header=T, stringsAsFactors=F)
	b=read.table(meanQualityAlignedFile, header=T, stringsAsFactors=F)
	maxY=max(a$MEAN_QUALITY, b$MEAN_QUALITY)
	plot(a$CYCLE, a$MEAN_QUALITY, col="green", type='l', ylim=c(0, maxY), xlab="Read Cycle", ylab="mean base quality")
	points(b$CYCLE, b$MEAN_QUALITY, col="blue", type='l', ylim=c(0, maxY))
	legend("bottomright", c("all reads", "aligned reads"), fill=c("green", "blue"), cex=0.75)
}

plotExonIntronFraction<-function (exonIntronFile, stringTitle="") {
	a=read.table(exonIntronFile, header=T, stringsAsFactors=F, nrow=1, fill=T, sep="\t")
	r=c(genic=sum(a$PCT_CODING_BASES, a$PCT_INTRONIC_BASES, a$PCT_UTR_BASES), exonic=sum(a$PCT_CODING_BASES+a$PCT_UTR_BASES),intronic=a$PCT_INTRONIC_BASES, intergenic=a$PCT_INTERGENIC_BASES, ribsomal=a$PCT_RIBOSOMAL_BASES)
	r=round(r, 3)
	r=r*100
	barplot(r, col="light blue")
	text(seq(0.65,5.5,length.out=5), max(r)/2, r)	
	title(stringTitle)
}

plotStartTagTrimming<-function (startTagTrimFile, totalNumReads) {
	
	a=read.table(startTagTrimFile, header=F, stringsAsFactors=F)[,1:2]
	a=a[order(as.numeric(a[,2])),]
	pctReadsTrimmed=round(sum (a$V1)/totalNumReads*100,2)
	
	plot (a$V2, a$V1, xlab="num bases trimmed 5'", ylab="number of reads", type='l')
	title(paste("% Reads trimmed by 5' trimmer: ", pctReadsTrimmed, sep=""))
}

plotPolyATagTrimming<-function (polyATagTrimFile, totalNumReads) {
	a=read.table(polyATagTrimFile, header=F, stringsAsFactors=F)[,1:2]
	a=a[order(as.numeric(a[,2])),]
	pctReadsTrimmed=round(sum (a$V1)/totalNumReads*100,2)
	
	plot (a$V2, a$V1, xlab="first base of PolyA tail trimmed", ylab="number of reads", type='l')
	title(paste("% Reads trimmed by 3' PolyA trimmer: ", pctReadsTrimmed, sep=""))
}

plotNumCellBarcodes<-function (cellBCCountsFile, cellBCCollapsedCountsFile, xlimit=NULL, selectedCellsFile=NULL) {
	
	a=read.table(cellBCCountsFile, header=F, stringsAsFactors=F)[,1:2]
	cell_barcodes=a[order(a$V1,decreasing=T),]
	
	b=read.table(cellBCCollapsedCountsFile, header=F, stringsAsFactors=F)[,1:2]
	cell_barcodes_collapsed=b[order(b$V1,decreasing=T),]

	# Convert to 64-bit ints to avoid overflow.
	cum=cumsum(as.numeric(cell_barcodes$V1))
	cum=cum/cum[length(cum)]
	
	cumB=cumsum(as.numeric(cell_barcodes_collapsed $V1))
	cumB=cumB/cumB[length(cumB)]
	
	
	numCellsFound=NULL
	
	if (is.null(xlimit)==F) {
		plot(1:length(cum), cum, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", xlim=range(xlimit))
		points(1:length(cumB), cumB, type='l', col="green")
	} else {
		plot(1:length(cum), cum, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads")
		points(1:length(cumB), cumB, type='l', col="green")
	}
	
	if (!is.null(selectedCellsFile)) {
		z=read.table(selectedCellsFile, header=F, stringsAsFactors=F)
		abline(v=dim(z)[1], lwd=2)
		numCellsFound=dim(z)[1]
	}
	
	if (!is.null(numCellsFound)) {
		title(paste("Cumulative fraction of reads per cell barcode, cells found:", numCellsFound))
	} else {
		title("Cumulative fraction of reads per cell barcode")
	}
	
	legend('bottomright', legend=c("cell barcodes", "collapsed cell barcodes"), fill=c("blue", "green"))
	
}


plotGenesPerCellBarcode<-function (digitalExpressionSummaryFile, xlimit, point.cex=0.75, organism=NULL) {
	a=read.table(digitalExpressionSummaryFile, header=T, stringsAsFactors=F)[,'NUM_GENES']
	a=sort(a,decreasing=T)
	if (is.null(xlimit)) {
		plot (a, xlab="cell barcodes sorted by #of genes", ylab="number of genes", pch=16, col="blue", cex=point.cex)
	} else {
		plot (a, xlab="cell barcodes sorted by #of genes", ylab="number of genes", pch=16, col="blue", xlim=range(xlimit), cex=point.cex)
	}
    if (is.null(organism)) {
        suffix = ""
    } else {
        suffix = paste0(" (", organism, ")")
    }
    title(paste0("Number of genes per cell barcode", suffix))
}

plotGenesPerCellBarcodeTwoRanges<-function (digitalExpressionSummaryFile, xlimits, point.cex=0.75, organism=NULL) {
    for (xlimit in xlimits) {
        plotGenesPerCellBarcode(digitalExpressionSummaryFile, xlimit=c(1, xlimit), point.cex=point.cex, organism=organism)
    }
}

plotNumReadsPerCellBarcode<-function (readsPerCellBarcodeFile, estimatedNumCells, plotRemainder=F, point.cex=0.75) {
	
	#readsPerCellBarcodeFile =getNumReadsPerCellBarcode(bamFile, cellTagCollapsed, readMapQuality=10, organism=NULL, tempDir, gapAnalysisDir)
	a=read.table(readsPerCellBarcodeFile, header=F, stringsAsFactors=F)	
	df=data.frame(barcode=a[,2], count=a[,1], stringsAsFactors=F)

	dd=restrictDF(df, estimatedNumCells, plotRemainder)
	if (is.null(dd) | dim(dd)[1]==0) return (NULL)
	if (plotRemainder==F) {
		strTitle="Number of reads per cell barcode"
	}  else {
		strTitle="Number of reads per cell barcode - non-core cell barcodes only"
	}
	
	plot(1:dim(dd)[1], dd$count, type='p', xlab=paste("cell barcodes"), ylab=paste("Number of reads"), pch=16, col="blue", cex=point.cex)
	
	title(strTitle)
}

#read summary instead of DE file (Modified)
getGenesAndTranscriptsPerCellBarcode<-function (digitalExpressionFile) {
	a=read.table(digitalExpressionFile, header=T, stringsAsFactors=F)
	colnames(a)=c("cellBC", "numReads", "numTranscripts", "numGenes")
	#colnames(a)=c("cellBC", "numGenes", "numTranscripts")
	#colnames(a)=c("cellBC", "numTranscripts")
	return (a)
}


plotNumTranscriptsPerCellBarcodeByOrganismPair<-function (df, estimatedNumCells=100, organismOne="HUMAN", organismTwo="MOUSE", plotRemainder=F, reportDir=NULL, point.cex=0.75) {
	if (is.null(organismOne) || is.null(organismTwo)) return(NULL)
		
	dd=restrictDF(df, estimatedNumCells, plotRemainder)
	#stop if there's no data.
	if (dim(dd)[1]==0) return (NULL)
	
	if (plotRemainder==F) {
		strTitle="Number of transcripts per cell barcode"
	}  else {
		strTitle="Number of transcripts per cell barcode - non-core cell barcodes only"
	}
	
	maxR=range(dd[,2], dd[,3], na.rm=T)
	plot(maxR, maxR, type='n', xlab=paste("Number of transcripts", organismOne), ylab=paste("Number of transcripts", organismTwo))
	points(dd[[organismOne]], dd[[organismTwo]], pch=16, col="blue", cex=point.cex)
	
	title(strTitle)
	if (!is.null(reportDir)) {
		colnames(df)[2]= organismOne
		colnames(df)[3]= organismTwo
		outFile=paste(reportDir, "/transcripts_per_cell_barcode.txt", sep="")
		write.table(df, outFile, row.names=F, col.names=T, quote=F, sep="\t")
	}
	
	#plot(df[1: estimatedNumCells,]$o1Count, df[1: estimatedNumCells,]$o2Count, pch=16, xlab=paste("Number of transcripts", organismOne), ylab=paste("Number of transcripts", organismTwo), col="blue")
	#title("Number of transcripts per cell barcode")
}

##########
# SINGLE ORGANISM PLOTS
##########
plotNumTranscriptsPerCellBarcode<-function (df, estimatedNumCells=100, plotRemainder=F, point.cex=0.75) {
	
	df=df[order(df$numTranscripts,decreasing=T),]
		
	dd=restrictDF(df, estimatedNumCells, plotRemainder)
	if (is.null(dd) | dim(dd)[1]==0) return (NULL)
	if (plotRemainder==F) {
		strTitle="Number of transcripts per cell barcode"
	}  else {
		strTitle="Number of transcripts per cell barcode - non-core cell barcodes only"
	}
	
	plot(1:dim(dd)[1], dd$numTranscripts, type='p', xlab=paste("cell barcodes"), ylab=paste("Number of transcripts"), pch=16, col="blue", cex=point.cex)
	
	title(strTitle)
}

plotNumGenesPerCellBarcode<-function (df, estimatedNumCells=100, plotRemainder=F, point.cex=0.75) {
	
	df=df[order(df$numTranscripts,decreasing=T),]
		
	dd=restrictDF(df, estimatedNumCells, plotRemainder)
	if (is.null(dd) | dim(dd)[1]==0) return (NULL)
	if (plotRemainder==F) {
		strTitle="Number of genes per cell barcode"
	}  else {
		strTitle="Number of genes per cell barcode - non-core cell barcodes only"
	}
	
	plot(1:dim(dd)[1], dd$numGenes, type='p', xlab=paste("cell barcodes"), ylab=paste("Number of genes"), pch=16, col="blue", cex=point.cex)
	
	title(strTitle)
}

plotCumulativeReadsPerUMI<-function (umiFile, organismName="", maxReadsPerUMI=NA) {
	transcripts=read.table(umiFile, header=T, stringsAsFactors=F, sep="\t")
	
	numObsTable=table(transcripts$Num_Obs)
	#multiply the number of reads per UMI by the number of observations.
	numObsTable=as.numeric(names (numObsTable))*numObsTable
	cs=cumsum(numObsTable)
	cs=cs/max(cs)
	
	#how many UMIs should we plot?
	maxX= maxReadsPerUMI
	if (is.na(maxReadsPerUMI)) {
		threshold=0.80
                idx=which(cs<=threshold)
                if (length(idx)==0) {
					maxX=length(cs)
				} else {
					maxX=min(10000, max(which(cs<=threshold)))
				}
	}	
	
	
	
	numObsTable2=table(transcripts$Num_Obs)
	cs2=cumsum(numObsTable2)
	cs2=cs2/max(cs2)
	old.par <- par(no.readonly = TRUE)
	on.exit(par(old.par))
	par(mfrow=c(2,1))
	par(mar=c(4,5,2,2))
	
	plot (as.integer(names (cs)), cs, type='l', xlim=c(1,maxX), xlab="Reads per UMI", ylab="Cumulative fraction of reads")
	title (organismName)
	
	plot (as.integer(names (cs2)), cs2, type='l', xlim=c(1,maxX), xlab="Reads per UMI", ylab="Fraction of UMIs")
	invisible()
	
}

plotBarcodeDegeneratePlots<-function (synthesisMetricsFile) {
    if (file.info(synthesisMetricsFile)$size > 0) {
        a=read.table(synthesisMetricsFile, header=T, stringsAsFactors=F)

        #errorBase=8
        plotOne<-function (errorBase, a) {
            cellsDegen=a[a$SYNTH_MISSING_BASE == errorBase,]$CELL_BARCODE
            plotTrimPlot(cellsDegen, experimentName=paste("Degenerate Base", errorBase))
        }
        errorBases=intersect(7:8,unique (a$SYNTH_MISSING_BASE))
        z=lapply(errorBases, plotOne, a)
        cellsNormal=a[a$SYNTH_MISSING_BASE==-1,]$CELL_BARCODE
        plotTrimPlot(cellsNormal, experimentName=paste("Normal Barcodes"))
	}
}

plotTrimPlot<-function (cellBarcodes, experimentName) {
	pardefault <- par(no.readonly = T)
	par(mar=c(3,5,3,5))
	numCells=length(cellBarcodes)
	d= downsampledBarcodes(cellTypesFile=NULL, cellBarcodes=cellBarcodes, numCells)
	par(mfrow=c(2,1))
	plotD(d$startTrimmed, numCells=numCells)
	title(paste(experimentName, "Start Trimmed"))
	plotD(d$endTrimmed, numCells=numCells)
	title(paste("End Trimmed"))
	on.exit(par(pardefault))
}




plotNumReadsPerCellBarcodeByOrganismPair<-function (df, estimatedNumCells, organismOne="HUMAN", organismTwo="MOUSE", plotRemainder=F, reportDir=NULL, point.cex=0.75) {
		
	dd=restrictDF(df, estimatedNumCells, plotRemainder)
	#stop if there's no data.
	if (dim(dd)[1]==0) return (NULL)
	
	if (plotRemainder==F) {
		strTitle="Number of reads per cell barcode"
	}  else {
		strTitle="Number of reads per cell barcode - non-core cell barcodes only"
	}
	maxR=range(dd[,2], dd[,3], na.rm=T)
	plot(maxR, maxR, type='n', xlab=paste("Number of reads", organismOne), ylab=paste("Number of reads", organismTwo), cex=point.cex)
	points(dd[[organismOne]], dd[[organismTwo]], pch=16, col="blue", cex=point.cex)
	
	title(strTitle)
	if (!is.null(reportDir)) {
		colnames(df)[2]= organismOne
		colnames(df)[3]= organismTwo
		outFile=paste(reportDir, "/reads_per_cell_barcode.txt", sep="")
		write.table(df, outFile, row.names=F, col.names=T, quote=F, sep="\t")
	}
}

getNumReadsPerCellBarcodeByOrganismPair<-function(readsPerCellBCOrganism1File, readsPerCellBCOrganism2File, organismOne, organismTwo) {
	
	o1= read.table(readsPerCellBCOrganism1File,header=F, stringsAsFactors=F)
	colnames(o1)=c("count", "barcode")
	o2= read.table(readsPerCellBCOrganism2File,header=F, stringsAsFactors=F)
	colnames(o2)=c("count", "barcode")
	
	commonBC=union(o1$barcode, o2$barcode)
	o1p=o1[match(commonBC, o1$barcode),]
	o2p=o2[match(commonBC, o2$barcode),]
	df=data.frame(tag=commonBC, o1Count=o1p$count, o2Count=o2p$count, stringsAsFactors=F)
	
	idx=which(is.na(df$o1Count))
	if (length(idx)>0) df[idx,]$o1Count=0
	
	idx=which(is.na(df$o2Count))
	if (length(idx)>0) df[idx,]$o2Count=0
	
	
	df$total=apply(df[,2:3], 1, sum, na.rm=T)
	df=df[order(df$total, decreasing=T),]
	colnames(df)[2]= organismOne
	colnames(df)[3]= organismTwo
	return (df)
	
}

getNumGenesPerCellBarcodeByOrganismPair<-function (digitalExpressionFileO1, digitalExpressionFileO2, organismOne, organismTwo) {
	if (is.null(organismOne) || is.null(organismTwo)) return(NULL)
	
	o1=getGenesAndTranscriptsPerCellBarcode(digitalExpressionFileO1)
	o2=getGenesAndTranscriptsPerCellBarcode(digitalExpressionFileO2)

	commonBC=union(o1$cellBC, o2$cellBC)
	o1p=o1[match(commonBC, o1$cellBC),]
	o2p=o2[match(commonBC, o2$cellBC),]
	df=data.frame(tag= commonBC, o1Count=o1p$numGenes, o2Count=o2p$numGenes, stringsAsFactors=F)
	
	idx=which(is.na(df$o1Count))
	if (length(idx)>0) df[idx,]$o1Count=0
	
	idx=which(is.na(df$o2Count))
	if (length(idx)>0) df[idx,]$o2Count=0
	
	
	df$total=apply(df[,2:3], 1, sum, na.rm=T)
	df=df[order(df$total, decreasing=T),]
	return (df)
}

plotNumGenesPerCellBarcodeByOrganismPair<-function (df, estimatedNumCells=100, organismOne="HUMAN", organismTwo="MOUSE", plotRemainder=F, point.cex=0.75) {
	if (is.null(organismOne) || is.null(organismTwo)) return(NULL)
		
	dd=restrictDF(df, estimatedNumCells, plotRemainder)
	if (dim(dd)[1]==0) return (NULL)
	if (plotRemainder==F) {
		strTitle="Number of genes per cell barcode"
	}  else {
		strTitle="Number of genes per cell barcode - non-core cell barcodes only"
	}
	
	maxR=range(dd[,2], dd[,3], na.rm=T)
	plot(maxR, maxR, type='n', xlab=paste("Number of genes", organismOne), ylab=paste("Number of genes", organismTwo))
	points(dd$o1Count, dd$o2Count, pch=16, col="blue", cex=point.cex)
	title(strTitle)
}


# Modified colors
plotCategorizedCellTypes<-function (cellTypesFile, point.cex, organismOne, organismTwo) {
	
	df=read.table(cellTypesFile, header=T, stringsAsFactors=F, check.names = FALSE)
	maxRange=max(df[[organismOne]], df[[organismTwo]])
	#maxRange = 20000
	
	dfOrganismOne=df[df$organism==organismOne,]
	dfOrganismTwo=df[df$organism==organismTwo,]
	dfMixed=df[df$organism=="Mixed",]
	
	plot(dfOrganismOne[[organismOne]], dfOrganismOne[[organismTwo]], col="cornflowerblue", pch=16, xlim=c(0,maxRange), ylim=c(0,maxRange),
	    xlab=paste(organismOne, "transcripts"), ylab=paste(organismTwo, "transcripts"), cex=point.cex,
	    cex.lab=1.1, cex.axis=1.25, cex.main=1.25, cex.sub=1.25)
	points(dfOrganismTwo[[organismOne]], dfOrganismTwo[[organismTwo]], col="firebrick3", pch=16, cex=point.cex)
	points(dfMixed[[organismOne]], dfMixed[[organismTwo]], col="purple", pch=16, cex=point.cex)
	l=c(paste(organismOne, dim(dfOrganismOne)[1]), paste(organismTwo, dim(dfOrganismTwo)[1]), paste("Mixed", dim(dfMixed)[1]))
	legend("topright", legend=l, fill= c("cornflowerblue", "firebrick3", "purple", "grey")) #c("blue", "red", "purple", "grey"))
		
}

plotBarcodeBaseDist<-function (baseDistFile, bases=c("A", "C", "G", "T"), strTitle="") {
	a=read.table(baseDistFile, header=T, stringsAsFactors=F)

    # Convert from integer to numeric to handle number > 32-bit integer
	basePctMatrix =a[,-1]/apply(a[,-1], 1, function(x) sum(as.numeric(x)))*100
	maxY=max(basePctMatrix)
	
	plot(c(1,dim(basePctMatrix)[1]), c(0,maxY), type='n', xlab="base position", ylab="fraction of reads")
	colors=brewer.pal(length(bases), "Set1")
	for (i in 1:length(bases)) {
		b= bases[i]
		points(1:dim(basePctMatrix)[1], basePctMatrix[,b], col=colors[i], type='p', pch=16)
	}
	legend(0.7, maxY*1.15, xpd=NA, legend=bases, fill=colors, ncol=length(bases), cex=0.75)
	title(strTitle)
	
}

plotPerCellAnnos<-function (exonIntronPerCellFile, readsPerCellBarcodeFile, titleStr="", cellTypesFile=NULL, organism=NULL, numCellsToPlot=NULL) {
	
	a=read.table(exonIntronPerCellFile, header=T, stringsAsFactors=F, sep='\t', fill=T)
	a=a[a$SAMPLE!="",]
	r=data.frame(cellBarcode=a$SAMPLE, genic=(a$PCT_CODING_BASES+ a$PCT_INTRONIC_BASES+a$PCT_UTR_BASES), exonic=a$PCT_CODING_BASES+a$PCT_UTR_BASES,intronic=a$PCT_INTRONIC_BASES, intergenic=a$PCT_INTERGENIC_BASES, ribsomal=a$PCT_RIBOSOMAL_BASES, stringsAsFactors=F)
	
	b=read.table(readsPerCellBarcodeFile, stringsAsFactors=F)
	idx=match(r$cellBarcode, b$V2)
	
	r$num_reads =b[idx,]$V1
	
	if (!is.null(cellTypesFile) & !is.null(organism)) {
		z=read.table(cellTypesFile, header=T, stringsAsFactors=F)
		z=z[z$organism==organism,]
                if (nrow(z) == 0) {
                   return(NULL)
                }	        

		barcodes=intersect(r$cellBarcode, z$tag)
		idx=match(barcodes, r$cellBarcode)
		r=r[idx,]
		idx=match(r$cellBarcode, z$tag)
		r$pctHuman=z[idx,]$ratio
		#titleStr=paste(titleStr, organism, sep=" ")
	}
		
	r=r[order(r$num_reads, decreasing=F),]
	if (!is.null(numCellsToPlot)) r=tail(r, n= numCellsToPlot)
	idx=match(c("genic", "exonic", "intronic", "intergenic", "ribsomal"), colnames(r))
	r[,idx]=r[,idx]*100
	colors=brewer.pal(5, "Set2")
	
	#plot(r$num_reads, r$genic, col=colors[1], xlab="# reads", ylab="% annotation", pch=16, cex=0.75)
	plot(1:dim(r)[1], r$genic, col=colors[1], xlab="# cell [sorted smallest to largest]", ylab="% annotation", ylim=c(0, 100), pch=16, cex=0.75)
	points(1:dim(r)[1], r$exonic, col=colors[2], pch=16, cex=0.75)
	points(1:dim(r)[1], r$intronic, col=colors[3], pch=16, cex=0.75)
	points(1:dim(r)[1], r$intergenic, col=colors[4], pch=16, cex=0.75)
	points(1:dim(r)[1], r$ribsomal, col=colors[5], pch=16, cex=0.75)
	legend(x=0, y=110, xpd=NA, legend=c("genic", "exonic", "intronic", "intergenic", "ribosomal"), fill=colors, ncol=5, box.col="white")
	title(paste(titleStr, organism))
	
	#options(scipen=20)
	#plot(r$num_reads, r$genic, col=colors[1], xlab="# reads", ylab="% annotation", ylim=c(0, 100), pch=16, cex=0.75)
	#points(r$num_reads, r$exonic, col=colors[2], pch=16, cex=0.75)
	#points(r$num_reads, r$intronic, col=colors[3], pch=16, cex=0.75)
	#points(r$num_reads, r$intergenic, col=colors[4], pch=16, cex=0.75)
	#points(r$num_reads, r$ribsomal, col=colors[5], pch=16, cex=0.75)
	#legend(x=0, y=110, xpd=NA, legend=c("genic", "exonic", "intronic", "intergenic", "ribosomal"), fill=colors, ncol=5, box.col="white")
	#title(paste(titleStr, organism))
	
	#plot(r$intronic, r$pctHuman, pch=16, cex=0.75, col="blue", xlab="fraction intronic", ylab="fraction human")	
}

restrictDF<-function (df, numRows, remainder=F) {
	n=dim(df)[1]
	if (remainder==T) {
		# empty result
		if (numRows>=n) return(df[0,])
		# have some rows remaining, return them.
		return (df[-1:-numRows,])
	}
	
	if (remainder==F && numRows<n) {
		# have some rows remaining, return them.
		return(df[1:numRows,])
	}
	return (df)
}


#if you're single organism, plot the estimated # cells.
#if you're multi organism, plot the categorized cells for an organism.
#alignmentQualityByCellFile="P5Retina_star.bam_ReadQualityMetricsByCell.txt.gz"
plotPerCellAlignment<-function (alignmentQualityByCellFile, readsPerCellBarcodeFile, estimatedNumCells=NULL, titleStr="", cellTypesFile=NULL, organism=NULL) {
    
    a=read.table(alignmentQualityByCellFile, header=T, stringsAsFactors=F, sep='\t', fill=T)
    idx=which(a$aggregate=="read quality")
    a=a[1:(idx-1),]
    a=a[-which(a$aggregate=="all"),]
    
    aa=data.frame(cellBarcode=a$aggregate, totalReads=as.numeric(a$totalReads), mappedReads=as.numeric(a$mappedReads), hqMappedReads=as.numeric(a$hqMappedReads), stringsAsFactors=F)
    aa=aa[order(aa$totalReads, decreasing=T),]
    aa$pctMapped=aa$mappedReads/aa$totalReads*100
    aa$pctHQ=aa$hqMappedReads/aa$totalReads*100
    
    if (!is.null(estimatedNumCells)) {
        aa=aa[1: estimatedNumCells,]
    }
    
    if (!is.null(cellTypesFile)) {
        z=read.table(cellTypesFile, header=T, stringsAsFactors=F)
        z=z[z$organism==organism,]
        barcodes=intersect(aa$cellBarcode, z$tag)
        idx=match(barcodes, aa$cellBarcode)
        aa=aa[idx,]
    }
    
    
    colors=brewer.pal(5, "Set2")
    maxY=max(aa$pctMapped, aa$pctHQ)
    
    plot (c(1,dim(aa)[1]), c(0,maxY*1.1), xlab="cells ordered by size", ylab="pct aligned", type='n')
    points(1:dim(aa)[1], aa$pctMapped, col=colors[1], pch=16, cex=0.75)
    points(1:dim(aa)[1], aa$pctHQ, col=colors[2], cex=0.75)
    legend("topleft", legend=c("%mapped", "%HQ"), fill=colors[1:2], cex=0.75)
    
    plot (range(aa$totalReads), c(0,maxY*1.1), xlab="total number reads", ylab="pct aligned", type='n')
    points(aa$totalReads, aa$pctMapped, col=colors[1], pch=16, cex=0.75)
    points(aa$totalReads, aa$pctHQ, col=colors[2], cex=0.75)
    legend("topleft", legend=c("%mapped", "%HQ"), fill=colors[1:2], cex=0.75)
    
    
}

# One of estimatedNumCells and selectedCellsFile must be specified
downsampleTranscriptsAndQuantiles<-function(
        outDownsamplingFile,
        outQuantileFile,
        molecularBarcodeDistributionByGeneFile,
        estimatedNumCells=NULL,
        selectedCellsFile=NULL,
        random.seed=1) {
    set.seed(random.seed)
    if (is.null(estimatedNumCells)) {
        estimatedNumCells = length(readLines(selectedCellsFile))
    }
	#I'm capturing returning of these functions so they don't spill in certain logs.
	z=runTranscriptDownsampling(molecularBarcodeDistributionByGeneFile, estimatedNumCells=estimatedNumCells, outFile=outDownsamplingFile)
	z=generateTranscriptQuantileTableSimple(outDownsamplingFile, outQuantileFile)
}

plotSensitivity<-function (downsamplingFile, quantileSequence=seq(0.1, 1, by=0.1), maxExpansion=2.5, readMultiplier=c(0.5, 2, 10), organism=NULL) {
    #load up data and run plots!
    d=read.table(downsamplingFile, header=T, stringsAsFactors=F, check.names=F)
    predictTranscriptsGainedPerRead (downsampledData=d, quantileSequence=quantileSequence, maxExpansion=maxExpansion, excludeLargestSeries=F, organism=organism)
    predictTranscriptsGainedPerRead (downsampledData=d, quantileSequence=quantileSequence, maxExpansion=maxExpansion, excludeLargestSeries=T, organism=organism)
    predictTranscriptsGainedPerReadAsPercent(downsampledData=d, quantileSequence=quantileSequence, maxExpansion=maxExpansion, organism=organism)
    dotChartTranscriptsGain(downsampledData=d, quantileSequence=quantileSequence, readMultiplier=readMultiplier, organism=organism)
}


plotStandardAnalysisSingleOrganism<-function(outPlot, digitalExpressionSummaryFile, molecularBarcodeDistributionByGeneFile,
											transcriptDownsamplingFile, transcriptQuantileFile, numCells=NULL, point.cex=1) {
    prevScipen = options(scipen=10000)[[1]]
    if (!is.null(outPlot)) {
        pdf(outPlot, compress=T)
    }
    df=getGenesAndTranscriptsPerCellBarcode(digitalExpressionSummaryFile)

    if (is.null(numCells)) {
        numCells = nrow(df)
    }

    plotNumTranscriptsPerCellBarcode(df, numCells, plotRemainder=F, point.cex=point.cex)
    plotNumGenesPerCellBarcode(df, numCells, plotRemainder=F, point.cex=point.cex)

    plotCumulativeReadsPerUMI(molecularBarcodeDistributionByGeneFile)

	plotSensitivity(transcriptDownsamplingFile)
	plotDecilesTable(transcriptQuantileFile)

	if (!is.null(outPlot)) {
        dev.off()
    }
    options(scipen=prevScipen)
}

plotStandardAnalysisPairOrganism<-function(outPlot, digitalExpressionSummaryFile1, digitalExpressionSummaryFile2,
                        molecularBarcodeDistributionByGeneFile1, molecularBarcodeDistributionByGeneFile2,
                        organism1, organism2, cellTypesFile,
                        estimatedNumCells, estimatedNumBeads,
                        transcriptDownsamplingFile1=NULL, transcriptDownsamplingFile2=NULL,
                        transcriptQuantileFile1=NULL, transcriptQuantileFile2=NULL,
                        point.cex=1) {
    prevScipen = options(scipen=10000)[[1]]
    if (!is.null(outPlot)) {
        pdf(outPlot, compress=T)
    }
    df=getNumTranscriptsPerCellBarcodeByOrganismPair(digitalExpressionSummaryFile1, digitalExpressionSummaryFile2, organism1, organism2)

    plotNumTranscriptsPerCellBarcodeByOrganismPair(df, estimatedNumCells, organism1, organism2, plotRemainder=F, point.cex=point.cex)
    plotNumTranscriptsPerCellBarcodeByOrganismPair(df, estimatedNumBeads, organism1, organism2, plotRemainder=F, point.cex=point.cex)
    #plotNumTranscriptsPerCellBarcodeByOrganismPair(df, estimatedNumBeads, organism1, organism2, plotRemainder=T)
    rm (df)

    #I should probably move this function into this file and split out the plotting from the categorization.
    plotCategorizedCellTypes(cellTypesFile, point.cex, organism1, organism2)

    df=getNumGenesPerCellBarcodeByOrganismPair(digitalExpressionSummaryFile1, digitalExpressionSummaryFile2, organism1, organism2)

    plotNumGenesPerCellBarcodeByOrganismPair(df, estimatedNumCells= estimatedNumCells, organism1, organism2, plotRemainder=F, point.cex=point.cex)
    plotNumGenesPerCellBarcodeByOrganismPair(df, estimatedNumCells= estimatedNumBeads, organism1, organism2, plotRemainder=F, point.cex=point.cex)
    #plotNumGenesPerCellBarcodeByOrganismPair(df, estimatedNumCells= estimatedNumBeads, organism1, organism2, plotRemainder=T)

    plotCumulativeReadsPerUMI(molecularBarcodeDistributionByGeneFile1, organism1)
    plotCumulativeReadsPerUMI(molecularBarcodeDistributionByGeneFile2, organism2)

    if (!is.null(transcriptDownsamplingFile1)) { plotSensitivity(transcriptDownsamplingFile1, organism=organism1) }
    if (!is.null(transcriptQuantileFile1)) { plotDecilesTable(transcriptQuantileFile1, organism=organism1) }
    if (!is.null(transcriptDownsamplingFile2)) { plotSensitivity(transcriptDownsamplingFile2, organism=organism2) }
    if (!is.null(transcriptQuantileFile2)) { plotDecilesTable(transcriptQuantileFile2, organism=organism2) }

    if (!is.null(outPlot)) {
        dev.off()
    }
    options(scipen=prevScipen)
}

makeTranscriptDownsamplingPath<-function(outPlot) {
	return(paste(sub(".pdf", "", outPlot), "_transcript_downsampling.txt", sep=""))
}

makeTranscriptQuantilePath<-function(outPlot) {
	return(paste(sub(".pdf", "", outPlot), "_transcript_downsampling_deciles.txt", sep=""))
}

