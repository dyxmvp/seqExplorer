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

##########################
# selectCellsByReadsSCM, categorizeCellsUsingKnee and dependent functions

selectCellsByReadsSCM<-function (bamFile, reportDir, readsPerCellBarcodeFile, estimatedNumBeads, outputFile=NULL) {
	a=read.table(readsPerCellBarcodeFile, header=F, stringsAsFactors=F)
	r= selectCellsByReads(readsPerCellBarcodeFile, estimatedNumBeads, plotFigure=T)
	a=a[1:r,"V2", drop=F]
	if (is.null(outputFile)) {
	    o=paste(reportDir, "/", basename(bamFile), "_selectedCellBarcodes.txt", sep="")
	} else {
	    o = outputFile
	}
	write.table(a, o, row.names=F, col.names=F, quote=F, sep="\t")
	return (o)
}


selectCellsByReads<-function (inFile, estimateNumBeads, all.knots=NA, plotFigure=F, expName="", xlim=NULL) {
	
    a=read.table(inFile, header=F, stringsAsFactors=F)
    v1=a$V1
     
    v1=cumsum(v1)
    v1=v1/v1[length(v1)]
    v1=head (v1, n=estimateNumBeads)
    
    xx=prepData(v1, all.knots=NA, normalizeX=F)
        
    cat (paste("Calling SelectCells on ", estimateNumBeads, "barcodes using inFile", inFile, "\n"))
    
    angleResult=selectByAngle(v1, plotData=F, all.knots=NA)
    #curvatureResult=selectByMaximumCurvature(v1, plotData=F, all.knots=all.knots)
    curvatureResult=NULL
    if (plotFigure) {
    	#par(mar=c(5, 4, 4, 6) + 0.1)
    	plotFullFigure(xx, angleResult, curvatureResult, isDataFitted=T, expName=expName)  	
    }
    cat (paste("Minimum Angle", angleResult$stat[angleResult$idx], "\n"))
    return (angleResult$idx)
}

#' all.knots means to use a knot at every point when smoothing.  Otherwise, use 1 knot per 10 points, or 1000 knots, whatever is smaller.
#' normalizeX means to change output X to be from 0-1 instead of 1:length(d)
prepData<-function (d, all.knots=T, normalizeX=T) {
	if (normalizeX) {
		xx=data.frame(x=seq(from=min(d), to=max(d), length.out=length(d)), y=d)	
	} else {
		xx=data.frame(x=1:length(d), y=d)
	}
	if (is.na(all.knots)) {
		xx$fittedX=xx$x
		xx$fittedY=xx$y
		return (xx)
	}
	
	if (all.knots) {
		z=smooth.spline(xx, all.knots=T)	    	
    } else {
    	nknots=min(length(d)/10, 1000)
    	cat (paste("Number of knots selected", nknots, "\n"))
    	z=smooth.spline(xx, all.knots=F, nknots= nknots)
    }
    
    xx=xx[1:length(z$x),]
    xx$fittedX=z$x
    xx$fittedY=z$y
	return (xx)
}

selectByAngle<-function (d, all.knots=T, plotData =F) {
    
    xx=prepData(d, all.knots)
    
    #'@param a is the current row
    #'@param s is the first row of xx
    #'@param e is the last row of xx
    #a=xx[1000,];s=s=xx[1,];e=xx[dim(xx)[1],];colX="x";colY="y"
    getAngleForPointFaster<-function(a, s, e, colX="fittedX", colY="fittedY") {
        l1=c(as.numeric(s[[colX]]-a[[colX]]), as.numeric(s[[colY]]-a[[colY]]))
        l2=c(as.numeric(e[[colX]]-a[[colX]]), as.numeric(e[[colY]]-a[[colY]]))
        #segments(x0=s[[colX]], x1=a[[colX]], y0=s[[colY]], y1=a[[colY]], col="blue", lwd=3)
        #segments(x0=e[[colX]], x1=a[[colX]], y0=e[[colY]], y1=a[[colY]], col="blue", lwd=3)
        dotProd <- l1%*%l2 
        normL1 <- norm(l1,type="2")
        normL2 <- norm(l2,type="2")
        theta <- acos(dotProd / (normL1 * normL2))
        as.numeric(theta)*(180/pi)
    }
    anglesRawPointsFaster=apply(xx, 1, getAngleForPointFaster, xx[1,], xx[dim(xx)[1],], colX="x", colY="y")
    
    idx=which.min(anglesRawPointsFaster)
    
    if (plotData) {
        plot (xx$x, xx$y, xlab="cells", ylab="cumulative fraction of reads")
        points(xx$fittedX, xx$fittedY, col='red', type='l')
        abline(v=xx[idx,"x"])
        title(paste("Number of cells", idx, " by angle"))
    }
    result=list(idx=idx, stat= anglesRawPointsFaster)
    return (result)
}

plotFullFigure<-function (xx, angleResult=NULL, curvatureResult=NULL, isDataFitted=F, expName="", xlim=NULL) {
	options(scipen=1000)
    
    ltext=c()
    lColors=c()
    plot (xx$x, xx$y, xlab="cells", ylab="cumulative fraction of reads", type='l', lwd=7, col="grey", xlim=xlim, axes=T)
    if (isDataFitted) {
    	points(xx$fittedX, xx$fittedY, type='l', lwd=1, col='red')
    }
    if (!is.null(angleResult)) {
    	points(xx[angleResult$idx,]$x, xx[angleResult$idx,]$y, col="blue", pch=20, cex=1, xlim=xlim)
    }
    if (!is.null(curvatureResult)) {
    	points(xx[curvatureResult$idx,]$x, xx[curvatureResult$idx,]$y, col="red", pch=20, cex=1, xlim=xlim)
    }
    title(paste("Cell Selection", expName))
    
    
    if (!is.null(curvatureResult)) {
    	par(new=TRUE)
	    plot (curvatureResult$stat, type='l', col="blue", axes=F, xlab="", ylab="", ylim=c(min(curvatureResult$stat)*1.5, (max(curvatureResult$stat)*1.5)), xlim=xlim)
    	points(curvatureResult$idx, curvatureResult$stat[curvatureResult$idx], pch=2, col="red", xlim=xlim)
    	
 		ltext=c(ltext, paste("Curvature", curvatureResult$idx))
 		lColors=c(lColors, "red")
    }	
    
 	if (!is.null(angleResult)) {
 		par(new=TRUE)
 		plot (angleResult$stat, type='l', col="blue", axes=F, xlab="", ylab="", xlim=xlim)
    	points(angleResult$idx, angleResult$stat[angleResult$idx], pch=6, col="blue", xlim=xlim)
    	
    	ltext=c(ltext, paste("Angle", angleResult$idx))
    	lColors=c(lColors, "blue")
 	}
	
    legend("bottomright", legend=ltext, fill=lColors)
    
     
}

################################
# categorizeCellsUsingKnee and dependent functions

categorizeCellsUsingKnee<-function (digitalExpressionFile1, digitalExpressionFile2, organismOne, organismTwo,
            minMixedRatio=0.2, pureRatio=0.2, selectedCellsFile, reportDir, bamFile=NULL, outFile=NULL,
            selectedCellsFile1=NULL, selectedCellsFile2=NULL) {
	a=read.table(selectedCellsFile, header=F)
	numCells=dim(a)[1]
	
	dfFull=getNumTranscriptsPerCellBarcodeByOrganismPair(digitalExpressionFile1, digitalExpressionFile2, organismOne, organismTwo)
	dfFull=dfFull[order(dfFull$total, decreasing=T),]
	dfFull$ratio=dfFull[,2]/dfFull[,4]
	df=head (dfFull, n= numCells)
	
	#dfNoCall=dfFull[-1:-numCells,]
	#dfNoCall$organism="No Call"
	df$organism="Mixed"

	#filter out cells with no reads.  This shouldn't happen, but it seems like people do silly things.
	idxBad=which(df[organismOne]==0 & df[organismTwo]==0)
	if (length(idxBad)>0) df=df[-idxBad,]
	
	idx=which(df$ratio> (1-pureRatio))
	#again, protection for no cells being part of that organism.
	if (length(idx)>0) df[idx,]$organism=organismOne
	
	idx=which(df$ratio< (pureRatio))
	if (length(idx)>0) df[idx,]$organism=organismTwo
	
	#dm=trimData(df[df$organism=="Mixed",]$total, trim=0.1)
	#threshold=mean(dm)
	
	#result=rbind(df, dfNoCall)
	result=df
	
	if (is.null(outFile)) {
	    outFile=paste(reportDir, "/", basename(bamFile), "_categorized_cellTypes.txt", sep="")
	}
	write.table(result, outFile, row.names=F, col.names=T, quote=F, sep="\t")

    if (!is.null(selectedCellsFile1)) {
        write.table(result[result$organism==organismOne,]$tag, file=selectedCellsFile1, quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
    if (!is.null(selectedCellsFile2)) {
        write.table(result[result$organism==organismTwo,]$tag, file=selectedCellsFile2, quote=FALSE, row.names=FALSE, col.names=FALSE)
    }

    return (outFile)
}

# Modified to add number of genes
getNumTranscriptsPerCellBarcodeByOrganismPair<-function (digitalExpressionFileO1, digitalExpressionFileO2, organismOne, organismTwo) {
	if (is.null(organismOne) || is.null(organismTwo)) return(NULL)
	
	o1=getGenesAndTranscriptsPerCellBarcode(digitalExpressionFileO1)
	o2=getGenesAndTranscriptsPerCellBarcode(digitalExpressionFileO2)

	commonBC=union(o1$cellBC, o2$cellBC)
	o1p=o1[match(commonBC, o1$cellBC),]
	o2p=o2[match(commonBC, o2$cellBC),]
	df=data.frame(tag=commonBC, o1Count=o1p$numTranscripts, o2Count=o2p$numTranscripts, stringsAsFactors=F) 
	
	idx1=which(is.na(df$o1Count))
	idx2=which(is.na(df$o2Count))
	if (length(idx1)>0) df[idx1,]$o1Count=0
	if (length(idx2)>0) df[idx2,]$o2Count=0

	df$total=apply(df[,2:3], 1, sum, na.rm=T)
	
	
	#df=df[order(df$total, decreasing=T),]         #V0
	colnames(df)[2]= organismOne
	colnames(df)[3]= organismTwo
	
	# added to include nGene
	df$o1Count_gene <- o1p$numGenes
	df$o2Count_gene <- o2p$numGenes
	
	idx1_gene=which(is.na(df$o1Count_gene))                    # added 
	idx2_gene=which(is.na(df$o2Count_gene))                    # added
	if (length(idx1_gene)>0) df[idx1_gene,]$o1Count_gene=0          # added
	if (length(idx2_gene)>0) df[idx2_gene,]$o2Count_gene=0          # added
	df$total_gene=apply(df[,5:6], 1, sum, na.rm=T)                 # added
	colnames(df)[5]= paste0(organismOne, "_gene")           # added
	colnames(df)[6]= paste0(organismTwo, "_gene")           # added
	
	
	df=df[order(df$total, decreasing=T),]
	return (df)
}

