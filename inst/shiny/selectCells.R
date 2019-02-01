#plot knee curve and select cells

source("cell_selection.R")

input <- commandArgs(trailingOnly = TRUE)
estimatedNumCells <- as.numeric(input[1])
samplename <- input[2]
scriptPath <- input[3]

dir <- dirname(scriptPath)
bamFile <- file.path(dir, samplename)
reportDir <- dir
readsPerCellBarcodeFile <- file.path(dir, paste(samplename, ".readcounts.txt.gz", sep = ""))
estimatedNumBeads = 20 * estimatedNumCells


tiff(file.path(dir, "selectCells.tif"), type="cairo")

selectCellsByReadsSCM(bamFile, reportDir, readsPerCellBarcodeFile, estimatedNumBeads)

dev.off()
