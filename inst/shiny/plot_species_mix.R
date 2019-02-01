#plot species mixing results


# Plot with summary files
plot_mix_scatter_sum <- function(){

  organismOne <- samplename_human
  organismTwo <- samplename_mouse
  Experiment_name <- paste0(organismOne, "_", organismTwo)
  reportDir <- dirname(selectedBCPath_mix)

  digitalExpressionFile1 <- sum_humanPath
  digitalExpressionFile2 <- sum_mousePath

  selectedCellsFile <- selectedBCPath_mix

  cellTypesFile <<- categorizeCellsUsingKnee(digitalExpressionFile1, digitalExpressionFile2, organismOne, organismTwo,
                           minMixedRatio=0.2, pureRatio=0.2,
                           selectedCellsFile, reportDir, bamFile=Experiment_name, outFile=NULL,
                           selectedCellsFile1=NULL, selectedCellsFile2=NULL)

  plotCategorizedCellTypes(cellTypesFile, point.cex = 0.5, organismOne, organismTwo)

}

plot_mix_hist_sum <- function(){

  df <- read.table(cellTypesFile, header=T, stringsAsFactors=F, check.names = FALSE)
  df <- df$ratio

  bins <- seq(min(df), max(df), length.out =  50)
  
  h <- hist(df, breaks = bins)
  ylim <- range( h$counts )
  hist(df, breaks = bins, col = "dodgerblue4", border = "white",
       xlab = "Fraction of transcripts uniquely mapped to the human",
       ylab = "Number of cells", main = NULL,
       cex.lab = 1.1, cex.axis = 1.25, cex.main = 1.25, cex.sub = 1.25, ylim = ylim*1.1)

}

plot_nTranscripts_sum <- function(){

  df <- read.table(cellTypesFile, header=T, stringsAsFactors=F, check.names = FALSE)

  dp <- ggplot(df[df$organism != "Mixed",], aes(x=organism, y=total, fill=organism)) +
    geom_violin(adjust=3, lwd=0.6) +
    geom_boxplot(width=0.1, fill="white") +
    labs(title=NULL,x=NULL, y = "Number of transcripts") +
    scale_x_discrete(limits=c(samplename_human, samplename_mouse)) +
    scale_fill_manual(values=c("brown", "dodgerblue3"))

  dp + guides(fill=FALSE) +
    theme(axis.text=element_text(size = 16), axis.title=element_text(size = 18, face="bold"))

}

plot_nGene_sum <- function(){

  df <- read.table(cellTypesFile, header=T, stringsAsFactors=F, check.names = FALSE)

  dp <- ggplot(df[df$organism != "Mixed",], aes(x=organism, y=total_gene, fill=organism)) +
    geom_violin(adjust=3, lwd=0.6)+
    geom_boxplot(width=0.1, fill="white") +
    labs(title=NULL,x=NULL, y = "Number of genes") +
    scale_x_discrete(limits=c(samplename_human, samplename_mouse)) +
    scale_fill_manual(values=c("brown", "dodgerblue3"))

  dp + guides(fill=FALSE) +
    theme(axis.text=element_text(size = 16), axis.title=element_text(size = 18, face="bold"))

}
# Plot with summary files END

# Plot with two dge files
two_dge_to_sum <- function(){

  human_2dge <- read.table(dge_humanPath, header=T, row.names = 1, stringsAsFactors=F)
  mouse_2dge <- read.table(dge_mousePath, header=T, row.names = 1, stringsAsFactors=F)

  # generate summary files

  #sampleName_human <- "Human"
  #sampleName_mouse <- "Mouse"

  human_BC <- colnames(human_2dge)
  mouse_BC <- colnames(mouse_2dge)

  human_nUMI <- colSums(human_2dge)
  mouse_nUMI <- colSums(mouse_2dge)

  human_genes <- colSums(human_2dge != 0)
  mouse_genes <- colSums(mouse_2dge != 0)


  sum_humanPath <<- file.path(dirname(dge_humanPath), paste0(samplename_human, ".dge.summary.txt"))
  sum_mousePath <<- file.path(dirname(dge_mousePath), paste0(samplename_mouse, ".dge.summary.txt"))

  human_sum <- data.frame(CELL_BARCODE=human_BC, NUM_GENIC_READS=0,
                          NUM_TRANSCRIPTS=human_nUMI, NUM_GENES=human_genes, stringsAsFactors=F, row.names = NULL)
  mouse_sum <- data.frame(CELL_BARCODE=mouse_BC, NUM_GENIC_READS=0,
                          NUM_TRANSCRIPTS=mouse_nUMI, NUM_GENES=mouse_genes, stringsAsFactors=F, row.names = NULL)

  write.table(human_sum, sum_humanPath, sep="\t", row.names = FALSE, quote = FALSE)
  write.table(mouse_sum, sum_mousePath, sep="\t", row.names = FALSE, quote = FALSE)

}
# Plot with two dge files END

# Plot with one merged dge files
merged_dge_to_2dge <- function(){

  # load human reference
  hg_ref <- read.table("useful/GSM1629193_hg19_ERCC.refFlat")
  hg <- data.frame(hg_ref$V1)

  # load mouse reference
  mm_ref <- read.table("useful/mm10.refFlat")
  mm <- data.frame(mm_ref$V1)

  # add a empty column for organism
  d1 <- read.table(dge_mergedPath, header=T, stringsAsFactors=F)
  d1$Organism <- NA

  # search and label human
  idx_h <- which(d1$GENE %in% hg$hg_ref.V1)
  len_h <- length(idx_h)
  d1[idx_h,]$Organism <- samplename_human

  # search and label mouse
  idx_m <- which(d1$GENE %in% mm$mm_ref.V1)
  len_m <- length(idx_m)
  d1[idx_m,]$Organism <- samplename_mouse

  # check the label
  #len_sum <- len_m + len_h
  #num_NA <- sum(is.na(d1$Organism))

  # generate separeted dge files by organism
  human_dge <- subset(d1, Organism == samplename_human, select = -c(Organism))
  mouse_dge <- subset(d1, Organism == samplename_mouse, select = -c(Organism))

  dge_humanPath <<- file.path(dirname(dge_mergedPath), paste0(samplename_human, ".dge.txt.gz"))
  dge_mousePath <<- file.path(dirname(dge_mergedPath), paste0(samplename_mouse, ".dge.txt.gz"))

  write.table(human_dge, dge_humanPath, sep="\t", row.names = FALSE)
  write.table(mouse_dge, dge_mousePath, sep="\t", row.names = FALSE)

}
# Plot with one merged dge files END



