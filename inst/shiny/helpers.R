
# Button indicators
withBusyIndicatorUI <- function(icon_name, button) {
  id <- button[['attribs']][['id']]
  n_s <- strsplit(id, split = "-")[[1]][1]
  div(
    `data-for-btn` = id,
    button,
    span(
      class = "btn-loading-container",
      shinyjs::hidden(
        img(id = paste0(n_s, "-load_", icon_name), src = "ajax-loader-bar.gif", class = "btn-loading-indicator"),
        img(id = paste0(n_s, "-check_", icon_name), src = "checkmark.jpg", class = "btn-loading-indicator")
      )
    )
  )
} # Button indicators END

appCSS <- "
.btn-loading-container {
margin-left: 10px;
font-size: 1.2em;
}
.btn-done-indicator {
color: green;
}
"

# Identify differential expressed genes across conditions
LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0,
                       adj.y.s = 0, text.size = 3.5, segment.size = 0.1) {
  for (i in genes) {
    x1 <- exp.mat[i, 1]
    y1 <- exp.mat[i, 2]
    plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t,
                            label = i, size = text.size)
    plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 +
                              adj.y.s, yend = y1, size = segment.size)
  }
  return(plot)
}

LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05,
                    adj.r.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t,
                    adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
}

LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05,
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t,
                    adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
}

# Capitalizing every first letter of a word
capwords <- function(s, strict = TRUE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

# Plot annotation confidence in SingleR
annotation.confidence.SR <- function (SingleR, dot.size = 1.5, plot.feature = "MaxScore", title = NULL) 
{

  df = data.frame(row.names = rownames(SingleR$cell.names), singler$meta.data$xy)
  
  if (plot.feature == "MaxScore") {
    df$Feature = apply(SingleR$scores, 1, max)
    tit = "Max Score"
  }
  else {
    df$Feature = -log10(SingleR$pval)
    tit = "-log10(p-value)"
  }

  if (is.null(title)) {
    title = tit
  }
  
  ggplot(df, aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point(aes(color = Feature), size = dot.size) + 
    scale_colour_gradient(low = "gray", high = "blue") + ggtitle(title) + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}


# import Seurat V3 object to monocle 2
importCDS_new <- function (otherCDS, import_all = FALSE) 
{
    data <- GetAssayData(object = otherCDS, slot = "counts")
    if (class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }
    pd <- tryCatch({
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, error = function(e) {
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)
      message("This Seurat object doesn't provide any meta data")
      pd
    })
    if (length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]
    }
    fData <- data.frame(gene_short_name = row.names(data), 
                        row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    lowerDetectionLimit <- 0
    if (all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    }
    else if (any(data < 0)) {
      expressionFamily <- uninormal()
    }
    else {
      expressionFamily <- tobit()
    }
    valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                                  lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily)
    if (import_all) {
      if ("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      }
      else {
        mist_list <- otherCDS
      }
    }
    else {
      mist_list <- list()
    }
    if (length(VariableFeatures(object = pbmc)) != 0) {
      var.genes <- setOrderingFilter(monocle_cds, VariableFeatures(object = otherCDS))
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
    return(monocle_cds)
  
}