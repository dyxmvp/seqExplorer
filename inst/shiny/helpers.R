
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
