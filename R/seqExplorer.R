#' @export
#' seqExplorer
seqExplorer <- function() {
  appDir <- system.file("shiny", package = "seqExplorer")

  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}
