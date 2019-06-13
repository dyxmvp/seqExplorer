introUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "intro_tab",
          fluidRow(
            
            box(title = "User Guide", width = 12, solidHeader = T, status = "primary",
                       includeMarkdown("intro.Rmd")
            )
          )
  )
}
