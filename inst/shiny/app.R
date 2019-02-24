require(shiny)
require(shinyFiles)
require(shinydashboard)
require(shinyjs)
require(shinyBS)
require(shinycssloaders)  # devtools::install_github("andrewsali/shinycssloaders")
require(DT)
require(Seurat)
require(dplyr)
require(Matrix)
require(plotly)

#########################   UI Start   ###############################################

source("intro.R")
source("fastqc.R")
source("dge.R")
source("species_mix.R")
source("quick_look.R")
source("seurat_normal.R")
source("seurat_integrated.R")

ui <- tagList(
  dashboardPage(
    dashboardHeader(title = "seq-Explorer"),
    dashboardSidebar(
      sidebarMenu(
        id = "tabs",
        menuItem("User Guide", tabName = "intro_tab", icon = icon("info-circle")),
        menuItem("FastQC", tabName = "fastqc_tab", icon = icon("info-circle")),
        menuItem("Generate DGE", tabName = "dge_tab", icon = icon("th")),
        menuItem("Species mix", tabName = "mix_tab", icon = icon("bar-chart")),
        menuItem("Quick look", tabName = "look_tab", icon = icon("search")),
        menuItem("Seurat", tabName = "seurat_tab", icon = icon("bar-chart"),
                 menuSubItem('Normal workflow', tabName = 'seurat_normal_tab', icon = icon('th')),
                 menuSubItem('Integrated analyses', tabName = 'seurat_integrated_tab', icon = icon('th')))
      )
    ),
    dashboardBody(

      shinyjs::useShinyjs(),
      tabItems(
        #source("ui-tab-intro.R", local = TRUE)$value,
        introUI("Intro"),
        fastqcUI("QC"),
        dgeUI("DGE"),
        mixUI("MIX"),
        lookUI("LOOK"),
        seurat_normal_UI("seuratNormal"),
        seurat_integrated_UI("seuratIntegrated")
      )

    )

  ),
  tags$footer(
    wellPanel(
      HTML(
        '
      <p align="center" width="4">Developed and maintained by: Yanxiang Deng@Rong Fan Lab, BME, Yale University</p>
      <p align="center" width="4">Email: yanxiang.deng@yale.edu</p>      
      <p align="center" width="4">The software is licensed under the GPL-3.0</p>
      '
      )
    )
  )
)



###############################################   UI END   ##################################################



###############################################   Server Start   ############################################

#max upload 300mb
options(shiny.maxRequestSize = 300*1024^2)


server <- function(input, output, session) {

  df_fastqc <- callModule(fastqcServer, "QC")

  df_dge <- callModule(dgeServer, "DGE")

  df_mix <- callModule(mixServer, "MIX")
  
  df_look <- callModule(lookServer, "LOOK")

  df_seurat_normal <- callModule(seurat_normal_Server, "seuratNormal")

  df_seurat_integrated <- callModule(seurat_integrated_Server, "seuratIntegrated")

}



############################################   Server END   ###############################################


shinyApp(ui, server) # Run App
