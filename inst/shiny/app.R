require(shiny) # install.packages("shiny")
require(shinyFiles) # install.packages('shinyFiles')
require(shinydashboard) # install.packages("shinydashboard")
require(shinyjs) # install.packages("shinyjs")
require(shinyBS) # install.packages("shinyBS")
require(shinycssloaders)  # devtools::install_github("andrewsali/shinycssloaders")
require(DT) # install.packages('DT')
require(Seurat) # install.packages('Seurat');install.packages("lmtest")
require(dplyr)
require(Matrix)
require(plotly) # install.packages("plotly")
require(cowplot) # install.packages("cowplot")
require(rmarkdown) # install.packages("rmarkdown")
require(SingleR) # install.packages("BiocManager"); devtools::install_github('dviraran/SingleR')
require(monocle) # if (!requireNamespace("BiocManager", quietly = TRUE));install.packages("BiocManager");BiocManager::install("monocle")

# For UMAP
#install.packages("reticulate")
#library(reticulate)
#conda_create("r-reticulate")
#reticulate::py_install(packages = 'umap-learn')
require(reticulate)
use_condaenv("r-reticulate")

#reticulate::py_install(packages = 'mlxtend')

#########################   UI Start   ###############################################

source("intro.R")
source("fastqc.R")
source("dge.R")
source("species_mix.R")
source("quick_look.R")
source("seurat_normal.R")
source("seurat_integrated.R")
source("monocle.R")

ui <- tagList(
  dashboardPage(
    dashboardHeader(title = "seq-Explorer"),
    dashboardSidebar(
      sidebarMenu(
        id = "tabs",
        menuItem("User Guide", tabName = "intro_tab", icon = icon("info-circle")),
        menuItem("FastQC", tabName = "fastqc_tab", icon = icon("info-circle")),
        menuItem("Generate DGE", tabName = "dge_tab", icon = icon("th")),
        menuItem("Species mix", tabName = "mix_tab", icon = icon("ruler-combined")),
        menuItem("Quick look", tabName = "look_tab", icon = icon("search")),
        menuItem("Seurat", tabName = "seurat_tab", icon = icon("bar-chart"),
                 menuSubItem('Normal workflow', tabName = 'seurat_normal_tab', icon = icon('th')),
                 menuSubItem('Integrated analyses', tabName = 'seurat_integrated_tab', icon = icon('th'))),
        menuItem("Monocle", tabName = "monocle_tab", icon = icon("pagelines"))
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
        seurat_integrated_UI("seuratIntegrated"),
        monocle_UI("MONOCLE")
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
  
  df_monocle <- callModule(monocle_Server, "MONOCLE")

}



############################################   Server END   ###############################################


shinyApp(ui, server) # Run App
