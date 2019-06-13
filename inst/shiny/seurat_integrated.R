

source("helpers.R")
# Data input and QC interface
seurat_integrated_UI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "seurat_integrated_tab",
          navbarPage(
            title = NULL,
            tabPanel("1. Data preprocessing",
                     fluidRow(

                       box(title = "1. Upload Data", width = 3, solidHeader = TRUE, status = "info",
                           collapsible = T,
                           tags$style(appCSS),
                           radioButtons(ns("data_type_integrated"), NULL,
                                        c(
                                          "DGE Data (*.dge.*)" = "dge_file_integrated",
                                          "10X Data (*.mtx and *.tsv)" = "10X_file_integrated",
                                          "Example Data (PBMCs)" = "demo_file_integrated"
                                        ), selected = "dge_file_integrated"),

                           conditionalPanel(condition  = "input.data_type_integrated == 'dge_file_integrated'",
                                            tags$b("Select DGE file:"),
                                            tags$p(),
                                            verbatimTextOutput(ns("input_dge_con"), placeholder = TRUE),
                                            shinyFilesButton(ns("con_dge"), "File from Control Sample",
                                                             "Please select a file", multiple = FALSE),

                                            tags$p(),
                                            verbatimTextOutput(ns("input_dge_stim"), placeholder = TRUE),
                                            shinyFilesButton(ns("stim_dge"), "File from Stimulated Sample",
                                                             "Please select a file", multiple = FALSE),
                                            ns = NS(id)
                           ),

                           conditionalPanel(condition  = "input.data_type_integrated == '10X_file_integrated'",
                                            tags$b("Select 10x files folder:"),
                                            tags$p(),
                                            verbatimTextOutput(ns("input_10x_con"), placeholder = TRUE),
                                            shinyDirButton(ns("con_10x"), "Folder from Control Sample", "Please select a folder"),

                                            tags$p(),
                                            verbatimTextOutput(ns("input_10x_stim"), placeholder = TRUE),
                                            shinyDirButton(ns("stim_10x"), "Folder from Stimulated Sample", "Please select a folder"),
                                            ns = NS(id)
                           )
                       ),

                       box(title = "2. Initialize Seurat object", width = 4, height = 400, solidHeader = TRUE, status = "warning",
                           collapsible = T,

                           splitLayout(
                             textInput(ns("project_name_ctrl"), value = "IMMUNE_CTRL", label = "Project Name (control)"),
                             textInput(ns("project_name_stim"), value = "IMMUNE_STIM", label = "Project Name (stimulation)")
                           ),

                           numericInput(ns("min_cells_integrated"),
                                        label="Keep genes expressed at least in (?) cells", min = 1, max = Inf, value = 5),
                           numericInput(ns("min_genes_integrated"),
                                        label="Keep cells with at least (?) genes",min = 1, max = Inf, value = 1),
                           withBusyIndicatorUI(icon_name = "sinit_integrated",
                                               actionButton(ns("seurat_init_integrated"),
                                                            "Process raw data", style = "width: 80%")
                           )

                       ),

                       box(title = "3. Cell QC", width = 2, height = 400, solidHeader = TRUE, status = "danger",
                           collapsible = T,

                           selectizeInput(ns("subsetNames_integrated"), label="Filter Names",
                                          choices = c("nGene", "percent.mito"),
                                          selected = "nGene",
                                          multiple=TRUE),
                           withBusyIndicatorUI(icon_name = "cfilter_integrated",
                                               actionButton(ns("cell_filter_integrated"),"Select cells", style = "width: 80%")
                           )
                       ),

                       box(
                         title = "4. Normalize data", width = 3, height = 400, solidHeader = TRUE, status = "success",
                         selectizeInput(ns("norm_methods_integrated"), label="Normalization method",
                                        choices = c("LogNormalize"),
                                        selected = c("LogNormalize"),
                                        multiple=FALSE),
                         numericInput(ns("scale_integrated"),
                                      label="Scale factor", min = 1, max = Inf, value = 10000),

                         selectizeInput(ns("vars_regress_integrated"), label="Vars.to.regress",
                                        choices = c("nUMI", "percent.mito"),
                                        multiple = TRUE),

                         withBusyIndicatorUI(icon_name = "norm_integrated",
                                             actionButton(ns("data_norm_integrated"),"Normalize and scale", style = "width: 80%")
                         )
                       )


                     ),

                     fluidRow(
                       box(
                         title = "3.1 VlnPlot (Cell QC)",  width = 8, height = 600, status = "danger",

                         sidebarPanel(
                           tableOutput(ns("nGene_range_integrated")),
                           uiOutput(ns("gene_low_high_integrated")),
                           hr(),
                           tableOutput(ns("mito_range_integrated")),
                           uiOutput(ns("mito_low_high_integrated"))
                         ),
                         mainPanel(
                           uiOutput(ns("filter_plots_integrated"))
                         )

                       ),
                       box(
                         title = "3.2 GenePlot (Cell QC)",  width = 4, height = 600, status = "danger",
                         uiOutput(ns("gene_plots_integrated"))
                       )
                     )
            ), # Tab 1 END

            tabPanel('2. Perform integration',
                     fluidRow(
                       column(width = 5,
                              box(title = "Perform integration", width = NULL, solidHeader = TRUE, status = "warning",
                                  collapsible = T,

                                  splitLayout(
                                    numericInput(ns("align_dims"),
                                                 label="Dimensions to align: ", min = 1, max = Inf, value = 20),
                                    numericInput(ns("ccs_res"),
                                                 label="Resolution (TSNE): ", min = 0, max = Inf, value = 0.5, step = 0.01)
                                  ),

                                  withBusyIndicatorUI(icon_name = "cca_align",
                                                      actionButton(ns("viz_cca_align"),
                                                                   "Perform integration",
                                                                   style = "width: 90%")
                                  ),
                                  br(),
                                  uiOutput(ns("old_name_cca_align")),
                                  textInput(ns("new_name_cca_align"), "Assign new names (use ',' to separate): "),

                                  withBusyIndicatorUI(icon_name = "ccs_names",
                                                      actionButton(ns("assign_ccs_names"),
                                                                   "Assign new names",
                                                                   style = "width: 90%")
                                  )
                                  )

                       ), #

                       column(width = 7,
                              box(
                                title = "Plots",  width = NULL, height = 850, status = "primary",
                                tabsetPanel(type = "tabs",
                                            tabPanel("UMAP",
                                                     br(),
                                                     uiOutput(ns("ccs_tsne_plot"))),
                                            tabPanel("UMAP_side_by_side",
                                                     br(),
                                                     uiOutput(ns("ccs_side_plot")))
                                )
                              )
                       )#
                     )
            ), # Tab 2 END

            tabPanel('3. Identify conserved cell type markers',

                     fluidRow(

                             box(title = "Identify conserved markers", width = 3, height = 500, solidHeader = TRUE, status = "warning",
                                 collapsible = T,

                                 radioButtons(ns("con_genes_radio"), label = NULL,
                                              choices = list("Single cluster" = "single_con",
                                                             "Multiple clusters" = "multiple_con",
                                                             "All" = "all_con"),
                                              selected = "single_con", inline = TRUE),

                                 uiOutput(ns("id1_con_genes")),
                                 uiOutput(ns("id2_con_genes")),

                                 splitLayout(
                                   numericInput(ns("min_pct_con"),
                                                label="Minimum percentage:", min = 0, max = 1, value = 0.25, step = 0.01),
                                   numericInput(ns("top_markers_con"),
                                                label="Top markers to show:", min = 1, max = Inf, value = 10)
                                 ),

                                 withBusyIndicatorUI(icon_name = "con_genes",
                                                     actionButton(ns("viz_con_genes"),
                                                                  "Cluster biomarkers",
                                                                  style = "width: 80%")
                                 ),
                                 
                                 br(),
                                 uiOutput(ns("id_con_genes")),

                                 withBusyIndicatorUI(icon_name = "con_genes_con",
                                                     actionButton(ns("viz_con_genes_con"),
                                                                  "Conserved markers",
                                                                  style = "width: 80%")
                                 )

                             ),


                              box(
                                title = "Conserved markers table",  width = 9, height = 500, status = "danger",
                                uiOutput(ns("con_genes_table"))
                              ),

                             box(
                               title = "Cluster biomarkers table",  width = 5, height = 900, solidHeader = TRUE, status = "success",
                               uiOutput(ns("ccs_genes_table"))
                             ),

                              box(
                                title = "Biomarkers visulization",  width = 7, height = 900, status = "primary",
                                tabsetPanel(type = "tabs",
                                            tabPanel("Heatmap",
                                                     br(),
                                                     uiOutput(ns("con_genes_heatmap"))),
                                            tabPanel("Probability distributions",
                                                     br(),
                                                     uiOutput(ns("con_genes_vlnplot"))),
                                            tabPanel("Features",
                                                     br(),
                                                     uiOutput(ns("con_genes_featureplot"))),
                                            tabPanel("RidgePlot",
                                                     br(),
                                                     uiOutput(ns("con_genes_ridgeplot"))),
                                            tabPanel("Features_Conserved markers",
                                                     br(),
                                                     uiOutput(ns("con_genes_featureplot_con"))),
                                            tabPanel("SplitDotPlot",
                                                     br(),
                                                     uiOutput(ns("con_genes_splitdotplot")))
                               )
                             )

                     )#
            ), # Tab 3 END

            tabPanel('4. Differential genes across conditions',
                     fluidRow(
                       column(width = 5,
                              box(title = "Differential genes", width = NULL, solidHeader = TRUE, status = "warning",
                                  collapsible = T,

                                  uiOutput(ns("id_con_genes_diff")),

                                  withBusyIndicatorUI(icon_name = "ccs_diff",
                                                      actionButton(ns("viz_ccs_diff"),
                                                                   "Find differential genes",
                                                                   style = "width: 90%")
                                  ),
                                  br(),
                                  uiOutput(ns("acs_genes_table"))
                                  ),

                              box(title = "Average expression", width = NULL, height = 580, status = "success",
                                  uiOutput(ns("acs_plot"))
                              )
                       ), #

                       column(width = 7,
                              box(
                                title = "Differential genes heatmap",  width = NULL, status = "primary",
                                tabsetPanel(type = "tabs",
                                            tabPanel("FeaturePlots",
                                                     br(),
                                                     uiOutput(ns("ccs_heatmap_diff"))),
                                            tabPanel("VlnPlot",
                                                     br(),
                                                     uiOutput(ns("ccs_vlnplot_diff")))
                                )
                              )
                       )#
                     )
            ), # Tab 4 END
            
            tabPanel('* Save and load Seurat object',
                     fluidRow(
                       
                       box(
                         title = "Save Seurat object (after RunCCA)", width = 6, solidHeader = TRUE, status = "primary",
                         textInput(ns("saveName_integrated"), "Sample name:", width = "100%"),
                         tags$b("Save to directory:"),
                         tags$p(),
                         shinyDirButton(ns("dir_save_integrated"), "Select a folder to save", "Please select a folder"),
                         tags$p(),
                         verbatimTextOutput(ns("save_dir_integrated"), placeholder = TRUE),
                         tags$p(),
                         
                         withBusyIndicatorUI(icon_name = "save_integrated",
                                             actionButton(ns("do_save_integrated"),
                                                          "Save (after RunCCA)",
                                                          style = "width: 90%")
                         )
                       ),
                       
                       box(
                         title = "Load Seurat object", width = 6, solidHeader = TRUE, status = "warning",
                         
                         tags$b("Load from directory:"),
                         tags$p(),
                         shinyFilesButton(ns("integrated_obj"), "Select Seurat object file", "Please select a file", multiple = FALSE),
                         tags$p(),
                         verbatimTextOutput(ns("input_obj_integrated"), placeholder = TRUE),
                         tags$p(),
                         
                         withBusyIndicatorUI(icon_name = "read_integrated",
                                             actionButton(ns("do_load_integrated"),
                                                          "Load",
                                                          style = "width: 90%")
                         )
                       )
                       
                     )#
            ), # Tab 5 END
            
            tabPanel('* Save figures',
                     fluidRow(
                       
                       box(
                         title = "Save figure", width = 3, solidHeader = TRUE, status = "primary",
                         
                         selectizeInput(ns("psave_select_integrated"), label="Choose a figure to save",
                                        choices = c("QC", "GenePlot", "CCA_plot", "CCA_side", "VlnPlot_markers", 
                                                    "FeaturePlot", "RidgePlot", "Features_conserved", "SplitDotPlotGG", 
                                                    "Scatter_plot", "FeatureHeatmap", "VlnPlot_side"),
                                        selected = c("QC"),
                                        multiple=FALSE
                         ),
                         
                         numericInput(ns("pWidth_integrated"), "Width (px)", value = 1200),   # Width is divided by 2 for display
                         numericInput(ns("pHeight_integrated"), "Height (px)", value = 1200), # Height is divided by 2 for display
                         numericInput(ns("pRes_integrated"), "Resolution (ppi)", value = 300),
                         numericInput(ns("pPsize_integrated"), "Pointsize", value = 12),
                         
                         downloadButton(ns("do_psave_integrated"), "Save figure", style = "width: 100%")
                       ),
                       
                       box(
                         title = "Figure to save", width = 9, status = "success",
                         
                         uiOutput(ns("figure_save_integrated"))
                       )
                       
                     )#
            )# Tab 6 END
          )
  )

}

seurat_integrated_Server <- function(input, output, session, fileRoot = NULL) {

  ns <- session$ns

  #volumes <- (c(Home = fs::path_home()))
  volumes <- (c(Home = fs::path_home(), getVolumes()()))

  shinyFileChoose(input, "con_dge", roots = volumes, session = session)
  shinyFileChoose(input, "stim_dge", roots = volumes, session = session)

  shinyDirChoose(input, "con_10x", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyDirChoose(input, "stim_10x", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  shinyDirChoose(input, "dir_save_integrated", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyFileChoose(input, "integrated_obj", roots = volumes, session = session, filetypes=c('', 'rds'))

  output$input_dge_con <- renderPrint({
    as.character(parseFilePaths(volumes, input$con_dge)[4])
  })

  output$input_dge_stim <- renderPrint({
    as.character(parseFilePaths(volumes, input$stim_dge)[4])
  })

  output$input_10x_con <- renderPrint({
    parseDirPath(volumes, input$con_10x)
  })

  output$input_10x_stim <- renderPrint({
    parseDirPath(volumes, input$stim_10x)
  })
  
  output$save_dir_integrated <- renderPrint({
    parseDirPath(volumes, input$dir_save_integrated)
  })
  
  output$input_obj_integrated <- renderPrint({
    as.character(parseFilePaths(volumes, input$integrated_obj)[4])
  })

  # Object initialization
  observeEvent(input$seurat_init_integrated, {

    withProgress(message = "Initializing Seurat Object ...",{

      shinyjs::show("load_sinit_integrated")

      if(input$data_type_integrated == "dge_file_integrated"){
        File_dge_con <- isolate(as.character(parseFilePaths(volumes, input$con_dge)[4]))
        File_dge_stim <- isolate(as.character(parseFilePaths(volumes, input$stim_dge)[4]))

        ctrl.data <- read.table(File_dge_con, header = T, row.names = 1, stringsAsFactors = F)
        stim.data <- read.table(File_dge_stim, header = T, row.names = 1, stringsAsFactors = F)
      }
      else if(input$data_type_integrated == "10X_file_integrated"){
        Path_10x_con <- isolate(as.character(parseDirPath(volumes, input$con_10x)))
        Path_10x_stim <- isolate(as.character(parseDirPath(volumes, input$stim_10x)))

        ctrl.data <- Read10X(data.dir = Path_10x_con)
        stim.data <- Read10X(data.dir = Path_10x_stim)
      }
      else{
        ctrl.data <- read.table("useful//integrated_analysis/immune_control_expression_matrix.txt.gz",
                                sep = "\t")
        stim.data <- read.table("useful//integrated_analysis/immune_stimulated_expression_matrix.txt.gz",
                                sep = "\t")
      }

      setProgress(value = 0.5)

      # Set up control object
      ctrl <<- CreateSeuratObject(counts = ctrl.data, min.cells = input$min_cells_integrated,
                                  min.features = input$min_genes_integrated,
                                  project = input$project_name_ctrl)
      ctrl$stim <- "CTRL"
      ctrl <<- ctrl

      ctrl[["percent.mito"]] <- PercentageFeatureSet(object = ctrl, pattern = "(?i)^MT-") #(?i) for case insensitive
      ctrl <<- ctrl
      

      # Set up stimulated object
      stim <<- CreateSeuratObject(counts = stim.data, min.cells = input$min_cells_integrated,
                                  min.features = input$min_genes_integrated,
                                  project = input$project_name_stim)
      stim$stim <- "STIM"
      stim <<- stim
      
      stim[["percent.mito"]] <- PercentageFeatureSet(object = stim, pattern = "(?i)^MT-") #(?i) for case insensitive
      stim <<- stim

      setProgress(value = 0.7)

      # Set up combined object
      pbmc.combined <<- merge(x = ctrl, y = stim, merge.data = FALSE,
                             add.cell.ids = c("ctrl", "stim"), project = "Integrated_analysis")
      
      pbmc.combined[["percent.mito"]] <- PercentageFeatureSet(object = pbmc.combined, pattern = "(?i)^MT-") #(?i) for case insensitive
      pbmc.combined <<- pbmc.combined

      setProgress(value = 1)

      shinyjs::hide("load_sinit_integrated")
      shinyjs::show("check_sinit_integrated")
    })

  }

  )# Object initialization END

  #Cell QC
  observeEvent(input$cell_filter_integrated, {

    shinyjs::show("load_cfilter_integrated")

    output$filter_plots_integrated <- renderUI({
      splitLayout(cellWidths = c("33%", "33%", "34%"),
                  withSpinner(plotOutput(ns("nGenePlot_integrated")), type = getOption("spinner.type", default = 4)),
                  withSpinner(plotOutput(ns("nUMIPlot_integrated")), type = getOption("spinner.type", default = 4)),
                  withSpinner(plotOutput(ns("mitoPlot_integrated")), type = getOption("spinner.type", default = 4))
      )
    }
    )

    output$gene_plots_integrated <- renderUI({
      splitLayout(cellWidths = c("49%", "49%"),
                  withSpinner(plotOutput(ns("gene_nUMI_integrated")), type = getOption("spinner.type", default = 4)),
                  withSpinner(plotOutput(ns("mito_nUMI_integrated")), type = getOption("spinner.type", default = 4))
      )
    }
    )

    min_gene_ctrl <- min(ctrl@meta.data$nFeature_RNA)
    max_gene_ctrl <- max(ctrl@meta.data$nFeature_RNA)

    min_gene_stim <- min(stim@meta.data$nFeature_RNA)
    max_gene_stim <- max(stim@meta.data$nFeature_RNA)

    min_mito_ctrl <- round(min(ctrl@meta.data$percent.mito), 2)
    max_mito_ctrl <- round(max(ctrl@meta.data$percent.mito), 2)

    min_mito_stim <- round(min(stim@meta.data$percent.mito), 2)
    max_mito_stim <- round(max(stim@meta.data$percent.mito), 2)

    output$nGenePlot_integrated <- renderPlot({

      if("nGene" %in% input$subsetNames_integrated){

        VlnPlot(object = pbmc.combined, features = c("nFeature_RNA"), nCol = 1) +
          geom_hline(yintercept = input$gene_low_integrated, color = 'darkgreen', linetype = "dashed", size = 1.5) +
          #geom_text(data = data.frame(x=1.2,y=input$gene_low_integrated), aes(x, y), label="low", vjust=1.5, hjust=0,color = "darkgreen",size = 5,fontface = "bold") +
          geom_hline(yintercept = input$gene_high_integrated, color = 'blue', linetype = "dashed", size = 1.5) +
          #geom_text(data = data.frame(x=1.2,y=input$gene_high_integrated), aes(x, y), label="high", vjust=-1, hjust=0,color = "blue",size = 5,fontface = "bold") +
          scale_y_continuous(limits=c(0, 1.1*max_gene_ctrl)) + NoLegend()
      }
      else{
        VlnPlot(object = pbmc.combined, features = c("nFeature_RNA"), nCol = 1) + NoLegend()
      }
    })
    output$nUMIPlot_integrated <- renderPlot({
      VlnPlot(object = pbmc.combined, features = c("nCount_RNA"), nCol = 1) + NoLegend()

    })
    output$mitoPlot_integrated <- renderPlot({

      if("percent.mito" %in% input$subsetNames_integrated){

        VlnPlot(object = pbmc.combined, features = c("percent.mito"), nCol = 1) +
          geom_hline(yintercept = input$mito_low_integrated, color = 'darkgreen', linetype = "dashed", size = 1.5) +
          #geom_text(data = data.frame(x=1.2,y=input$mito_low_integrated), aes(x, y), label="low", vjust=1.5, hjust=0,color = "darkgreen",size = 5,fontface = "bold") +
          geom_hline(yintercept = input$mito_high_integrated, color = 'blue', linetype = "dashed", size = 1.5) +
          #geom_text(data = data.frame(x=1.2,y=input$mito_high_integrated), aes(x, y), label="high", vjust=-1, hjust=0,color = "blue",size = 5,fontface = "bold") +
          scale_y_continuous(limits=c(min_mito_ctrl-0.01, max_mito_ctrl+0.01)) + NoLegend()
      }
      else{
        VlnPlot(object = pbmc.combined, features = c("percent.mito"), nCol = 1) + NoLegend()
      }
    })
    
    p1_integrated <<- VlnPlot(object = pbmc.combined, 
                              features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), 
                              ncol = 3) + NoLegend()
    
    median_gene_ctrl <- median(ctrl@meta.data$nFeature_RNA)
    median_gene_stim <- median(stim@meta.data$nFeature_RNA)

    min_UMI_ctrl <- min(ctrl@meta.data$nCount_RNA)
    max_UMI_ctrl <- max(ctrl@meta.data$nCount_RNA)
    median_UMI_ctrl <- median(ctrl@meta.data$nCount_RNA)

    min_UMI_stim <- min(stim@meta.data$nCount_RNA)
    max_UMI_stim <- max(stim@meta.data$nCount_RNA)
    median_UMI_stim <- median(stim@meta.data$nCount_RNA)

    nGene_summary <- data.frame(c(min_gene_ctrl, median_gene_ctrl, max_gene_ctrl),
                                c(min_UMI_ctrl, median_UMI_ctrl, max_UMI_ctrl),
                                c(min_gene_stim, median_gene_stim, max_gene_stim),
                                c(min_UMI_stim, median_UMI_stim, max_UMI_stim),
                                row.names = c("Min", "Median", "Max"))
    nGene_summary <- t(nGene_summary)
    row.names(nGene_summary) <- c("nGene_ctrl", "nUMI_ctrl", "nGene_stim", "nUMI_stim")

    output$nGene_range_integrated <- renderTable(spacing = "xs", align = "c", rownames = TRUE,
                                      striped = TRUE, hover = TRUE, bordered = TRUE, digits = 0,
                                      { nGene_summary }
    )

    output$gene_low_high_integrated <- renderUI({
      splitLayout(
        numericInput(ns("gene_low_integrated"), "nGene_lowest", value = min_gene_ctrl),
        numericInput(ns("gene_high_integrated"), "nGene_highest", value = max_gene_ctrl)
      )
    }
    )

    median_mito_ctrl <- median(ctrl@meta.data$percent.mito)
    median_mito_stim <- median(stim@meta.data$percent.mito)
    mito_summary <- data.frame(c(min_mito_ctrl, median_mito_ctrl, max_mito_ctrl),
                               c(min_mito_stim, median_mito_stim, max_mito_stim),
                               row.names = c("Min", "Median", "Max"))
    mito_summary <- t(mito_summary)
    row.names(mito_summary) <- c("p.mito_ctrl", "p.mito_stim")

    output$mito_range_integrated <- renderTable(spacing = "xs", align = "c", rownames = TRUE,
                                     striped = TRUE, hover = TRUE, bordered = TRUE,
                                     { mito_summary }
    )

    output$mito_low_high_integrated <- renderUI({
      splitLayout(
        numericInput(ns("mito_low_integrated"), "mito_lowest", value = min_mito_ctrl - 0.01, step = 0.01),
        numericInput(ns("mito_high_integrated"), "mito_highest", value = max_mito_ctrl, step = 0.01)
      )
    }
    )

    output$mito_nUMI_integrated <- renderPlot({

      if("percent.mito" %in% input$subsetNames_integrated){

        FeatureScatter(object = pbmc.combined, feature1 = "nCount_RNA", feature2 = "percent.mito") + 
        geom_hline(yintercept = input$mito_low_integrated, color = 'darkgreen', linetype = "dashed", size = 1.5) +
        geom_hline(yintercept = input$mito_high_integrated, color = 'blue', linetype = "dashed", size = 1.5) +
        NoLegend()
      }
      else{
        FeatureScatter(object = pbmc.combined, feature1 = "nCount_RNA", feature2 = "percent.mito") + NoLegend()
      }
    })

    output$gene_nUMI_integrated <- renderPlot({

      if("nGene" %in% input$subsetNames_integrated){

        FeatureScatter(object = pbmc.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
          geom_hline(yintercept = input$gene_low_integrated, color = 'darkgreen', linetype = "dashed", size = 1.5) +
          geom_hline(yintercept = input$gene_high_integrated, color = 'blue', linetype = "dashed", size = 1.5) +
          NoLegend()
      }
      else{
        FeatureScatter(object = pbmc.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
      }
    })

    shinyjs::hide("load_cfilter_integrated")
    shinyjs::show("check_cfilter_integrated")

  }

  )# Cell QC END


  # Normalization
  observeEvent(input$data_norm_integrated, {

    withProgress(message = "Normalizing data...",{

      shinyjs::show("load_norm_integrated")
      setProgress(value = 0.35)
      
      gene_low_sel <<- isolate(as.numeric(input$gene_low_integrated))
      gene_high_sel <<- isolate(as.numeric(input$gene_high_integrated))
      mito_high_sel <<- isolate(as.numeric(input$mito_high_integrated))

      if("nGene" %in% input$subsetNames_integrated && "percent.mito" %in% input$subsetNames_integrated){

        #print("c1")
        
        ctrl <<- subset(x = ctrl, subset = nFeature_RNA > gene_low_sel & nFeature_RNA < gene_high_sel & 
                          percent.mito < mito_high_sel)
        
        stim <<- subset(x = stim, subset = nFeature_RNA > gene_low_sel & nFeature_RNA < gene_high_sel & 
                          percent.mito < mito_high_sel)

      }
      else if("nGene" %in% input$subsetNames_integrated){

        #print("c2")
        
        ctrl <<- subset(x = ctrl, subset = nFeature_RNA > gene_low_sel & nFeature_RNA < gene_high_sel)
        
        stim <<- subset(x = stim, subset = nFeature_RNA > gene_low_sel & nFeature_RNA < gene_high_sel)

      }
      else{

        #print("c3")
        
        ctrl <<- subset(x = ctrl, subset = percent.mito < mito_high_sel)
        
        stim <<- subset(x = stim, subset = percent.mito < mito_high_sel)

      }

      ctrl <<- NormalizeData(ctrl, verbose = FALSE)
      stim <<- NormalizeData(stim, verbose = FALSE)
      
      ctrl <<- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
      stim <<- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)
      
      setProgress(value = 1)
      shinyjs::hide("load_norm_integrated")
      shinyjs::show("check_norm_integrated")
    })

  }

  )# Normalization END
  

  # Perform integration
  observeEvent(input$viz_cca_align, {

    output$align_cca_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_align_cca"), height = "400px"), type = getOption("spinner.type", default = 4))
    }
    )

    output$ccs_tsne_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_cca_tsne"), height = "650px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    output$ccs_side_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_cca_side"), height = "650px"), type = getOption("spinner.type", default = 4))
    }
    )

    withProgress(message = "Plotting UMAP ...",{

      shinyjs::show("load_cca_align")

      setProgress(value = 0.35)
      
      dims <- isolate(as.numeric(input$align_dims))
      resolution <- isolate(as.numeric(input$ccs_res))

      immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:dims)
      immune.combined <<- IntegrateData(anchorset = immune.anchors, dims = 1:dims)
      
      DefaultAssay(immune.combined) <<- "integrated"
      immune.combined <<- immune.combined
      
      # Run the standard workflow for visualization and clustering
      #immune.combined <<- ScaleData(immune.combined, verbose = FALSE)
      immune.combined <<- ScaleData(immune.combined, vars.to.regress = c("percent.mito"), verbose = FALSE)
      #immune.combined <<- ScaleData(immune.combined, vars.to.regress = c("CC.Difference", "percent.mito"), verbose = FALSE)
      immune.combined <<- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
      
      setProgress(value = 0.65)
      # UMAP and Clustering
      immune.combined <<- RunUMAP(immune.combined, reduction = "pca", dims = 1:dims)
      immune.combined <<- FindNeighbors(immune.combined, reduction = "pca", dims = 1:dims)
      immune.combined <<- FindClusters(immune.combined, resolution = resolution)
      
      setProgress(value = 0.8)

      output$plot_cca_tsne <- renderPlot({

        # Visualization
        p1_dim <<- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
        p2_dim <<- DimPlot(immune.combined, reduction = "umap", label = TRUE)
        plot_grid(p1_dim, p2_dim)

      })
      
      output$plot_cca_side <- renderPlot({
        # Visualization
        p3_dim <<- DimPlot(immune.combined, reduction = "umap", split.by = "stim")
        p3_dim
      })

      setProgress(value = 1)
      shinyjs::hide("load_cca_align")
      shinyjs::show("check_cca_align")
    })

    output$old_name_cca_align <- renderUI({

      selectizeInput(ns("ccs_old_name"), label = "Select names to change: ",
                     choices = levels(x = immune.combined),
                     selected = levels(x = immune.combined)[1],
                     multiple = TRUE)
    }
    )

    output$id_con_genes <- renderUI({

      selectizeInput(ns("con_genes_id"), label="Find conserved markers of cluster:",
                     choices = levels(x = immune.combined),
                     selected = levels(x = immune.combined)[1],
                     multiple = FALSE)
    }
    )


    output$id1_con_genes <- renderUI({

      selectizeInput(ns("con_genes_id1"), label="Find all markers of cluster:",
                     choices = levels(x = immune.combined),
                     selected = levels(x = immune.combined)[1],
                     multiple = FALSE)
    }
    )

    output$id2_con_genes <- renderUI({

      selectizeInput(ns("con_genes_id2"), label="Distinguishing from clusters:",
                     choices = levels(x = immune.combined)[levels(x = immune.combined) != input$con_genes_id1],
                     selected = levels(x = immune.combined)[levels(x = immune.combined) != input$con_genes_id1][1],
                     multiple = TRUE)
    }
    )

    output$id_con_genes_diff <- renderUI({

      selectizeInput(ns("con_genes_id_diff"), label="Find average expression of cluster:",
                     choices = levels(x = immune.combined),
                     selected = levels(x = immune.combined)[1],
                     multiple = FALSE)
    }
    )



  }
  )# Perform integration END

  # Assign new names
  observeEvent(input$assign_ccs_names, {

    output$ccs_tsne_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_cca_tsne"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    withProgress(message = "Plotting UMAP ...",{

      shinyjs::show("load_ccs_names")
      setProgress(value = 0.35)

      current.cluster.ids <- input$ccs_old_name
      new.cluster.ids <- input$new_name_cca_align
      new.cluster.ids <- unlist(strsplit(new.cluster.ids, split = ",|, "))

      validate(need(length(current.cluster.ids) == length(new.cluster.ids), "Error: Number of names should match!"))

      Idents(object = immune.combined) <- plyr::mapvalues(x = Idents(object = immune.combined), 
                                                         from = current.cluster.ids, 
                                                         to = new.cluster.ids)
      immune.combined <<- immune.combined

      output$plot_cca_tsne <- renderPlot({

        DimPlot(object = immune.combined, reduction = "umap", label = TRUE, pt.size = 1.5)
        
      })

      setProgress(value = 1)
      shinyjs::hide("load_ccs_names")
      shinyjs::show("check_ccs_names")
    })

    output$old_name_cca_align <- renderUI({

      selectizeInput(ns("ccs_old_name"), label = "Select names to change: ",
                     choices = levels(x = immune.combined),
                     selected = levels(x = immune.combined)[1],
                     multiple = TRUE)
    }
    )

    output$id_con_genes <- renderUI({

      selectizeInput(ns("con_genes_id"), label="Find conserved markers of cluster:",
                     choices = levels(x = immune.combined),
                     selected = levels(x = immune.combined)[1],
                     multiple = FALSE)
    }
    )

    output$id1_con_genes <- renderUI({

      selectizeInput(ns("con_genes_id1"), label="Find all markers of cluster:",
                     choices = levels(x = immune.combined),
                     selected = levels(x = immune.combined)[1],
                     multiple = FALSE)
    }
    )

    output$id2_con_genes <- renderUI({

      selectizeInput(ns("con_genes_id2"), label="Distinguishing from clusters:",
                     choices = levels(x = immune.combined)[levels(x = immune.combined) != input$con_genes_id1],
                     selected = levels(x = immune.combined)[levels(x = immune.combined) != input$con_genes_id1][1],
                     multiple = TRUE)
    }
    )

    output$id_con_genes_diff <- renderUI({

      selectizeInput(ns("con_genes_id_diff"), label="Find average expression of cluster:",
                     choices = levels(x = immune.combined),
                     selected = levels(x = immune.combined)[1],
                     multiple = FALSE)
    }
    )

  }
  )# Assign new names END

  # Cluster biomarkers
  observeEvent(input$viz_con_genes, {



    output$ccs_genes_table <- renderUI({
      withSpinner(DTOutput(ns("table_ccs_genes")), type = getOption("spinner.type", default = 4))
    }
    )

    output$con_genes_vlnplot <- renderUI({
      withSpinner(plotOutput(ns("vlnplot_con_genes"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    output$con_genes_featureplot <- renderUI({
      withSpinner(plotOutput(ns("featureplot_con_genes"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    output$con_genes_ridgeplot <- renderUI({
      withSpinner(plotOutput(ns("ridgeplot_con_genes"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )


    withProgress(message = "Clustering biomarkers ...",{

      shinyjs::show("load_con_genes")
      setProgress(value = 0.35)

      if(input$con_genes_radio == "single_con"){
        cluster.markers <- FindMarkers(object = immune.combined, ident.1 = input$con_genes_id1, min.pct = 0.25)
        is.num <- sapply(cluster.markers, is.numeric)
        cluster.markers[is.num] <- lapply(cluster.markers[is.num], round, 3)

        output$table_ccs_genes <- renderDT(datatable(cluster.markers, rownames = TRUE,
                                                      selection = list(selected = 1, target = 'row')
        ))

      }
      else if(input$con_genes_radio == "multiple_con"){
        cluster.markers <- FindMarkers(object = immune.combined, ident.1 = input$con_genes_id1, ident.2 = input$con_genes_id2,
                                       min.pct = 0.25)
        is.num <- sapply(cluster.markers, is.numeric)
        cluster.markers[is.num] <- lapply(cluster.markers[is.num], round, 3)

        output$table_ccs_genes <- renderDT(datatable(cluster.markers, rownames = TRUE,
                                                      selection = list(selected = 1, target = 'row')
        ))
      }
      else{
        cluster.markers <- FindAllMarkers(object = immune.combined, only.pos = TRUE, min.pct = 0.25,
                                          logfc.threshold = 0.25)

        cluster.markers <- cluster.markers %>% group_by(cluster) %>% top_n(input$top_markers_con, avg_logFC)

        output$con_genes_heatmap <- renderUI({
          withSpinner(plotOutput(ns("heatmap_con_genes"), height = "750px"), type = getOption("spinner.type", default = 4))
        }
        )

        output$heatmap_con_genes <- renderPlot({
          DoHeatmap(object = immune.combined, features = cluster.markers$gene) + NoLegend()
        })

        is.num <- sapply(cluster.markers, is.numeric)
        cluster.markers[is.num] <- lapply(cluster.markers[is.num], round, 3)

        output$table_ccs_genes <- renderDT(datatable(cluster.markers, rownames = FALSE,
                                                      selection =  list(selected = 1, target = 'row')
        ))

      }

      setProgress(value = 1)
      shinyjs::hide("load_con_genes")
      shinyjs::show("check_con_genes")
    })

    if(input$con_genes_radio == "all_con"){

      output$vlnplot_con_genes <- renderPlot({

        validate(need(length(input$table_ccs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        p7_integrated <<- VlnPlot(object = immune.combined,
                features = cluster.markers$gene[input$table_ccs_genes_rows_selected])
        p7_integrated
      })

      output$featureplot_con_genes <- renderPlot({

        features_plot_integrated <<- cluster.markers$gene[input$table_ccs_genes_rows_selected]
        
        validate(need(length(input$table_ccs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        FeaturePlot(object = immune.combined,
                    features = cluster.markers$gene[input$table_ccs_genes_rows_selected])
      })

      output$ridgeplot_con_genes <- renderPlot({
        
        features_plot_integrated <<- cluster.markers$gene[input$table_ccs_genes_rows_selected]

        validate(need(length(input$table_ccs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        RidgePlot(object = immune.combined,
                  features = cluster.markers$gene[input$table_ccs_genes_rows_selected], ncol = 2)
      })


    }
    else{

      output$vlnplot_con_genes <- renderPlot({

        validate(need(length(input$table_ccs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        p7_integrated <<- VlnPlot(object = immune.combined,
                features = rownames(cluster.markers)[input$table_ccs_genes_rows_selected])
        p7_integrated
      })

      output$featureplot_con_genes <- renderPlot({

        features_plot_integrated <<- rownames(cluster.markers)[input$table_ccs_genes_rows_selected]
        
        validate(need(length(input$table_ccs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        FeaturePlot(object = immune.combined,
                    features = rownames(cluster.markers)[input$table_ccs_genes_rows_selected])
      })

      output$ridgeplot_con_genes <- renderPlot({
        
        features_plot_integrated <<- rownames(cluster.markers)[input$table_ccs_genes_rows_selected]

        validate(need(length(input$table_ccs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        RidgePlot(object = immune.combined,
                  features = rownames(cluster.markers)[input$table_ccs_genes_rows_selected], ncol = 2)
      })

    }

  }
  )# Cluster biomarkers END

  # Conserved biomarkers
  observeEvent(input$viz_con_genes_con, {



    output$con_genes_table <- renderUI({
      withSpinner(DTOutput(ns("table_con_genes")), type = getOption("spinner.type", default = 4))
    }
    )

    output$con_genes_featureplot_con <- renderUI({
      withSpinner(plotOutput(ns("featureplot_con_genes_con"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    output$con_genes_splitdotplot <- renderUI({
      withSpinner(plotOutput(ns("splitdotplot_con_genes"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )


    withProgress(message = "Clustering biomarkers ...",{

      shinyjs::show("load_con_genes_con")
      setProgress(value = 0.35)

      # Feature plot (conserved markers)
      
      DefaultAssay(immune.combined) <<- "RNA"
      Idents(immune.combined) <<- "seurat_clusters"
      immune.combined <<- immune.combined
      
      nk.markers <- FindConservedMarkers(immune.combined, ident.1 = input$con_genes_id, grouping.var = "stim")

      nk.markers <- nk.markers[, -c(5, 10, 11, 12)]

      is.num <- sapply(nk.markers, is.numeric)
      nk.markers[is.num] <- lapply(nk.markers[is.num], round, 3)

      output$table_con_genes <- renderDT(datatable(nk.markers, rownames = TRUE,
                                                   selection = list(selected = 1, target = 'row'),
                                                   options = list(lengthMenu = c(5, 6, 7), pageLength = 5))

                                         )

      output$featureplot_con_genes_con <- renderPlot({

        features_plot_conserved <<- rownames(nk.markers)[input$table_con_genes_rows_selected]
        
        validate(need(length(input$table_con_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        FeaturePlot(object = immune.combined,
                    features = rownames(nk.markers)[input$table_con_genes_rows_selected],
                    min.cutoff = "q9", pt.size = 1)
      })

      # SplitDotPlot
      output$splitdotplot_con_genes <- renderPlot({

        validate(need(length(input$table_con_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        
        sdp <<- DotPlot(immune.combined, features = rev(rownames(nk.markers)[input$table_con_genes_rows_selected]), 
                cols = c("blue", "red"), dot.scale = 8, 
                split.by = "stim") + RotatedAxis()
        sdp
      })

      setProgress(value = 1)
      shinyjs::hide("load_con_genes_con")
      shinyjs::show("check_con_genes_con")
    })

  }
  )# Conserved biomarkers END

  # Identify differential expressed genes
  observeEvent(input$viz_ccs_diff, {



    output$acs_genes_table <- renderUI({
      withSpinner(DTOutput(ns("table_acs_genes")), type = getOption("spinner.type", default = 4))
    }
    )

    output$acs_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_acs"), height = "500px", click = ns("plot_acs_click")), type = getOption("spinner.type", default = 4))
    }
    )

    output$ccs_heatmap_diff <- renderUI({
      withSpinner(plotOutput(ns("heatmap_ccs_diff"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    output$ccs_vlnplot_diff <- renderUI({
      withSpinner(plotOutput(ns("vlnplot_ccs_diff"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    withProgress(message = "Generating plots ...",{

      shinyjs::show("load_ccs_diff")
      setProgress(value = 0.35)

      t.cells <- subset(immune.combined, idents = input$con_genes_id_diff)
      Idents(t.cells) <- "stim"
      avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
      avg.t.cells$gene <- rownames(avg.t.cells)

      immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
      immune.combined$celltype <- Idents(immune.combined)
      Idents(immune.combined) <- "celltype.stim"
      immune.combined <<- immune.combined
      
      b.interferon.response <- FindMarkers(immune.combined, ident.1 = paste0(input$con_genes_id_diff, "_STIM"),
                                           ident.2 = paste0(input$con_genes_id_diff, "_CTRL"),
                                           verbose = FALSE)

      is.num <- sapply(b.interferon.response, is.numeric)
      b.interferon.response[is.num] <- lapply(b.interferon.response[is.num], round, 3)

      output$table_acs_genes <- renderDT(datatable(b.interferon.response, rownames = TRUE,
                                                        selection = list(selected = 1, target = 'row'),
                                                        options = list(lengthMenu = c(5, 6, 7), pageLength = 5))

      )

      # Scatter plot
      output$plot_acs <- renderPlot({

        #validate(need(length(input$table_ccs_genes_diff_rows_selected) != "0", "Error: Please select at least one feature!"))
        p8_integrated <<- LabelUL(ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle(input$con_genes_id_diff),
                genes = rownames(b.interferon.response)[input$table_acs_genes_rows_selected],
                avg.t.cells, adj.u.t = 0.5, adj.u.s = 0.4,
                adj.l.t = 0.25, adj.l.s = 0.25)
        p8_integrated
      })

      observeEvent(input$plot_acs_click, {
        g <- nearPoints(avg.t.cells, input$plot_acs_click, addDist = TRUE)$gene
        g_n <- which(rownames(b.interferon.response) == g)
        s_n <- isolate(as.numeric(input$table_acs_genes_rows_selected))
        output$table_acs_genes <- renderDT(datatable(b.interferon.response, rownames = TRUE,
                                                          selection = list(selected = c(s_n ,g_n), target = 'row'),
                                                          options = list(lengthMenu = c(5, 6, 7), pageLength = 5)))
      })

      # FeaturePlots
      output$heatmap_ccs_diff <- renderPlot({

        features_heatmap <<- rownames(b.interferon.response)[input$table_acs_genes_rows_selected]
        
        validate(need(length(input$table_acs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        FeaturePlot(immune.combined,
                    features = rownames(b.interferon.response)[input$table_acs_genes_rows_selected],
                    split.by = "stim", pt.size = 1.5, max.cutoff = 3, 
                    cols = c("grey", "red"))
      })
      
      # VlnPlot
      output$vlnplot_ccs_diff <- renderPlot({
        
        features_vlnplot <<- rownames(b.interferon.response)[input$table_acs_genes_rows_selected]
        
        validate(need(length(input$table_acs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        
        plots <<- VlnPlot(immune.combined, features = rownames(b.interferon.response)[input$table_acs_genes_rows_selected], 
                         split.by = "stim", group.by = "celltype", 
                         pt.size = 0, combine = FALSE)
        CombinePlots(plots = plots, ncol = 1)
      })

      setProgress(value = 1)
      shinyjs::hide("load_ccs_diff")
      shinyjs::show("check_ccs_diff")
    })

  }
  )# Identify differential expressed genes END
  
  
  # Save
  observeEvent(input$do_save_integrated, {
    
    
    withProgress(message = "Saving ...",{
      
      shinyjs::show("load_save_integrated")
      setProgress(value = 0.35)
      
      save_Path <- isolate(as.character(parseDirPath(volumes, input$dir_save_integrated)))
      
      saveRDS(immune.combined, file = file.path(save_Path, paste0(input$saveName_integrated, ".rds"))) 
      
      setProgress(value = 1)
      shinyjs::hide("load_save_integrated")
      shinyjs::show("check_save_integrated")
    })
    
  }
  
  )# Save END
  
  
  # Load
  observeEvent(input$do_load_integrated, {
    
    
    withProgress(message = "Loading ...",{
      
      shinyjs::show("load_read_integrated")
      setProgress(value = 0.35)
      
      load_File <- isolate(as.character(parseFilePaths(volumes, input$integrated_obj)[4]))
      
      immune.combined <<- readRDS(load_File) 
      
      setProgress(value = 1)
      shinyjs::hide("load_read_integrated")
      shinyjs::show("check_read_integrated")
    })
    
  }
  
  )# Load END
  
  
  output$figure_save_integrated <- renderUI({
    withSpinner(plotOutput(ns("save_figure_integrated"), width = paste0(input$pWidth_integrated/2, "px"), 
                           height = paste0(input$pHeight_integrated/2, "px")),
                type = getOption("spinner.type", default = 4))
  }
  )
  
  output$save_figure_integrated <- renderPlot({
    
    switch(input$psave_select_integrated,
           "QC" = p1_integrated,
           "GenePlot" = {
             gp1 <- FeatureScatter(object = pbmc.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
             gp2 <- FeatureScatter(object = pbmc.combined, feature1 = "nCount_RNA", feature2 = "percent.mito") + NoLegend()
             CombinePlots(plots = list(gp1, gp2))
           },
           "CCA_plot" = plot_grid(p1_dim, p2_dim),
           "CCA_side" = p3_dim,
           "VlnPlot_markers" = p7_integrated,
           "FeaturePlot" = FeaturePlot(object = immune.combined,
                                       features = features_plot_integrated),
           "RidgePlot" = RidgePlot(object = immune.combined,
                                   features = features_plot_integrated, ncol = 2),
           "Features_conserved" = FeaturePlot(object = immune.combined,
                                       features = features_plot_conserved),
           "SplitDotPlotGG" = sdp,
           "Scatter_plot" = p8_integrated,
           "FeatureHeatmap" = FeaturePlot(immune.combined,
                                          features = features_heatmap,
                                          split.by = "stim", pt.size = 1.5, max.cutoff = 3, 
                                          cols = c("grey", "red")),
           "VlnPlot_side" = CombinePlots(plots = plots, ncol = 1)
    )
  })
  
  plotInput = function() {
    
    switch(input$psave_select_integrated,
           "QC" = p1_integrated,
           "GenePlot" = {
             gp1 <- FeatureScatter(object = pbmc.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
             gp2 <- FeatureScatter(object = pbmc.combined, feature1 = "nCount_RNA", feature2 = "percent.mito") + NoLegend()
             CombinePlots(plots = list(gp1, gp2))
           },
           "CCA_plot" = plot_grid(p1_dim, p2_dim),
           "CCA_side" = p3_dim,
           "VlnPlot_markers" = p7_integrated,
           "FeaturePlot" = FeaturePlot(object = immune.combined,
                                       features = features_plot_integrated),
           "RidgePlot" = RidgePlot(object = immune.combined,
                                   features = features_plot_integrated, ncol = 2),
           "Features_conserved" = FeaturePlot(object = immune.combined,
                                              features = features_plot_conserved),
           "SplitDotPlotGG" = sdp,
           "Scatter_plot" = p8_integrated,
           "FeatureHeatmap" = FeaturePlot(immune.combined,
                                          features = features_heatmap,
                                          split.by = "stim", pt.size = 1.5, max.cutoff = 3, 
                                          cols = c("grey", "red")),
           "VlnPlot_side" = CombinePlots(plots = plots, ncol = 1)
    )
  }
  
  output$do_psave_integrated = downloadHandler(
    filename = "figure.png",
    content = function(file) {
      
      png(file, width = input$pWidth_integrated, height = input$pHeight_integrated,
          res = input$pRes_integrated, pointsize = input$pPsize_integrated, type = "cairo")
      print(plotInput())
      dev.off()
    }
  )
}
