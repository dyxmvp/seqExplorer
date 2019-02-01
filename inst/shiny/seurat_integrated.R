

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

                       box(title = "3. Cell selection", width = 2, height = 400, solidHeader = TRUE, status = "danger",
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
                         title = "4. Normalize and scale data", width = 3, height = 400, solidHeader = TRUE, status = "success",
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
                         title = "3.1 VlnPlot (Cell selection)",  width = 8, height = 600, status = "danger",

                         sidebarPanel(
                           tableOutput(ns("nGene_range_integrated")),
                           uiOutput(ns("gene_low_high_integrated")),
                           hr(),
                           tableOutput(ns("mito_range_integrated")),
                           uiOutput(ns("mito_low_high_integrated"))
                           #uiOutput(ns("nGene_range")),
                           #uiOutput(ns("mito_range"))
                         ),
                         mainPanel(
                           uiOutput(ns("filter_plots_integrated"))
                         )

                       ),
                       box(
                         title = "3.2 GenePlot (Cell selection)",  width = 4, height = 600, status = "danger",
                         uiOutput(ns("gene_plots_integrated"))
                       )
                     )
            ), # Tab 1 END

            tabPanel('2. Canonical correlation analysis (CCA)',
                     fluidRow(

                       box(title = "CCA", width = 3, height = 600, solidHeader = TRUE, status = "warning",
                           collapsible = T,


                           selectizeInput(ns("mean_fun_integrated"), label="Mean function",
                                          choices = c("ExpMean"),
                                          selected = c("ExpMean"),
                                          multiple = FALSE),

                           selectizeInput(ns("disper_fun_integrated"), label="Dispersion function",
                                          choices = c("LogVMR"),
                                          selected = c("LogVMR"),
                                          multiple = FALSE),

                           splitLayout(

                           numericInput(ns("x_low_integrated"),
                                        label="X low_cutoff", min = 0, max = Inf, value = 0.1, step = 0.0001),
                           numericInput(ns("x_high_integrated"),
                                        label="X high_cutoff", min = 0, max = Inf, value = 8, step = 0.1)
                           ),

                           splitLayout(
                           numericInput(ns("y_cut_integrated"),
                                        label="Y cutoff", min = 0, max = Inf, value = 1, step = 0.1),

                           numericInput(ns("union_genes_integrated"),
                                        label="Take union of (?) genes", min = 0, max = Inf, value = 1000)
                           ),

                           splitLayout(
                             numericInput(ns("cca_dim1_print"),
                                          label="Print from CC: ", min = 1, max = Inf, value = 1),
                             numericInput(ns("cca_dim2_print"),
                                          label="to CC: ", min = 1, max = Inf, value = 2)
                           ),

                           splitLayout(
                             numericInput(ns("cca_genes_print"),
                                          label="Number of genes to print", min = 1, max = Inf, value = 10),
                             numericInput(ns("cc_num"),
                                          label="Number of cc to calculate", min = 1, max = Inf, value = 30)
                           ),


                           withBusyIndicatorUI(icon_name = "cca_integrated",
                                               actionButton(ns("find_cca_integrated"),
                                                            "CCA",
                                                            style = "width: 80%")
                           )
                       ),
                       box(
                         title = "CCs",  width = 4, height = 600, status = "success",
                         uiOutput(ns("cca_table_integrated"))
                       ),
                       box(
                         title = "CCA plots",  width = 5, height = 600, status = "primary",
                         uiOutput(ns("cca_plots_integrated"))
                       )
                     )
            ), # Tab 2 END

            tabPanel('3. Choose CCs',
                     fluidRow(
                       column(width = 5,
                              box(title = "Choose CCs", width = NULL, solidHeader = TRUE, status = "warning",
                                  collapsible = T,

                                  numericInput(ns("dims_eval"),
                                               label="Dimensions to evalutate: ", min = 1, max = Inf, value = 30),

                                  splitLayout(
                                    numericInput(ns("ccs_dims_heatmap"),
                                                 label="Dimensions to show in Heatmap:", min = 1, max = Inf, value = 9),
                                    numericInput(ns("ccs_cells_heatmap"),
                                                 label="Cells to use in heatmap:", min = 1, max = Inf, value = 500)
                                  ),

                                  withBusyIndicatorUI(icon_name = "ccs",
                                                      actionButton(ns("viz_ccs"),
                                                                   "CCs Plots",
                                                                   style = "width: 90%")
                                  )),

                              box(title = "CCs plot", width = NULL, height = 580, status = "success",
                                  uiOutput(ns("ccs_plot"))
                              )
                       ), #

                       column(width = 7,
                              box(
                                title = "CCs Heatmap",  width = NULL, height = 850, status = "primary",
                                uiOutput(ns("ccs_heatmap"))
                              )
                       )#
                     )
            ), # Tab 3 END

            tabPanel('4. Align CCA subspaces',
                     fluidRow(
                       column(width = 5,
                              box(title = "Align CCA", width = NULL, solidHeader = TRUE, status = "warning",
                                  collapsible = T,

                                  splitLayout(
                                    numericInput(ns("align_dims"),
                                                 label="Dimensions to align: ", min = 1, max = Inf, value = 20),
                                    numericInput(ns("ccs_res"),
                                                 label="Resolution (TSNE): ", min = 0, max = Inf, value = 0.6, step = 0.01),
                                    numericInput(ns("ccs_perp"),
                                                 label="Perplexity (TSNE): ", min = 1, max = Inf, value = 30)
                                  ),

                                  withBusyIndicatorUI(icon_name = "cca_align",
                                                      actionButton(ns("viz_cca_align"),
                                                                   "Align CCA",
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
                                  ),

                              box(title = "Aligned CCA plot", width = NULL, height = 500, status = "success",

                                  uiOutput(ns("align_cca_plot"))

                                  )

                       ), #

                       column(width = 7,
                              box(
                                title = "TSNE Plot",  width = NULL, height = 850, status = "primary",
                                uiOutput(ns("ccs_tsne_plot"))
                              )
                       )#
                     )
            ), # Tab 4 END

            tabPanel('5. Identify conserved cell type markers',

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
            ), # Tab 5 END

            tabPanel('6. Differential genes across conditions',
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
                                title = "Differential genes heatmap",  width = NULL, height = 850, status = "primary",
                                uiOutput(ns("ccs_heatmap_diff"))
                              )
                       )#
                     )
            ), # Tab 6 END
            
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
            ), # Tab 7 END
            
            tabPanel('* Save figures',
                     fluidRow(
                       
                       box(
                         title = "Save figure", width = 3, solidHeader = TRUE, status = "primary",
                         
                         selectizeInput(ns("psave_select_integrated"), label="Choose a figure to save",
                                        choices = c("QC", "GenePlot", "CCA_plot", "MetageneBicorPlot", "DimHeatmap",
                                                    "CCA_aligned", "TSNEPlot", "VlnPlot_markers", "FeaturePlot", 
                                                    "RidgePlot", "Features_conserved", "SplitDotPlotGG", 
                                                    "Scatter_plot", "FeatureHeatmap"),
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
            )# Tab 8 END
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
      ctrl <<- CreateSeuratObject(raw.data = ctrl.data, min.cells = input$min_cells_integrated,
                                  min.genes = input$min_genes_integrated,
                                  project = input$project_name_ctrl)
      ctrl@meta.data$stim <- "CTRL"
      ctrl <<- ctrl

      mito.genes.ctrl <- grep(pattern = "^MT-", x = rownames(x = ctrl@data), value = TRUE, ignore.case = TRUE)
      percent.mito.ctrl <- Matrix::colSums(ctrl@raw.data[mito.genes.ctrl, ])/Matrix::colSums(ctrl@raw.data)
      ctrl <<- AddMetaData(object = ctrl, metadata = percent.mito.ctrl, col.name = "percent.mito.ctrl")

      # Set up stimulated object
      stim <<- CreateSeuratObject(raw.data = stim.data, min.cells = input$min_cells_integrated,
                                  min.genes = input$min_genes_integrated,
                                  project = input$project_name_stim)
      stim@meta.data$stim <- "STIM"
      stim <<- stim

      mito.genes.stim <- grep(pattern = "^MT-", x = rownames(x = stim@data), value = TRUE, ignore.case = TRUE)
      percent.mito.stim <- Matrix::colSums(stim@raw.data[mito.genes.stim, ])/Matrix::colSums(stim@raw.data)
      stim <<- AddMetaData(object = stim, metadata = percent.mito.stim, col.name = "percent.mito.stim")

      setProgress(value = 0.7)

      # Set up combined object
      pbmc.combined <<- MergeSeurat(object1 = ctrl, object2 = stim, add.cell.id1 = "ctrl",
                                    add.cell.id2 = "stim", project = "Integrated_analysis")

      mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc.combined@data), value = TRUE, ignore.case = TRUE)
      percent.mito <- Matrix::colSums(pbmc.combined@raw.data[mito.genes, ])/Matrix::colSums(pbmc.combined@raw.data)
      pbmc.combined <<- AddMetaData(object = pbmc.combined, metadata = percent.mito, col.name = "percent.mito")

      setProgress(value = 1)

      shinyjs::hide("load_sinit_integrated")
      shinyjs::show("check_sinit_integrated")
    })

  }

  )# Object initialization END

  #Cell selection
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

    min_gene_ctrl <- min(ctrl@meta.data$nGene)
    max_gene_ctrl <- max(ctrl@meta.data$nGene)

    min_gene_stim <- min(stim@meta.data$nGene)
    max_gene_stim <- max(stim@meta.data$nGene)

    min_mito_ctrl <- round(min(ctrl@meta.data$percent.mito.ctrl), 2)
    max_mito_ctrl <- round(max(ctrl@meta.data$percent.mito.ctrl), 2)

    min_mito_stim <- round(min(stim@meta.data$percent.mito.stim), 2)
    max_mito_stim <- round(max(stim@meta.data$percent.mito.stim), 2)

    output$nGenePlot_integrated <- renderPlot({

      if("nGene" %in% input$subsetNames_integrated){

        VlnPlot(object = pbmc.combined, features.plot = c("nGene"), nCol = 1) +
          geom_hline(yintercept = input$gene_low_integrated, color = 'darkgreen', linetype = "dashed", size = 1.5) +
          geom_text(data = data.frame(x=1.2,y=input$gene_low_integrated), aes(x, y), label="low", vjust=1.5, hjust=0,color = "darkgreen",size = 5,fontface = "bold") +
          geom_hline(yintercept = input$gene_high_integrated, color = 'blue', linetype = "dashed", size = 1.5) +
          geom_text(data = data.frame(x=1.2,y=input$gene_high_integrated), aes(x, y), label="high", vjust=-1, hjust=0,color = "blue",size = 5,fontface = "bold") +
          scale_y_continuous(limits=c(0, 1.1*max_gene_ctrl))
      }
      else{
        VlnPlot(object = pbmc.combined, features.plot = c("nGene"), nCol = 1)
      }
    })
    output$nUMIPlot_integrated <- renderPlot({
      VlnPlot(object = pbmc.combined, features.plot = c("nUMI"), nCol = 1)

    })
    output$mitoPlot_integrated <- renderPlot({

      if("percent.mito" %in% input$subsetNames_integrated){

        VlnPlot(object = pbmc.combined, features.plot = c("percent.mito"), nCol = 1) +
          geom_hline(yintercept = input$mito_low_integrated, color = 'darkgreen', linetype = "dashed", size = 1.5) +
          geom_text(data = data.frame(x=1.2,y=input$mito_low_integrated), aes(x, y), label="low", vjust=1.5, hjust=0,color = "darkgreen",size = 5,fontface = "bold") +
          geom_hline(yintercept = input$mito_high_integrated, color = 'blue', linetype = "dashed", size = 1.5) +
          geom_text(data = data.frame(x=1.2,y=input$mito_high_integrated), aes(x, y), label="high", vjust=-1, hjust=0,color = "blue",size = 5,fontface = "bold") +
          scale_y_continuous(limits=c(min_mito_ctrl-0.01, max_mito_ctrl+0.01))
      }
      else{
        VlnPlot(object = pbmc.combined, features.plot = c("percent.mito"), nCol = 1)
      }
    })
    
    p1_integrated <<- VlnPlot(object = pbmc.combined, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
    
    median_gene_ctrl <- median(ctrl@meta.data$nGene)
    median_gene_stim <- median(stim@meta.data$nGene)

    min_UMI_ctrl <- min(ctrl@meta.data$nUMI)
    max_UMI_ctrl <- max(ctrl@meta.data$nUMI)
    median_UMI_ctrl <- median(ctrl@meta.data$nUMI)

    min_UMI_stim <- min(stim@meta.data$nUMI)
    max_UMI_stim <- max(stim@meta.data$nUMI)
    median_UMI_stim <- median(stim@meta.data$nUMI)

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

    median_mito_ctrl <- median(ctrl@meta.data$percent.mito.ctrl)
    median_mito_stim <- median(stim@meta.data$percent.mito.stim)
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

        GenePlot(object = pbmc.combined, gene1 = "nUMI", gene2 = "percent.mito")
        abline(h = input$mito_low_integrated, col="darkgreen", lwd=2, lty="dashed")
        abline(h = input$mito_high_integrated, col="blue", lwd=2, lty="dashed")
      }
      else{
        GenePlot(object = pbmc.combined, gene1 = "nUMI", gene2 = "percent.mito")
      }
    })

    output$gene_nUMI_integrated <- renderPlot({

      if("nGene" %in% input$subsetNames_integrated){

        GenePlot(object = pbmc.combined, gene1 = "nUMI", gene2 = "nGene")
        abline(h = input$gene_low_integrated, col="darkgreen", lwd=2, lty="dashed")
        abline(h = input$gene_high_integrated, col="blue", lwd=2, lty="dashed")
      }
      else{
        GenePlot(object = pbmc.combined, gene1 = "nUMI", gene2 = "nGene")
      }
    })

    shinyjs::hide("load_cfilter_integrated")
    shinyjs::show("check_cfilter_integrated")

  }

  )# Cell selection END


  # Normalization
  observeEvent(input$data_norm_integrated, {

    withProgress(message = "Normalization and scaling...",{

      shinyjs::show("load_norm_integrated")
      setProgress(value = 0.35)

      if("nGene" %in% input$subsetNames_integrated && "percent.mito" %in% input$subsetNames_integrated){

        #print("c1")

         ctrl <<- FilterCells(object = ctrl, subset.names = c("nGene", "percent.mito.ctrl"),
                              low.thresholds = c(input$gene_low_integrated, input$mito_low_integrated),
                              high.thresholds = c(input$gene_high_integrated, input$mito_high_integrated))

         stim <<- FilterCells(object = stim, subset.names = c("nGene", "percent.mito.stim"),
                              low.thresholds = c(input$gene_low_integrated, input$mito_low_integrated),
                              high.thresholds = c(input$gene_high_integrated, input$mito_high_integrated))
      }
      else if("nGene" %in% input$subsetNames_integrated){

        #print("c2")

        ctrl <<- FilterCells(object = ctrl, subset.names = "nGene",
                             low.thresholds = input$gene_low_integrated,
                             high.thresholds = input$gene_high_integrated)

        stim <<- FilterCells(object = stim, subset.names = "nGene",
                             low.thresholds = input$gene_low_integrated,
                             high.thresholds = input$gene_high_integrated)
      }
      else{

        #print("c3")

        ctrl <<- FilterCells(object = ctrl, subset.names = "percent.mito.ctrl",
                             low.thresholds = input$mito_low_integrated,
                             high.thresholds = input$mito_high_integrated)

        stim <<- FilterCells(object = stim, subset.names = "percent.mito.stim",
                             low.thresholds = input$mito_low_integrated,
                             high.thresholds = input$mito_high_integrated)

      }

      ctrl <<- NormalizeData(object = ctrl, normalization.method = input$norm_methods_integrated,
                             scale.factor = input$scale_integrated)
      stim <<- NormalizeData(object = stim, normalization.method = input$norm_methods_integrated,
                             scale.factor = input$scale_integrated)

      if(length(input$vars_regress_integrated) > 0){
        if("nUMI" %in% input$vars_regress_integrated && "percent.mito" %in% input$vars_regress_integrated){
          ctrl <<- ScaleData(object = ctrl, vars.to.regress = c("nUMI", "percent.mito.ctrl"), display.progress = F)
          stim <<- ScaleData(object = stim, vars.to.regress = c("nUMI", "percent.mito.stim"), display.progress = F)
        }
        else if("percent.mito" %in% input$vars_regress_integrated){
          ctrl <<- ScaleData(object = ctrl, vars.to.regress = c("percent.mito.ctrl"), display.progress = F)
          stim <<- ScaleData(object = stim, vars.to.regress = c("percent.mito.stim"), display.progress = F)
        }
        else{
          ctrl <<- ScaleData(object = ctrl, vars.to.regress = input$vars_regress_integrated, display.progress = F)
          stim <<- ScaleData(object = stim, vars.to.regress = input$vars_regress_integrated, display.progress = F)
        }
      }
      else{
        ctrl <<- ScaleData(ctrl, display.progress = F)
        stim <<- ScaleData(stim, display.progress = F)
      }

      setProgress(value = 1)
      shinyjs::hide("load_norm_integrated")
      shinyjs::show("check_norm_integrated")
    })

  }

  )# Normalization END

  # CCA
  observeEvent(input$find_cca_integrated, {

    shinyjs::show("load_cca_integrated")

    output$cca_plots_integrated <- renderUI({

      withSpinner(plotOutput(ns("plot_cca_integrated"), height = "450px"), type = getOption("spinner.type", default = 4))

    }
    )

    output$cca_table_integrated <- renderUI({
      withSpinner(DTOutput(ns("cca_genes_integrated")), type = getOption("spinner.type", default = 4))
    }
    )

    withProgress(message = "Performing CCA ...",{

      setProgress(value = 0.35)

      ctrl <<- FindVariableGenes(object = ctrl, mean.function = input$mean_fun_integrated,
                                 dispersion.function = input$disper_fun_integrated, x.low.cutoff = input$x_low_integrated,
                                 x.high.cutoff = input$x_high_integrated, y.cutoff = input$y_cut_integrated,
                                 do.plot = F)

      stim <<- FindVariableGenes(object = stim, mean.function = input$mean_fun_integrated,
                                 dispersion.function = input$disper_fun_integrated, x.low.cutoff = input$x_low_integrated,
                                 x.high.cutoff = input$x_high_integrated, y.cutoff = input$y_cut_integrated,
                                 do.plot = F)

      g.1 <- head(rownames(ctrl@hvg.info), input$union_genes_integrated)
      g.2 <- head(rownames(stim@hvg.info), input$union_genes_integrated)

      genes.use <- unique(c(g.1, g.2))
      genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
      genes.use <- intersect(genes.use, rownames(stim@scale.data))

      setProgress(value = 0.65)

      immune.combined <<- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = input$cc_num)


      setProgress(value = 1)

      output$plot_cca_integrated <- renderPlot({
        p2_integrated <<- DimPlot(object = immune.combined, reduction.use = "cca", group.by = "stim",
                      pt.size = 0.5, do.return = TRUE)
        p3_integrated <<- VlnPlot(object = immune.combined, features.plot = "CC1", group.by = "stim",
                      do.return = TRUE)
        plot_grid(p2_integrated, p3_integrated)
      })

      shinyjs::hide("load_cca_integrated")
      shinyjs::show("check_cca_integrated")


    })

    ccaGenes <- capture.output(PrintDim(object = immune.combined, reduction.type = "cca", dims.print = input$cca_dim1_print:input$cca_dim2_print,
                                         genes.print = input$cca_genes_print), type = "output")
    ccaGenes <- gsub("\\D[[1-9]]", "", ccaGenes)
    ccaGenes <- matrix(ccaGenes)
    ccaGenes <- ccaGenes[!apply(ccaGenes == " \"\"", 1, all),]
    ccaGenes <- unlist(strsplit(ccaGenes, " "))
    ccaGenes <- matrix(ccaGenes)
    ccaGenes <- matrix(ccaGenes[!apply(ccaGenes == "", 1, all),], nrow = 2*input$cca_genes_print+1)
    ccaGenes <- data.frame(ccaGenes)
    ccaGenes <- as.data.frame(sapply(ccaGenes, function(x) gsub("\"", "", x)))
    colnames(ccaGenes) <- as.character(unlist(ccaGenes[1,]))
    ccaGenes <- ccaGenes[-1,]

    output$cca_genes_integrated <- renderDT(datatable(ccaGenes, rownames = FALSE, selection = "none"))

  }

  )# CCA END

  # CCs
  observeEvent(input$viz_ccs, {

    output$ccs_heatmap <- renderUI({
      withSpinner(plotOutput(ns("plot_ccs_heatmap"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    output$ccs_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_ccs"), height = "500px"), type = getOption("spinner.type", default = 4))
    }
    )

    withProgress(message = "Plotting CCs ...",{

      shinyjs::show("load_ccs")
      setProgress(value = 0.35)

      output$plot_ccs_heatmap <- renderPlot({
        DimHeatmap(object = immune.combined, reduction.type = "cca", cells.use = input$ccs_cells_heatmap,
                   dim.use = 1:input$ccs_dims_heatmap, do.balanced = TRUE)
      })

      output$plot_ccs <- renderPlot({
        p4_integrated <<- MetageneBicorPlot(immune.combined, grouping.var = "stim", dims.eval = 1:input$dims_eval,
                                display.progress = FALSE)
      })

      setProgress(value = 1)
      shinyjs::hide("load_ccs")
      shinyjs::show("check_ccs")
    })

  }

  )# CCs END


  # Align the CCA subspaces
  observeEvent(input$viz_cca_align, {

    output$align_cca_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_align_cca"), height = "400px"), type = getOption("spinner.type", default = 4))
    }
    )

    output$ccs_tsne_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_cca_tsne"), height = "650px"), type = getOption("spinner.type", default = 4))
    }
    )

    withProgress(message = "Plotting aligned CCA ...",{

      shinyjs::show("load_cca_align")

      setProgress(value = 0.35)
      immune.combined <<- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim",
                                       dims.align = 1:input$align_dims)

      output$plot_align_cca <- renderPlot({

        p5_integrated <<- VlnPlot(object = immune.combined, features.plot = "ACC1", group.by = "stim",
                      do.return = TRUE)
        p6_integrated <<- VlnPlot(object = immune.combined, features.plot = "ACC2", group.by = "stim",
                      do.return = TRUE)
        plot_grid(p5_integrated, p6_integrated)

      })

      setProgress(value = 0.6)

      immune.combined <<- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:input$align_dims,
                                  perplexity = input$ccs_perp, do.fast = T)

      immune.combined <<- FindClusters(immune.combined, reduction.type = "cca.aligned",
                                      resolution = input$ccs_res, dims.use = 1:input$align_dims)

      setProgress(value = 0.8)

      output$plot_cca_tsne <- renderPlot({

        p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 1, group.by = "stim")
        p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 1)
        plot_grid(p1, p2)

      })

      setProgress(value = 1)
      shinyjs::hide("load_cca_align")
      shinyjs::show("check_cca_align")
    })

    output$old_name_cca_align <- renderUI({

      selectizeInput(ns("ccs_old_name"), label = "Select names to change: ",
                     choices = levels(immune.combined@ident),
                     selected = levels(immune.combined@ident)[1],
                     multiple = TRUE)
    }
    )

    output$id_con_genes <- renderUI({

      selectizeInput(ns("con_genes_id"), label="Find conserved markers of cluster:",
                     choices = levels(immune.combined@ident),
                     selected = levels(immune.combined@ident)[1],
                     multiple = FALSE)
    }
    )


    output$id1_con_genes <- renderUI({

      selectizeInput(ns("con_genes_id1"), label="Find all markers of cluster:",
                     choices = levels(immune.combined@ident),
                     selected = levels(immune.combined@ident)[1],
                     multiple = FALSE)
    }
    )

    output$id2_con_genes <- renderUI({

      selectizeInput(ns("con_genes_id2"), label="Distinguishing from clusters:",
                     choices = levels(immune.combined@ident)[levels(immune.combined@ident) != input$con_genes_id1],
                     selected = levels(immune.combined@ident)[levels(immune.combined@ident) != input$con_genes_id1][1],
                     multiple = TRUE)
    }
    )

    output$id_con_genes_diff <- renderUI({

      selectizeInput(ns("con_genes_id_diff"), label="Find average expression of cluster:",
                     choices = levels(immune.combined@ident),
                     selected = levels(immune.combined@ident)[1],
                     multiple = FALSE)
    }
    )



  }
  )# Align the CCA subspaces END

  # Assign new names
  observeEvent(input$assign_ccs_names, {

    output$ccs_tsne_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_cca_tsne"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    withProgress(message = "Plotting TSNE ...",{

      shinyjs::show("load_ccs_names")
      setProgress(value = 0.35)

      current.cluster.ids <- input$ccs_old_name
      new.cluster.ids <- input$new_name_cca_align
      new.cluster.ids <- unlist(strsplit(new.cluster.ids, split = ",|, "))

      validate(need(length(current.cluster.ids) == length(new.cluster.ids), "Error: Number of names should match!"))

      immune.combined@ident <- plyr::mapvalues(x = immune.combined@ident, from = current.cluster.ids, to = new.cluster.ids)

      immune.combined <<- immune.combined

      output$plot_cca_tsne <- renderPlot({

        TSNEPlot(object = immune.combined, do.label = T, do.hover = FALSE, pt.size = 2)

      })

      setProgress(value = 1)
      shinyjs::hide("load_ccs_names")
      shinyjs::show("check_ccs_names")
    })

    output$old_name_cca_align <- renderUI({

      selectizeInput(ns("ccs_old_name"), label = "Select names to change: ",
                     choices = levels(immune.combined@ident),
                     selected = levels(immune.combined@ident)[1],
                     multiple = TRUE)
    }
    )

    output$id_con_genes <- renderUI({

      selectizeInput(ns("con_genes_id"), label="Find conserved markers of cluster:",
                     choices = levels(immune.combined@ident),
                     selected = levels(immune.combined@ident)[1],
                     multiple = FALSE)
    }
    )

    output$id1_con_genes <- renderUI({

      selectizeInput(ns("con_genes_id1"), label="Find all markers of cluster:",
                     choices = levels(immune.combined@ident),
                     selected = levels(immune.combined@ident)[1],
                     multiple = FALSE)
    }
    )

    output$id2_con_genes <- renderUI({

      selectizeInput(ns("con_genes_id2"), label="Distinguishing from clusters:",
                     choices = levels(immune.combined@ident)[levels(immune.combined@ident) != input$con_genes_id1],
                     selected = levels(immune.combined@ident)[levels(immune.combined@ident) != input$con_genes_id1][1],
                     multiple = TRUE)
    }
    )

    output$id_con_genes_diff <- renderUI({

      selectizeInput(ns("con_genes_id_diff"), label="Find average expression of cluster:",
                     choices = levels(immune.combined@ident),
                     selected = levels(immune.combined@ident)[1],
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
                                          thresh.use = 0.25)

        cluster.markers <- cluster.markers %>% group_by(cluster) %>% top_n(input$top_markers_con, avg_logFC)

        output$con_genes_heatmap <- renderUI({
          withSpinner(plotOutput(ns("heatmap_con_genes"), height = "750px"), type = getOption("spinner.type", default = 4))
        }
        )

        output$heatmap_con_genes <- renderPlot({
          DoHeatmap(object = immune.combined, genes.use = cluster.markers$gene, slim.col.label = TRUE, remove.key = TRUE)
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
        p7_integrated <<- VlnPlot(object = immune.combined, x.lab.rot = TRUE,
                features.plot = cluster.markers$gene[input$table_ccs_genes_rows_selected])
        p7_integrated
      })

      output$featureplot_con_genes <- renderPlot({

        features_plot_integrated <<- cluster.markers$gene[input$table_ccs_genes_rows_selected]
        
        validate(need(length(input$table_ccs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        FeaturePlot(object = immune.combined,
                    features.plot = cluster.markers$gene[input$table_ccs_genes_rows_selected],
                    cols.use = c("grey", "blue"),
                    reduction.use = "tsne")
      })

      output$ridgeplot_con_genes <- renderPlot({
        
        features_plot_integrated <<- cluster.markers$gene[input$table_ccs_genes_rows_selected]

        validate(need(length(input$table_ccs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        RidgePlot(object = immune.combined,
                  features.plot = cluster.markers$gene[input$table_ccs_genes_rows_selected], nCol = 2)
      })


    }
    else{

      output$vlnplot_con_genes <- renderPlot({

        validate(need(length(input$table_ccs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        p7_integrated <<- VlnPlot(object = immune.combined, x.lab.rot = TRUE,
                features.plot = rownames(cluster.markers)[input$table_ccs_genes_rows_selected])
        p7_integrated
      })

      output$featureplot_con_genes <- renderPlot({

        features_plot_integrated <<- rownames(cluster.markers)[input$table_ccs_genes_rows_selected]
        
        validate(need(length(input$table_ccs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        FeaturePlot(object = immune.combined,
                    features.plot = rownames(cluster.markers)[input$table_ccs_genes_rows_selected],
                    cols.use = c("grey", "blue"),
                    reduction.use = "tsne")
      })

      output$ridgeplot_con_genes <- renderPlot({
        
        features_plot_integrated <<- rownames(cluster.markers)[input$table_ccs_genes_rows_selected]

        validate(need(length(input$table_ccs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        RidgePlot(object = immune.combined,
                  features.plot = rownames(cluster.markers)[input$table_ccs_genes_rows_selected], nCol = 2)
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
      nk.markers <- FindConservedMarkers(immune.combined, ident.1 = input$con_genes_id, grouping.var = "stim",
                                         print.bar = FALSE)

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
                    features.plot = rownames(nk.markers)[input$table_con_genes_rows_selected],
                    min.cutoff = "q9",
                    cols.use = c("lightgrey", "blue"),
                    pt.size = 0.5)
      })

      # SplitDotPlot
      output$splitdotplot_con_genes <- renderPlot({

        validate(need(length(input$table_con_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        sdp <<- SplitDotPlotGG(immune.combined, genes.plot = rev(rownames(nk.markers)[input$table_con_genes_rows_selected]),
                              cols.use = c("blue", "red"),
                              x.lab.rot = T, plot.legend = T, dot.scale = 8, do.return = T, grouping.var = "stim")
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

    withProgress(message = "Generating plots ...",{

      shinyjs::show("load_ccs_diff")
      setProgress(value = 0.35)

      t.cells <- SubsetData(immune.combined, ident.use = input$con_genes_id_diff, subset.raw = T)
      t.cells <- SetAllIdent(t.cells, id = "stim")
      avg.t.cells <- log1p(AverageExpression(t.cells, show.progress = FALSE))
      avg.t.cells$gene <- rownames(avg.t.cells)

      immune.combined@meta.data$celltype.stim <- paste0(immune.combined@ident, "_",
                                                          immune.combined@meta.data$stim)
      immune.combined <- StashIdent(immune.combined, save.name = "celltype")
      immune.combined <- SetAllIdent(immune.combined, id = "celltype.stim")
      b.interferon.response <- FindMarkers(immune.combined, ident.1 = paste0(input$con_genes_id_diff, "_STIM"),
                                           ident.2 = paste0(input$con_genes_id_diff, "_CTRL"),
                                           print.bar = FALSE)

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

      # FeatureHeatmap
      output$heatmap_ccs_diff <- renderPlot({

        features_heatmap <<- rownames(b.interferon.response)[input$table_acs_genes_rows_selected]
        
        validate(need(length(input$table_acs_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        FeatureHeatmap(immune.combined,
                       features.plot = rownames(b.interferon.response)[input$table_acs_genes_rows_selected],
                       group.by = "stim", pt.size = 0.5, key.position = "top",
                       max.exp = 3)
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
      
      load_File <- isolate(as.character(parseFilePaths(volumes, input$normal_obj_integrated)[4]))
      
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
             par(mfrow = c(1, 2))
             GenePlot(object = pbmc.combined, gene1 = "nUMI", gene2 = "nGene")
             GenePlot(object = pbmc.combined, gene1 = "nUMI", gene2 = "percent.mito")
           },
           "CCA_plot" = plot_grid(p2_integrated, p3_integrated),
           "MetageneBicorPlot" = p4_integrated,
           "DimHeatmap" = DimHeatmap(object = immune.combined, reduction.type = "cca", cells.use = input$ccs_cells_heatmap,
                                     dim.use = 1:input$ccs_dims_heatmap, do.balanced = TRUE),
           "CCA_aligned" = plot_grid(p5_integrated, p6_integrated),
           "TSNEPlot" = TSNEPlot(immune.combined, do.label = T, do.hover = FALSE, pt.size = 1.5),
           "VlnPlot_markers" = p7_integrated,
           "FeaturePlot" = FeaturePlot(object = immune.combined,
                                       features.plot = features_plot_integrated,
                                       cols.use = c("grey", "blue"),
                                       reduction.use = "tsne"),
           "RidgePlot" = RidgePlot(object = immune.combined,
                                   features.plot = features_plot_integrated, nCol = 2),
           "Features_conserved" = FeaturePlot(object = immune.combined,
                                       features.plot = features_plot_conserved,
                                       cols.use = c("grey", "blue"),
                                       reduction.use = "tsne"),
           "SplitDotPlotGG" = sdp,
           "Scatter_plot" = p8_integrated,
           "FeatureHeatmap" = FeatureHeatmap(immune.combined, features.plot = features_heatmap, 
                                             group.by = "stim", pt.size = 0.25, key.position = "top", 
                                             max.exp = 3)
    )
  })
  
  plotInput = function() {
    
    switch(input$psave_select_integrated,
           "QC" = p1_integrated,
           "GenePlot" = {
             par(mfrow = c(1, 2))
             GenePlot(object = pbmc.combined, gene1 = "nUMI", gene2 = "nGene")
             GenePlot(object = pbmc.combined, gene1 = "nUMI", gene2 = "percent.mito")
           },
           "CCA_plot" = plot_grid(p2_integrated, p3_integrated),
           "MetageneBicorPlot" = p4_integrated,
           "DimHeatmap" = DimHeatmap(object = immune.combined, reduction.type = "cca", cells.use = input$ccs_cells_heatmap,
                                     dim.use = 1:input$ccs_dims_heatmap, do.balanced = TRUE),
           "CCA_aligned" = plot_grid(p5_integrated, p6_integrated),
           "TSNEPlot" = TSNEPlot(immune.combined, do.label = T, do.hover = FALSE, pt.size = 1.5),
           "VlnPlot_markers" = p7_integrated,
           "FeaturePlot" = FeaturePlot(object = immune.combined,
                                       features.plot = features_plot_integrated,
                                       cols.use = c("grey", "blue"),
                                       reduction.use = "tsne"),
           "RidgePlot" = RidgePlot(object = immune.combined,
                                   features.plot = features_plot_integrated, nCol = 2),
           "Features_conserved" = FeaturePlot(object = immune.combined,
                                              features.plot = features_plot_conserved,
                                              cols.use = c("grey", "blue"),
                                              reduction.use = "tsne"),
           "SplitDotPlotGG" = sdp,
           "Scatter_plot" = p8_integrated,
           "FeatureHeatmap" = FeatureHeatmap(immune.combined, features.plot = features_heatmap, 
                                             group.by = "stim", pt.size = 0.25, key.position = "top", 
                                             max.exp = 3)
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
