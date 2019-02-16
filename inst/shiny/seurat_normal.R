

source("helpers.R")
# Data input and QC interface
seurat_normal_UI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "seurat_normal_tab",
          navbarPage(
            title = NULL,
            tabPanel("1. Data preprocessing",
                     fluidRow(

                         box(title = "1. Upload Data", width = 3, height = 350, solidHeader = TRUE, status = "info",
                           collapsible = T, id = "data_upload",
                           tags$style(appCSS),
                           radioButtons(ns("data_type"), NULL,
                                        c(
                                          "DGE Data (*.dge.*)" = "dge_file",
                                          "10X Data (*.mtx and *.tsv)" = "10X_file",
                                          "Example Data (2700 PBMCs)" = "demo_file"
                                        ), selected = "dge_file"),

                           conditionalPanel(condition  = "input.data_type == 'dge_file'",
                                            tags$b("Select DGE file:"),
                                            tags$p(),
                                            verbatimTextOutput(ns("input_dge_normal"), placeholder = TRUE),
                                            shinyFilesButton(ns("normal_dge"), "Select DGE file",
                                                             "Please select file", multiple = FALSE),
                                            ns = NS(id)
                           ),

                           conditionalPanel(condition  = "input.data_type == '10X_file'",
                                            tags$b("Select 10x files folder:"),
                                            tags$p(),
                                            verbatimTextOutput(ns("input_10x_normal"), placeholder = TRUE),
                                            shinyDirButton(ns("normal_10x"), "Select 10x file folder", "Please select a folder"),
                                            ns = NS(id)
                           )

                         ),

                         box(title = "2. Initialize Seurat object", width = 3, height = 350, solidHeader = TRUE, status = "warning",
                             collapsible = T, id = "object_init",

                             textInput(ns("project_name"), value = "Project1", label = "Project Name"),
                             numericInput(ns("min_cells"),
                                          label="Keep genes expressed at least in (?) cells", min = 1, max = Inf, value = 3),
                             numericInput(ns("min_genes"),
                                          label="Keep cells with at least (?) genes",min = 1, max = Inf, value = 200),
                             withBusyIndicatorUI(icon_name = "sinit",
                               actionButton(ns("seurat_init"),"Process raw data", style = "width: 80%")
                             )

                         ),

                         box(title = "3. Cell selection for further analysis", width = 3, height = 350, solidHeader = TRUE, status = "danger",
                             collapsible = T, id = "cell_qc",

                             selectizeInput(ns("subsetNames"), label="Filter Names",
                                            choices = c("nGene", "percent.mito"),
                                            selected = c("nGene", "percent.mito"),
                                            multiple=TRUE),
                             withBusyIndicatorUI(icon_name = "cfilter",
                               actionButton(ns("cell_filter"),"Select cells", style = "width: 80%")
                             )
                         ),

                         box(
                           title = "4. Normalizing the data", width = 3, height = 350, solidHeader = TRUE, status = "success",
                           selectizeInput(ns("norm_methods"), label="Normalization method",
                                          choices = c("LogNormalize"),
                                          selected = c("LogNormalize"),
                                          multiple=FALSE),
                           numericInput(ns("scale"),
                                        label="Scale factor", min = 1, max = Inf, value = 10000),

                           withBusyIndicatorUI(icon_name = "norm",
                             actionButton(ns("data_norm"),"Normalize data", style = "width: 80%")
                           )
                         )


                       ),

                     fluidRow(
                              box(
                                title = "3.1 VlnPlot",  width = 8, height = 500, status = "danger",

                                sidebarPanel(
                                  tableOutput(ns("nGene_range")),
                                  uiOutput(ns("gene_low_high")),
                                  hr(),
                                  tableOutput(ns("mito_range")),
                                  uiOutput(ns("mito_low_high"))
                                  #uiOutput(ns("nGene_range")),
                                  #uiOutput(ns("mito_range"))
                                ),
                                mainPanel(
                                  uiOutput(ns("filter_plots"))
                                )

                              ),
                              box(
                                title = "3.2 GenePlot",  width = 4, height = 500, status = "danger",
                                uiOutput(ns("gene_plots"))
                              )
                     )
            ), # Tab 1 END

            tabPanel('2. Detection of variable genes',
                     fluidRow(

                       box(title = "Variable genes detection", width = 4, height = 600, solidHeader = TRUE, status = "warning",
                           collapsible = T, id = "vgene_detect",

                           selectizeInput(ns("mean_fun"), label="Mean function",
                                          choices = c("ExpMean"),
                                          selected = c("ExpMean"),
                                          multiple = FALSE),

                           selectizeInput(ns("disper_fun"), label="Dispersion function",
                                          choices = c("LogVMR"),
                                          selected = c("LogVMR"),
                                          multiple = FALSE),

                           numericInput(ns("x_low"),
                                        label="X low_cutoff", min = 0, max = Inf, value = 0.0125, step = 0.0001),
                           numericInput(ns("x_high"),
                                        label="X high_cutoff", min = 0, max = Inf, value = 3, step = 0.1),
                           numericInput(ns("y_cut"),
                                        label="Y cutoff", min = 0, max = Inf, value = 0.5, step = 0.1),

                           selectizeInput(ns("vars_regress"), label="Vars.to.regress",
                                          choices = c("nUMI", "percent.mito"),
                                          selected = c("nUMI", "percent.mito"),
                                          multiple = TRUE),

                           withBusyIndicatorUI(icon_name = "vgenes",
                                               actionButton(ns("find_vgenes"),
                                                            "Detect variable genes",
                                                            style = "width: 90%")
                           ),
                           br(),
                           htmlOutput(ns("selected_vgenes"))
                           #textOutput(ns("selected_vgenes"))

                       ),
                       box(
                         title = "Variable gene plot",  width = 7, height = 600, status = "primary",
                         uiOutput(ns("vgenes_plots"))
                       )
                     )
            ), # Tab 2 END

            tabPanel('3. PCA-1',
                     fluidRow(

                       box(title = "PCA-1", width = 6, height = 850, solidHeader = TRUE, status = "warning",
                           collapsible = T, id = "pca1",
                           splitLayout(
                           numericInput(ns("pca1_pcs1_print"),
                                        label="Print from PC: ", min = 1, max = Inf, value = 1),
                           numericInput(ns("pca1_pcs2_print"),
                                        label="to PC: ", min = 1, max = Inf, value = 5),
                           numericInput(ns("pca1_genes_print"),
                                        label="Number of genes to print", min = 1, max = Inf, value = 5)
                           ),
                           splitLayout(
                           numericInput(ns("pca1_pcs1_plot"),
                                        label="Plot from PC:", min = 1, max = Inf, value = 1),
                           numericInput(ns("pca1_pcs2_plot"),
                                        label="To PC:", min = 2, max = Inf, value = 2)
                           ),

                           withBusyIndicatorUI(icon_name = "pca1",
                                               actionButton(ns("viz_pca1"),
                                                            "Run PCA",
                                                            style = "width: 90%")
                           ),
                           br(),
                           uiOutput(ns("pca1_table"))
                           #DTOutput(ns("pca1_genes"))
                       ),

                       box(
                         title = "Variable gene plot",  width = 6, height = 850, status = "primary",
                         uiOutput(ns("pca1_plot"))
                       )

                     )#
            ), # Tab 3 END

            tabPanel('4. PCA-2',
                     fluidRow(
                       column(width = 5,
                         box(title = "PCA-2", width = NULL, solidHeader = TRUE, status = "warning",
                             collapsible = T, id = "pca2",
                             splitLayout(
                               numericInput(ns("pca2_pcs1_plot"),
                                            label="Plot PC: ", min = 1, max = Inf, value = 1),
                               numericInput(ns("pca2_pcs2_plot"),
                                            label="and PC: ", min = 1, max = Inf, value = 2)
                             ),
                             splitLayout(
                               numericInput(ns("pca2_pcs1_heatmap"),
                                            label="Heatmap uses from PC:", min = 1, max = Inf, value = 1),
                               numericInput(ns("pca2_pcs2_heatmap"),
                                            label="To PC:", min = 1, max = Inf, value = 1),
                               numericInput(ns("pca2_cells_heatmap"),
                                            label="Cells to use in heatmap:", min = 1, max = Inf, value = 500)
                             ),

                             withBusyIndicatorUI(icon_name = "pca2",
                                                 actionButton(ns("viz_pca2"),
                                                              "PCA Plots",
                                                              style = "width: 90%")
                             )),

                             box(title = "PCA plot", width = NULL, height = 580, status = "success",
                               uiOutput(ns("pca2_plot"))
                             )
                       ), #

                       column(width = 7,
                         box(
                           title = "PCHeatmap",  width = NULL, height = 850, status = "primary",
                           uiOutput(ns("pca2_heatmap"))
                         )
                       )#
                    )
            ), # Tab 4 END

            tabPanel('5. Determine significant PCs',
                     fluidRow(
                       column(width = 5,
                              box(title = "Significant principal components", width = NULL, solidHeader = TRUE, status = "warning",
                                  collapsible = T, id = "SPCs_parameters",
                                  splitLayout(
                                    checkboxInput(ns("plot_js"), "Generate JackStrawPlot ?", TRUE),
                                    numericInput(ns("num_replicate"),
                                                 label="Num.replicate: ", min = 1, max = Inf, value = 100)
                                  ),
                                  splitLayout(
                                    numericInput(ns("JS_pcs1"),
                                                 label="JackStrawPlot PC from: ", min = 1, max = Inf, value = 1),
                                    numericInput(ns("JS_pcs2"),
                                                 label="to PC: ", min = 1, max = Inf, value = 12)
                                  ),

                                  withBusyIndicatorUI(icon_name = "spcs",
                                                      actionButton(ns("viz_spcs"),
                                                                   "Significant PCs",
                                                                   style = "width: 90%")
                                  )),

                              box(title = "PC Elbow Plot", width = NULL, height = 580, status = "success",
                                  uiOutput(ns("pc_elbow_plot"))
                              )
                       ), #

                       column(width = 7,
                              box(
                                title = "JackStraw Plot",  width = NULL, height = 850, status = "primary",
                                uiOutput(ns("jackstraw_plot"))
                              )
                       )#
                     )
                     ), # Tab 5 END

            tabPanel('6. Cluster the cells',
                     fluidRow(
                       column(width = 5,
                              box(title = "Cluster the cells", width = NULL, solidHeader = TRUE, status = "warning",
                                  collapsible = T, id = "cluster_parameters",

                                  splitLayout(
                                    numericInput(ns("cluster_dim1"),
                                                 label="Find cluster from dimension: ", min = 1, max = Inf, value = 1),
                                    numericInput(ns("cluster_dim2"),
                                                 label="to dimension: ", min = 1, max = Inf, value = 10)
                                  ),
                                  splitLayout(
                                    numericInput(ns("res"),
                                                 label="Resolution (TSNE)", min = 0, max = Inf, value = 0.6, step = 0.01),
                                    numericInput(ns("perp"),
                                                 label="Perplexity (TSNE): ", min = 1, max = Inf, value = 30)
                                  ),

                                  withBusyIndicatorUI(icon_name = "clusters",
                                                      actionButton(ns("viz_clusters"),
                                                                   "Run TSNE",
                                                                   style = "width: 90%")
                                  )),

                              box(title = "Assigning cell type identity to clusters", width = NULL, solidHeader = TRUE,
                                  height = 580, status = "success",

                                  uiOutput(ns("old_name_clusters")),
                                  textInput(ns("new_name_clusters"), "Assign new names (use ',' to separate): "),

                                  withBusyIndicatorUI(icon_name = "clusters_names",
                                                      actionButton(ns("assign_names"),
                                                                   "Assign new names",
                                                                   style = "width: 90%")
                                  )
                              )
                       ), #

                       column(width = 7,
                              box(
                                title = "TSNE Plot",  width = NULL, height = 850, status = "primary",
                                uiOutput(ns("tsne_plot1"))
                              )
                       )#
                     )
                     ), # Tab 6 END

            tabPanel('7. Differentially expressed genes',

                     fluidRow(

                       box(title = "Finding differentially expressed genes", width = 5, height = 950, solidHeader = TRUE, status = "warning",
                           collapsible = T, id = "find_diff_genes",

                           radioButtons(ns("diff_genes_radio"), label = NULL,
                                        choices = list("Single cluster" = "single",
                                                       "Multiple clusters" = "multiple",
                                                       "All clusters" = "all"),
                                        selected = "single", inline = TRUE),

                           uiOutput(ns("id1_diff_genes")),
                           uiOutput(ns("id2_diff_genes")),
                           #uiOutput(ns("id2"))

                           splitLayout(
                             numericInput(ns("min_pct"),
                                          label="Minimum percentage:", min = 0, max = 1, value = 0.25, step = 0.01),
                             numericInput(ns("top_markers"),
                                          label="Top markers to show:", min = 1, max = Inf, value = 10)
                           ),

                           withBusyIndicatorUI(icon_name = "diff_genes",
                                               actionButton(ns("viz_diff_genes"),
                                                            "Cluster biomarkers",
                                                            style = "width: 90%")
                           ),
                           br(),
                           uiOutput(ns("diff_genes_table"))
                           #DTOutput(ns("pca1_genes"))
                       ),

                       box(
                         title = "Cluster biomarkers visulization",  width = 7, height = 950, solidHeader = TRUE, status = "primary",
                         tabsetPanel(type = "tabs",
                                     tabPanel("Expression heatmap",
                                              br(),
                                              uiOutput(ns("diff_genes_heatmap"))),
                                     tabPanel("Expression probability distributions",
                                              br(),
                                              #verbatimTextOutput(ns('x4')),
                                              uiOutput(ns("diff_genes_vlnplot"))),
                                     tabPanel("Feature plot",
                                              br(),
                                              uiOutput(ns("diff_genes_featureplot"))),
                                     tabPanel("Ridge plot",
                                              br(),
                                              uiOutput(ns("diff_genes_ridgeplot")))
                         )
                       )
                     )#
                     ), # Tab 7 END
            
            tabPanel('* Save and load Seurat object',
                     fluidRow(
                       
                       box(
                         title = "Save Seurat object", width = 4, solidHeader = TRUE, status = "primary",
                         textInput(ns("saveName"), "Sample name:", width = "100%"),
                         tags$b("Save to directory:"),
                         tags$p(),
                         shinyDirButton(ns("dir_save"), "Select a folder to save ...", "Please select a folder"),
                         tags$p(),
                         verbatimTextOutput(ns("save_dir"), placeholder = TRUE),
                         tags$p(),
                         
                         withBusyIndicatorUI(icon_name = "save",
                                             actionButton(ns("do_save"),
                                                          "Save",
                                                          style = "width: 90%")
                         )
                       ),
                       
                       box(
                         title = "Load Seurat object", width = 4, solidHeader = TRUE, status = "warning",
    
                         tags$b("Load from directory:"),
                         tags$p(),
                         shinyFilesButton(ns("normal_obj"), "Select Seurat object file", "Please select a file", multiple = FALSE),
                         tags$p(),
                         verbatimTextOutput(ns("input_obj_normal"), placeholder = TRUE),
                         tags$p(),
                         
                         withBusyIndicatorUI(icon_name = "read",
                                             actionButton(ns("do_load"),
                                                          "Load",
                                                          style = "width: 90%")
                         )
                       ),
                       
                       box(title = "Add samples into existing Seurat object", width = 4, solidHeader = TRUE, status = "success",
                           
                           textInput(ns("addName"), "New sample name:", width = "100%"),
                           radioButtons(ns("data_type_add"), NULL,
                                        c(
                                          "DGE Data (*.dge.*)" = "dge_file_add",
                                          "10X Data (*.mtx and *.tsv)" = "10X_file_add"
                                        ), selected = "dge_file_add"),
                           
                           conditionalPanel(condition  = "input.data_type_add == 'dge_file_add'",
                                            tags$b("Select DGE file:"),
                                            tags$p(),
                                            shinyFilesButton(ns("normal_dge_add"), "Select DGE file",
                                                             "Please select file", multiple = FALSE),
                                            tags$p(),
                                            verbatimTextOutput(ns("dge_normal_add"), placeholder = TRUE),
                                            ns = NS(id)
                           ),
                           
                           conditionalPanel(condition  = "input.data_type_add == '10X_file_add'",
                                            tags$b("Select 10x files folder:"),
                                            tags$p(),
                                            shinyDirButton(ns("normal_10x_add"), "Select 10x file folder", "Please select a folder"),
                                            tags$p(),
                                            verbatimTextOutput(ns("input_10x_add"), placeholder = TRUE),
                                            ns = NS(id)
                           ),
                           
                           tags$p(),
                           
                           withBusyIndicatorUI(icon_name = "add",
                                               actionButton(ns("do_add"),
                                                            "Add new sample",
                                                            style = "width: 90%")
                           )
                           
                       )
                       
                     )#
            ),# Tab 8 END
            
            tabPanel('* Save figures',
                     fluidRow(
                       
                       box(
                         title = "Save figure", width = 3, solidHeader = TRUE, status = "primary",
                         
                         selectizeInput(ns("psave_select"), label="Choose a figure to save",
                                        choices = c("QC", "GenePlot", "Dispersion_plot", "VizPCA", "PCAPlot",
                                                    "PCHeatmap", "JackStrawPlot", "PCElbowPlot", "TSNEPlot",
                                                    "VlnPlot_markers", "FeaturePlot", "RidgePlot", "DoHeatmap"),
                                        selected = c("QC"),
                                        multiple=FALSE
                                        ),

                         numericInput(ns("pWidth"), "Width (px)", value = 1200),   # Width is divided by 2 for display
                         numericInput(ns("pHeight"), "Height (px)", value = 1200), # Height is divided by 2 for display
                         numericInput(ns("pRes"), "Resolution (ppi)", value = 300),
                         numericInput(ns("pPsize"), "Pointsize", value = 12),
                         
                         downloadButton(ns("do_psave"), "Save figure", style = "width: 100%")
                       ),
                       
                       box(
                         title = "Figure to save", width = 9, status = "success",
                         
                         uiOutput(ns("figure_save"))
                       )
                       
                     )#
            )# Tab 9 END
          )
  )

}

seurat_normal_Server <- function(input, output, session, fileRoot = NULL) {

  ns <- session$ns

  volumes <- (c(Home = fs::path_home(), getVolumes()()))

  shinyFileChoose(input, "normal_dge", roots = volumes, session = session)
  shinyDirChoose(input, "normal_10x", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  shinyDirChoose(input, "dir_save", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyFileChoose(input, "normal_obj", roots = volumes, session = session, filetypes=c('', 'rds'))
  
  shinyFileChoose(input, "normal_dge_add", roots = volumes, session = session)
  shinyDirChoose(input, "normal_10x_add", roots = volumes, session = session, restrictions = system.file(package = "base"))

  output$input_dge_normal <- renderPrint({
    as.character(parseFilePaths(volumes, input$normal_dge)[4])
  })

  output$input_10x_normal <- renderPrint({
    parseDirPath(volumes, input$normal_10x)
  })
  
  output$save_dir <- renderPrint({
    parseDirPath(volumes, input$dir_save)
  })
  
  output$input_obj_normal <- renderPrint({
    as.character(parseFilePaths(volumes, input$normal_obj)[4])
  })
  
  output$dge_normal_add <- renderPrint({
    as.character(parseFilePaths(volumes, input$normal_dge_add)[4])
  })
  
  output$input_10x_add <- renderPrint({
    parseDirPath(volumes, input$normal_10x_add)
  })

  # Object initialization
  observeEvent(input$seurat_init, {

    withProgress(message = "Initializing Seurat Object ...",{
    shinyjs::show("load_sinit")

    if(input$data_type == "dge_file"){
      File_dge_normal <- isolate(as.character(parseFilePaths(volumes, input$normal_dge)[4]))
      rawdata <- read.table(File_dge_normal, header = T, row.names = 1, stringsAsFactors = F)
    }
    else if(input$data_type == "10X_file"){
      Path_10x_normal <- isolate(as.character(parseDirPath(volumes, input$normal_10x)))
      rawdata <- Read10X(data.dir = Path_10x_normal)
    }
    else{
      rawdata <- Read10X(data.dir = "useful/normal_analysis/")
    }

    setProgress(value = 0.5)

    pbmc <<- CreateSeuratObject(raw.data = rawdata, min.cells = input$min_cells, min.genes = input$min_genes,
                                project = input$project_name)

    setProgress(value = 0.7)

    mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE, ignore.case = TRUE)
    percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
    pbmc <<- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

    setProgress(value = 1)

    shinyjs::hide("load_sinit")
    shinyjs::show("check_sinit")
    })
  }

  )# Object initialization END

  #Cell selection
  observeEvent(input$cell_filter, {

    shinyjs::show("load_cfilter")

    output$filter_plots <- renderUI({
      splitLayout(cellWidths = c("33%", "33%", "34%"),
                  withSpinner(plotOutput(ns("nGenePlot")), type = getOption("spinner.type", default = 4)),
                  withSpinner(plotOutput(ns("nUMIPlot")), type = getOption("spinner.type", default = 4)),
                  withSpinner(plotOutput(ns("mitoPlot")), type = getOption("spinner.type", default = 4))
      )
    }
    )

    output$gene_plots <- renderUI({
      splitLayout(cellWidths = c("49%", "49%"),
                  withSpinner(plotOutput(ns("gene_nUMI")), type = getOption("spinner.type", default = 4)),
                  withSpinner(plotOutput(ns("mito_nUMI")), type = getOption("spinner.type", default = 4))
      )
    }
    )

    min_gene <- min(pbmc@meta.data$nGene)
    max_gene <- max(pbmc@meta.data$nGene)

    min_mito <- round(min(pbmc@meta.data$percent.mito), 2)
    max_mito <- round(max(pbmc@meta.data$percent.mito), 2)

    output$nGenePlot <- renderPlot({

      if("nGene" %in% input$subsetNames){

        VlnPlot(object = pbmc, features.plot = c("nGene"), nCol = 1) +
        geom_hline(yintercept = input$gene_low, color = 'darkgreen', linetype = "dashed", size = 1.5) +
        geom_text(data = data.frame(x=1.2,y=input$gene_low), aes(x, y), label="low", vjust=1.5, hjust=0,color = "darkgreen",size = 5,fontface = "bold") +
        geom_hline(yintercept = input$gene_high, color = 'blue', linetype = "dashed", size = 1.5) +
        geom_text(data = data.frame(x=1.2,y=input$gene_high), aes(x, y), label="high", vjust=-1, hjust=0,color = "blue",size = 5,fontface = "bold") +
        scale_y_continuous(limits=c(0, 1.1*max_gene))
      }
      else{
        VlnPlot(object = pbmc, features.plot = c("nGene"), nCol = 1)
      }
    })
    output$nUMIPlot <- renderPlot({
      VlnPlot(object = pbmc, features.plot = c("nUMI"), nCol = 1)

    })
    output$mitoPlot <- renderPlot({

      if("percent.mito" %in% input$subsetNames){

        VlnPlot(object = pbmc, features.plot = c("percent.mito"), nCol = 1) +
        geom_hline(yintercept = input$mito_low, color = 'darkgreen', linetype = "dashed", size = 1.5) +
        geom_text(data = data.frame(x=1.2,y=input$mito_low), aes(x, y), label="low", vjust=1.5, hjust=0,color = "darkgreen",size = 5,fontface = "bold") +
        geom_hline(yintercept = input$mito_high, color = 'blue', linetype = "dashed", size = 1.5) +
        geom_text(data = data.frame(x=1.2,y=input$mito_high), aes(x, y), label="high", vjust=-1, hjust=0,color = "blue",size = 5,fontface = "bold") +
        scale_y_continuous(limits=c(min_mito-0.01, max_mito+0.01))
      }
      else{
        VlnPlot(object = pbmc, features.plot = c("percent.mito"), nCol = 1)
      }
    })
    
    p1 <<- VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

    median_gene <- median(pbmc@meta.data$nGene)

    min_UMI <- min(pbmc@meta.data$nUMI)
    max_UMI <- max(pbmc@meta.data$nUMI)
    median_UMI <- median(pbmc@meta.data$nUMI)

    nGene_summary <- data.frame(c(min_gene, median_gene, max_gene),
                                c(min_UMI, median_UMI, max_UMI),
                                row.names = c("Min", "Median", "Max"))
    nGene_summary <- t(nGene_summary)
    row.names(nGene_summary) <- c("nGene", "nUMI")

    #output$nGene_range <- renderUI({
    #  sliderInput(ns("nGene_slide"), "Choose nGene Range:",
    #              min = min_gene, max = max_gene, value = c(min_gene, max_gene))
    #}
    #)

    output$nGene_range <- renderTable(spacing = "xs", align = "c", rownames = TRUE,
                                     striped = TRUE, hover = TRUE, bordered = TRUE, digits = 0,
                                     { nGene_summary }
                                     )

    output$gene_low_high <- renderUI({
      splitLayout(
        numericInput(ns("gene_low"), "nGene_lowest", value = min_gene),
        numericInput(ns("gene_high"), "nGene_highest", value = max_gene)
      )
    }
    )

    median_mito <- median(pbmc@meta.data$percent.mito)
    mito_summary <- data.frame(c(min_mito, median_mito, max_mito),
                                row.names = c("Min", "Median", "Max"))
    mito_summary <- t(mito_summary )
    row.names(mito_summary) <- "p.mito"

    #output$mito_range <- renderUI({
    #  sliderInput(ns("mito_slide"), "Choose percent.mito Range:",
    #              min = min_mito, max = max_mito, value = c(min_mito, max_mito))
    #}
    #)

    output$mito_range <- renderTable(spacing = "xs", align = "c", rownames = TRUE,
                                      striped = TRUE, hover = TRUE, bordered = TRUE,
                                      { mito_summary }
    )

    output$mito_low_high <- renderUI({
      splitLayout(
        numericInput(ns("mito_low"), "mito_lowest", value = min_mito - 0.01, step = 0.01),
        numericInput(ns("mito_high"), "mito_highest", value = max_mito, step = 0.01)
      )
    }
    )

    output$mito_nUMI <- renderPlot({

      if("percent.mito" %in% input$subsetNames){

        GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
        abline(h = input$mito_low, col="darkgreen", lwd=2, lty="dashed")
        abline(h = input$mito_high, col="blue", lwd=2, lty="dashed")
      }
      else{
        GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
      }
    })

    output$gene_nUMI <- renderPlot({

      if("nGene" %in% input$subsetNames){

        GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
        abline(h = input$gene_low, col="darkgreen", lwd=2, lty="dashed")
        abline(h = input$gene_high, col="blue", lwd=2, lty="dashed")
      }
      else{
        GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
      }
    })
    
    pbmc_o <<- pbmc
    
    shinyjs::hide("load_cfilter")
    shinyjs::show("check_cfilter")

  }

  )# Cell selection END


  # Normalization
  observeEvent(input$data_norm, {

    withProgress(message = "Performing log-normalization ...",{

      shinyjs::show("load_norm")
      setProgress(value = 0.35)

      if("nGene" %in% input$subsetNames && "percent.mito" %in% input$subsetNames){

        pbmc <<- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"),
                             low.thresholds = c(input$gene_low, input$mito_low),
                             high.thresholds = c(input$gene_high, input$mito_high))
      }
      else if("nGene" %in% input$subsetNames){

        pbmc <<- FilterCells(object = pbmc, subset.names = "nGene",
                             low.thresholds = input$gene_low,
                             high.thresholds = input$gene_high)
      }
      else{

        pbmc <<- FilterCells(object = pbmc, subset.names = "percent.mito",
                             low.thresholds = input$mito_low,
                             high.thresholds = input$mito_high)
      }

      pbmc <<- NormalizeData(object = pbmc, normalization.method = input$norm_methods,
                             scale.factor = input$scale)

      setProgress(value = 1)
      shinyjs::hide("load_norm")
      shinyjs::show("check_norm")
    })

  }

  )# Normalization END

  # Detect variable genes
  observeEvent(input$find_vgenes, {

    shinyjs::show("load_vgenes")

    output$vgenes_plots <- renderUI({

      withSpinner(plotOutput(ns("vgenePlot"), height = "500px"), type = getOption("spinner.type", default = 4))

    }
    )

    withProgress(message = "Detecting variable genes ...",{

      setProgress(value = 0.35)

      pbmc <<- FindVariableGenes(object = pbmc, mean.function = input$mean_fun, dispersion.function = input$disper_fun,
                                 x.low.cutoff = input$x_low, x.high.cutoff = input$x_high, y.cutoff = input$y_cut)

      pbmc <<- ScaleData(object = pbmc, vars.to.regress = input$vars_regress)

      setProgress(value = 1)

      output$selected_vgenes <- renderUI({

        HTML(paste("<b> ", length(x = pbmc@var.genes), "variable genes have been selected.</b>"))

      })

      output$vgenePlot <- renderPlot({
        VariableGenePlot(object = pbmc)
      })

      shinyjs::hide("load_vgenes")
      shinyjs::show("check_vgenes")


    })

  }

  )# Detect variable genes END

  # PCA1
  observeEvent(input$viz_pca1, {

    output$pca1_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_pca1"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    output$pca1_table <- renderUI({
      withSpinner(DTOutput(ns("pca1_genes")), type = getOption("spinner.type", default = 4))
    }
    )

    withProgress(message = "Running PCA ...",{

      shinyjs::show("load_pca1")
      setProgress(value = 0.35)

      pbmc <<- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = FALSE)
                      #pcs.print = input$pca1_pcs1_print:input$pca1_pcs2_print,
                      #genes.print = input$pca1_genes_print)

      output$plot_pca1 <- renderPlot({
        VizPCA(object = pbmc, pcs.use = input$pca1_pcs1_plot:input$pca1_pcs2_plot)
      })

      setProgress(value = 1)
      shinyjs::hide("load_pca1")
      shinyjs::show("check_pca1")
    })

    pca1Genes <- capture.output(PrintPCA(object = pbmc, pcs.print = input$pca1_pcs1_print:input$pca1_pcs2_print,
                                         genes.print = input$pca1_genes_print, use.full = FALSE), type = "output")
    pca1Genes <- gsub("\\D[[1-9]]", "", pca1Genes)
    pca1Genes <- matrix(pca1Genes) #, nrow = input$pca1_genes_print+1)
    pca1Genes <- pca1Genes[!apply(pca1Genes == " \"\"", 1, all),]
    pca1Genes <- unlist(strsplit(pca1Genes, " "))
    pca1Genes <- matrix(pca1Genes)
    pca1Genes <- matrix(pca1Genes[!apply(pca1Genes == "", 1, all),], nrow = 2*input$pca1_genes_print+1)
    pca1Genes <- data.frame(pca1Genes)
    pca1Genes <- as.data.frame(sapply(pca1Genes, function(x) gsub("\"", "", x)))
    colnames(pca1Genes) <- as.character(unlist(pca1Genes[1,]))
    pca1Genes <- pca1Genes[-1,]

    output$pca1_genes <- renderDT(datatable(pca1Genes, rownames = FALSE, selection = "none"))

  }

  )# PCA1 END


  # PCA2
  observeEvent(input$viz_pca2, {

    output$pca2_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_pca2"), height = "500px"), type = getOption("spinner.type", default = 4))
    }
    )

    output$pca2_heatmap <- renderUI({
      withSpinner(plotOutput(ns("plot_heatmap"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    withProgress(message = "Plotting PCA ...",{

      shinyjs::show("load_pca2")
      setProgress(value = 0.35)

      pbmc <<- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = FALSE)

      output$plot_pca2 <- renderPlot({
        p5 <<- PCAPlot(object = pbmc, dim.1 = input$pca2_pcs1_plot, dim.2 = input$pca2_pcs2_plot)
        p5
      })

      pbmc <<- ProjectPCA(object = pbmc, do.print = FALSE)

      output$plot_heatmap <- renderPlot({
        PCHeatmap(object = pbmc, pc.use = input$pca2_pcs1_heatmap:input$pca2_pcs2_heatmap,
                  cells.use = input$pca2_cells_heatmap,
                  do.balanced = TRUE, label.columns = FALSE)
      })

      setProgress(value = 1)
      shinyjs::hide("load_pca2")
      shinyjs::show("check_pca2")
    })

  }

  )# PCA2 END


  # SPCs
  observeEvent(input$viz_spcs, {

    output$pc_elbow_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_elbow"), height = "500px"), type = getOption("spinner.type", default = 4))
    }
    )

    if(input$plot_js == TRUE){
      output$jackstraw_plot <- renderUI({
        withSpinner(plotOutput(ns("plot_jackstraw"), height = "750px"), type = getOption("spinner.type", default = 4))
      }
      )
    }

    withProgress(message = "Plotting significant PCs ...",{

      shinyjs::show("load_spcs")
      setProgress(value = 0.35)

      if(input$plot_js == TRUE){
        pbmc <- JackStraw(object = pbmc, num.replicate = input$num_replicate, display.progress = FALSE)
        output$plot_jackstraw <- renderPlot({
          p7 <<- JackStrawPlot(object = pbmc, PCs = input$JS_pcs1:input$JS_pcs2)
          p7
        })
      }

      output$plot_elbow <- renderPlot({
        p8 <<- PCElbowPlot(object = pbmc)
        p8
      })

      setProgress(value = 1)
      shinyjs::hide("load_spcs")
      shinyjs::show("check_spcs")
    })

  }

  )# SPCs END


  # Clusters
  observeEvent(input$viz_clusters, {

    output$tsne_plot1 <- renderUI({
      withSpinner(plotOutput(ns("plot_tsne1"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    withProgress(message = "Plotting TSNE ...",{

      shinyjs::show("load_clusters")
      setProgress(value = 0.35)

      pbmc <<- FindClusters(object = pbmc, reduction.type = "pca", dims.use = input$cluster_dim1:input$cluster_dim2,
                           resolution = input$res, print.output = 0, save.SNN = TRUE)

      pbmc <<- RunTSNE(object = pbmc, dims.use = input$cluster_dim1:input$cluster_dim2,
                       do.fast = TRUE, perplexity = input$perp)

      output$plot_tsne1 <- renderPlot({
        TSNEPlot(object = pbmc, do.hover = FALSE, pt.size = 1.5)

      })

      setProgress(value = 1)
      shinyjs::hide("load_clusters")
      shinyjs::show("check_clusters")
    })

    output$old_name_clusters <- renderUI({

      selectizeInput(ns("clusters_old_name"), label="Current names for each cluster: ",
                     choices = levels(pbmc@ident),
                     selected = levels(pbmc@ident)[1],
                     multiple = TRUE)
    }
    )

    output$id1_diff_genes <- renderUI({

      selectizeInput(ns("diff_genes_id1"), label="Find all markers of cluster:",
                     choices = levels(pbmc@ident),
                     selected = levels(pbmc@ident)[1],
                     multiple = FALSE)
    }
    )

    output$id2_diff_genes <- renderUI({

      selectizeInput(ns("diff_genes_id2"), label="Distinguishing from clusters:",
                     choices = levels(pbmc@ident)[levels(pbmc@ident) != input$diff_genes_id1],
                     selected = levels(pbmc@ident)[levels(pbmc@ident) != input$diff_genes_id1][1],
                     multiple = TRUE)
    }
    )

  }
  )# Clusters END

  # Assign new names
  observeEvent(input$assign_names, {

    output$tsne_plot1 <- renderUI({
      withSpinner(plotOutput(ns("plot_tsne1"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    withProgress(message = "Plotting TSNE ...",{

      shinyjs::show("load_clusters_names")
      setProgress(value = 0.35)

      current.cluster.ids <- input$clusters_old_name
      new.cluster.ids <- input$new_name_clusters
      new.cluster.ids <- unlist(strsplit(new.cluster.ids, split = ",|, "))

      pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)

      pbmc <<- pbmc

      output$plot_tsne1 <- renderPlot({

        TSNEPlot(object = pbmc, do.hover = FALSE, pt.size = 1.5)

      })

      setProgress(value = 1)
      shinyjs::hide("load_clusters_names")
      shinyjs::show("check_clusters_names")
    })

    output$old_name_clusters <- renderUI({

      selectizeInput(ns("clusters_old_name"), label="Select names to change: ",
                     choices = levels(pbmc@ident),
                     selected = levels(pbmc@ident)[1],
                     multiple = TRUE)
    }
    )

    output$id1_diff_genes <- renderUI({

      selectizeInput(ns("diff_genes_id1"), label="Find all markers of cluster:",
                     choices = levels(pbmc@ident),
                     selected = levels(pbmc@ident)[1],
                     multiple = FALSE)
    }
    )

    output$id2_diff_genes <- renderUI({

      selectizeInput(ns("diff_genes_id2"), label="Distinguishing from clusters:",
                     choices = levels(pbmc@ident)[levels(pbmc@ident) != input$diff_genes_id1],
                     selected = levels(pbmc@ident)[levels(pbmc@ident) != input$diff_genes_id1][1],
                     multiple = TRUE)
    }
    )

  }
  )# Assign new names END

  # Cluster biomarkers
  observeEvent(input$viz_diff_genes, {



    output$diff_genes_table <- renderUI({
      withSpinner(DTOutput(ns("table_diff_genes")), type = getOption("spinner.type", default = 4))
    }
    )

    output$diff_genes_vlnplot <- renderUI({
      withSpinner(plotOutput(ns("vlnplot_diff_genes"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    output$diff_genes_featureplot <- renderUI({
      withSpinner(plotOutput(ns("featureplot_diff_genes"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )

    output$diff_genes_ridgeplot <- renderUI({
      withSpinner(plotOutput(ns("ridgeplot_diff_genes"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )


    withProgress(message = "Clustering biomarkers ...",{

      shinyjs::show("load_diff_genes")
      setProgress(value = 0.35)

      if(input$diff_genes_radio == "single"){
        cluster.markers <- FindMarkers(object = pbmc, ident.1 = input$diff_genes_id1, min.pct = 0.25)
        is.num <- sapply(cluster.markers, is.numeric)
        cluster.markers[is.num] <- lapply(cluster.markers[is.num], round, 3)

        output$table_diff_genes <- renderDT(datatable(cluster.markers, rownames = TRUE,
                                                      selection = list(selected = 1, target = 'row')
                                                      ))


      }
      else if(input$diff_genes_radio == "multiple"){
        cluster.markers <- FindMarkers(object = pbmc, ident.1 = input$diff_genes_id1, ident.2 = input$diff_genes_id2,
                                       min.pct = 0.25)
        is.num <- sapply(cluster.markers, is.numeric)
        cluster.markers[is.num] <- lapply(cluster.markers[is.num], round, 3)

        output$table_diff_genes <- renderDT(datatable(cluster.markers, rownames = TRUE,
                                                      selection = list(selected = 1, target = 'row')
                                                      ))
      }
      else{
        cluster.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25,
                                          thresh.use = 0.25)

        cluster.markers <- cluster.markers %>% group_by(cluster) %>% top_n(input$top_markers, avg_logFC)

        output$diff_genes_heatmap <- renderUI({
          withSpinner(plotOutput(ns("heatmap_diff_genes"), height = "750px"), type = getOption("spinner.type", default = 4))
        }
        )

        output$heatmap_diff_genes <- renderPlot({
          p13 <<- DoHeatmap(object = pbmc, genes.use = cluster.markers$gene, slim.col.label = TRUE, remove.key = TRUE)
          p13
        })

        is.num <- sapply(cluster.markers, is.numeric)
        cluster.markers[is.num] <- lapply(cluster.markers[is.num], round, 3)

        output$table_diff_genes <- renderDT(datatable(cluster.markers, rownames = FALSE,
                                                      selection =  list(selected = 1, target = 'row')
                                                      ))


      }

      setProgress(value = 1)
      shinyjs::hide("load_diff_genes")
      shinyjs::show("check_diff_genes")
    })

    if(input$diff_genes_radio == "all"){

      output$vlnplot_diff_genes <- renderPlot({

        validate(need(length(input$table_diff_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        p10 <<- VlnPlot(object = pbmc, x.lab.rot = TRUE,
                features.plot = cluster.markers$gene[input$table_diff_genes_rows_selected])
        p10
        })

        output$featureplot_diff_genes <- renderPlot({

          features_plot <<- cluster.markers$gene[input$table_diff_genes_rows_selected]
          
          validate(need(length(input$table_diff_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
          FeaturePlot(object = pbmc,
                      features.plot = cluster.markers$gene[input$table_diff_genes_rows_selected],
                      cols.use = c("grey", "blue"),
                      reduction.use = "tsne")
          
        })

        output$ridgeplot_diff_genes <- renderPlot({

          features_plot <<- cluster.markers$gene[input$table_diff_genes_rows_selected]
          
          validate(need(length(input$table_diff_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
          RidgePlot(object = pbmc,
                      features.plot = cluster.markers$gene[input$table_diff_genes_rows_selected], nCol = 2)
          
        })


    }
    else{

      output$vlnplot_diff_genes <- renderPlot({

        validate(need(length(input$table_diff_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        p10 <<- VlnPlot(object = pbmc, x.lab.rot = TRUE,
                features.plot = rownames(cluster.markers)[input$table_diff_genes_rows_selected])
        p10
      })

      output$featureplot_diff_genes <- renderPlot({

        validate(need(length(input$table_diff_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        FeaturePlot(object = pbmc,
                    features.plot = rownames(cluster.markers)[input$table_diff_genes_rows_selected],
                    cols.use = c("grey", "blue"),
                    reduction.use = "tsne")
        
        features_plot <<- rownames(cluster.markers)[input$table_diff_genes_rows_selected]
      })

      output$ridgeplot_diff_genes <- renderPlot({

        features_plot <<- rownames(cluster.markers)[input$table_diff_genes_rows_selected]
        
        validate(need(length(input$table_diff_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        RidgePlot(object = pbmc,
                  features.plot = rownames(cluster.markers)[input$table_diff_genes_rows_selected], nCol = 2)
      })

    }

  }
  )# Cluster biomarkers END
  
  
  # Save
  observeEvent(input$do_save, {

    
    withProgress(message = "Saving ...",{
      
      shinyjs::show("load_save")
      setProgress(value = 0.35)
      
      save_Path <- isolate(as.character(parseDirPath(volumes, input$dir_save)))
      
      saveRDS(pbmc, file = file.path(save_Path, paste0(input$saveName, ".rds"))) 
      
      setProgress(value = 1)
      shinyjs::hide("load_save")
      shinyjs::show("check_save")
    })
    
  }
  
  )# Save END
  
  # Load
  observeEvent(input$do_load, {
    
    
    withProgress(message = "Loading ...",{
      
      shinyjs::show("load_read")
      setProgress(value = 0.35)
      
      load_File <- isolate(as.character(parseFilePaths(volumes, input$normal_obj)[4]))
      
      pbmc <<- readRDS(load_File) 
      
      setProgress(value = 1)
      shinyjs::hide("load_read")
      shinyjs::show("check_read")
    })
    
  }
  
  )# Load END
  
  
  # Add samples to object
  observeEvent(input$do_add, {
    
    withProgress(message = "Adding data to object ...",{
      shinyjs::show("load_add")
      
      if(input$data_type_add == "dge_file_add"){
        File_dge_add <- isolate(as.character(parseFilePaths(volumes, input$normal_dge_add)[4]))
        rawdata <- read.table(File_dge_add, header = T, row.names = 1, stringsAsFactors = F)
      }
      else{
        Path_10x_add <- isolate(as.character(parseDirPath(volumes, input$normal_10x_add)))
        rawdata <- Read10X(data.dir = Path_10x_add)
      }
      
      setProgress(value = 0.5)
      
      new_obj <- CreateSeuratObject(raw.data = rawdata, min.cells = input$min_cells, min.genes = input$min_genes, 
                                    project = input$addName)
        
      if(length(unique(pbmc@meta.data$orig.ident)) == 1) {
        pbmc.combined <- MergeSeurat(object1 = pbmc, object2 = new_obj, do.normalize = FALSE,
                                     add.cell.id1 = pbmc@meta.data$orig.ident, add.cell.id2 = input$addName)
      }
      else{
        pbmc.combined <- MergeSeurat(object1 = pbmc, object2 = new_obj, 
                                     add.cell.id2 = input$addName, do.normalize = FALSE)
      }
      setProgress(value = 0.7)
      
      mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc.combined@data), value = TRUE, ignore.case = TRUE)
      percent.mito <- Matrix::colSums(pbmc.combined@raw.data[mito.genes, ])/Matrix::colSums(pbmc.combined@raw.data)
      pbmc.combined <- AddMetaData(object = pbmc.combined, metadata = percent.mito, col.name = "percent.mito")
      
      pbmc <<- pbmc.combined
      
      setProgress(value = 1)
      
      shinyjs::hide("load_add")
      shinyjs::show("check_add")
    })
  }
  
  )# Add samples to object END
  
  output$figure_save <- renderUI({
    withSpinner(plotOutput(ns("save_figure"), width = paste0(input$pWidth/2, "px"), height = paste0(input$pHeight/2, "px")),
                type = getOption("spinner.type", default = 4))
  }
  )
  
  output$save_figure <- renderPlot({
    
    switch(input$psave_select,
           "QC" = p1,
           "GenePlot" = {
             par(mfrow = c(1, 2))
             GenePlot(object = pbmc_o, gene1 = "nUMI", gene2 = "nGene")
             GenePlot(object = pbmc_o, gene1 = "nUMI", gene2 = "percent.mito")
             },
           "Dispersion_plot" = VariableGenePlot(object = pbmc),
           "VizPCA" = VizPCA(object = pbmc, pcs.use = input$pca1_pcs1_plot:input$pca1_pcs2_plot),
           "PCAPlot" = p5,
           "PCHeatmap" = PCHeatmap(object = pbmc, pc.use = input$pca2_pcs1_heatmap:input$pca2_pcs2_heatmap,
                                   cells.use = input$pca2_cells_heatmap,
                                   do.balanced = TRUE, label.columns = FALSE),
           "JackStrawPlot" = p7@dr$pca@misc$jackstraw.plot,
           "PCElbowPlot" = p8,
           "TSNEPlot" = TSNEPlot(object = pbmc, do.hover = FALSE, pt.size = 1.5),
           "VlnPlot_markers" = p10,
           "FeaturePlot" = FeaturePlot(object = pbmc,
                                       features.plot = features_plot,
                                       cols.use = c("grey", "blue"),
                                       reduction.use = "tsne"),
           "RidgePlot" = RidgePlot(object = pbmc,
                                   features.plot = features_plot, nCol = 2),
           "DoHeatmap" = p13
    )
  })
  
  plotInput = function() {
    
    switch(input$psave_select,
           "QC" = p1,
           "GenePlot" = {
             par(mfrow = c(1, 2))
             GenePlot(object = pbmc_o, gene1 = "nUMI", gene2 = "nGene")
             GenePlot(object = pbmc_o, gene1 = "nUMI", gene2 = "percent.mito")
           },
           "Dispersion_plot" = VariableGenePlot(object = pbmc),
           "VizPCA" = VizPCA(object = pbmc, pcs.use = input$pca1_pcs1_plot:input$pca1_pcs2_plot),
           "PCAPlot" = p5,
           "PCHeatmap" = PCHeatmap(object = pbmc, pc.use = input$pca2_pcs1_heatmap:input$pca2_pcs2_heatmap,
                                   cells.use = input$pca2_cells_heatmap,
                                   do.balanced = TRUE, label.columns = FALSE),
           "JackStrawPlot" = p7@dr$pca@misc$jackstraw.plot,
           "PCElbowPlot" = p8,
           "TSNEPlot" = TSNEPlot(object = pbmc, do.hover = FALSE, pt.size = 1.5),
           "VlnPlot_markers" = p10,
           "FeaturePlot" = FeaturePlot(object = pbmc,
                                       features.plot = features_plot,
                                       cols.use = c("grey", "blue"),
                                       reduction.use = "tsne"),
           "RidgePlot" = RidgePlot(object = pbmc,
                                   features.plot = features_plot, nCol = 2),
           "DoHeatmap" = p13
    )
  }
  
  output$do_psave = downloadHandler(
    filename = "figure.png",
    content = function(file) {
      
      png(file, width = input$pWidth, height = input$pHeight,
          res = input$pRes, pointsize = input$pPsize, type = "cairo")
      print(plotInput())
      dev.off()
      
      #device <- function(..., width, height) {
      #  grDevices::png(..., width = 300, height = 300,
      #                 res = input$pRes, pointsize = input$pPsize, type = "cairo", units = "in")
      #}
      #ggsave(file, plot = plotInput(), device = device)
    }
    )

} ##### END
