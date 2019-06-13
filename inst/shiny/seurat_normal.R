

source("helpers.R")
# Data input and QC interface
seurat_normal_UI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "seurat_normal_tab",
          navbarPage(
            title = NULL,
            tabPanel("1. Data preprocessing",
                     fluidRow(

                         box(title = "1. Upload Data", width = 4, height = 350, solidHeader = TRUE, status = "info",
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

                         box(title = "2. Initialize Seurat object", width = 4, height = 350, solidHeader = TRUE, status = "warning",
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

                         box(title = "3. Cell selection for further analysis", width = 4, height = 350, solidHeader = TRUE, status = "danger",
                             collapsible = T, id = "cell_qc",

                             selectizeInput(ns("subsetNames"), label="Filter Names",
                                            choices = c("nGene", "percent.mito"),
                                            selected = c("nGene", "percent.mito"),
                                            multiple=TRUE),
                             withBusyIndicatorUI(icon_name = "cQC",
                               actionButton(ns("cell_QC"),"Check data QC", style = "width: 80%")
                             ),
                             br(),
                             htmlOutput(ns("raw_ncells")),
                             br(),
                             tags$b("-> Set cell selection criteria from 3.1 and 3.2"),
                             br(),
                             br(),
                             withBusyIndicatorUI(icon_name = "cfilter",
                                                 actionButton(ns("cell_filter"),"Select cells", style = "width: 80%")),
                             br(),
                             htmlOutput(ns("filtered_ncells"))
    
                         )

                       ),

                     fluidRow(
                              box(
                                title = "3.1 Data QC (VlnPlot)",  width = 8, height = 500, status = "danger",

                                sidebarPanel(
                                  tableOutput(ns("nGene_range")),
                                  uiOutput(ns("gene_low_high")),
                                  hr(),
                                  tableOutput(ns("mito_range")),
                                  uiOutput(ns("mito_low_high"))
                                ),
                                mainPanel(
                                  uiOutput(ns("filter_plots"))
                                )

                              ),
                              box(
                                title = "3.2 Data QC (FeatureScatter)",  width = 4, height = 500, status = "danger",
                                uiOutput(ns("gene_plots"))
                              )
                     )
            ), # Tab 1 END

            tabPanel('2. Normalization and scaling',
                     fluidRow(
                       column(width = 5,
                              
                              box(title = "Check cell cycle effects (Optional)", width = NULL, solidHeader = TRUE, status = "warning",
                                  collapsible = T,
                                  
                                  tags$b("-> Select the species of the sample"),
                                  br(),
                                  radioButtons(ns("species_cc"), NULL,
                                               c(
                                                 "Human" = "human_cc",
                                                 "Mouse" = "mouse_cc"
                                               ), selected = "human_cc", inline = TRUE),
                                  
                                  tags$b("-> Select the method to remove cell cycle effects"),
                                  br(),
                                  radioButtons(ns("methods_cc"), NULL,
                                               c(
                                                 "Normal workflow" = "normal_cc",
                                                 "Alternate workflow" = "alt_cc"
                                               ), selected = "normal_cc", inline = TRUE),
                                  
                                  withBusyIndicatorUI(icon_name = "cc",
                                                      actionButton(ns("cell_cycle"),
                                                                   "Check cell cycle effects",
                                                                   style = "width: 80%"))
                              ),
                              
                              box(title = "Normalization and scaling", width = NULL, solidHeader = TRUE, status = "success",
                                  collapsible = T,
                                  
                                  checkboxInput(ns("use_cc"), "Remove cell cycle effects? (select species and method above)", FALSE),
                                  br(),
                                  
                                  tags$b("-> Select the method for normalization and scaling"),
                                  br(),
                                  radioButtons(ns("methods_norm"), NULL,
                                               c(
                                                 "SCTransform" = "SCT_norm",
                                                 "Normal workflow" = "normal_norm"
                                               ), selected = "SCT_norm", inline = TRUE),

                                  numericInput(ns("n_vars"),
                                               label="Number of variable features to use", 
                                               min = 1, max = Inf, value = 3000, step = 1),
                                  
                                  withBusyIndicatorUI(icon_name = "norm",
                                                      actionButton(ns("data_norm"),
                                                                   "Normalize data", style = "width: 80%"))
                              )
                              ),
                     
                     column(width = 7,
                            box(title = "Cell cycle effects PCA plot", width = NULL, status = "primary",
                                uiOutput(ns("cc_plot")))
                     )
                     )
            ), # Tab 2 END

            tabPanel('3. RunPCA',
                     fluidRow(

                       box(title = "RunPCA", width = 6, solidHeader = TRUE, status = "warning",
                           collapsible = T, id = "pca1",
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
                           )
                       ),

                       box(
                         title = "Variable gene plot",  width = 6, height = 850, status = "primary",
                         uiOutput(ns("pca1_plot"))
                       )

                     )#
            ), # Tab 3 END

            tabPanel('4. PCA plots',
                     fluidRow(
                       column(width = 5,
                         box(title = "PCA plots", width = NULL, solidHeader = TRUE, status = "warning",
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

            tabPanel('5. Cluster the cells',
                     fluidRow(
                       column(width = 5,
                              box(title = "Cluster the cells", width = NULL, solidHeader = TRUE, status = "primary",
                                  collapsible = T, id = "cluster_parameters",

                                  splitLayout(
                                    numericInput(ns("cluster_dim1"),
                                                 label="Find cluster from dimension: ", min = 1, max = Inf, value = 1),
                                    numericInput(ns("cluster_dim2"),
                                                 label="to dimension: ", min = 1, max = Inf, value = 30)
                                  ),
                                  splitLayout(
                                    numericInput(ns("res"),
                                                 label="Resolution (UMAP)", min = 0, max = Inf, value = 0.5, step = 0.01)
                                  ),

                                  withBusyIndicatorUI(icon_name = "clusters",
                                                      actionButton(ns("viz_clusters"),
                                                                   "Run UMAP",
                                                                   style = "width: 90%")
                                  )),

                              box(title = "Assigning cell type identity to clusters", width = NULL, 
                                  solidHeader = TRUE, status = "success",

                                  uiOutput(ns("old_name_clusters")),
                                  textInput(ns("new_name_clusters"), "Assign new names (use ',' to separate): "),

                                  withBusyIndicatorUI(icon_name = "clusters_names",
                                                      actionButton(ns("assign_names"),
                                                                   "Assign new names",
                                                                   style = "width: 90%")
                                  )
                              ),
                              
                              box(title = "Group clusters", width = NULL, 
                                  solidHeader = TRUE, status = "warning",
                                  
                                  uiOutput(ns("group_clusters")),
                                  
                                  withBusyIndicatorUI(icon_name = "group",
                                                      actionButton(ns("group_by"),
                                                                   "Group clusters",
                                                                   style = "width: 90%")
                                  )
                              ),
                              
                              box(title = "Compare with previous identities", width = NULL, 
                                  solidHeader = TRUE, status = "info",

                                  withBusyIndicatorUI(icon_name = "compare",
                                                      actionButton(ns("compare_ids"),
                                                                   "Compare with previous identities",
                                                                   style = "width: 90%")
                                  )
                              )
                              

                       ), #

                       column(width = 7,
                              box(
                                title = "UMAP Plot",  width = NULL, status = "primary",
                                
                                tabsetPanel(type = "tabs",
                                            tabPanel("UMAP",
                                                     br(),
                                                     uiOutput(ns("tsne_plot1"))),
                                            tabPanel("Tabulate cells by cluster ID",
                                                     br(),
                                                     uiOutput(ns("umap_fracs")),
                                                     uiOutput(ns("umap_fracs_orig")))
                                )
                                
                              )
                       )#
                     )
                     ), # Tab 5 END

            tabPanel('6. Differentially expressed genes',

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
                     ), # Tab 6 END
            
            tabPanel('7. Cell type recognition (SingleR)',
                     
                     fluidRow(
                       column(width = 3,
                         box(title = "Cell type recognition (SingleR)", width = NULL, solidHeader = TRUE, status = "warning",
                             collapsible = T,
                             
                             radioButtons(ns("SR_species"), label = NULL,
                                          choices = list("Human" = "Human",
                                                         "Mouse" = "Mouse"),
                                          selected = "Human", inline = TRUE),
                             
                             checkboxInput(ns("use_fine_SR"), "Use fine-tuning? (More accurate but take a long time)", FALSE),
                             
                             withBusyIndicatorUI(icon_name = "SR",
                                                 actionButton(ns("run_sR"),
                                                              "Run SingleR",
                                                              style = "width: 80%"))
                         ),
                         
                         box(title = "Cell type visualization", width = NULL, solidHeader = TRUE, status = "success",
                             collapsible = T,
                             
                             conditionalPanel(condition  = "input.SR_species == 'Human'",
                                              selectizeInput(ns("ref_SR_human"), label="Reference set",
                                                             choices = c("Human Primary Cell Atlas (HPCA)", 
                                                                         "Blueprint+Encode"),
                                                             selected = c("Human Primary Cell Atlas (HPCA)"),
                                                             multiple=FALSE),
                                              ns = NS(id)
                             ),
                             
                             conditionalPanel(condition  = "input.SR_species == 'Mouse'",
                                              selectizeInput(ns("ref_SR_mouse"), label="Reference set",
                                                             choices = c("Immunological Genome Project (ImmGen)", 
                                                                         "Mouse-RNAseq"),
                                                             selected = c("Immunological Genome Project (ImmGen)"),
                                                             multiple=FALSE),
                                              ns = NS(id)
                             ),
                             
                             checkboxInput(ns("labels_more_SR"), "Show more cell types", FALSE),
                             
                             withBusyIndicatorUI(icon_name = "SR_viz",
                                                 actionButton(ns("viz_SR"),
                                                              "Visualize cell types",
                                                              style = "width: 80%"))
                         ),
                         
                         box(title = "Select cell types", width = NULL, solidHeader = TRUE, status = "info",
                             collapsible = T,
                             
                             uiOutput(ns("cell_types_list")),
                             br(),
                             
                             checkboxInput(ns("use_subset"), 
                                           "Use cell subset in Seurat?", 
                                           FALSE),
                             
                             withBusyIndicatorUI(icon_name = "SR_subset",
                                                 actionButton(ns("subset_SR"),
                                                              "Select cells",
                                                              style = "width: 80%"))
                         )
                       ),
                       
                       column(width = 9,
                       box(
                         title = "Cell types visulization",  width = NULL, solidHeader = TRUE, status = "primary",
                         tabsetPanel(type = "tabs",
                                     
                                     tabPanel("Cell types clusters",
                                              br(),
                                              uiOutput(ns("cluster_cell_types"))),
                                     tabPanel("Cell-type fractions",
                                              br(),
                                              uiOutput(ns("fracs_SR")),
                                              uiOutput(ns("fracs_hm_SR"))),
                                     tabPanel("Heat map",
                                              br(),
                                              uiOutput(ns("heatmap_SR"))),
                                     tabPanel("Score heatmap",
                                              br(),
                                              plotOutput(ns("score_heatmap"), height = "750px")#uiOutput(ns("heatmap_score"))
                                              ),
                                     tabPanel("Confidence of annotations",
                                              br(),
                                              splitLayout(
                                                uiOutput(ns("max_score")),
                                                uiOutput(ns("pval_SR"))
                                              )
                                     )
                         )
                       ))
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
                                        choices = c("QC", "FeatureScatter", "VizDimLoadings", "PCAPlot",
                                                    "PCHeatmap", "UMAP", "VlnPlot_markers", "FeaturePlot", 
                                                    "RidgePlot", "DoHeatmap", "UMAP_SingleR", "Heatmap_cellType"),
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

    pbmc <<- CreateSeuratObject(counts = rawdata, min.cells = input$min_cells, min.features = input$min_genes,
                                project = input$project_name)

    setProgress(value = 0.7)

    pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "(?i)^MT-") #(?i) for case insensitive
    pbmc <<- pbmc

    setProgress(value = 1)

    shinyjs::hide("load_sinit")
    shinyjs::show("check_sinit")
    })
  }

  )# Object initialization END

  #Cell QC
  observeEvent(input$cell_QC, {

    shinyjs::show("load_cQC")

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

    min_gene <- min(pbmc@meta.data$nFeature_RNA)
    max_gene <- max(pbmc@meta.data$nFeature_RNA)

    min_mito <- round(min(pbmc@meta.data$percent.mt), 2)
    max_mito <- round(max(pbmc@meta.data$percent.mt), 2)

    output$nGenePlot <- renderPlot({

      if("nGene" %in% input$subsetNames){

        VlnPlot(object = pbmc, features = c("nFeature_RNA"), nCol = 1) +
        geom_hline(yintercept = input$gene_low, color = 'darkgreen', linetype = "dashed", size = 1.5) +
        #geom_text(data = data.frame(x=1.2,y=input$gene_low), aes(x, y), label="low", vjust=1.5, hjust=0,color = "darkgreen",size = 5,fontface = "bold") +
        geom_hline(yintercept = input$gene_high, color = 'blue', linetype = "dashed", size = 1.5) +
        #geom_text(data = data.frame(x=1.2,y=input$gene_high), aes(x, y), label="high", vjust=-1, hjust=0,color = "blue",size = 5,fontface = "bold") +
        scale_y_continuous(limits=c(0, 1.1*max_gene)) + NoLegend()
      }
      else{
        VlnPlot(object = pbmc, features = c("nFeature_RNA"), nCol = 1) + NoLegend()
      }
    })
    output$nUMIPlot <- renderPlot({
      VlnPlot(object = pbmc, features = c("nCount_RNA"), nCol = 1) + NoLegend()

    })
    output$mitoPlot <- renderPlot({

      if("percent.mito" %in% input$subsetNames){

        VlnPlot(object = pbmc, features = c("percent.mt"), nCol = 1) +
        geom_hline(yintercept = input$mito_low, color = 'darkgreen', linetype = "dashed", size = 1.5) +
        #geom_text(data = data.frame(x=1.2,y=input$mito_low), aes(x, y), label="low", vjust=1.5, hjust=0,color = "darkgreen",size = 5,fontface = "bold") +
        geom_hline(yintercept = input$mito_high, color = 'blue', linetype = "dashed", size = 1.5) +
        #geom_text(data = data.frame(x=1.2,y=input$mito_high), aes(x, y), label="high", vjust=-1, hjust=0,color = "blue",size = 5,fontface = "bold") +
        scale_y_continuous(limits=c(min_mito-0.01, max_mito+0.01)) + NoLegend()
      }
      else{
        VlnPlot(object = pbmc, features = c("percent.mt"), nCol = 1) + NoLegend()
      }
    })
    
    p1 <<- VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + NoLegend()

    median_gene <- median(pbmc@meta.data$nFeature_RNA)

    min_UMI <- min(pbmc@meta.data$nCount_RNA)
    max_UMI <- max(pbmc@meta.data$nCount_RNA)
    median_UMI <- median(pbmc@meta.data$nCount_RNA)

    nGene_summary <- data.frame(c(min_gene, median_gene, max_gene),
                                c(min_UMI, median_UMI, max_UMI),
                                row.names = c("Min", "Median", "Max"))
    nGene_summary <- t(nGene_summary)
    row.names(nGene_summary) <- c("nGene", "nUMI")

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

    median_mito <- median(pbmc@meta.data$percent.mt)
    mito_summary <- data.frame(c(min_mito, median_mito, max_mito),
                                row.names = c("Min", "Median", "Max"))
    mito_summary <- t(mito_summary)
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

        FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
        geom_hline(yintercept = input$mito_low, color = 'darkgreen', linetype = "dashed", size = 1.5) +
        geom_hline(yintercept = input$mito_high, color = 'blue', linetype = "dashed", size = 1.5) + 
        NoLegend()
        #abline(h = input$mito_low, col="darkgreen", lwd=2, lty="dashed")
        #abline(h = input$mito_high, col="blue", lwd=2, lty="dashed")
      }
      else{
        FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
      }
    })

    output$gene_nUMI <- renderPlot({

      if("nGene" %in% input$subsetNames){

        FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
        geom_hline(yintercept = input$gene_low, color = 'darkgreen', linetype = "dashed", size = 1.5) +
        geom_hline(yintercept = input$gene_high, color = 'blue', linetype = "dashed", size = 1.5) +
        NoLegend()
      }
      else{
        FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
      }
    })
    
    pbmc_o <<- pbmc
    
    output$raw_ncells <- renderUI({
      
      HTML(paste("<b> ", dim(GetAssayData(object = pbmc, slot = "counts"))[2], "cells in raw data.</b>"))
      
    })
    
    shinyjs::hide("load_cQC")
    shinyjs::show("check_cQC")

  }

  )# Cell QC END

  # Cell selection
  observeEvent(input$cell_filter, {
    
    withProgress(message = "Selecting cells ...",{
      
      shinyjs::show("load_cfilter")
      setProgress(value = 0.35)
      
      gene_low_sel <<- isolate(as.numeric(input$gene_low))
      gene_high_sel <<- isolate(as.numeric(input$gene_high))
      mito_high_sel <<- isolate(as.numeric(input$mito_high))
      
      if("nGene" %in% input$subsetNames && "percent.mito" %in% input$subsetNames){
        
        pbmc <<- subset(x = pbmc, subset = nFeature_RNA > gene_low_sel & nFeature_RNA < gene_high_sel & 
                          percent.mt < mito_high_sel)
        
      }
      else if("nGene" %in% input$subsetNames){
        
        pbmc <<- subset(x = pbmc, subset = nFeature_RNA > gene_low_sel & nFeature_RNA < gene_high_sel)
        
      }
      else{
        
        pbmc <<- subset(x = pbmc, subset = percent.mt < mito_high_sel)
        
      }
      
      filtered_data <<- GetAssayData(object = pbmc) # saved data for SingleR
      
      output$filtered_ncells <- renderUI({
        
        HTML(paste("<b> ", dim(GetAssayData(object = pbmc))[2], "cells are selected.</b>"))
        
      })

      
      setProgress(value = 1)
      shinyjs::hide("load_cfilter")
      shinyjs::show("check_cfilter")
    })
    
  }
  
  )# Cell selection END
  
  # Removce cell cycle effects
  observeEvent(input$cell_cycle, {
    
    output$cc_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_cc")), type = getOption("spinner.type", default = 4))
    }
    )
    
    withProgress(message = "Anaylzing cell cycle effects ...",{
      
      shinyjs::show("load_cc")
      setProgress(value = 0.35)
      
      # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
      # segregate this list into markers of G2/M phase and markers of S phase
      s.genes <- cc.genes$s.genes
      g2m.genes <- cc.genes$g2m.genes
      
      if(input$species_cc == "mouse_cc"){
        s.genes <<- capwords(s.genes, strict = TRUE)
        g2m.genes <<- capwords(g2m.genes, strict = TRUE)
      }
      
      pbmc_cc <- NormalizeData(object = pbmc)
      pbmc_cc <- FindVariableFeatures(object = pbmc_cc, selection.method = "vst")
      pbmc_cc <- ScaleData(object = pbmc_cc)
      
      # Before cell cycle effects removing
      pbmc_cc <- CellCycleScoring(object = pbmc_cc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
      pbmc_cc <- RunPCA(object = pbmc_cc, features = c(s.genes, g2m.genes), verbose = FALSE)
      p1_cc <- DimPlot(object = pbmc_cc, reduction = "pca")
      
      # After cell cycle effects removing
      if(input$methods_cc == "normal_cc") {
        pbmc_cc <- ScaleData(object = pbmc_cc, vars.to.regress = c("S.Score", "G2M.Score"))
      }
      else{
        pbmc_cc$CC.Difference <- pbmc_cc$S.Score - pbmc_cc$G2M.Score
        pbmc_cc <- ScaleData(object = pbmc_cc, vars.to.regress = "CC.Difference")
      }
      pbmc_cc <- RunPCA(object = pbmc_cc, features = c(s.genes, g2m.genes), verbose = FALSE)
      p2_cc <-DimPlot(object = pbmc_cc, reduction = "pca")
      
      output$plot_cc <- renderPlot({
        CombinePlots(plots = list(p1_cc, p2_cc))

      })
      
      setProgress(value = 1)
      shinyjs::hide("load_cc")
      shinyjs::show("check_cc")
    })
    
  }
  
  )# Removce cell cycle effects END

  # Normalization and scaling
  observeEvent(input$data_norm, {

    withProgress(message = "Normalizing and scaling data ...",{

      shinyjs::show("load_norm")
      setProgress(value = 0.35)
      
      if(input$methods_norm == "SCT_norm") {
        pbmc <<- SCTransform(object = pbmc, vars.to.regress = "percent.mt", 
                             verbose = FALSE, variable.features.n = input$n_vars)
      }
      else{
        pbmc <<- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
        pbmc <<- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = input$n_vars)
        pbmc <<- ScaleData(object = pbmc)
        
        if(input$use_cc == TRUE){
          
          s.genes <- cc.genes$s.genes
          g2m.genes <- cc.genes$g2m.genes
          
          if(input$species_cc == "mouse_cc"){
            s.genes <<- capwords(s.genes, strict = TRUE)
            g2m.genes <<- capwords(g2m.genes, strict = TRUE)
          }
          
          pbmc <<- CellCycleScoring(object = pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
          
          if(input$methods_cc == "normal_cc") {
            pbmc <<- ScaleData(object = pbmc, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
          }
          else{
            pbmc$CC.Difference <- pbmc$S.Score - pbmc$G2M.Score
            pbmc <<- pbmc
            pbmc <<- ScaleData(object = pbmc, vars.to.regress = c("percent.mt", "CC.Difference"))
          }
        }
        else{
          pbmc <<- ScaleData(object = pbmc, vars.to.regress = "percent.mt")
        }
      }

      setProgress(value = 1)
      shinyjs::hide("load_norm")
      shinyjs::show("check_norm")
    })

  }

  )# Normalization and scaling END

  # Detect variable genes
  observeEvent(input$find_vgenes, {

    shinyjs::show("load_vgenes")

    output$vgenes_plots <- renderUI({

      withSpinner(plotOutput(ns("vgenePlot"), height = "500px"), type = getOption("spinner.type", default = 4))

    }
    )

    withProgress(message = "Detecting variable genes ...",{

      setProgress(value = 0.35)

      pbmc <<- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)

      pbmc <<- ScaleData(object = pbmc, vars.to.regress = input$vars_regress)
      
      all.genes <- rownames(x = pbmc)
      pbmc <<- ScaleData(object = pbmc, features = all.genes)

      setProgress(value = 1)

      output$selected_vgenes <- renderUI({

        HTML(paste("<b> ", length(x = VariableFeatures(object = pbmc)), "variable genes have been selected.</b>"))

      })

      output$vgenePlot <- renderPlot({
        
        top10_vst <- head(x = VariableFeatures(object = pbmc), 10)
        LabelPoints(plot = VariableFeaturePlot(object = pbmc), points = top10_vst, repel = TRUE)
      })

      shinyjs::hide("load_vgenes")
      shinyjs::show("check_vgenes")


    })

  }

  )# Detect variable genes END

  # RunPCA
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

      pbmc <<- RunPCA(object = pbmc, verbose = FALSE, features = VariableFeatures(object = pbmc))

      output$plot_pca1 <- renderPlot({
        VizDimLoadings(object = pbmc, dims = input$pca1_pcs1_plot:input$pca1_pcs2_plot, reduction = "pca")
      })

      setProgress(value = 1)
      shinyjs::hide("load_pca1")
      shinyjs::show("check_pca1")
    })

  }

  )# RunPCA END


  # PCA plots
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

      #pbmc <<- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = FALSE)

      output$plot_pca2 <- renderPlot({
        p5 <<- DimPlot(object = pbmc, reduction = "pca")
        p5
      })

      output$plot_heatmap <- renderPlot({
        DimHeatmap(object = pbmc, dims = input$pca2_pcs1_heatmap:input$pca2_pcs2_heatmap,
                  cells = input$pca2_cells_heatmap, balanced = TRUE)
      })

      setProgress(value = 1)
      shinyjs::hide("load_pca2")
      shinyjs::show("check_pca2")
    })

  }

  )# PCA plots END

  # Clusters
  observeEvent(input$viz_clusters, {

    output$tsne_plot1 <- renderUI({
      withSpinner(plotOutput(ns("plot_tsne1"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    output$umap_fracs <- renderUI({
      withSpinner(plotlyOutput(ns("fracs_umap"), height = "350px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    output$umap_fracs_orig <- renderUI({
      withSpinner(plotlyOutput(ns("fracs_umap_orig"), height = "450px"), type = getOption("spinner.type", default = 4))
    }
    )

    withProgress(message = "Plotting UMAP ...",{

      shinyjs::show("load_clusters")
      setProgress(value = 0.35)

      # Stash identities for later
      #pbmc$ClusterNames_old <- Idents(object = pbmc)
      #pbmc <<- pbmc
      
      pbmc <<- FindNeighbors(object = pbmc, dims = input$cluster_dim1:input$cluster_dim2, verbose = FALSE)
      pbmc <<- FindClusters(object = pbmc, resolution = input$res)
      pbmc <<- RunUMAP(object = pbmc, dims = input$cluster_dim1:input$cluster_dim2, verbose = FALSE)

      output$plot_tsne1 <- renderPlot({
        DimPlot(object = pbmc, reduction = "umap", label = TRUE)

      })
      
      output$fracs_umap <- renderPlotly({
        
        plot_ly(as.data.frame(table(Idents(pbmc))), labels = ~Var1, values = ~Freq, type = 'pie') %>%
          layout(title = 'Cell fractions by clusters',
                 xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                 yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
      })
      
      output$fracs_umap_orig <- renderPlotly({
        
        sx <- list(
          autotick = FALSE,
          ticklen = 0,
          dtick = 1,
          title = "Clusters"
        )
        
        sy <- list(
          autotick = FALSE,
          ticklen = 0,
          dtick = 1,
          autorange = "reversed",
          title = "Orignal ID"
        )
        
        plot_ly(
          x = row.names(prop.table(table(Idents(pbmc), pbmc$orig.ident))), 
          y = colnames(prop.table(table(Idents(pbmc), pbmc$orig.ident))),
          z = t(prop.table(table(Idents(pbmc), pbmc$orig.ident))), type = "heatmap") %>% 
          layout(xaxis = sx, yaxis = sy)
        
      })

      setProgress(value = 1)
      shinyjs::hide("load_clusters")
      shinyjs::show("check_clusters")
    })

    output$old_name_clusters <- renderUI({

      selectizeInput(ns("clusters_old_name"), label="Current names for each cluster: ",
                     choices = levels(x = pbmc),
                     selected = levels(x = pbmc)[1],
                     multiple = TRUE)
    }
    )
    
    output$group_clusters <- renderUI({
      
      selectizeInput(ns("clusters_group"), label="Current groups for clusters: ",
                     choices = colnames(pbmc@meta.data),
                     selected = colnames(pbmc@meta.data)[1],
                     multiple = FALSE)
    }
    )

    output$id1_diff_genes <- renderUI({

      selectizeInput(ns("diff_genes_id1"), label="Find all markers of cluster:",
                     choices = levels(x = pbmc),
                     selected = levels(x = pbmc)[1],
                     multiple = FALSE)
    }
    )

    output$id2_diff_genes <- renderUI({

      selectizeInput(ns("diff_genes_id2"), label="Distinguishing from clusters:",
                     choices = levels(x = pbmc)[levels(x = pbmc) != input$diff_genes_id1],
                     selected = levels(x = pbmc)[levels(x = pbmc) != input$diff_genes_id1][1],
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
    
    output$umap_fracs <- renderUI({
      withSpinner(plotlyOutput(ns("fracs_umap"), height = "350px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    output$umap_fracs_orig <- renderUI({
      withSpinner(plotlyOutput(ns("fracs_umap_orig"), height = "450px"), type = getOption("spinner.type", default = 4))
    }
    )

    withProgress(message = "Plotting UMAP ...",{

      shinyjs::show("load_clusters_names")
      setProgress(value = 0.35)
      
      current.cluster.ids <- input$clusters_old_name
      new.cluster.ids <- input$new_name_clusters
      new.cluster.ids <- unlist(strsplit(new.cluster.ids, split = ",|, "))
      
      # Stash identities for later
      pbmc$ClusterNames_old <- Idents(object = pbmc)
      Idents(object = pbmc) <- plyr::mapvalues(x = Idents(object = pbmc), from = current.cluster.ids, to = new.cluster.ids)
      pbmc <<- pbmc
      
      #names(x = new.cluster.ids) <- levels(x = pbmc)
      #pbmc <<- RenameIdents(object = pbmc, new.cluster.ids)

      output$plot_tsne1 <- renderPlot({

        DimPlot(object = pbmc, reduction = "umap", label = TRUE, pt.size = 1.5)

      })
      
      output$fracs_umap <- renderPlotly({
        
        plot_ly(as.data.frame(table(Idents(pbmc))), labels = ~Var1, values = ~Freq, type = 'pie') %>%
          layout(title = 'Cell fractions by clusters',
                 xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                 yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
      })
      
      output$fracs_umap_orig <- renderPlotly({
        
        sx <- list(
          autotick = FALSE,
          ticklen = 0,
          dtick = 1,
          title = "Clusters"
        )
        
        sy <- list(
          autotick = FALSE,
          ticklen = 0,
          dtick = 1,
          autorange = "reversed",
          title = "Orignal ID"
        )
        
        plot_ly(
          x = row.names(prop.table(table(Idents(pbmc), pbmc$orig.ident))), 
          y = colnames(prop.table(table(Idents(pbmc), pbmc$orig.ident))),
          z = t(prop.table(table(Idents(pbmc), pbmc$orig.ident))), type = "heatmap") %>% 
          layout(xaxis = sx, yaxis = sy)
        
      })
      
      setProgress(value = 1)
      shinyjs::hide("load_clusters_names")
      shinyjs::show("check_clusters_names")
    })

    output$old_name_clusters <- renderUI({

      selectizeInput(ns("clusters_old_name"), label="Select names to change: ",
                     choices = levels(x = pbmc),
                     selected = levels(x = pbmc)[1],
                     multiple = TRUE)
    }
    )
    
    output$group_clusters <- renderUI({
      
      selectizeInput(ns("clusters_group"), label="Current groups for clusters: ",
                     choices = colnames(pbmc@meta.data),
                     selected = colnames(pbmc@meta.data)[1],
                     multiple = FALSE)
    }
    )

    output$id1_diff_genes <- renderUI({

      selectizeInput(ns("diff_genes_id1"), label="Find all markers of cluster:",
                     choices = levels(x = pbmc),
                     selected = levels(x = pbmc)[1],
                     multiple = FALSE)
    }
    )

    output$id2_diff_genes <- renderUI({

      selectizeInput(ns("diff_genes_id2"), label="Distinguishing from clusters:",
                     choices = levels(x = pbmc)[levels(x = pbmc) != input$diff_genes_id1],
                     selected = levels(x = pbmc)[levels(x = pbmc) != input$diff_genes_id1][1],
                     multiple = TRUE)
    }
    )

  }
  )# Assign new names END
  
  # Compare with old identities
  observeEvent(input$compare_ids, {
    
    output$tsne_plot1 <- renderUI({
      withSpinner(plotOutput(ns("plot_tsne1"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    withProgress(message = "Plotting UMAP ...",{
      
      shinyjs::show("load_compare")
      setProgress(value = 0.35)
      
      output$plot_tsne1 <- renderPlot({
        
        plot1_compare <- DimPlot(object = pbmc, reduction = "umap", label = TRUE, pt.size = 1.5) + NoLegend()
        plot2_compare <- DimPlot(object = pbmc, reduction = "umap", 
                                 label = TRUE, pt.size = 1.5, group.by = "orig.ident") + NoLegend()
        CombinePlots(plots = list(plot1_compare, plot2_compare))
        
      })
      
      setProgress(value = 1)
      shinyjs::hide("load_compare")
      shinyjs::show("check_compare")
    })
  })# Compare with old identities END
  
  # Group clusters
  observeEvent(input$group_by, {
    
    output$tsne_plot1 <- renderUI({
      withSpinner(plotOutput(ns("plot_tsne1"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    withProgress(message = "Plotting UMAP ...",{
      
      shinyjs::show("load_group")
      setProgress(value = 0.35)
      
      output$plot_tsne1 <- renderPlot({
        
       DimPlot(object = pbmc, reduction = "umap", 
               label = TRUE, pt.size = 1.5, group.by = input$clusters_group)
      })
      
      setProgress(value = 1)
      shinyjs::hide("load_group")
      shinyjs::show("check_group")
    })
  })# Group clusters END

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
        cluster.markers <<- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25,
                                          logfc.threshold = 0.25)

        cluster.markers <<- cluster.markers %>% group_by(cluster) %>% top_n(input$top_markers, avg_logFC)

        output$diff_genes_heatmap <- renderUI({
          withSpinner(plotOutput(ns("heatmap_diff_genes"), height = "750px"), type = getOption("spinner.type", default = 4))
        }
        )

        output$heatmap_diff_genes <- renderPlot({
          p13 <<- DoHeatmap(object = pbmc, features = cluster.markers$gene)
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
        p10 <<- VlnPlot(object = pbmc, features = cluster.markers$gene[input$table_diff_genes_rows_selected])
        p10
        })

        output$featureplot_diff_genes <- renderPlot({

          features_plot <<- cluster.markers$gene[input$table_diff_genes_rows_selected]
          
          validate(need(length(input$table_diff_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
          FeaturePlot(object = pbmc,
                      features = cluster.markers$gene[input$table_diff_genes_rows_selected])
          
        })

        output$ridgeplot_diff_genes <- renderPlot({

          features_plot <<- cluster.markers$gene[input$table_diff_genes_rows_selected]
          
          validate(need(length(input$table_diff_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
          RidgePlot(object = pbmc,
                      features = cluster.markers$gene[input$table_diff_genes_rows_selected], ncol = 2)
          
        })


    }
    else{

      output$vlnplot_diff_genes <- renderPlot({

        validate(need(length(input$table_diff_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        p10 <<- VlnPlot(object = pbmc,
                features = rownames(cluster.markers)[input$table_diff_genes_rows_selected])
        p10
      })

      output$featureplot_diff_genes <- renderPlot({

        features_plot <<- rownames(cluster.markers)[input$table_diff_genes_rows_selected]
        validate(need(length(input$table_diff_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        FeaturePlot(object = pbmc,
                    features = features_plot)
        
        
      })

      output$ridgeplot_diff_genes <- renderPlot({

        features_plot <<- rownames(cluster.markers)[input$table_diff_genes_rows_selected]
        
        validate(need(length(input$table_diff_genes_rows_selected) != "0", "Error: Please select at least one feature!"))
        RidgePlot(object = pbmc,
                  features = features_plot, ncol = 2)
      })

    }

  }
  )# Cluster biomarkers END
  
  # Run SingleR
  observeEvent(input$run_sR, {
    
    withProgress(message = "Running SingleR ...",{
      
      shinyjs::show("load_SR")
      setProgress(value = 0.35)
      
      #filtered_data <<- GetAssayData(object = pbmc, slot = "counts")
      singler <<- CreateSinglerObject(filtered_data, annot = NULL, project.name = "SingleR", min.genes = 0,
                                      technology = "RNAseq", species = input$SR_species, citation = "",
                                      normalize.gene.length = F, variable.genes = "de",
                                      fine.tune = input$use_fine_SR, do.signatures = T, clusters = NULL, do.main.types = T, 
                                      reduce.file.size = T, numCores = SingleR.numCores)
      
      seurat.object <- pbmc
      
      singler$seurat <- seurat.object
      singler <<- singler
      
      singler$meta.data$orig.ident <- seurat.object@meta.data$orig.ident
      singler <<- singler
      
      singler$meta.data$xy = seurat.object@reductions$umap@cell.embeddings
      singler <<- singler
      
      singler$meta.data$clusters = seurat.object@active.ident
      singler <<- singler
      
      setProgress(value = 1)
      shinyjs::hide("load_SR")
      shinyjs::show("check_SR")
    })
    

  }
  )# Run SingleR END
  
  
  # Cell types visualization
  observeEvent(input$viz_SR, {
    
    output$cluster_cell_types <- renderUI({
      withSpinner(plotOutput(ns("cell_types_cluster"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    output$fracs_SR <- renderUI({
      withSpinner(plotlyOutput(ns("SR_fracs"), height = "450px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    output$fracs_hm_SR <- renderUI({
      withSpinner(plotlyOutput(ns("SR_fracs_hm"), height = "450px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    output$heatmap_SR <- renderUI({
      withSpinner(plotOutput(ns("SR_heatmap"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    #output$heatmap_score <- renderUI({
    #  withSpinner(plotOutput(ns("score_heatmap")), type = getOption("spinner.type", default = 4))
    #}
    #)
    
    output$max_score <- renderUI({
      withSpinner(plotOutput(ns("score_max"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    output$pval_SR <- renderUI({
      withSpinner(plotOutput(ns("SR_pval"), height = "750px"), type = getOption("spinner.type", default = 4))
    }
    )
    
    
    withProgress(message = "Plotting ...",{
      
      shinyjs::show("load_SR_viz")
      setProgress(value = 0.35)
      
      if(input$SR_species == "Human" && input$ref_SR_human == "Human Primary Cell Atlas (HPCA)" || 
         input$SR_species == "Mouse" && input$ref_SR_mouse == "Immunological Genome Project (ImmGen)"){
        
        labels_obj <- singler$singler[[1]]$SingleR.single.main
        labels_more_obj <- singler$singler[[1]]$SingleR.single
        
        labels = singler$singler[[1]]$SingleR.single.main$labels
        labels_more = singler$singler[[1]]$SingleR.single$labels
        
      }
      else{
        
        labels_obj <- singler$singler[[2]]$SingleR.single.main
        labels_more_obj <- singler$singler[[2]]$SingleR.single
        
        labels = singler$singler[[2]]$SingleR.single.main$labels
        labels_more = singler$singler[[2]]$SingleR.single$labels
      }
      
      #tbl = table(labels)
      #labels[labels %in% names(tbl)[tbl<20]] = 'X'

      #tbl_more = table(labels_more)
      #labels_more[labels_more %in% names(tbl_more)[tbl_more<20]] = 'X'
      
      pbmc@meta.data$cell_types <- labels
      pbmc <<- pbmc
      
      pbmc@meta.data$cell_types_more <- labels_more
      pbmc <<- pbmc
      
      output$cell_types_cluster <- renderPlot({
        
        if(input$labels_more_SR == FALSE){
          DimPlot(object = pbmc, label = TRUE, group.by = "cell_types")
        }
        else{
          DimPlot(object = pbmc, label = TRUE, group.by = "cell_types_more")
        }

      })
      
      output$SR_fracs <- renderPlotly({
        
        if(input$labels_more_SR == FALSE){
          plot_ly(as.data.frame(table(labels)), labels = ~labels, values = ~Freq, type = 'pie') %>%
            layout(title = 'Cell type fractions',
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
        }
        else{
          plot_ly(as.data.frame(table(labels_more)), labels = ~labels_more, values = ~Freq, type = 'pie') %>%
            layout(title = 'Cell type fractions',
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
        }
        
      })
      
      output$SR_fracs_hm <- renderPlotly({
        
        s <- list(
               autotick = FALSE,
               ticklen = 0,
               dtick = 1,
               autorange = "reversed",
               title = "Cluster in UMAP1"
          )
        
        if(input$labels_more_SR == FALSE){
          fracs_clusters <- as.matrix(table(pbmc$seurat_clusters, pbmc$cell_types))
          fracs_clusters <- fracs_clusters/rowSums(fracs_clusters)
        }
        else{
          fracs_clusters <- as.matrix(table(pbmc$seurat_clusters, pbmc$cell_types_more))
          fracs_clusters <- fracs_clusters/rowSums(fracs_clusters)
        }
        
        plot_ly(
          x = colnames(fracs_clusters), y = row.names(fracs_clusters),
          z = fracs_clusters, type = "heatmap") %>% layout(xaxis = list(ticklen = 0), yaxis = s)
        
      })
      
      output$SR_heatmap <- renderPlot({
        DoHeatmap(object = pbmc, features = cluster.markers$gene, group.by = "cell_types")
      })
      
      output$score_heatmap <- renderPlot({
        
        
        if(input$labels_more_SR == FALSE){
          p_s_SR <- SingleR.DrawHeatmap(labels_obj, top.n = Inf, clusters = pbmc$seurat_clusters)
          p_s_SR
          dev.off()
          p_s_SR
        
        }
        else{
          p_s_SR <- SingleR.DrawHeatmap(labels_more_obj, top.n = 50, clusters = pbmc$seurat_clusters)
          p_s_SR
          dev.off()
          p_s_SR
        }
        
      })
      
      output$score_max <- renderPlot({
        
        
        if(input$labels_more_SR == FALSE){
          annotation.confidence.SR(labels_obj, plot.feature = "MaxScore")
          
        }
        else{
          annotation.confidence.SR(labels_more_obj, plot.feature = "MaxScore")
          
        }
        
      })
      
      output$SR_pval <- renderPlot({
        
        
        if(input$labels_more_SR == FALSE){
          annotation.confidence.SR(labels_obj, plot.feature = "p-value")
          
        }
        else{
          annotation.confidence.SR(labels_more_obj, plot.feature = "p-value")
          
        }
        
      })
      

      setProgress(value = 1)
      shinyjs::hide("load_SR_viz")
      shinyjs::show("check_SR_viz")
    })
    
    output$cell_types_list <- renderUI({
      
      if(input$labels_more_SR == FALSE){
        choices_cell_types = unique(pbmc$cell_types)
      }
      else{
        choices_cell_types = unique(pbmc$cell_types_more)
      }
      
      selectizeInput(ns("list_cell_types"), label="Select cell types: ",
                     choices = choices_cell_types,
                     selected = choices_cell_types[1],
                     multiple = TRUE)
    }
    )
    
  }
  )# Cell types visualization END
  
  
  # Select cell types
  observeEvent(input$subset_SR, {
    
    
    withProgress(message = "Selecting cells ...",{
      
      shinyjs::show("load_SR_subset")
      setProgress(value = 0.35)
      
      if(input$labels_more_SR == FALSE){
        Idents(pbmc) <- "cell_types"
      }
      else{
        Idents(pbmc) <- "cell_types_more"
      }
      
      pbmc_subset <<- subset(pbmc, idents = input$list_cell_types)
      
      pbmc_subset <<- SCTransform(object = pbmc_subset, vars.to.regress = "percent.mt", 
                                  verbose = FALSE, variable.features.n = input$n_vars)
      setProgress(value = 0.6)
      pbmc_subset <<- RunPCA(object = pbmc_subset, verbose = FALSE, 
                             features = VariableFeatures(object = pbmc_subset))
      
      if(input$use_subset == FALSE){
        pbmc_subset <<- FindNeighbors(object = pbmc_subset, dims = input$cluster_dim1:input$cluster_dim2, 
                                      verbose = FALSE)
        pbmc_subset <<- FindClusters(object = pbmc_subset, resolution = input$res)
        
        setProgress(value = 0.8)
        pbmc_subset <<- RunUMAP(object = pbmc_subset, dims = input$cluster_dim1:input$cluster_dim2, 
                                verbose = FALSE)
      }
      else{
        pbmc <<- pbmc_subset
      }
      
      setProgress(value = 1)
      shinyjs::hide("load_SR_subset")
      shinyjs::show("check_SR_subset")
    })
    
  }
  
  )# Select cell types END
  
  
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
      
      new_obj <- CreateSeuratObject(counts = rawdata, min.cells = input$min_cells, min.features = input$min_genes, 
                                    project = input$addName)
        
      if(length(unique(pbmc@meta.data$orig.ident)) == 1) {
        pbmc.combined <- merge(x = pbmc, y = new_obj, merge.data = FALSE,
                               add.cell.ids = c(as.character(unique(pbmc@meta.data$orig.ident)), input$addName))
      }
      else{
        pbmc.combined <- merge(x = pbmc, y = new_obj, merge.data = FALSE,
                               add.cell.ids = c("c", input$addName))
      }
      setProgress(value = 0.7)
      
      pbmc.combined[["percent.mt"]] <- PercentageFeatureSet(object = pbmc.combined, pattern = "(?i)^MT-") #(?i) for case insensitive
      pbmc.combined <<- pbmc.combined
      
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
           "FeatureScatter" = {
             plot1 <- FeatureScatter(object = pbmc_o, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
             plot2 <- FeatureScatter(object = pbmc_o, feature1 = "nCount_RNA", feature2 = "percent.mt")
             CombinePlots(plots = list(plot1, plot2))
             },
           "VizDimLoadings" = VizDimLoadings(object = pbmc, pcs.use = input$pca1_pcs1_plot:input$pca1_pcs2_plot, 
                                             reduction = "pca"),
           "PCAPlot" = p5,
           "PCHeatmap" = DimHeatmap(object = pbmc, dims = input$pca2_pcs1_heatmap:input$pca2_pcs2_heatmap,
                                    cells = input$pca2_cells_heatmap, balanced = TRUE),
           "UMAP" = DimPlot(object = pbmc, reduction = "umap", pt.size = 1.5),
           "VlnPlot_markers" = p10,
           "FeaturePlot" = FeaturePlot(object = pbmc,
                                       features = features_plot),
           "RidgePlot" = RidgePlot(object = pbmc, features = features_plot, ncol = 2),
           "DoHeatmap" = p13,
           "UMAP_SingleR" = DimPlot(object = pbmc, label = TRUE, group.by = "cell_types", pt.size = 1.5),
           "Heatmap_cellType" = DoHeatmap(object = pbmc, features = cluster.markers$gene, group.by = "cell_types", label = FALSE)
    )
  })
  
  plotInput = function() {
    
    switch(input$psave_select,
           "QC" = p1,
           "FeatureScatter" = {
             plot1 <- FeatureScatter(object = pbmc_o, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
             plot2 <- FeatureScatter(object = pbmc_o, feature1 = "nCount_RNA", feature2 = "percent.mt")
             CombinePlots(plots = list(plot1, plot2))
           },
           "VizDimLoadings" = VizDimLoadings(object = pbmc, pcs.use = input$pca1_pcs1_plot:input$pca1_pcs2_plot,
                                             reduction = "pca"),
           "PCAPlot" = p5,
           "PCHeatmap" = DimHeatmap(object = pbmc, dims = input$pca2_pcs1_heatmap:input$pca2_pcs2_heatmap,
                                    cells = input$pca2_cells_heatmap, balanced = TRUE),
           "UMAP" = DimPlot(object = pbmc, reduction = "umap", pt.size = 1.5),
           "VlnPlot_markers" = p10,
           "FeaturePlot" = FeaturePlot(object = pbmc,
                                       features = features_plot),
           "RidgePlot" = RidgePlot(object = pbmc, features = features_plot, ncol = 2),
           "DoHeatmap" = p13,
           "UMAP_SingleR" = DimPlot(object = pbmc, label = TRUE, group.by = "cell_types", pt.size = 1.5),
           "Heatmap_cellType" = DoHeatmap(object = pbmc, features = cluster.markers$gene, group.by = "cell_types", label = FALSE)
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
