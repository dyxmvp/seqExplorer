
source("helpers.R")

monocle_UI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "monocle_tab",
          navbarPage(
            title = NULL,
            tabPanel("1. Data selection",
                     fluidRow(
                       column(width = 3,
                         box(title = "1. Select Data", width = NULL, solidHeader = TRUE, status = "info",
                           collapsible = T,
                           tags$style(appCSS),
                           radioButtons(ns("obj_mono"), NULL,
                                        c(
                                          "All seurat object" = "all_seurat",
                                          "Subset seurat object" = "sub_seurat"
                                        ), selected = "all_seurat"),
                           withBusyIndicatorUI(icon_name = "minit",
                                               actionButton(ns("mono_init"),"Load data", style = "width: 80%")
                           )

                         ),

                         box(title = "Group clusters", width = NULL, 
                             solidHeader = TRUE, status = "warning",
                             
                             uiOutput(ns("group_clusters_mono")),
                             
                             withBusyIndicatorUI(icon_name = "group_mono",
                                                 actionButton(ns("group_by_mono"),
                                                              "Visualize",
                                                              style = "width: 80%")
                             )
                         )
                        ),
                       
                       column(width = 9,
                              box(
                                title = "Single-cell trajectories",  width = NULL, status = "primary",
                                
                                tabsetPanel(type = "tabs",
                                            tabPanel("Single-cell trajectories",
                                                     br(),
                                                     uiOutput(ns("mono_plot"))),
                                            tabPanel("Genes change over pseudotime",
                                                     br(),
                                                     uiOutput(ns("mono_gene"))),
                                            tabPanel("Pseudotemporal expression pattern",
                                                     br(),
                                                     uiOutput(ns("mono_heatmap")))
                                )
                                
                              )
                       )
                     )
                         
            ), # Tab 1 END

            tabPanel('* Save figures',
                     fluidRow(
                       
                       box(
                         title = "Save figure", width = 3, solidHeader = TRUE, status = "primary",
                         
                         selectizeInput(ns("psave_select"), label="Choose a figure to save",
                                        choices = c("QC", "FeatureScatter", "VizDimLoadings", "PCAPlot",
                                                    "PCHeatmap", "UMAP", "VlnPlot_markers", "FeaturePlot", 
                                                    "RidgePlot", "DoHeatmap"),
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
            )# Tab 2 END
          )
  )

}

monocle_Server <- function(input, output, session, fileRoot = NULL) {

  ns <- session$ns

  # Object initialization
  observeEvent(input$mono_init, {

    withProgress(message = "Initializing Monocle Object ...",{
    shinyjs::show("load_minit")

    setProgress(value = 0.35)
    
    if(input$obj_mono == "all_seurat"){
      
      lung <- importCDS_new(pbmc, import_all = TRUE)
      dCellData <- estimateSizeFactors(lung)
      dCellData <- estimateDispersions(dCellData)
      
      setProgress(value = 0.5)
      
      dCellData <- setOrderingFilter(dCellData, cluster.markers$gene)
      
      setProgress(value = 0.7)
      dCellDataSet <- reduceDimension(dCellData, norm_method = "log")
      dCellDataSet <<- orderCells(dCellDataSet, reverse = FALSE)
    }
    else{
      
      lung <- importCDS_new(pbmc_subset, import_all = TRUE)
      dCellData <- estimateSizeFactors(lung)
      dCellData <- estimateDispersions(dCellData)
      cluster.markers_sub <- FindAllMarkers(object = pbmc_subset, only.pos = TRUE, min.pct = 0.25,
                                         logfc.threshold = 0.25)

      setProgress(value = 0.6)
      dCellData <- setOrderingFilter(dCellData, cluster.markers_sub$gene)
      
      setProgress(value = 0.8)
      dCellDataSet <<- reduceDimension(dCellData, norm_method = "log")
      dCellDataSet <<- orderCells(dCellDataSet, reverse = FALSE)
    }

    setProgress(value = 1)

    shinyjs::hide("load_minit")
    shinyjs::show("check_minit")
    })
    
    output$group_clusters_mono <- renderUI({
      
      if(input$obj_mono == "all_seurat"){
        selectizeInput(ns("clusters_group_mono"), label="Current groups for clusters: ",
                       choices = c("State", colnames(pbmc@meta.data)),
                       selected = colnames(pbmc@meta.data)[1],
                       multiple = FALSE)
      }
      else{
        selectizeInput(ns("clusters_group_mono"), label="Current groups for clusters: ",
                       choices = c("State", colnames(pbmc_subset@meta.data)),
                       selected = colnames(pbmc_subset@meta.data)[1],
                       multiple = FALSE)
      }
    }
    )
  }

  )# Object initialization END

  # Visualize
  observeEvent(input$group_by_mono, {

    output$mono_plot <- renderUI({
      withSpinner(plotOutput(ns("plot_mono"), height = "750px"), type = getOption("spinner.type", default = 4))
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

    withProgress(message = "Plotting figure ...",{

      shinyjs::show("load_group_mono")
      setProgress(value = 0.35)

      output$plot_mono <- renderPlot({
        plot_cell_trajectory(dCellDataSet, color_by = input$clusters_group_mono)

      })
      

      setProgress(value = 1)
      shinyjs::hide("load_group_mono")
      shinyjs::show("check_group_mono")
    })

  }
  )# Visualize END

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
           "RidgePlot" = RidgePlot(object = pbmc, features = features_plot, nCol = 2),
           "DoHeatmap" = p13
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
           "RidgePlot" = RidgePlot(object = pbmc, features = features_plot, nCol = 2),
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
    }
    )

} ##### END
