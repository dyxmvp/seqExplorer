
# Global varibles


# Quick interface
lookUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "look_tab",
          navbarPage(title = NULL,
                     tabPanel("Data preprocessing",
                              fluidRow(
                                box(
                                  title = "Optional: Filter cells based on genes", width = 3, solidHeader = TRUE, status = "danger", collapsible = T,
                                  
                                  textInput(ns("genes_filter"), "Type gene list (use ',' to separate): "),
                                  
                                  tags$b("Select dge files:"),
                                  tags$p(),
                                  verbatimTextOutput(ns("samples_files_filter"), placeholder = TRUE),
                                  shinyFilesButton(ns("file_samples_filter"), "Select dge files", "Please select files", multiple = TRUE),
                                  tags$p(),
                                  
                                  radioButtons(ns("keep_del"), label = "Keep or delete selected cells?",
                                               choices = list("Keep" = "1", 
                                                              "Delete" = "2"), 
                                               selected = "1"),

                                  withBusyIndicatorUI(icon_name = "filter",
                                                      actionButton(ns("do_filter"),
                                                                   "Filter",
                                                                   style = "width: 70%")
                                  )
                                ),       
                                
                                box(
                                         title = "Step 1: Select samples", width = 4, solidHeader = TRUE, status = "warning", collapsible = T,
                                         
                                         selectizeInput(ns("look_filetype"), label="Select a file type",
                                                        choices = c("DGE files (*.dge.*)" = "dge_look",
                                                                    "Combined files (*.txt*)" = "combined_look"
                                                        ),
                                                        selected = c("DGE files (*.dge.*)"), multiple = FALSE),
                                         
                                         conditionalPanel(condition  = "input.look_filetype == 'dge_look'",
                                             tags$b("Select dge files:"),
                                             tags$p(),
                                             verbatimTextOutput(ns("samples_files_look"), placeholder = TRUE),
                                             shinyFilesButton(ns("file_samples_look"), "Select dge files", "Please select files", multiple = TRUE),
                                             tags$p(),
                                             ns = NS(id)
                                         ),
                                         
                                         conditionalPanel(condition  = "input.look_filetype == 'combined_look'",
                                                          tags$b("Select a cell fraction file:"),
                                                          tags$p(),
                                                          verbatimTextOutput(ns("combined_files_look"), placeholder = TRUE),
                                                          shinyFilesButton(ns("file_combined_look"), "Select a cell fraction file", "Please select a file", multiple = FALSE),
                                                          tags$p(),
                                                          tags$b("Select an average transcripts file:"),
                                                          tags$p(),
                                                          verbatimTextOutput(ns("avg_files_look"), placeholder = TRUE),
                                                          shinyFilesButton(ns("file_avg_look"), "Select an average transcripts file", "Please select a file", multiple = FALSE),
                                                          tags$p(),
                                                          ns = NS(id)
                                         ),
                                         
                                         withBusyIndicatorUI(icon_name = "upload",
                                                             actionButton(ns("do_upload"),
                                                                          "Upload data",
                                                                          style = "width: 90%")
                                         )
                                       ),
                                       
                                       box(
                                         title = "Step 2: Input gene list", width = 5, solidHeader = TRUE, status = "success", collapsible = T,
                                         
                                         textInput(ns("gene_list"), "Type gene list (use ',' to separate): "),
                                         
                                         numericInput(ns("p_ncols"), label="Number of columns in the plot",
                                                      min = 1, max = Inf, value = 5),
                                         numericInput(ns("bar_width"), label="Bar width in the plot",
                                                      min = 0, max = 1, value = 0.9, step = 0.01),

                                         withBusyIndicatorUI(icon_name = "quick_look",
                                                             actionButton(ns("do_look"),
                                                                          "Plot",
                                                                          style = "width: 90%")
                                         )
                                       ),

                                       box(
                                         title = "Percentage of cells expressing selected genes",  width = 12, status = "primary", collapsible = T,
                                           uiOutput(ns("gene_plots_look"))
                                       ),
                                       
                                       box(
                                         title = "Average number of transcripts each gene expressed",  width = 12, status = "info", collapsible = T,
                                         uiOutput(ns("gene_table_look"))
                                       )
                                       
                              )
                     ), # Tab 1 END
                     
                     
                     tabPanel('* Save figures',
                              fluidRow(
                                
                                box(
                                  title = "Save figure", width = 2, solidHeader = TRUE, status = "primary",
                                  
                                  numericInput(ns("pWidth_look"), "Width (px)", value = 1200),   # Width is divided by 2 for display
                                  numericInput(ns("pHeight_look"), "Height (px)", value = 1200), # Height is divided by 2 for display
                                  numericInput(ns("pRes_look"), "Resolution (ppi)", value = 300),
                                  numericInput(ns("pPsize_look"), "Pointsize", value = 12),
                                  
                                  downloadButton(ns("do_psave_look"), "Save figure", style = "width: 100%")
                                ),
                                
                                box(
                                  title = "Figure to save", width = 10, status = "success",
                                  
                                  uiOutput(ns("figure_save_look"))
                                )
                                
                              )#
                     )# Tab 2 END
                     
          ) # navbarPage END
  ) # TabItem END
  
}

lookServer <- function(input, output, session, fileRoot = NULL) {
  
  ns <- session$ns
  
  volumes <- (c(Home = fs::path_home(), getVolumes()()))
  
  shinyFileChoose(input, "file_samples_filter", roots = volumes, session = session, filetypes=c('', 'txt', 'gz'))
  shinyFileChoose(input, "file_samples_look", roots = volumes, session = session, filetypes=c('', 'txt', 'gz'))
  shinyFileChoose(input, "file_combined_look", roots = volumes, session = session, filetypes=c('', 'txt'))
  shinyFileChoose(input, "file_avg_look", roots = volumes, session = session, filetypes=c('', 'txt'))

  output$samples_files_filter <- renderPrint({
    as.character(parseFilePaths(volumes, input$file_samples_filter)[4])
  })
  output$samples_files_look <- renderPrint({
    as.character(parseFilePaths(volumes, input$file_samples_look)[4])
  })
  output$combined_files_look <- renderPrint({
    as.character(parseFilePaths(volumes, input$file_combined_look)[4])
  })
  output$avg_files_look <- renderPrint({
    as.character(parseFilePaths(volumes, input$file_avg_look)[4])
  })
  
  
  # Filter cells 
  observeEvent(input$do_filter, {
    
    withProgress(message = "Please wait ...",{
      
      shinyjs::show("load_filter")

      setProgress(value = 0.35)
      
      keep_filter <<- isolate(as.character(input$keep_del))

      sample_list_filter <<- isolate(as.character(parseFilePaths(volumes, input$file_samples_filter)$datapath))
      gene_list_filter <<- input$genes_filter
      gene_list_filter <<- unlist(strsplit(gene_list_filter, split = ",|, "))
      
      lapply(sample_list_filter, upload_dges_filter)
      
      setProgress(value = 1)
      shinyjs::hide("load_filter")
      shinyjs::show("check_filter")
    })
  })# Filter cells END
  
  # Upload data 
  observeEvent(input$do_upload, {
    
    output$gene_table_look <- renderUI({
        withSpinner(DTOutput(ns("look_gene_table")), type = getOption("spinner.type", default = 4))
    })
    
    withProgress(message = "Please wait ...",{
      
      shinyjs::show("load_upload")
      
      #sample_list <<- isolate(as.character(parseFilePaths(volumes, input$file_samples_look)$datapath))
      #sample_list <<- isolate(as.character(parseFilePaths(volumes, input$file_samples_look)$name))
      #sample_list <<- sub(pattern = "(.*?)\\..*$", replacement = "\\1", sample_list)

      setProgress(value = 0.35)
      
      #genes_look <- input$gene_list
      #genes_look <<- unlist(strsplit(genes_look, split = ",|, "))
      
      #gene_frac_list <<- lapply(sample_list, gene_frac_df)
      
      if(input$look_filetype == 'dge_look'){
        
        sample_list <<- isolate(as.character(parseFilePaths(volumes, input$file_samples_look)$datapath))
        lapply(sample_list, upload_dges)
        
        if(length(sample_list) > 1){
          
          df_merge <<- df_merge[, -1]
          df_merge_avg <<- df_merge_avg[, -1]
        }
        
        df_merge[is.na(df_merge)] <- 0
        df_merge_avg[is.na(df_merge_avg)] <- 0
        
        write.table(df_merge, file.path(dirname(sample_list[1]), "cell_frac.txt"), append = FALSE, sep = "\t", dec = ".",
                    row.names = TRUE, col.names = TRUE)
        write.table(df_merge_avg, file.path(dirname(sample_list[1]), "genes_avg.txt"), append = FALSE, sep = "\t", dec = ".",
                    row.names = TRUE, col.names = TRUE)
      }
      else{
        
        sample_list <<- isolate(as.character(parseFilePaths(volumes, input$file_combined_look)$datapath))
        df_merge <<- read.table(sample_list, header = T, row.names = 1, stringsAsFactors = F)
        
        sample_list_avg <- isolate(as.character(parseFilePaths(volumes, input$file_avg_look)$datapath))
        df_merge_avg <<- read.table(sample_list_avg, header = T, row.names = 1, stringsAsFactors = F)
      }
      
      setProgress(value = 0.6)

      output$look_gene_table <- renderDT(datatable(df_merge_avg, rownames = TRUE, selection = "none"))
      
      setProgress(value = 1)
      shinyjs::hide("load_upload")
      shinyjs::show("check_upload")
    })
  })# Upload data END
  

  # plots
  observeEvent(input$do_look, {
    
    withProgress(message = "Please wait ...",{
      
      shinyjs::show("load_quick_look")
      
      
      setProgress(value = 0.35)
      
      width_bar <<- input$bar_width
      
      genes_look <<- input$gene_list
      genes_look <<- unlist(strsplit(genes_look, split = ",|, "))
      
      if(length(colnames(df_merge)) == 1){
        
        df_look <<- as.data.frame(df_merge[match(genes_look, rownames(df_merge)),])
        colnames(df_look) <<- colnames(df_merge)
        rownames(df_look) <<- genes_look
      }
      else{
        
        df_look <<- df_merge[match(genes_look, rownames(df_merge)),]
      }
      
      #sample_list <- colnames(df_look)
      
      pl <<- lapply(colnames(df_look), gene_frac_plot)
      
      output$gene_plots_look <- renderUI({
          withSpinner(plotOutput(ns("look_plot")), type = getOption("spinner.type", default = 4))
      })

      output$look_plot <- renderPlot({
          plot_grid(plotlist = pl, ncol = input$p_ncols)
      })
      
      setProgress(value = 1)
      shinyjs::hide("load_quick_look")
      shinyjs::show("check_quick_look")
    })
  })# plots END
  
  
  output$figure_save_look <- renderUI({
    withSpinner(plotOutput(ns("save_figure_look"), width = paste0(input$pWidth_look/2, "px"), 
                           height = paste0(input$pHeight_look/2, "px")),
                type = getOption("spinner.type", default = 4))
  }
  )
  
  output$save_figure_look <- renderPlot({
    plot_grid(plotlist = pl, ncol = input$p_ncols)
  })
  
  plotInput = function() {
    plot_grid(plotlist = pl, ncol = input$p_ncols)
  }
  
  output$do_psave_look = downloadHandler(
    filename = "figure.png",
    content = function(file) {
      png(file, width = input$pWidth_look, height = input$pHeight_look,
          res = input$pRes_look, pointsize = input$pPsize_look, type = "cairo")
      print(plotInput())
      dev.off()
    }
  )
  
} # END



#########################   Gene fraction Start   ##########################################
df_merge <<- data.frame(empty_name=character(), stringsAsFactors=FALSE)
df_merge_avg <<- data.frame(empty_name=character(), stringsAsFactors=FALSE)


upload_dges_filter <- function(sample_file){
  
  sample_name <- sub(pattern = "(.*?)\\..*$", replacement = "\\1", basename(sample_file))
  
  df <- read.table(sample_file, header = T, row.names = 1, stringsAsFactors = F)
  
  for (g in gene_list_filter){
    
    acols <- (df[which(row.names(df) == g),] > 0)
    bcols <- !acols
    
    a <- as.data.frame(df[,acols])
    colnames(a) <- colnames(df)[acols]
    b <- as.data.frame(df[,bcols])
    colnames(b) <- colnames(df)[bcols]
    
    if(keep_filter == 1){
      df <- a
    }
    else{
      df <- b
    }
    
    
  }
  
  output_dir <- file.path(dirname(sample_list_filter[1]), "dges_filtered")
  dir.create(output_dir)
  output_file <- file.path(output_dir, paste0(sample_name, "_filtered", ".dge.txt"))
    
  write.table(df, output_file, append = FALSE, sep = "\t", dec = ".", 
              row.names = TRUE, col.names = TRUE)

  R.utils::gzip(output_file)
}

upload_dges <- function(sample_file){
  
  # Input separated dge files. Returns combined df.
  sample_name <- sub(pattern = "(.*?)\\..*$", replacement = "\\1", basename(sample_file))
  
  df <- read.table(sample_file, header = T, row.names = 1, stringsAsFactors = F)
  
  # Calculate percentage of cells expressing all genes
  df_frac <- data.frame(lapply(rowMeans(df!=0), round, 2))
  df_frac <- t(df_frac)
  colnames(df_frac) <- sample_name
  rownames(df_frac) <- rownames(df)
  
  if(length(sample_list) == 1){
    df_merge <<- df_frac
  }
  else{
    df_merge <<- transform(merge(df_merge, df_frac, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
  }
  
  # Calculate the average transcripts of each gene expressed
  df_avg <- data.frame(lapply(rowMeans(df), round, 2))
  df_avg <- t(df_avg)
  colnames(df_avg) <- paste0(sample_name, "_avg")
  rownames(df_avg) <- rownames(df)
  
  if(length(sample_list) == 1){
    df_merge_avg <<- df_avg
  }
  else{
    df_merge_avg <<- transform(merge(df_merge_avg, df_avg, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
  }
}

gene_frac_plot <- function(sample_name){
  
  #ind <- which(sample_list == sample_file)
  #sample_name <- sub(pattern = "(.*?)\\..*$", replacement = "\\1", basename(sample_file))
  #df <- data.frame(genes_look, lapply(gene_frac_list[ind], round, 2))

  df <- data.frame(genes_look, df_look[, sample_name])
  colnames(df) <- c("gene_list", "gene_frac")
  
  df$gene_list <-factor(df$gene_list, levels = genes_look)
  
  p <- ggplot(data = df, aes(x = gene_list, y = gene_frac)) + # scale_x_discrete(limits = genes_look) +
       geom_bar(stat = "identity", fill = "steelblue", width = width_bar) + coord_flip() + theme_minimal() +
       geom_text(aes(label = gene_frac), vjust = -0.5, angle = -90, size = 3.5) +
       labs(title = sample_name) + theme_bw() + # labs(title = sub("_avg$", "", sample_name))
       theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
             panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
             axis.text = element_text(size = 10)) +
       scale_y_continuous(limits = c(0, 1))
  return (p)
}

#########################   Gene fraction End   ############################################