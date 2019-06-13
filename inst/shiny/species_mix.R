
# Global varibles

source("cell_selection.R")
source("plot.R")
source("plot_species_mix.R")

# Species mix interface
mixUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "mix_tab",
    navbarPage(title = NULL,
      tabPanel("Data preprocessing",
          fluidRow(
            column(width = 3,
                   box(
                     title = "Step 1: Experiment infomation", width = NULL, solidHeader = TRUE, status = "warning",

                     selectizeInput(ns("mix_filetype"), label="Select a file type",
                                    choices = c("BAM file (*.bam)" = "bam_file_mix",
                                                "Summary Files (*.dge.summary.*)" = "sum_file_mix",
                                                "Two filtered DGE files (*.dge.*)" = "2_dge_mix",
                                                "Merged DGE file (*.dge.*)" = "1_dge_mix"
                                    ),
                                    selected = c("BAM file (*.bam)"), multiple = FALSE),
                     
                     splitLayout(
                       textInput(ns("name_human"), "Human sample name:", width = "100%"),
                       textInput(ns("name_mouse"), "Mouse sample name:", width = "100%")
                     ),
                     
                     tags$b("Selected cell barcodes file:"),
                     tags$p(),
                     verbatimTextOutput(ns("cell_bC_mix"), placeholder = TRUE),
                     shinyFilesButton(ns("file_cell_bc_mix"), "Select cell barcodes file", "Please select a BC file", multiple = FALSE),
                     tags$p(),
                     
                     
                     conditionalPanel(condition  = "input.mix_filetype == 'bam_file_mix'",
                                      tags$b("BAM file:"),
                                      tags$p(),
                                      verbatimTextOutput(ns("file_bam"), placeholder = TRUE),
                                      shinyFilesButton(ns("file_mix_bam"), "Select a bam file", "Please select a bam file", multiple = FALSE),
                                      tags$p(),
                                      
                                      tags$b("Drop-seq tool folder:"),
                                      tags$p(),
                                      verbatimTextOutput(ns("dropseq_tool_mix"), placeholder = TRUE),
                                      shinyDirButton(ns("dropseq_mix"), "Select drop-seq tool folder", "Please select a folder"),
                                      ns = NS(id)
                     ),
                     
                     conditionalPanel(condition  = "input.mix_filetype == 'sum_file_mix'",
                                      tags$b("Summary files (two files):"),
                                      tags$p(),
                                      verbatimTextOutput(ns("file_sum_human"), placeholder = TRUE),
                                      shinyFilesButton(ns("sum_mix_human"), "Select a summary file (Human)", "Please select a summary file", multiple = FALSE),
                                      tags$p(),
                                      verbatimTextOutput(ns("file_sum_mouse"), placeholder = TRUE),
                                      shinyFilesButton(ns("sum_mix_mouse"), "Select a summary file (Mouse)", "Please select a summary file", multiple = FALSE),
                                      ns = NS(id)
                     ),
                     
                     conditionalPanel(condition  = "input.mix_filetype == '2_dge_mix'",
                                      tags$b("DGE files (two files):"),
                                      tags$p(),
                                      verbatimTextOutput(ns("file_dge_human"), placeholder = TRUE),
                                      shinyFilesButton(ns("dge_mix_human"), "Select a dge file (Human)", "Please select a dge file", multiple = FALSE),
                                      tags$p(),
                                      verbatimTextOutput(ns("file_dge_mouse"), placeholder = TRUE),
                                      shinyFilesButton(ns("dge_mix_mouse"), "Select a dge file (Mouse)", "Please select a dge file", multiple = FALSE),
                                      ns = NS(id)
                     ),
                     
                     conditionalPanel(condition  = "input.mix_filetype == '1_dge_mix'",
                                      tags$b("Merged DGE file:"),
                                      tags$p(),
                                      verbatimTextOutput(ns("file_dge_merged"), placeholder = TRUE),
                                      shinyFilesButton(ns("dge_merged"), "Select a dge file", "Please select a dge file", multiple = FALSE),
                                      tags$p(),
                                      radioButtons(ns("NA_name"), label = "Replace NA with",
                                                   choices = list("Human sample name" = "1", 
                                                                  "Mouse sample name" = "2"), 
                                                   selected = "1"),
                                      checkboxInput(ns("dge_only"), "Generate dge files only? (No cell barcodes file)", FALSE),
                                      ns = NS(id)
                     )
                   ),
                   
                   box(
                     title = "Step 2: System settings (BAM file mode)", width = NULL, solidHeader = TRUE, status = "success",
                     sliderInput(ns("CPUs_mix"), "CPU number:", 1, 20, 20, step = 1),
                     sliderInput(ns("RAMs_mix"), "Memory size (G):", 10, 100, 64, step = 1),
                     sliderInput(ns("Times_mix"), "Time limit (hours):", 10, 120, 120, step = 1),
                     textInput(ns("email_mix"), "Email:", width = "100%"),
                     selectInput(ns("mtype_mix"), "Send an email when:",
                                 c("Finished" = "END",
                                   "Begin and Finished" = "ALL",
                                   "None" = "NONE"))
                   )
            ),
            
            column(width = 3,
                   box(
                     tags$p(),
                     title = "Step 3: Submit the job", width = NULL, solidHeader = TRUE, status = "danger",
                     
                     withBusyIndicatorUI(icon_name = "mix_genes",
                                         actionButton(ns("do_mix"), tags$b("Submit job"), width = "85%")
                     )
                   ),
                   
                   box(
                     title = "Step 4: Check job status (BAM file mode)", width = NULL, solidHeader = TRUE, status = "primary",
                     tags$b("BAM file folder:"),
                     tags$p(),
                     shinyFilesButton(ns("output_mix_bam"), "Select a BAM file", "Please select a BAM file", multiple = FALSE),
                     tags$p(),
                     verbatimTextOutput(ns("output_mix"), placeholder = TRUE),
                     #HTML('&nbsp&nbsp&nbsp'),
                     actionButton(ns("st_mix"), tags$b("Job status"), width = "45%"),
                     HTML('&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp'),
                     actionButton(ns("cancel_mix"), tags$b("Cancel job"), width = "45%")
                   ),
                   
                   box(
                     title = "Status", width = NULL, height = 400, solidHeader = TRUE, status = "info",
                     htmlOutput(ns("status_mix"))
                   )
            ),
            
            column(width = 6,
                   box(
                     title = "Species mixing plots",  width = NULL, height = 900, status = "primary",
                     splitLayout(
                       uiOutput(ns("pmix1")),
                       uiOutput(ns("pmix2"))
                       #plotOutput(ns("mix_plot1")),
                       #plotOutput(ns("mix_plot2"))
                     ),
                     
                     splitLayout(
                       uiOutput(ns("pmix3")),
                       uiOutput(ns("pmix4"))
                       #plotOutput(ns("mix_plot3")),
                       #plotOutput(ns("mix_plot4"))
                     )
                   )
            )
            
          )
      ), # Tab 1 END
      
      
      tabPanel('* Save figures',
               fluidRow(
                 
                 box(
                   title = "Save figure", width = 3, solidHeader = TRUE, status = "primary",
                   
                   selectizeInput(ns("psave_select_mix"), label="Choose a figure to save",
                                  choices = c("Figure 1", "Figure 2", "Figure 3", "Figure 4"),
                                  selected = c("Figure 1"),
                                  multiple=FALSE
                   ),
                   
                   numericInput(ns("pWidth_mix"), "Width (px)", value = 1200),   # Width is divided by 2 for display
                   numericInput(ns("pHeight_mix"), "Height (px)", value = 1200), # Height is divided by 2 for display
                   numericInput(ns("pRes_mix"), "Resolution (ppi)", value = 300),
                   numericInput(ns("pPsize_mix"), "Pointsize", value = 12),
                   
                   downloadButton(ns("do_psave_mix"), "Save figure", style = "width: 100%")
                 ),
                 
                 box(
                   title = "Figure to save", width = 9, status = "success",
                   
                   uiOutput(ns("figure_save_mix"))
                 )
                 
               )#
      )# Tab 2 END
      
    ) # navbarPage END
  ) # TabItem END
  
}

mixServer <- function(input, output, session, fileRoot = NULL) {
  
  ns <- session$ns

  volumes <- (c(Home = fs::path_home(), getVolumes()()))#volumes <- (c(Home = fs::path_home())) #(c("Current" = getwd())) #(c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()()))
  
  shinyFileChoose(input, "file_cell_bc_mix", roots = volumes, session = session, filetypes=c('', 'txt'))
  
  shinyFileChoose(input, "file_mix_bam", roots = volumes, session = session, filetypes=c('', 'bam'))
  shinyDirChoose(input, "dropseq_mix", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  shinyFileChoose(input, "sum_mix_human", roots = volumes, session = session)
  shinyFileChoose(input, "sum_mix_mouse", roots = volumes, session = session)
  
  shinyFileChoose(input, "dge_mix_human", roots = volumes, session = session)
  shinyFileChoose(input, "dge_mix_mouse", roots = volumes, session = session)
  
  shinyFileChoose(input, "dge_merged", roots = volumes, session = session)
  
  shinyFileChoose(input, "output_mix_bam", roots = volumes, session = session, filetypes=c('', 'bam'))

  output$cell_bC_mix <- renderPrint({
      as.character(parseFilePaths(volumes, input$file_cell_bc_mix)[4])
  })
  
  output$file_bam <- renderPrint({
      as.character(parseFilePaths(volumes, input$file_mix_bam)[4])
  })
  
  output$dropseq_tool_mix <- renderPrint({
    parseDirPath(volumes, input$dropseq_mix)
  })
  
  output$file_sum_human <- renderPrint({
      as.character(parseFilePaths(volumes, input$sum_mix_human)[4])
  })
  
  output$file_sum_mouse <- renderPrint({
      as.character(parseFilePaths(volumes, input$sum_mix_mouse)[4])
  })
  
  output$file_dge_human <- renderPrint({
      as.character(parseFilePaths(volumes, input$dge_mix_human)[4])
  })
  
  output$file_dge_mouse <- renderPrint({
      as.character(parseFilePaths(volumes, input$dge_mix_mouse)[4])
  })
  
  output$file_dge_merged <- renderPrint({
      as.character(parseFilePaths(volumes, input$dge_merged)[4])
  })
  
  output$output_mix <- renderPrint({
    as.character(parseFilePaths(volumes, input$output_mix_bam)[4])
  })

  # Submit jobs 
  observeEvent(input$do_mix, {
    
    withProgress(message = "Please wait ...",{
      
      shinyjs::show("load_mix_genes")
      
      samplename_human <<- isolate(as.character(input$name_human))
      samplename_mouse <<- isolate(as.character(input$name_mouse))
      
      choose_NA_name <- isolate(as.character(input$NA_name))
      
      if(choose_NA_name == "1"){
        samplename_NA <<- samplename_human
      }
      else{
        samplename_NA <<- samplename_mouse
      }
      
      selectedBCPath_mix <<- isolate(as.character(parseFilePaths(volumes, input$file_cell_bc_mix)[4]))
      
      bamPath_mix <<- isolate(as.character(parseFilePaths(volumes, input$file_mix_bam)[4]))
      dropseq_mixPath <<- isolate(as.character(parseDirPath(volumes, input$dropseq_mix)))
      
      sum_humanPath <<- isolate(as.character(parseFilePaths(volumes, input$sum_mix_human)[4]))
      sum_mousePath <<- isolate(as.character(parseFilePaths(volumes, input$sum_mix_mouse)[4]))
  
      dge_humanPath <<- isolate(as.character(parseFilePaths(volumes, input$dge_mix_human)[4]))
      dge_mousePath <<- isolate(as.character(parseFilePaths(volumes, input$dge_mix_mouse)[4]))
      
      dge_mergedPath <<- isolate(as.character(parseFilePaths(volumes, input$dge_merged)[4]))
  
      ncores_mix <<- isolate(as.character(input$CPUs_mix))
      mem_mix <<- isolate(as.character(input$RAMs_mix))
      timelimit_mix <<- isolate(as.character(input$Times_mix))
      mailType_mix <<- isolate(as.character(input$mtype_mix))
      Email_mix <<- isolate(as.character(input$email_mix))
      
      setProgress(value = 0.35)
      
      if(input$mix_filetype == "bam_file_mix"){
        generate_sh_mix()
        job_ids_mix <<- submit_sh_mix()     # real run
        job_sts_mix <- job_status_mix()    # real run
      
        #job_sts_mix <- "PENDING" # for debug only
      
        output$status_mix <- renderUI({
          list(
            HTML("Generating bash scripts ...<br/>"),
            HTML("-> filterBAM.sh is generated.<br/><br/>"),
            
            HTML("Submitting bash scripts ...<br/>"),
            HTML(paste("-> filterBAM.sh is ", "<b>", job_sts_mix, "</b>.<br/>", sep = ""))
          )
        })
        
        file.create(file.path(dirname(bamPath_mix), "job_ids_mix.txt"))   # save jobIDs to file
        cat(job_ids_mix, sep = ",", file = file.path(dirname(bamPath_mix), "job_ids_mix.txt"))         # real run
      }
      else{
        
        if(input$mix_filetype == "2_dge_mix"){
          two_dge_to_sum()
        }
        
        if(input$mix_filetype == "1_dge_mix"){
          merged_dge_to_2dge()
          if(input$dge_only == FALSE){
            two_dge_to_sum()
          }
        }
        
        output$pmix1 <- renderUI({
          withSpinner(plotOutput(ns("mix_plot1")), type = getOption("spinner.type", default = 4))
        })
        output$pmix2 <- renderUI({
          withSpinner(plotOutput(ns("mix_plot2")), type = getOption("spinner.type", default = 4))
        })
        output$pmix3 <- renderUI({
          withSpinner(plotOutput(ns("mix_plot3")), type = getOption("spinner.type", default = 4))
        })
        output$pmix4 <- renderUI({
          withSpinner(plotOutput(ns("mix_plot4")), type = getOption("spinner.type", default = 4))
        })
        
        output$mix_plot1 <- renderPlot({
          plot_mix_scatter_sum()
        })
        
        output$mix_plot2 <- renderPlot({
          plot_mix_hist_sum()
        })
        
        output$mix_plot3 <- renderPlot({
          plot_nTranscripts_sum()
        })
        output$mix_plot4 <- renderPlot({
          plot_nGene_sum()
        })
      }
      
      setProgress(value = 1)
      shinyjs::hide("load_mix_genes")
      shinyjs::show("check_mix_genes")
    })
  })# Submit jobs END
  
  
  # Check job status
  observeEvent(input$st_mix, {
    
    outputPath_mix <<- isolate(as.character(parseFilePaths(volumes, input$output_mix_bam)[4]))
    
    job_ids_mix <<- scan(file.path(dirname(outputPath_mix), "job_ids_mix.txt"), sep = ",", character())
    
    job_sts_mix <- job_status_mix()     # real run
    #job_sts_mix <- "RUNNING" # for debug only
    
    output$status_mix <- renderUI({
      list(
        HTML("Checking job status ...<br/><br/>"),
        HTML(paste("-> filterBAM.sh is ", "<b>", job_sts_mix, "</b>.<br/>", sep = ""))
      )
    })
    
  })# Check job status END
  
  
  # Cancel jobs
  observeEvent(input$cancel_mix, {
    
    outputPath_mix <<- isolate(as.character(parseFilePaths(volumes, input$output_mix_bam)[4]))
    
    job_ids_mix <<- scan(file.path(dirname(outputPath_mix), "job_ids_mix.txt"), sep = ",", character())
    job_cancel_mix()     # real run
    
    output$status_mix <- renderUI({
      list(
        HTML("Cancelling the job ...<br/><br/>"),
        HTML("-> filterBAM.sh has been cancelled.<br/>")
      )
    })
    
  })# Check job status END
  
  
  output$figure_save_mix <- renderUI({
    withSpinner(plotOutput(ns("save_figure_mix"), width = paste0(input$pWidth_mix/2, "px"), 
                           height = paste0(input$pHeight_mix/2, "px")),
                type = getOption("spinner.type", default = 4))
  }
  )
  
  output$save_figure_mix <- renderPlot({
    
    switch(input$psave_select_mix,
           "Figure 1" = plot_mix_scatter_sum(),
           "Figure 2" = plot_mix_hist_sum(),
           "Figure 3" = plot_nTranscripts_sum(),
           "Figure 4" = plot_nGene_sum()
    )
  })
  
  plotInput = function() {
    
    switch(input$psave_select_mix,
           "Figure 1" = {	par(mar=c(4,4,1,1))
             plot_mix_scatter_sum()},
           "Figure 2" = {par(mar=c(4,4,1,1))
             plot_mix_hist_sum()},
           "Figure 3" = plot_nTranscripts_sum(),
           "Figure 4" = plot_nGene_sum()
    )
  }
  
  output$do_psave_mix = downloadHandler(
    filename = "figure.png",
    content = function(file) {
      
      png(file, width = input$pWidth_mix, height = input$pHeight_mix,
          res = input$pRes_mix, pointsize = input$pPsize_mix, type = "cairo")
      print(plotInput())
      dev.off()
    }
  )
  
} # END


#########################   Genereate bash files Start   ##########################################
generate_sh_mix <- function(){
  
  scriptPath_mix <<- file.path(dirname(bamPath_mix), "scripts")
  dir.create(scriptPath_mix, showWarnings = FALSE)
  
  rfilePath1_mix <- file.path(getwd(), "plot.R")
  rfilePath2_mix <- file.path(getwd(), "cell_selection.R")
  
  file.copy(from = rfilePath1_mix, to = scriptPath_mix, overwrite = TRUE)
  file.copy(from = rfilePath2_mix, to = scriptPath_mix, overwrite = TRUE)
  
  #final_filename_mix <- "/error_detected.bam" # file name after Drop_seq.sh (drop-seq_tools-1.13)
  final_filename_mix <- "/final.bam"    # file name after Drop_seq.sh (drop-seq_tools-2.0.0)
  
  #filter_filename <- "/FilterBAM " # file name for bamtaghistogram.sh (drop-seq_tools-1.13)
  filter_filename <- "/FilterBam " # file name for bamtaghistogram.sh (drop-seq_tools-2.0.0)
  
  # Function to configure the cluster
  sh_head_mix <- function(filename, ncores, mem, timelimit, mailType, email){
    
    file.create(file.path(scriptPath_mix, paste(filename, ".sh", sep = "")))
    cat("#!/bin/bash \n",
        "#SBATCH --partition=general \n",
        paste("#SBATCH --job-name=", filename, "\n", sep=""), 
        paste("#SBATCH --ntasks=1 --cpus-per-task=", ncores, "\n", sep=""),
        paste("#SBATCH --mem=", mem, "g", "\n", sep=""),
        paste("#SBATCH --time=", timelimit, ":00:00", "\n", sep=""),
        paste("#SBATCH --mail-type=", mailType, "\n", sep=""),
        paste("#SBATCH --mail-user=", email, "\n\n\n\n", sep=""),
        sep = "", file = file.path(scriptPath_mix, paste(filename, ".sh", sep = "")))
  }
  
  # filterbam.sh
  sh_head_mix("filterbam", ncores_mix, mem_mix, timelimit_mix, mailType_mix, Email_mix)
  cat(paste(dropseq_mixPath, filter_filename, " \\", "\n", sep = ""),
      paste("I= ", bamPath_mix, " \\", "\n", sep = ""),
      paste("O= ", dirname(bamPath_mix), "/", samplename_human, ".bam", " \\", "\n", sep = ""),
      paste("REF_SOFT_MATCHED_RETAINED=HUMAN", "\n\n", sep = ""),
      sep = "", append = TRUE, file = file.path(scriptPath_mix, "filterbam.sh"))
  
  cat(paste(dropseq_mixPath, filter_filename, " \\", "\n", sep = ""),
      paste("I= ", bamPath_mix, " \\", "\n", sep = ""),
      paste("O= ", dirname(bamPath_mix), "/", samplename_mouse, ".bam", " \\", "\n", sep = ""),
      paste("REF_SOFT_MATCHED_RETAINED=MOUSE", "\n\n", sep = ""),
      sep = "", append = TRUE, file = file.path(scriptPath_mix, "filterbam.sh"))
  
  cat(paste(dropseq_mixPath, "/DigitalExpression ", " \\", "\n", sep = ""),
      paste("I= ", dirname(bamPath_mix), "/", samplename_human, ".bam", " \\", "\n", sep = ""),
      paste("O= ", dirname(bamPath_mix), "/", samplename_human, ".dge.txt.gz", " \\", "\n", sep = ""),
      paste("SUMMARY= ", dirname(bamPath_mix), "/", samplename_human, ".dge.summary.txt", " \\", "\n", sep = ""),
      paste("CELL_BC_FILE= ", selectedBCPath_mix, "\n\n", sep = ""),
      sep = "", append = TRUE, file = file.path(scriptPath_mix, "filterbam.sh"))
  
  cat(paste(dropseq_mixPath, "/DigitalExpression ", " \\", "\n", sep = ""),
      paste("I= ", dirname(bamPath_mix), "/", samplename_mouse, ".bam", " \\", "\n", sep = ""),
      paste("O= ", dirname(bamPath_mix), "/", samplename_mouse, ".dge.txt.gz", " \\", "\n", sep = ""),
      paste("SUMMARY= ", dirname(bamPath_mix), "/", samplename_mouse, ".dge.summary.txt", " \\", "\n", sep = ""),
      paste("CELL_BC_FILE= ", selectedBCPath_mix, "\n\n", sep = ""),
      sep = "", append = TRUE, file = file.path(scriptPath_mix, "filterbam.sh"))
  
  # 
  system(paste("dos2unix", file.path(scriptPath_mix, "filterbam.sh")))
  system(paste("chmod 777", file.path(scriptPath_mix, "filterbam.sh")))
  
}

#########################   genereate bash files End  ##########################################


#########################   Submit bash files Start   ##########################################

submit_sh_mix <- function(){
  job_mix <- system(paste("sbatch", file.path(scriptPath_mix, "filterbam.sh")), intern = TRUE)
  jobID_mix <- gsub("\\D", "", job_mix)

  return(jobID_mix)
}

#########################   Submit bash files End   ##########################################


#########################   Job status Start   ##########################################

job_status_mix <- function(){
  
  jobID_mix <- job_ids_mix
  
  jobST_mix <- system(paste("sacct -n -X -o STATE -j", jobID_mix), intern = TRUE)
  
  jobST_mix <- gsub("\\s", "", jobST_mix)

  return(jobST_mix)
}

#########################   Job status End   ##########################################


#########################   Cancel jobs Start   ##########################################

job_cancel_mix <- function(){
  
  jobID_mix <- job_ids_mix
  
  jobST_mix <- system(paste("scancel", jobID_mix), intern = TRUE)
  
}

#########################   Cancel jobs End   ##########################################
