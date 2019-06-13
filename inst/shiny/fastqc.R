
# Global varibles

# fastqc function
qc_report <- function(fastqc_dir) {
  
  qc_file <- file(fastqc_dir,open="rt")
  
  while (TRUE) {
    if(grepl(">>", readLines(qc_file, n = 1))){
      qc_sum <- read.table(qc_file, sep = '\t', nrows = 7)
      close(qc_file)
      break
    }
  }
  
  qc_sum <- qc_sum[-c(2,3),]
  qc_sum <- t(qc_sum)
  
  qc_sum <- as.data.frame(qc_sum)
  
  colnames(qc_sum) <- as.character(unlist(qc_sum[1,]))
  qc_sum <- qc_sum[-1,]
  colnames(qc_sum)[3] <- "Poor quality reads"
  
  qc_file <- file(fastqc_dir,open="rt")
  n <- 0
  x <- TRUE
  while (x) {
    if(grepl(">>Per sequence quality scores", readLines(qc_file, n = 1))){
      while(TRUE){
        n <- n + 1
        if(grepl(">>", readLines(qc_file, n = 1))){
          x <- FALSE
          close(qc_file)
          break
        }
      }
      
    }
  }
  
  qc_file <- file(fastqc_dir,open="rt")
  
  while (TRUE) {
    if(grepl(">>Per sequence quality scores", readLines(qc_file, n = 1))){
      qc <- read.table(qc_file, sep = '\t', nrows = n-2)
      close(qc_file)
      break
    }
  }
  
  total <- sum(qc$V2)
  Q10 <- sum(qc$V2[which(qc$V1 >= 10)]) / total * 100
  Q20 <- sum(qc$V2[which(qc$V1 >= 20)]) / total * 100
  Q30 <- sum(qc$V2[which(qc$V1 >= 30)]) / total * 100
  
  qc_sum$"%Q10" <- Q10
  qc_sum$"%Q20" <- Q20
  qc_sum$"%Q30" <- Q30
  
  return(qc_sum)
}
# fastqc function END


# Fastqc interface
fastqcUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "fastqc_tab",
          fluidRow(
            column(width = 2,
                   box(
                     title = "Step 1: Experiment info", width = NULL, solidHeader = TRUE, status = "warning",
                     
                     radioButtons(ns("fastqc_filetype"), NULL,
                                  c("Single-Read" = "single_read",
                                    "Paired-End" = "paired_read"
                                  ), selected = "paired_read"),
                
                     tags$b("Select sequence files:"),
                     tags$p(),

                     conditionalPanel(condition  = "input.fastqc_filetype == 'single_read'",
                                      verbatimTextOutput(ns("R1"), placeholder = TRUE),
                                      shinyFilesButton(ns("file_R1"), "Select read 1 file", "Please select a file", multiple = FALSE),
                                      tags$p(),
                                      ns = NS(id)
                     ),
                     
                     conditionalPanel(condition  = "input.fastqc_filetype == 'paired_read'",
                                      verbatimTextOutput(ns("R1_p"), placeholder = TRUE),
                                      shinyFilesButton(ns("file_R1_p"), "Select read 1 file", "Please select a file", multiple = FALSE),
                                      tags$p(),
                                      verbatimTextOutput(ns("R2"), placeholder = TRUE),
                                      shinyFilesButton(ns("file_R2"), "Select read 2 file", "Please select a file", multiple = FALSE),
                                      ns = NS(id)
                     ),
                     tags$p(),
                     tags$b("Output folder:"),
                     tags$p(),
                     verbatimTextOutput(ns("output_fastqc"), placeholder = TRUE),
                     shinyDirButton(ns("dir_fastqc"), "Select a output folder", "Please select a folder")
                   
                     ),
                   
                   box(
                     title = "Step 2: System settings", width = NULL, solidHeader = TRUE, status = "success",
                     sliderInput(ns("CPUs_fastqc"), "CPU number:", 1, 20, 20, step = 1),
                     sliderInput(ns("RAMs_fastqc"), "Memory size (G):", 10, 100, 64, step = 1),
                     sliderInput(ns("Times_fastqc"), "Time limit (hours):", 10, 120, 120, step = 1),
                     textInput(ns("email_fastqc"), "Email:", width = "100%"),
                     selectInput(ns("mtype_fastqc"), "Send an email when:",
                                 c("Finished" = "END",
                                   "Begin and Finished" = "ALL",
                                   "None" = "NONE"))
                   )
            ),
            
            column(width = 2,
                   box(
                     tags$p(),
                     title = "Step 3: Submit the job", width = NULL, solidHeader = TRUE, status = "danger",
                     
                     withBusyIndicatorUI(icon_name = "fastqc",
                                         actionButton(ns("do_fastqc"), tags$b("Submit job"), width = "75%")
                     )
                   ),
                   
                   box(
                     title = "Step 4: Check job status", width = NULL, solidHeader = TRUE, status = "primary",
                     tags$b("Output folder:"),
                     tags$p(),
                     shinyDirButton(ns("dir_fastqc_js"), "Select a output folder", "Please select a folder"),
                     tags$p(),
                     verbatimTextOutput(ns("fastqc_js"), placeholder = TRUE),
                     tags$p(),
                     actionButton(ns("st_fastqc"), tags$b("Job status"), width = "45%"),
                     HTML('&nbsp&nbsp&nbsp&nbsp'),
                     actionButton(ns("cancel_fastqc"), tags$b("Cancel job"), width = "45%")
                   ),
                   
                   box(
                     title = "Step 5: Check resutls", width = NULL, solidHeader = TRUE, status = "success",
                     
                     conditionalPanel(condition  = "input.fastqc_filetype == 'single_read'",
                                      
                                      shinyDirButton(ns("result_R1"), "Fastqc read 1 folder", "Please select a folder"),
                                      tags$p(),
                                      verbatimTextOutput(ns("R1_result"), placeholder = TRUE),
                                      tags$p(),
                                      ns = NS(id)
                     ),
                     
                     conditionalPanel(condition  = "input.fastqc_filetype == 'paired_read'",
                                      shinyDirButton(ns("result_R1_p"), "Fastqc read 1 folder", "Please select a folder"),
                                      tags$p(),
                                      verbatimTextOutput(ns("R1_result_p"), placeholder = TRUE),
                                      tags$p(),
                                      shinyDirButton(ns("result_R2"), "Fastqc read 2 folder", "Please select a folder"),
                                      tags$p(),
                                      verbatimTextOutput(ns("R2_result"), placeholder = TRUE),
                                      
                                      ns = NS(id)
                     ),
                     actionButton(ns("result_fastqc"), tags$b("Fastqc results"), width = "100%")
                   ),
                   
                   box(
                     title = "Status", width = NULL, solidHeader = TRUE, status = "info",
                     htmlOutput(ns("status_fastqc"))
                   )
            ),
            
            column(width = 8,
                   
                   box(
                     title = "Fastqc summary",  width = NULL, status = "warning",
                     tableOutput(ns("fastqc_sum"))
                   ),
                     
                   box(
                     title = "Fastqc plots",  width = NULL, height = 900, status = "success",
                     
                     tabsetPanel(type = "tabs",
                                 tabPanel("Base", 
                                          tabsetPanel(type = "pills",
                                                      tabPanel("R1", plotOutput(ns("fastqc_plot1_R1"))),
                                                      tabPanel("R2", plotOutput(ns("fastqc_plot1_R2")))
                                          )
                                 ),
                                 tabPanel("Tile",  
                                          tabsetPanel(type = "pills",
                                                      tabPanel("R1", plotOutput(ns("fastqc_plot2_R1"))),
                                                      tabPanel("R2", plotOutput(ns("fastqc_plot2_R2")))
                                          )
                                 ),
                                 tabPanel("Quality distribution",
                                          tabsetPanel(type = "pills",
                                                      tabPanel("R1", plotOutput(ns("fastqc_plot3_R1"))),
                                                      tabPanel("R2", plotOutput(ns("fastqc_plot3_R2")))
                                          )
                                 ),
                                 tabPanel("ATGC",
                                          tabsetPanel(type = "pills",
                                                      tabPanel("R1", plotOutput(ns("fastqc_plot4_R1"))),
                                                      tabPanel("R2", plotOutput(ns("fastqc_plot4_R2")))
                                          )
                                 ),
                                 tabPanel("GC",
                                          tabsetPanel(type = "pills",
                                                      tabPanel("R1", plotOutput(ns("fastqc_plot5_R1"))),
                                                      tabPanel("R2", plotOutput(ns("fastqc_plot5_R2")))
                                          )
                                 ),
                                 tabPanel("N", 
                                          tabsetPanel(type = "pills",
                                                      tabPanel("R1", plotOutput(ns("fastqc_plot6_R1"))),
                                                      tabPanel("R2", plotOutput(ns("fastqc_plot6_R2")))
                                          )
                                 ),
                                 tabPanel("bp", 
                                          tabsetPanel(type = "pills",
                                                      tabPanel("R1", plotOutput(ns("fastqc_plot7_R1"))),
                                                      tabPanel("R2", plotOutput(ns("fastqc_plot7_R2")))
                                          )
                                 ),
                                 tabPanel("Duplication ",
                                          tabsetPanel(type = "pills",
                                                      tabPanel("R1", plotOutput(ns("fastqc_plot8_R1"))),
                                                      tabPanel("R2", plotOutput(ns("fastqc_plot8_R2")))
                                          )
                                 ),
                                 tabPanel("Adapter", 
                                          tabsetPanel(type = "pills",
                                                      tabPanel("R1", plotOutput(ns("fastqc_plot9_R1"))),
                                                      tabPanel("R2", plotOutput(ns("fastqc_plot9_R2")))
                                          )
                                 ),
                                 tabPanel("Kmer", 
                                          tabsetPanel(type = "pills",
                                                      tabPanel("R1", plotOutput(ns("fastqc_plot10_R1"))),
                                                      tabPanel("R2", plotOutput(ns("fastqc_plot10_R2")))
                                          )
                                 )
                     )
                   )
            )
            
          )
  )
  
}

fastqcServer <- function(input, output, session, fileRoot = NULL) {
  
  ns <- session$ns
  
  volumes <- (c(Home = fs::path_home()))
  
  shinyFileChoose(input, "file_R1", roots = volumes, session = session)
  shinyFileChoose(input, "file_R1_p", roots = volumes, session = session)
  shinyFileChoose(input, "file_R2", roots = volumes, session = session) 
  
  shinyDirChoose(input, "dir_fastqc", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyDirChoose(input, "dir_fastqc_js", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  shinyDirChoose(input, "result_R1", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyDirChoose(input, "result_R1_p", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyDirChoose(input, "result_R2", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  output$R1 <- renderPrint({
    as.character(parseFilePaths(volumes, input$file_R1)[4])
  })
  
  output$R1_p <- renderPrint({
    as.character(parseFilePaths(volumes, input$file_R1_p)[4])
  })
  
  output$R2 <- renderPrint({
    as.character(parseFilePaths(volumes, input$file_R2)[4])
  })
  
  output$output_fastqc <- renderPrint({
    parseDirPath(volumes, input$dir_fastqc)
  })
  
  output$fastqc_js <- renderPrint({
    parseDirPath(volumes, input$dir_fastqc_js)
  })
  
  output$R1_result <- renderPrint({
    parseDirPath(volumes, input$result_R1)
  })
  
  output$R1_result_p <- renderPrint({
    parseDirPath(volumes, input$result_R1_p)
  })
  
  output$R2_result <- renderPrint({
    parseDirPath(volumes, input$result_R2)
  })
  
  # Submit jobs 
  observeEvent(input$do_fastqc, {
    
    withProgress(message = "Please wait ...",{
      
      shinyjs::show("load_fastqc")
      
      if(input$fastqc_filetype == "single_read"){
        R1Path <<- isolate(as.character(parseFilePaths(volumes, input$file_R1)[4]))
      }
      else{
        R1Path <<- isolate(as.character(parseFilePaths(volumes, input$file_R1_p)[4]))
        R2Path <<- isolate(as.character(parseFilePaths(volumes, input$file_R2)[4]))
      }
      
      qc_outPath <<- isolate(as.character(parseDirPath(volumes, input$dir_fastqc)))
      
      ncores_fastqc <<- isolate(as.character(input$CPUs_fastqc))
      mem_fastqc <<- isolate(as.character(input$RAMs_fastqc))
      timelimit_fastqc <<- isolate(as.character(input$Times_fastqc))
      mailType_fastqc <<- isolate(as.character(input$mtype_fastqc))
      Email_fastqc <<- isolate(as.character(input$email_fastqc))
      
      setProgress(value = 0.35)
      
      generate_sh_fastqc(input)
      job_ids_fastqc <<- submit_sh_fastqc()     # real run
      job_sts_fastqc <- job_status_fastqc()    # real run
      
      #job_sts_fastqc <- "PENDING" # for debug only
      
      output$status_fastqc<- renderUI({
        list(
          HTML("Generating bash scripts ...<br/>"),
          HTML("-> Fastqc.sh is generated.<br/><br/>"),
          
          HTML("Submitting bash scripts ...<br/>"),
          HTML(paste("-> Fastqc.sh is ", "<b>", job_sts_fastqc, "</b>.<br/>", sep = ""))
        )
      })
      
      file.create(file.path(qc_outPath, "job_ids_fastqc.txt"))   # save jobIDs to file
      cat(job_ids_fastqc, sep = ",", file = file.path(qc_outPath, "job_ids_fastqc.txt"))         # real run
      
      setProgress(value = 1)
      shinyjs::hide("load_fastqc")
      shinyjs::show("check_fastqc")
    })
  })# Submit jobs END
  
  
  # Check job status
  observeEvent(input$st_fastqc, {
    
    js_outPath <<- isolate(as.character(parseDirPath(volumes, input$dir_fastqc_js)))
    job_ids_fastqc <<- scan(file.path(qc_outPath, "job_ids_fastqc.txt"), sep = ",", character())
    
    job_sts_fastqc  <- job_status_fastqc()     # real run
    #job_sts_fastqc  <- "RUNNING" # for debug only
    
    output$status_fastqc  <- renderUI({
      list(
        HTML("Checking job status ...<br/><br/>"),
        HTML(paste("-> Fastqc.sh is ", "<b>", job_sts_fastqc, "</b>.<br/>", sep = ""))
      )
    })
    
  })# Check job status END
  
  
  # Cancel jobs
  observeEvent(input$cancel_fastqc, {
    
    js_outPath <<- isolate(as.character(parseDirPath(volumes, input$dir_fastqc_js)))
    job_ids_fastqc <<- scan(file.path(qc_outPath, "job_ids_fastqc.txt"), sep = ",", character())
    
    job_cancel_fastqc()     # real run
    
    output$status_fastqc <- renderUI({
      list(
        HTML("Cancelling the job ...<br/><br/>"),
        HTML("-> Fastqc.sh has been cancelled.<br/>")
      )
    })
    
  })# Check job status END
  
  # Show results
  observeEvent(input$result_fastqc, {
    
    if(input$fastqc_filetype == "single_read"){
      R1_resultPath <<- isolate(as.character(parseDirPath(volumes, input$result_R1)))
      
      shinyjs::hide("fastqc_plot1_R2")
      shinyjs::hide("fastqc_plot2_R2")
      shinyjs::hide("fastqc_plot3_R2")
      shinyjs::hide("fastqc_plot4_R2")
      shinyjs::hide("fastqc_plot5_R2")
      shinyjs::hide("fastqc_plot6_R2")
      shinyjs::hide("fastqc_plot7_R2")
      shinyjs::hide("fastqc_plot8_R2")
      shinyjs::hide("fastqc_plot9_R2")
      shinyjs::hide("fastqc_plot10_R2")
    }
    else{
      R1_resultPath <<- isolate(as.character(parseDirPath(volumes, input$result_R1_p)))
      R2_resultPath <<- isolate(as.character(parseDirPath(volumes, input$result_R2)))
      
      shinyjs::show("fastqc_plot1_R2")
      shinyjs::show("fastqc_plot2_R2")
      shinyjs::show("fastqc_plot3_R2")
      shinyjs::show("fastqc_plot4_R2")
      shinyjs::show("fastqc_plot5_R2")
      shinyjs::show("fastqc_plot6_R2")
      shinyjs::show("fastqc_plot7_R2")
      shinyjs::show("fastqc_plot8_R2")
      shinyjs::show("fastqc_plot9_R2")
      shinyjs::show("fastqc_plot10_R2")
    }
    
    if(input$fastqc_filetype == "single_read" || input$fastqc_filetype == "paired_read"){
      
      # Table
      
      R1_data <- file.path(R1_resultPath, "fastqc_data.txt") 
      qc_summary <- qc_report(R1_data)
      
      # Images
      imgPath_R1 <- R1_resultPath
      imgPath_R1 <- file.path(imgPath_R1, "Images")
      
      output$fastqc_plot1_R1 <- renderImage({
        validate(need(imgPath_R1 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R1, paste0("per_base_quality", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot2_R1 <- renderImage({
        validate(need(imgPath_R1 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R1, paste0("per_tile_quality", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot3_R1 <- renderImage({
        validate(need(imgPath_R1 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R1, paste0("per_sequence_quality", ".png")))
      }, deleteFile = FALSE)

      output$fastqc_plot4_R1 <- renderImage({
        validate(need(imgPath_R1 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R1, paste0("per_base_sequence_content", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot5_R1 <- renderImage({
        validate(need(imgPath_R1 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R1, paste0("per_sequence_gc_content", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot6_R1 <- renderImage({
        validate(need(imgPath_R1 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R1, paste0("per_base_n_content", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot7_R1 <- renderImage({
        validate(need(imgPath_R1 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R1, paste0("sequence_length_distribution", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot8_R1 <- renderImage({
        validate(need(imgPath_R1 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R1, paste0("duplication_levels", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot9_R1 <- renderImage({
        validate(need(imgPath_R1 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R1, paste0("adapter_content", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot10_R1 <- renderImage({
        validate(need(imgPath_R1 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R1, paste0("kmer_profiles", ".png")))
      }, deleteFile = FALSE)
    }
    if(input$fastqc_filetype == "paired_read"){

      # Table
      R2_data <- file.path(R2_resultPath, "fastqc_data.txt")
      qc_R2_sum <- qc_report(R2_data)
      qc_summary <- rbind(qc_summary, qc_R2_sum)
      
      # Images
      imgPath_R2 <- R2_resultPath
      imgPath_R2 <- file.path(imgPath_R2, "Images")
      
      output$fastqc_plot1_R2 <- renderImage({
        validate(need(imgPath_R2 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R2, paste0("per_base_quality", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot2_R2 <- renderImage({
        validate(need(imgPath_R2 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R2, paste0("per_tile_quality", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot3_R2 <- renderImage({
        validate(need(imgPath_R2 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R2, paste0("per_sequence_quality", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot4_R2 <- renderImage({
        validate(need(imgPath_R2 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R2, paste0("per_base_sequence_content", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot5_R2 <- renderImage({
        validate(need(imgPath_R2 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R2, paste0("per_sequence_gc_content", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot6_R2 <- renderImage({
        validate(need(imgPath_R2 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R2, paste0("per_base_n_content", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot7_R2 <- renderImage({
        validate(need(imgPath_R2 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R2, paste0("sequence_length_distribution", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot8_R2 <- renderImage({
        validate(need(imgPath_R2 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R2, paste0("duplication_levels", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot9_R2 <- renderImage({
        validate(need(imgPath_R2 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R2, paste0("adapter_content", ".png")))
      }, deleteFile = FALSE)
      
      output$fastqc_plot10_R2 <- renderImage({
        validate(need(imgPath_R2 != "0", "Error: Please select the output folder!"))
        list(src = file.path(imgPath_R2, paste0("kmer_profiles", ".png")))
      }, deleteFile = FALSE)
      
    }
    
    output$fastqc_sum <- renderTable(align = "c", rownames = FALSE,
                                     striped = TRUE, hover = TRUE, bordered = TRUE,
                                     { qc_summary }
    )
    
  })# Show results END
  
} # END


#########################   Genereate bash files Start   ##########################################
generate_sh_fastqc <- function(input){
  
  scriptPath_fastqc <<- file.path(qc_outPath, "scripts")
  dir.create(scriptPath_fastqc, showWarnings = FALSE)

  # Function to configure the cluster
  sh_head_fastqc <- function(filename, ncores, mem, timelimit, mailType, email){
    
    file.create(file.path(scriptPath_fastqc, paste(filename, ".sh", sep = "")))
    cat("#!/bin/bash \n",
        "#SBATCH --partition=general \n",
        paste("#SBATCH --job-name=", filename, "\n", sep=""), 
        paste("#SBATCH --ntasks=1 --cpus-per-task=", ncores, "\n", sep=""),
        paste("#SBATCH --mem=", mem, "g", "\n", sep=""),
        paste("#SBATCH --time=", timelimit, ":00:00", "\n", sep=""),
        paste("#SBATCH --mail-type=", mailType, "\n", sep=""),
        paste("#SBATCH --mail-user=", email, "\n\n\n\n", sep=""),
        sep = "", file = file.path(scriptPath_fastqc, paste(filename, ".sh", sep = "")))
  }
  
  # Fastqc.sh
  sh_head_fastqc("Fastqc", ncores_fastqc, mem_fastqc, timelimit_fastqc, mailType_fastqc, Email_fastqc)
  
  if(input$fastqc_filetype == "single_read"){
  
    cat(paste0("module load FastQC", "\n"),
        paste0("fastqc ", "-o ", qc_outPath, " --extract", " \\", "\n"),
        paste0(R1Path),
        sep = "", append = TRUE, file = file.path(scriptPath_fastqc, "Fastqc.sh"))
  }
  else{
    cat(paste0("module load FastQC", "\n"),
        paste0("fastqc ", "-o ", qc_outPath, " --extract", " \\", "\n"),
        paste(R1Path, R2Path),
        sep = "", append = TRUE, file = file.path(scriptPath_fastqc, "Fastqc.sh"))
  }
  
  # 
  system(paste("dos2unix", file.path(scriptPath_fastqc, "Fastqc.sh")))
  system(paste("chmod 777", file.path(scriptPath_fastqc, "Fastqc.sh")))
  
}

#########################   genereate bash files End  ##########################################


#########################   Submit bash files Start   ##########################################

submit_sh_fastqc <- function(){
  job_fastqc <- system(paste("sbatch", file.path(scriptPath_fastqc, "Fastqc.sh")), intern = TRUE)
  jobID_fastqc <- gsub("\\D", "", job_fastqc)
  
  return(jobID_fastqc)
}

#########################   Submit bash files End   ##########################################


#########################   Job status Start   ##########################################

job_status_fastqc <- function(){
  
  jobID_fastqc <- job_ids_fastqc
  
  jobST_fastqc <- system(paste("sacct -n -X -o STATE -j", jobID_fastqc), intern = TRUE)
  
  jobST_fastqc <- gsub("\\s", "", jobST_fastqc)
  
  return(jobST_fastqc)
}

#########################   Job status End   ##########################################


#########################   Cancel jobs Start   ##########################################

job_cancel_fastqc <- function(){
  
  jobID_fastqc <- job_ids_fastqc
  
  jobST_fastqc <- system(paste("scancel", jobID_fastqc), intern = TRUE)
  
}

#########################   Cancel jobs End   ##########################################