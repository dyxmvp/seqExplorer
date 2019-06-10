
# Global varibles

samplename = ""
cellnumber = 0
seqPath = ""
outputPath = "" 
refGenPath = ""

picardPath = ""
starPath = ""
drop_seqPath = ""

ncores = ""
mem = ""
timelimit = ""
mailType = ""
Email = ""

scriptPath = ""

job_ids = ""



# DGE interface
dgeUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "dge_tab",
          fluidRow(
            column(width = 4,
                   box(
                     title = "Step 1: Experiment infomation", width = NULL, solidHeader = TRUE, status = "info",
                     textInput(ns("sampleName"), "Sample name:", width = "100%"),
                     textInput(ns("cellNumber"), "Estimated cell number:", width = "100%")
                   ),
                   
                   box(
                     title = "Step 2: File directories", width = NULL, solidHeader = TRUE, status = "warning",
                     tags$b("Sequencing data folder:"),
                     tags$p(),
                     verbatimTextOutput(ns("seqData_path"), placeholder = TRUE),
                     shinyDirButton(ns("seqData"), "Select sequencing data folder", "Please select a folder"),
                     tags$p(),
                     #tags$hr(),
                     
                     tags$b("Output folder:"),
                     tags$p(),
                     verbatimTextOutput(ns("output_path"), placeholder = TRUE),
                     shinyDirButton(ns("output"), "Select output folder", "Please select a folder"),
                     tags$p(),
                     
                     tags$b("Reference genome folder:"),
                     tags$p(),
                     verbatimTextOutput(ns("refGen_path"), placeholder = TRUE),
                     shinyDirButton(ns("refGen"), "Select reference genome folder", "Please select a folder"),
                     tags$p()
                   ),
                   
                   box(
                     title = "Step 3: Software directories", width = NULL, solidHeader = TRUE, status = "primary",
                     #tags$p(),
                     tags$b("Picard folder:"),
                     tags$p(),
                     verbatimTextOutput(ns("picard_path"), placeholder = TRUE),
                     shinyDirButton(ns("picard"), "Select picard folder", "Please select a folder"),
                     tags$p(),
                     #tags$hr(),
                     
                     tags$b("Star folder:"),
                     tags$p(),
                     verbatimTextOutput(ns("star_path"), placeholder = TRUE),
                     shinyDirButton(ns("star"), "Select star folder", "Please select a folder"),
                     tags$p(),
                     
                     tags$b("Drop-seq tool folder:"),
                     tags$p(),
                     verbatimTextOutput(ns("drop_seq_path"), placeholder = TRUE),
                     shinyDirButton(ns("drop_seq"), "Select drop-seq tool folder", "Please select a folder"),
                     tags$p()
                   )
            ),
            
            column(width = 4,
                   box(
                     title = "Step 4: System settings", width = NULL, solidHeader = TRUE, status = "success",
                     sliderInput(ns("CPUs"), "CPU number:", 1, 20, 20, step = 1),
                     sliderInput(ns("RAMs"), "Memory size (G):", 10, 100, 64, step = 1),
                     sliderInput(ns("Times"), "Time limit (hours):", 10, 120, 120, step = 1),
                     textInput(ns("email"), "Email:", width = "100%"),
                     selectInput(ns("mtype"), "Send an email when:",
                                 c("Finished" = "END",
                                   "Begin and Finished" = "ALL",
                                   "None" = "NONE"),
                                 tableOutput("data"))
                   ),
                   
                   box(
                     tags$p(),
                     title = "Step 5: Submit the job", width = NULL, solidHeader = TRUE, status = "danger",
                     HTML('&nbsp&nbsp&nbsp&nbsp'),
                     actionButton(ns("do"), tags$b("Submit jobs")),
                     HTML('&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp'),
                     actionButton(ns("st"), tags$b("Job status")),
                     HTML('&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp'),
                     actionButton(ns("cancel"), tags$b("Cancel jobs")),
                     HTML('&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp'),
                     actionButton(ns("results"), tags$b("Show results"))
                   ),
                   
                   box(
                     title = "Status", width = NULL, height = 450, solidHeader = TRUE, status = "info",
                     htmlOutput(ns("status"))
                   )
            ),
            
            column(width = 4,
                   box(
                     title = "Cumulative distribution of reads",  width = NULL, height = 550, status = "primary",
                     plotOutput(ns("cell_select"))
                   ),
                   box(
                     title = "Generate dge with updated cell number", width = NULL, solidHeader = TRUE, status = "warning",
                     tags$b("1. Input sample name in Step 1;"),
                     tags$p(),
                     tags$b("2. Update cell number in Step 1;"),
                     tags$p(),
                     tags$b("3. Select the output folder in Step 2;"),
                     tags$p(),
                     tags$b("4. Select Drop-seq tool folder folder in Step 3;"),
                     tags$p(),
                     tags$b("5. Set system parameters in Step 4;"),
                     tags$p(),
                     actionButton(ns("do_new"), tags$b("6. Submit a new job"), style = "width: 100%")
                   )
            )
          )
  )
  
}

dgeServer <- function(input, output, session, fileRoot = NULL) {
  
  #ns <- session$ns (c("Root" = dirname(getwd())))
  volumes <- (c(Home = fs::path_home()))#(c("Root" = dirname(getwd()))) #(c("Current" = getwd()))
  shinyDirChoose(input, "seqData", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyDirChoose(input, "output", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyDirChoose(input, "refGen", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyDirChoose(input, "picard", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyDirChoose(input, "star", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyDirChoose(input, "drop_seq", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  
  output$seqData_path <- renderPrint({
    parseDirPath(volumes, input$seqData)
  })
  
  output$output_path <- renderPrint({
    parseDirPath(volumes, input$output)
  })
  
  output$refGen_path <- renderPrint({
    parseDirPath(volumes, input$refGen)
  })
  
  output$picard_path <- renderPrint({
    parseDirPath(volumes, input$picard)
  })
  
  output$star_path <- renderPrint({
    parseDirPath(volumes, input$star)
  })
  
  output$drop_seq_path <- renderPrint({
    parseDirPath(volumes, input$drop_seq)
  })
  
  # Submit jobs 
  observeEvent(input$do, {
    samplename <<- isolate(as.character(input$sampleName))
    cellnumber <<- isolate(as.character(input$cellNumber))
    seqPath <<- isolate(as.character(parseDirPath(volumes, input$seqData)))
    outputPath <<- isolate(as.character(parseDirPath(volumes, input$output)))
    refGenPath <<- isolate(as.character(parseDirPath(volumes, input$refGen)))
    
    picardPath <<- isolate(as.character(parseDirPath(volumes, input$picard)))
    starPath <<- isolate(as.character(parseDirPath(volumes, input$star)))
    drop_seqPath <<- isolate(as.character(parseDirPath(volumes, input$drop_seq)))
    
    ncores <<- isolate(as.character(input$CPUs))
    mem <<- isolate(as.character(input$RAMs))
    timelimit <<- isolate(as.character(input$Times))
    mailType <<- isolate(as.character(input$mtype))
    Email <<- isolate(as.character(input$email))
    
    generate_sh()
    job_ids <<- submit_sh()
    job_sts <- job_status()
    
    ###job_sts <- c("PENDING", "PENDING", "PENDING", "PENDING") # for debug only
    
    output$status <- renderUI({
      list(
        HTML("Generating bash scripts ...<br/>"),
        HTML("-> fastqtoSam.sh is generated.<br/>"),
        HTML("-> Drop_seq.sh is generated.<br/>"),
        HTML("-> bamtaghistogram.sh is generated.<br/>"),
        HTML("-> digitalexpression.sh is generated.<br/><br/>"),
        
        HTML("Submitting bash scripts ...<br/>"),
        HTML(paste("-> fastqtoSam.sh is ", job_sts[1], ".<br/>", sep = "")),
        HTML(paste("-> Drop_seq.sh is ", job_sts[2], ".<br/>", sep = "")),
        HTML(paste("-> bamtaghistogram.sh is ", job_sts[3], ".<br/>", sep = "")),
        HTML(paste("-> digitalexpression.sh is ", job_sts[4], ".<br/>", sep = ""))
      )
    })
    
    file.create(file.path(outputPath, "job_ids.txt"))   # save jobIDs to file
    cat(job_ids, sep = ",", file = file.path(outputPath, "job_ids.txt"))
    
  })# Submit jobs END
  
  
  # Check job status
  observeEvent(input$st, {
    
    outputPath <<- isolate(as.character(parseDirPath(volumes, input$output)))
    job_ids <<- scan(file.path(outputPath, "job_ids.txt"), sep = ",", character())
    
    job_sts <- job_status()
    ###job_sts <- c("RUNNING", "PENDING", "PENDING", "PENDING") # for debug only
    
    output$status <- renderUI({
      list(
        HTML("Checking job status ...<br/><br/>"),
        HTML(paste("-> fastqtoSam.sh is ", "<b>", job_sts[1], "</b>.<br/>", sep = "")),
        HTML(paste("-> Drop_seq.sh is ", "<b>", job_sts[2], "</b>.<br/>", sep = "")),
        HTML(paste("-> bamtaghistogram.sh is ", "<b>", job_sts[3], "</b>.<br/>", sep = "")),
        HTML(paste("-> digitalexpression.sh is ", "<b>", job_sts[4], "</b>.<br/>", sep = ""))
      )
    })
    
  })# Check job status END
  
  
  # Cancel jobs
  observeEvent(input$cancel, {
    
    outputPath <<- isolate(as.character(parseDirPath(volumes, input$output)))
    job_ids <<- scan(file.path(outputPath, "job_ids.txt"), sep = ",", character())
    
    job_cancel()
    
    output$status <- renderUI({
      list(
        HTML("Cancelling ALL jobs ...<br/><br/>"),
        HTML("-> fastqtoSam.sh has been cancelled.<br/>"),
        HTML("-> Drop_seq.sh has been cancelled.<br/>"),
        HTML("-> bamtaghistogram.sh has been cancelled.<br/>"),
        HTML("-> digitalexpression.sh has been cancelled.<br/>")
      )
    })
    
  })# Check job status END
  
  
  # Show results
  observeEvent(input$results, {
    
    imgPath <- as.character(input$output)
    
    output$cell_select <- renderImage({
      validate(need(imgPath != "0", "Error: Please select the output folder!"))
      outputPath <<- isolate(as.character(parseDirPath(volumes, input$output)))
      list(src = file.path(outputPath, "selectCells.tif"))
    }, deleteFile = FALSE)
    
  })# Show results END
  
  # New cell number
  observeEvent(input$do_new, {
    
    cellnumber <<- isolate(as.character(input$cellNumber))
    
    generate_sh_new()
    job_ids <<- submit_sh_new()
    job_sts <- job_status_new()
    
    ###job_sts <- c("PENDING", "PENDING", "PENDING", "PENDING") # for debug only
    
    output$status <- renderUI({
      list(
        HTML("Generating bash scripts ...<br/>"),
        HTML("-> bamtaghistogram.sh is generated.<br/>"),
        HTML("-> digitalexpression.sh is generated.<br/><br/>"),
        
        HTML("Submitting bash scripts ...<br/>"),
        HTML(paste("-> bamtaghistogram.sh is ", job_sts[1], ".<br/>", sep = "")),
        HTML(paste("-> digitalexpression.sh is ", job_sts[2], ".<br/>", sep = ""))
      )
    })
    
    file.create(file.path(outputPath, "job_ids.txt"))   # save jobIDs to file
    cat(job_ids, sep = ",", file = file.path(outputPath, "job_ids.txt"))
    
  })# New cell number END
  
} # END


#########################   Genereate bash files Start   ##########################################
generate_sh <- function(){

scriptPath <<- file.path(outputPath, "scripts")
dir.create(scriptPath, showWarnings = FALSE)

rfilePath1 <- file.path(getwd(), "selectCells.R")
rfilePath2 <- file.path(getwd(), "cell_selection.R")

file.copy(from = rfilePath1, to = scriptPath, overwrite = TRUE)
file.copy(from = rfilePath2, to = scriptPath, overwrite = TRUE)

#final_filename <- "/error_detected.bam" # file name after Drop_seq.sh (drop-seq_tools-1.13)
final_filename <- "/final.bam"    # file name after Drop_seq.sh (drop-seq_tools-2.0.0)
  
# Function to configure the cluster
sh_head <- function(filename, ncores, mem, timelimit, mailType, email){

file.create(file.path(scriptPath, paste(filename, ".sh", sep = "")))
cat("#!/bin/bash
#SBATCH --partition=general \n",
    paste("#SBATCH --job-name=", filename, "\n", sep=""), 
    paste("#SBATCH --ntasks=1 --cpus-per-task=", ncores, "\n", sep=""),
    paste("#SBATCH --mem=", mem, "g", "\n", sep=""),
    paste("#SBATCH --time=", timelimit, ":00:00", "\n", sep=""),
    paste("#SBATCH --mail-type=", mailType, "\n", sep=""),
    paste("#SBATCH --mail-user=", email, "\n\n\n\n", sep=""),
    sep = "", file = file.path(scriptPath, paste(filename, ".sh", sep = "")))
}
  
# fastqtoSam.sh
r1Path <- system(paste("find", seqPath, "-type f -print | grep '1.fastq\\|1.fq'"), intern = TRUE)
r2Path <- system(paste("find", seqPath, "-type f -print | grep '2.fastq\\|2.fq'"), intern = TRUE)

sh_head("fastqToSam", ncores, mem, timelimit, mailType, Email)
cat(paste("java -jar ", picardPath, "/picard.jar", " \\", "\n", sep = ""),
    paste("FastqToSam \\", "\n", sep = ""),
    paste("FASTQ=", r1Path, " \\", "\n", sep = ""),
    paste("FASTQ2=", r2Path, " \\", "\n", sep = ""),
    paste("OUTPUT=", outputPath, "/", samplename, ".bam", " \\", "\n", sep = ""),
    paste("SAMPLE_NAME=", samplename, sep = ""),
    sep = "", append = TRUE, file = file.path(scriptPath, "fastqToSam.sh"))

# Drop_seq.sh (Note: need to add "--runThreadN 20" in the drop-seq_tools-2.0.0)
sh_head("Drop_seq", ncores, mem, timelimit, mailType, Email)
cat(paste(drop_seqPath, "/Drop-seq_alignment.sh ", "-o ", outputPath, " \\", "\n", sep = ""),
    paste("-t ", outputPath, " \\", "\n", sep = ""),
    paste("-g ", refGenPath, " \\", "\n", sep = ""),
    paste("-r ", refGenPath, "/", list.files(path = refGenPath, pattern = "\\.fasta$"), " \\", "\n", sep = ""),
    paste("-d ", drop_seqPath, " \\", "\n", sep = ""),
    paste("-s ", starPath, "/bin/Linux_x86_64/STAR", " \\", "\n", sep = ""),
    #paste("-n ", cellnumber, " \\", "\n", sep = ""),      # drop-seq_tools-1.13, comment this line in drop-seq_tools-2.0.0
    paste(outputPath, "/", samplename, ".bam", sep = ""),
    sep = "", append = TRUE, file = file.path(scriptPath, "Drop_seq.sh"))

# bamtaghistogram.sh.
#hist_filename <- "/BAMTagHistogram " # file name for bamtaghistogram.sh (drop-seq_tools-1.13)
hist_filename <- "/BamTagHistogram " # file name for bamtaghistogram.sh (drop-seq_tools-2.0.0)

sh_head("bamtaghistogram", ncores, mem, timelimit, mailType, Email)
cat(paste(drop_seqPath, hist_filename, " \\", "\n", sep = ""),
    paste("I= ", outputPath, final_filename, " \\", "\n", sep = ""),
    paste("O= ", outputPath, "/", samplename, ".readcounts.txt.gz", " \\", "\n", sep = ""),
    paste("TAG=XC", "\n\n", sep = ""),
    paste("module load R", "\n", sep = ""),
    paste("Rscript", file.path(scriptPath, "selectCells.R"), cellnumber, samplename, scriptPath),
    sep = "", append = TRUE, file = file.path(scriptPath, "bamtaghistogram.sh"))

# digitalexpression.sh
sh_head("digitalexpression", ncores, mem, timelimit, mailType, Email)
cat(paste(drop_seqPath, "/DigitalExpression ", " \\", "\n", sep = ""),
    paste("I= ", outputPath, final_filename, " \\", "\n", sep = ""),
    paste("O= ", outputPath, "/", samplename, ".dge.txt.gz", " \\", "\n", sep = ""),
    paste("SUMMARY= ", outputPath, "/", samplename, ".dge.summary.txt", " \\", "\n", sep = ""),
    paste("CELL_BC_FILE= ", outputPath, "/", samplename, "_selectedCellBarcodes.txt", sep = ""),
    sep = "", append = TRUE, file = file.path(scriptPath, "digitalexpression.sh"))

# 
system(paste("dos2unix", file.path(scriptPath, "fastqToSam.sh")))
system(paste("chmod 777", file.path(scriptPath, "fastqToSam.sh")))

system(paste("dos2unix", file.path(scriptPath, "Drop_seq.sh")))
system(paste("chmod 777", file.path(scriptPath, "Drop_seq.sh")))

system(paste("dos2unix", file.path(scriptPath, "bamtaghistogram.sh")))
system(paste("chmod 777", file.path(scriptPath, "bamtaghistogram.sh")))

system(paste("dos2unix", file.path(scriptPath, "digitalexpression.sh")))
system(paste("chmod 777", file.path(scriptPath, "digitalexpression.sh")))

}

#########################   genereate bash files End  ##########################################


#########################   Submit bash files Start   ##########################################

submit_sh <- function(){
  job_1 <- system(paste("sbatch", file.path(scriptPath, "fastqToSam.sh")), intern = TRUE)
  jobID_1 <- gsub("\\D", "", job_1)
  
  job_2 <- system(paste(paste("sbatch --dependency=afterok:", jobID_1, sep=''),
                        file.path(scriptPath, "Drop_seq.sh")), intern = TRUE)
  jobID_2 <- gsub("\\D", "", job_2)
  
  job_3 <- system(paste(paste("sbatch --dependency=afterok:", jobID_2, sep=''),
                        file.path(scriptPath, "bamtaghistogram.sh")), intern = TRUE)
  jobID_3 <- gsub("\\D", "", job_3)
  
  job_4 <- system(paste(paste("sbatch --dependency=afterok:", jobID_3, sep=''),
                        file.path(scriptPath, "digitalexpression.sh")), intern = TRUE)
  jobID_4 <- gsub("\\D", "", job_4)
  
  return(c(jobID_1, jobID_2, jobID_3, jobID_4))
}

#########################   Submit bash files End   ##########################################


#########################   Job status Start   ##########################################

job_status <- function(){
  
  jobID_1 <- job_ids[1]
  jobID_2 <- job_ids[2]
  jobID_3 <- job_ids[3]
  jobID_4 <- job_ids[4]
  
  jobST_1 <- system(paste("sacct -n -X -o STATE -j", jobID_1), intern = TRUE)
  jobST_2 <- system(paste("sacct -n -X -o STATE -j", jobID_2), intern = TRUE)
  jobST_3 <- system(paste("sacct -n -X -o STATE -j", jobID_3), intern = TRUE)
  jobST_4 <- system(paste("sacct -n -X -o STATE -j", jobID_4), intern = TRUE)
  
  jobST_1 <- gsub("\\s", "", jobST_1)
  jobST_2 <- gsub("\\s", "", jobST_2)
  jobST_3 <- gsub("\\s", "", jobST_3)
  jobST_4 <- gsub("\\s", "", jobST_4)
  
  #jobST_1 <- system(paste("squeue -j", jobID_1, "-o \"%%T\""), intern = TRUE)[2]
  #jobST_2 <- system(paste("squeue -j", jobID_2, "-o \"%%T\""), intern = TRUE)[2]
  #jobST_3 <- system(paste("squeue -j", jobID_3, "-o \"%%T\""), intern = TRUE)[2]
  #jobST_4 <- system(paste("squeue -j", jobID_4, "-o \"%%T\""), intern = TRUE)[2]
  
  return(c(jobST_1, jobST_2, jobST_3, jobST_4))
}

#########################   Job status End   ##########################################


#########################   Cancel jobs Start   ##########################################

job_cancel <- function(){
  
  jobID_1 <- job_ids[1]
  jobID_2 <- job_ids[2]
  jobID_3 <- job_ids[3]
  jobID_4 <- job_ids[4]
  
  jobST_1 <- system(paste("scancel", jobID_1), intern = TRUE)
  jobST_2 <- system(paste("scancel", jobID_2), intern = TRUE)
  jobST_3 <- system(paste("scancel", jobID_3), intern = TRUE)
  jobST_4 <- system(paste("scancel", jobID_4), intern = TRUE)
  
}

#########################   Cancel jobs End   ##########################################


#########################   Genereate bash files with new cell number Start   ##########################################
generate_sh_new <- function(){
  
  scriptPath <<- file.path(outputPath, "scripts")
  dir.create(scriptPath, showWarnings = FALSE)
  
  #final_filename <- "/error_detected.bam" # file name after Drop_seq.sh (drop-seq_tools-1.13)
  final_filename <- "/final.bam"    # file name after Drop_seq.sh (drop-seq_tools-2.0.0)
  
  # Function to configure the cluster
  sh_head <- function(filename, ncores, mem, timelimit, mailType, email){
    
    file.create(file.path(scriptPath, paste(filename, ".sh", sep = "")))
    cat("#!/bin/bash
        #SBATCH --partition=general \n",
        paste("#SBATCH --job-name=", filename, "\n", sep=""), 
        paste("#SBATCH --ntasks=1 --cpus-per-task=", ncores, "\n", sep=""),
        paste("#SBATCH --mem=", mem, "g", "\n", sep=""),
        paste("#SBATCH --time=", timelimit, ":00:00", "\n", sep=""),
        paste("#SBATCH --mail-type=", mailType, "\n", sep=""),
        paste("#SBATCH --mail-user=", email, "\n", sep=""),
        paste("#SBATCH -o $s.%%j.out\n", sep=""),
        paste("#SBATCH -e $s.%%j.err", "\n\n\n\n", sep=""),
        sep = "", file = file.path(scriptPath, paste(filename, ".sh", sep = "")))
  }
  
  # bamtaghistogram.sh
  #hist_filename <- "/BAMTagHistogram " # file name for bamtaghistogram.sh (drop-seq_tools-1.13)
  hist_filename <- "/BamTagHistogram " # file name for bamtaghistogram.sh (drop-seq_tools-2.0.0)
  
  sh_head("bamtaghistogram", ncores, mem, timelimit, mailType, Email)
  cat(paste(drop_seqPath, hist_filename, " \\", "\n", sep = ""),
      paste("I= ", outputPath, final_filename, " \\", "\n", sep = ""),
      paste("O= ", outputPath, "/", samplename, ".readcounts.txt.gz", " \\", "\n", sep = ""),
      paste("TAG=XC", "\n\n", sep = ""),
      paste("module load R", "\n", sep = ""),
      paste("Rscript", file.path(scriptPath, "selectCells.R"), cellnumber, samplename, scriptPath),
      sep = "", append = TRUE, file = file.path(scriptPath, "bamtaghistogram.sh"))
  
  # digitalexpression.sh
  sh_head("digitalexpression", ncores, mem, timelimit, mailType, Email)
  cat(paste(drop_seqPath, "/DigitalExpression ", " \\", "\n", sep = ""),
      paste("I= ", outputPath, final_filename, " \\", "\n", sep = ""),
      paste("O= ", outputPath, "/", samplename, ".dge.txt.gz", " \\", "\n", sep = ""),
      paste("SUMMARY= ", outputPath, "/", samplename, ".dge.summary.txt", " \\", "\n", sep = ""),
      paste("CELL_BC_FILE= ", outputPath, "/", samplename, "_selectedCellBarcodes.txt", sep = ""),
      sep = "", append = TRUE, file = file.path(scriptPath, "digitalexpression.sh"))
  
  # 
  system(paste("dos2unix", file.path(scriptPath, "bamtaghistogram.sh")))
  system(paste("chmod 777", file.path(scriptPath, "bamtaghistogram.sh")))
  
  system(paste("dos2unix", file.path(scriptPath, "digitalexpression.sh")))
  system(paste("chmod 777", file.path(scriptPath, "digitalexpression.sh")))
  
}

#########################   genereate bash files with new cell number End  ##########################################


#########################   Submit new bash files Start   ##########################################

submit_sh_new <- function(){
  
  job_1 <- system(paste("sbatch", file.path(scriptPath, "bamtaghistogram.sh")), intern = TRUE)
  jobID_1 <- gsub("\\D", "", job_1)
  
  job_2 <- system(paste(paste("sbatch --dependency=afterok:", jobID_1, sep=''),
                        file.path(scriptPath, "digitalexpression.sh")), intern = TRUE)
  jobID_2 <- gsub("\\D", "", job_2)
  
  return(c(jobID_1, jobID_2))
}

#########################   Submit new bash files End   ##########################################


#########################   Job status new Start   ##########################################

job_status_new <- function(){
  
  jobID_1 <- job_ids[1]
  jobID_2 <- job_ids[2]
  
  jobST_1 <- system(paste("sacct -n -X -o STATE -j", jobID_1), intern = TRUE)
  jobST_2 <- system(paste("sacct -n -X -o STATE -j", jobID_2), intern = TRUE)
  
  jobST_1 <- gsub("\\s", "", jobST_1)
  jobST_2 <- gsub("\\s", "", jobST_2)
  
  return(c(jobST_1, jobST_2))
}

#########################   Job status new End   ##########################################


#########################   Cancel jobs new Start   ##########################################

job_cancel_new <- function(){
  
  jobID_1 <- job_ids[1]
  jobID_2 <- job_ids[2]
  
  jobST_1 <- system(paste("scancel", jobID_1), intern = TRUE)
  jobST_2 <- system(paste("scancel", jobID_2), intern = TRUE)
  
}

#########################   Cancel jobs new End   ##########################################
