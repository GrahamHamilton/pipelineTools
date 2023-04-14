#' Run the Bowtie2
#'
#' @description Runs the Bowtie2 tool
#'
#' @import parallel
#'
#' @param input1 List of the paths to files containing to the forward reads
#' @param input2 List of the paths to files containing to the reverse reads
#' @param index Path to the reference genome index
#' @param sample.name List of the sample names
#' @param out.dir Name of the directory from the Bowtie2 output
#' @param threads Number of threads for Bowtie2 to use, default set to 10
#' @param end2end Select presets for end to end alignments
#' @param sensitive Select the very sensitive preset
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param bowtie2 Path to the Bowtie2 program
#' @param version Returns the version number
#'
#' @return A list with the Bowtie2 commands
#'
#' @examples
#' \dontrun{

#' path <- "/software/bowtie_v2-2.3.5.1/bowtie2"
#'
#' bowtie2.version <- run_bowtie2(bowtie2 = path,
#'                                version = TRUE)
#' bowtie2.version[1]
#'
#' out <- "bowtie2_alignments"
#' trimmed_reads_dir <- "trimmed_reads"
#'
#' mate1 <- list.files(path = trimmed_reads_dir,
#'                     pattern = "*_R1_001.fastq$",
#'                     full.names = TRUE)
#' mate2 <- list.files(path = trimmed_reads_dir,
#'                    pattern = "*_R2_001.fastq$",
#'                    full.names = TRUE)
#' sample.names <- unlist(lapply(strsplit(list.files(path = trimmed_reads_dir,
#'                                                   pattern = "*_R1_001.fastq$",
#'                                                   full.names = FALSE),"_"), `[[`, 1))
#' index <- "/export/jessie3/gmh5n/BactoCap/JoHalliday/Project1543/MLSTReference/MLST"
#'
#' bowtie2.cmds <- run_bowtie2(input1 = mate1,
#'                             input2 = mate2,
#'                             index = index,
#'                             out.dir = out,
#'                             sample.name = sample.names,
#'                             bowtie2 = path,
#'                             version = FALSE)
#' bowtie2.cmds
#' }
#'
#' @export
#'

run_bowtie2 <- function(input1 = NULL,
                    input2 = NULL,
                    index = NULL,
                    sample.name = NULL,
                    out.dir = NULL,
                    threads = 10,
                    end2end = FALSE,
                    sensitive = FALSE,
                    parallel = FALSE,
                    cores = 4,
                    execute = TRUE,
                    bowtie2 = NULL,
                    version = FALSE){
  # Check bowtie2 program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", bowtie2)

  # Version
  if (isTRUE(version)){
    bowtie2.run <- sprintf('%s --version',
                       bowtie2)
    result <- system(bowtie2.run, intern = TRUE)
    return(result)
  }

  # Create the sample directories for the per sample bowtie2 results
  lapply(paste(out.dir,sample.name, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))

  # Set the additional arguments
  args <- ""
  # End to End Alignment
  if (isTRUE(end2end)){
    args <- paste(args,"--end-to-end",sep = " ")
    }
  # Sensitive
  if (isTRUE(sensitive)){
    args <- paste(args,"--very-sensitive",sep = " ")
  }
  # Threads
  if (!is.null(threads)){
    args <- paste(args,paste("--threads",threads,sep = " "),sep = " ")
  }

  # Set the names for the alignment and logfiles
  log.files <- paste(out.dir,sample.name,paste(sample.name,"log",sep = "."),sep = "/")
  sam.files <- paste(out.dir,sample.name,paste(sample.name,"sam",sep = "."),sep = "/")

  # Paired end
  if(!is.null(input2)){
    bowtie2.run <- sprintf('%s %s -x %s -S %s -1 %s -2 %s > %s 2>&1',
                           bowtie2,args,index,sam.files,input1,input2,log.files)
  }else{
    bowtie2.run <- sprintf('%s %s -x %s -S %s -U %s > %s 2>&1',
                           bowtie2,args,index,sam.files,input1,log.files)
  }

  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, bowtie2.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(bowtie2.run, function (cmd)  system(cmd))
    }
  }

  return(bowtie2.run)
}
