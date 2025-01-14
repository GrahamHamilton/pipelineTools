#' Run Seqtk
#'
#' @description Runs the seqtk tool
#'
#' @import parallel
#'
#' @param command Seqtk command to run, at present can choose from 'seq' and 'sample', required
#' @param input List of fastq files, can be gzipped
#' @param output list of output file names, fastq format
#' @param seed Random seed for read selection
#' @param reverse Reverse complement the fastq sequences
#' @param num_reads Number of reads to output for sub sampling
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param seqtk Path to the seqtk program, required
#'
#' @return A file with the Seqtk commands
#'
#' @examples
#' \dontrun{
#' num_reads <- 1000
#' cmd <- "sample"
#' seed <- 23
#' reads <- c("reads1.fq.gz","reads2.fq.gz")
#' out <- c("reads1_out.fq.","reads2_out.fq")
#' seqtk_path <- "/usr/bin/seqtk"
#'
#' run_seqtk(command = cmd,
#'           input = reads,
#'           output = out,
#'           seed = seed,
#'           num_reads = num_reads,
#'           parallel = TRUE,
#'           cores = 2,
#'           execute = FALSE,
#'           seqtk = seqtk_path)
#' }
#'
#' @export
run_seqtk <- function(command = NULL,
                   input = NULL,
                   output = NULL,
                   seed = NULL,
                   reverse = NULL,
                   num_reads = NULL,
                   parallel = FALSE,
                   cores = 4,
                   execute = TRUE,
                   seqtk = NULL){

  # Set the additional arguments
  args <- ""

  # Seed
  if (!is.null(seed)){
    args <- paste(args,"-s",seed,sep = "")
  }
  # Reverse
  if (!is.null(reverse)){
    args <- paste(args,"-r",reverse,sep = " ")
  }

  # Seq
  if (command == "seq"){
    seqtk.run <- sprintf('%s %s %s %s > %s',
                            seqtk,command,args,input,output)
  }
  #Subseq
  if (command == "sample"){
    seqtk.run <- sprintf('%s %s %s %s %s > %s',
                         seqtk,command,args,input,num_reads,output)
  }

  # Run the seqtk commands
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- parallel::makeCluster(cores)
      parLapply(cluster, seqtk.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(seqtk.run, function (cmd)  system(cmd))
    }
  }

  return(seqtk.run)
}

