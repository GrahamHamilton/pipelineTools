#' Run FastQC
#'
#' @description Runs the FastQC tool, a quality control tool for high throughput sequence data
#'
#' @param input List of all of the fastq files, forward and reverse, required
#' @param out.dir Name of the directory to write the FastQC results, required
#' @param threads Number of threads for FastQC
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param fastqc Path to the FastQC program, required
#' @param version Returns the version number
#'
#' @return A list with the FastQC commands
#'
#' @examples
#' \dontrun{
#' fastqc.path <- "/software/FastQC-v0.11.8/fastqc"
#'
#' # Version
#' fastqc.version <- ""
#' fastqc.version <- run_fastqc(fastqc = fastqc.path,
#'                              version = TRUE)
#' fastqc.version
#'
#' # Run fastqc
#' reads.path <- "raw_reads"
#' mate1 <- list.files(path = reads.path, pattern = "*_R1_001.fastq.gz$", full.names = TRUE)
#' mate2 <- list.files(path = reads.path, pattern = "*_R2_001.fastq.gz$", full.names = TRUE)
#' # Merge the reads file names
#' all.reads <- c(mate1,mate2)
#'
#' out.dir <- "fastqc"
#'
#' fastqc.cmds <- run_fastqc(input = all.reads,
#'                           out.dir = out.dir,
#'                           fastqc = fastqc.path)
#' fastqc.cmds
#' }
#' @export
#'
run_fastqc <- function(input = NULL,
                       out.dir = out.dir,
                       threads = NULL,
                       parallel = FALSE,
                       cores = 4,
                       execute = TRUE,
                       fastqc = NULL,
                       version = FALSE){
  # Check fastqc program can be found
  sprintf("type %s &>//dev//null && echo 'TRUE' || echo 'FALSE'", fastqc)

  # Version
  if (isTRUE(version)){
    fastqc.run <- sprintf('%s --version',
                          fastqc)
    result <- system(fastqc.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""
  # Threads
  if (!is.null(threads)){
    args <- paste(args,"--threads",threads,sep = " ")
  }

  fastqc.run <- sprintf('%s %s -o %s %s',
                        fastqc,args,out.dir,input)

  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, fastqc.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }
    else{
      lapply(fastqc.run, function (cmd)  system(cmd))
    }
  }

  return(fastqc.run)
}
