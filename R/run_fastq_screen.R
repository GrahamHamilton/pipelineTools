#' Run Fastq-Screen
#' @description Run the fastq-screen tool
#'
#' @param fq.files Vector of the full names and paths of the fastq files, usually only the forward reads set, required
#' @param out.dir Path to the directory to which to write the results. If NULL,
#'   which is the default, a directory named "fastq_screen" is created in the current
#'   working directory.
#' @param threads The number of threads to use.  The default is 2.
#' @param fastq_screen The path to the fastq_screen executable.
#' @param aligner The name of the program to perform the mapping, default bwa. Valid mappers are "bwa", "bowtie" and "bowtie2". Default set to bwa
#' @param conf The path to the configuration file, required
#' @param top Create a temporary datset by selecting the number of reads followed by how many to skip e.g deafaul is set to 500000,1000000 map 500 thousand reads and skip the first 1 million reads
#' @param version Returns the version number
#'
#' @return A file with the FastqScreen commands
#'
#' @examples
#'  \dontrun{
#' path <- "raw_reads"
#' mate1 <- list.files(path = path, pattern = "*_R1_001.fastq.gz$", full.names = TRUE)
#'
#' fastq_screen_cmds <- run_fastq_screen(fq.files = mate1,
#'                                       out.dir = fastq.screen.dir,
#'                                       conf = "/export/jessie3/gmh5n/PipelineTest/fastq_screen.conf",
#'                                       top = "100000,5000",
#'                                       fastq_screen = "/software/fastq_screen_v0.13.0/fastq_screen")
#' }
#'
#' @export

run_fastq_screen <- function(fq.files = NULL,
                             out.dir = NULL,
                             aligner = "bwa",
                             conf = NULL,
                             top = "500000,1000000",
                             threads = 2,
                             fastq_screen = NULL,
                             version = FALSE) {
  # Check fastq_screen program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", fastq_screen)

  # Version
  if (isTRUE(version)){
    fastqscreen.run <- sprintf('%s --version',
                               fastq_screen)
    result <- system(fastqscreen.run, intern = TRUE)
    return(result)
  }

  # Create the directory to write the data, if it is not present
  if (is.null(out.dir)){
    out.dir <- "fastq_screen"
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Set the additional arguments
  args <- ""
  # Force
  args <- paste(args,"--force",sep = " ")
  # Threads
  if (!is.null(threads)){
    args <- paste(args,"--threads",threads,sep = " ")
  }

  fastqscreen.run <- sprintf('%s %s --aligner %s --top %s --conf %s --outdir %s %s',
                             fastq_screen,args,aligner,top,conf,out.dir,fq.files)
  # Run Fastq Screen commands
  lapply(fastqscreen.run,function (cmd) system(cmd, intern = FALSE, wait = TRUE))

  return(fastqscreen.run)
}
