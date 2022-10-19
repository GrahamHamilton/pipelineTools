#' Run bowtie
#'
#' @description Runs the bowtie tool, currently set up primarily to get unaligned reads
#'
#' @import parallel
#'
#' @param input1 List of the paths to files containing to the forward reads, required
#' @param input2 List of the paths to files containing to the reverse reads, required for paired end sequence data
#' @param index Path to the reference transcriptome kallisto index, required
#' @param sample.name List of the sample names, required
#' @param out.dir Name of the directory from the Bowtie output. If NULL,
#'   which is the default, a directory named "bowtie_alignments" is created in the current working directory.
#' @param unaligned Name of the directory for the unaligned reads, fastq formatted
#' @param format Format for the aligned reads, currently only SAM is supported
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param bowtie Path to the bowtie proram, required
#' @param version Returns the version number
#'
#' @return A list with the bowtie commands
#'
#' @examples
#'  \dontrun{
#' bowtie.path <- "/software/bowtie-v1.2.3/bin/bowtie"
#'
#' bowtie.version <- run_bowtie(bowtie = bowtie.path,
#'                              version = TRUE)
#'
#' trimmed_reads_dir <- "trimmed_reads"
#' input1 <- list.files(path = trimmed_reads_dir, pattern = "*_R1_001.fastq$", full.names = TRUE)
#'
#' index <- "path/to/index/indes.ebwt"
#'
#' sample_names <- unlist(lapply(strsplit(list.files(path = trimmed_reads_dir,
#'                                                   pattern = "*_R1_001.fastq$",
#'                                                   full.names = FALSE),"_"), `[[`, 1))
#'
#' unaligned <- "rRNA_removed_reads"
#' out.dir <- "/dev/null" # Send the reads that align to be deleted
#'
#' bowtie.cmds <- run_bowtie(input1 = mate1,
#'                           index = index, sample.name = sample_names,
#'                           unaligned = unaligned,
#'                           out.dir = out.dir,
#'                           bowtie = bowtie.path)
#'}
#'
#' @export
#'

run_bowtie <- function(input1 = NULL,
                       input2 = NULL,
                       index = NULL,
                       sample.name = NULL,
                       out.dir = NULL,
                       unaligned = NULL,
                       format = "SAM",
                       parallel = FALSE,
                       cores = 4,
                       execute = TRUE,
                       bowtie = NULL,
                       version = FALSE){
  # Check bowtie program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", bowtie)

  # Version
  if (isTRUE(version)){
    bowtie_run <- sprintf('%s --version',
                          bowtie)
    result <- system(bowtie_run, intern = TRUE)
    return(result)
  }

  # Create the directory to write the data, if it is not present
  if (is.null(out.dir)) out.dir <- "bowtie_alignments"
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

  # Set the additional arguments
  args <- ""
  # Output format
  if (!is.null(format)){
    if (format == "SAM"){
      args <- paste(args,"--sam",sep = " ")
    }
  }

  # Create the output files for unaligned reads if set
  if (!is.null(unaligned)){
    dir.create(unaligned, showWarnings = FALSE, recursive = TRUE)
    unaligned.files <- paste(unaligned,paste(sample.name,"fastq",sep = "."),sep = "/")
    args <- paste(args,"--un",unaligned.files,sep = " ")
  }

  # Create output files
  if (out.dir != "/dev/null"){
    aligned.files <- paste(out.dir,sample.name,paste(sample.name,"sam",sep = "."),sep = "/")
    # Create the sample output directories
    lapply(paste(out.dir,sample.name, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))
  }

  # Create the log files
  log.files <- paste(out.dir,sample.name,paste(sample.name,"log",sep = "."),sep = "/")

  if (out.dir == "/dev/null"){
    bowtie.run <- sprintf('%s %s %s %s > %s 2> %s',
                          bowtie,args,index,input1,out.dir,log.files)
  }else{
    bowtie.run <- sprintf('%s %s %s %s > %s 2> %s',
                          bowtie,args,index,input1,aligned.files,log.files)
  }

  # Run the bowtie commands
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, bowtie.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(bowtie.run, function (cmd)  system(cmd))
    }
  }

  # Return the list of bowtie commands
  return(bowtie.run)
}
