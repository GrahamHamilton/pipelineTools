#' Run the BWA aligner
#'
#' @description Runs the BWA aligner
#'
#' @import parallel
#'
#' @param command BWA command to run, at present can choose 'mem',
#'                which is most useful fo Illumina, PacBio and Nanopore reads, required
#' @param mate1 List of the paths to files containing to the forward reads, required
#' @param mate2 List of the paths to files containing to the reverse reads, required for paired end sequence data
#' @param index Path to the indexed reference genome, required
#' @param sample.name List of sample names, required
#' @param out.dir Name of the directory from the BWA output
#' @param threads Number of threads for BWA to use, default set to 10
#' @param seed.length Minimum seed length
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param bwa Path to the bwa program, required
#' @param version Returns the version number
#'
#' @return A list with the BWA commands
#'
#' @examples
#'  \dontrun{
#'  path <- "/software/bwa/bwa"
#'
#' # Get version number
#' bwa.version <- run_bwa(bwa = path,
#'                        version = TRUE)
#' bwa.version[3] # print the line with version number
#'
#' reads.path <- "/path/to/fastqs"
#' reads.patt.1 <- "*_R1_001.fastq.gz$"
#' reads.patt.2 <- "*_R2_001.fastq.gz$"
#'
#' # Organise the samples
#' sample.dataframe <- prepare_samples(reads.path, c(reads.patt.1,reads.patt.2),trimmed.reads.dir)
#'
#' mate1 <- as.character(sample.dataframe$reads.path.1)
#' mate1.trim <- as.character(sample.dataframe$trimmed.reads.path.1)
#' # Paired end only
#'  mate2 <- as.character(sample.dataframe$reads.path.2)
#'  mate2.trim <- as.character(sample.dataframe$trimmed.reads.path.2)
#'
#' sample.names <- as.character(sample.dataframe$sample.names)
#'
#' genome <- "/path/to/reference/genome"
#'
#' command <- "mem"
#' out.dir <- "bwa_alignments"
#'
#' # Paired end
#' pe <- run_bwa(command = command,
#'               mate1 = mate1.trim,
#'               mate2 = mate2.trim,
#'               index = genome,
#'               threads = 20,
#'               out.dir = out.dir,
#'               sample.name = sample.names,
#'               seed.length = 21,
#'               bwa = path)
#'
#' # Single end
#' se <- run_bwa(command = command,
#'               mate1 = mate1.trim,
#'               index = genome,
#'               threads = 20,
#'               out.dir = out.dir,
#'               sample.name = sample.names,
#'               seed.length = 21,
#'               bwa = path)
#' }
#'
#' @export
#'
run_bwa <- function(command = NULL,
                    mate1 = NULL,
                    mate2 = NULL,
                    index = NULL,
                    sample.name = NULL,
                    out.dir = NULL,
                    threads = 10,
                    seed.length = NULL,
                    parallel = FALSE,
                    cores = 4,
                    bwa = NULL,
                    version = FALSE){
  # Check bwa program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", bwa)

  # Version
  if (isTRUE(version)){
    bwa.run <- sprintf('%s  2>&1',
                      bwa)
    result <- try(system(bwa.run, intern = TRUE))
    return(result)
  }

  # Create the sample directories for the per sample bwa results
  lapply(paste(out.dir,sample.name, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))

  # Set the additional arguments
  args <- ""
  # Threads
  if (!is.null(threads)){
    args <- paste(args,paste("-t",threads,sep = " "),sep = " ")
  }
  if (!is.null(seed.length)){
    args <- paste(args,paste("-k",seed.length,sep = " "),sep = " ")
  }

  # bwa mem -M -R $read_group_info -t $threads $genome $trimmed_reads_dir/$cut1 $trimmed_reads_dir/$cut2
  # Set the names for the alignment
  sam.files <- paste(out.dir,sample.name,paste(sample.name,"sam",sep = "."),sep = "/")
  # Paired end
  if(!is.null(mate2)){
    bwa.run <- sprintf('%s %s %s %s %s %s > %s',
                       bwa,command,args,index,mate1,mate2,sam.files)
  }else{
    bwa.run <- sprintf('%s %s %s %s %s > %s ',
                       bwa,command,args,index,mate1,sam.files)
  }

  if (isTRUE(parallel)){
    cluster <- makeCluster(cores)
    parLapply(cluster, bwa.run, function (cmd)  system(cmd))
    stopCluster(cluster)
  }else{
    lapply(bwa.run, function (cmd)  system(cmd))
  }

  return(bwa.run)
}
