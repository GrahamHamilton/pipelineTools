#' Run Cutadapt
#' @description Run the Cutadapt tool to remove sequencing adapters and low quality bases.
#'
#' @param mate1 List of the paths to files containing to the forward reads
#' @param mate2 List of the paths to files containing to the reverse reads
#' @param mate1.trim List of paths to the files to write the trimmed forward reads
#' @param mate2.trim List of paths to the files to write the trimmed reverse reads
#' @param quality The lower limit for the phred score
#' @param nextseq Was the sequence data generated on a NextSeq 500, trims dark cycle bases appearing as high-quality G bases
#' @param minimum The length at which a trimmed read will be discarded
#' @param trim.only Only keep reads that have had adapters trimmed
#' @param cut.for Remove the first 'n' bases form the 5' end of the forward read
#' @param cut.rev Remove the first 'n' bases form the 5' end of the reverse read
#' @param length Shorten each read down to a certain length
#' @param adapter1 Sequence for the adapter for the forward read
#' @param adapter2 Sequence for the adapter for the reverse read
#' @param polyA Number of A's
#' @param adapter.5.prime Sequence of te 5 prime adapter
#' @param maximum.error.rate Maximum number of errors tolerated, default 0.1 or 10 percent
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param cutadapt Path to the Cutadapt program, required
#' @param version Returns the version number
#'
#' @return A file with the Cutadapt commands and creates a directory of adapter and quality trimmed reads
#'
#' @examples
#'  \dontrun{
#'  # Version number
#'  run_cutadapt(cutadapt = cutadapt.path,
#'                          version = TRUE)
#'
#' # Trimmed reads directory
#' trimmed.reads.dir <- "trimmed_reads"
#' #Create the directory for the trimmed reads
#' dir.create(trimmed.reads.dir, showWarnings = FALSE)
#'
#' read1.pattern <- "*_R1_001.fastq.gz$"
#' read2.pattern <- "*_R2_001.fastq.gz$"
#'
#' reads.path <- "/export/buzz2/gpfrawdata/nextseq01/FastQ/2019/NeilBulleid/
#'                1466_NS205_0325_MarieAnnPringle_20190313"
#'
#' mate1 <- list.files(path = reads.path,
#'                     pattern = read1.pattern,
#'                     full.names = TRUE)
#' mate1.trim <- paste(trimmed.reads.dir,(list.files(path = reads.path,
#'                                                  pattern = read1.pattern,
#'                                                  full.names = FALSE)), sep = "/")
#'
#' mate2 <- list.files(path = reads.path,
#'                     pattern = read2.pattern,
#'                     full.names = TRUE)
#' mate2.trim <- paste(trimmed.reads.dir,(list.files(path = reads.path,
#'                                                  pattern = read2.pattern,
#'                                                  full.names = FALSE)), sep = "/")
#'
#' # Single end
#' run_cutadapt(mate1 = mate1,
#'              mate1.trim = mate1.trim,
#'              quality = 25,
#'              minimum = 17,
#'              trim.only = TRUE,
#'              cut = 1,
#'              adapter1 = "AGATCGGAAGAGCACACGTCT",
#'              cutadapt = cutadapt.path)
#'
#' # Paired end
#' run_cutadapt(mate1 = mate1,
#'              mate2 = mate2,
#'              mate1.trim = mate1.trim,
#'              mate2.trim = mate2.trim,
#'              quality = 25,
#'              minimum = 17,
#'              trim.only = TRUE,
#'              cut = 1,
#'              adapter1 = "AGATCGGAAGAGCACACGTCT",
#'              adapter1 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
#'              cutadapt = cutadapt.path)
#'  }
#'
#' @export
#'
run_cutadapt <- function(mate1 = NULL,
                         mate2 = NULL,
                         mate1.trim = NULL,
                         mate2.trim = NULL,
                         quality = NULL,
                         nextseq = FALSE,
                         minimum = NULL,
                         trim.only = FALSE,
                         cut.for = NULL,
                         cut.rev = NULL,
                         length = NULL,
                         adapter1 = NULL,
                         adapter2 = NULL,
                         polyA = NULL,
                         adapter.5.prime = NULL,
                         maximum.error.rate = NULL,
                         parallel = FALSE,
                         cores = 4,
                         cutadapt = NULL,
                         version = FALSE){
  # Check cutadapt program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", cutadapt)

  # Version
  if (isTRUE(version)){
    cutadapt.run <- sprintf('%s --version',
                             cutadapt)
    result <- system(cutadapt.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""
  # Quality
  if (!is.null(quality)){
    if (isTRUE(nextseq)){
      args <- paste(args,"--nextseq-trim", quality, sep = " ")
    }else{
      args <- paste(args,"-q", quality, sep = " ")
    }
  }
  # Minimum read length
  if (!is.null(minimum)){
    args <- paste(args,"-m", minimum, sep = " ")
  }
  # Keep the reads that are trimmed
  if (isTRUE(trim.only)){
    args <- paste(args,"--trimmed-only", sep = " ")
  }
  # Cut
  if (!is.null(cut.for)){
    args <- paste(args,"-u", cut.for, sep = " ")
  }
  if (!is.null(cut.rev)){
    args <- paste(args,"-U", cut.rev, sep = " ")
  }
  # Length
  if (!is.null(length)){
    args <- paste(args,"--length", length, sep = " ")
  }
  # Adapter1
  if (!is.null(adapter1)){
    args <- paste(args,"-a",adapter1, sep = " ")
  }
  # Adapter2
  if (!is.null(adapter2)){
    args <- paste(args,"-A",adapter2, sep = " ")
  }
  # Poly A trimming
  if (!is.null(polyA)){
    args <- paste(args,"-a",paste("A{",polyA,"}", sep = ""), sep = " ")
    # If trimming poly A tail as well as adapter need to set trimmimg number
    args <- paste(args,"--times=2", sep = " ")
  }
  # 5' adapter to trim
  if(!is.null(adapter.5.prime)){
    args <- paste(args,"-g",adapter.5.prime, sep = " ")
  }
  # Maximum error rate
  if(!is.null(maximum.error.rate)){
    args <- paste(args,"-e",maximum.error.rate,sep = " ")
  }


  # Single end
  if (is.null(mate2)){
    cutadapt.run <- sprintf('%s %s -o %s %s',
                            cutadapt,args,mate1.trim,mate1)
  }
  # Paired end
  if (!is.null(mate2) && !is.null(mate2.trim)){
    cutadapt.run <- sprintf('%s %s -o %s -p %s %s %s',
                            cutadapt,args,mate1.trim,mate2.trim,mate1, mate2)
  }

  # Run Cutadapts commands
  if (isTRUE(parallel)){
    cluster <- makeCluster(cores)
    parLapply(cluster, cutadapt.run, function (cmd) system(cmd, intern = FALSE, wait = TRUE))
    stopCluster(cluster)
  }else{
    lapply(cutadapt.run,function (cmd) system(cmd, intern = FALSE, wait = TRUE))
  }

  return(cutadapt.run)
}

