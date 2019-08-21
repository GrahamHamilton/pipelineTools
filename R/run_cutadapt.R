#' Run Cutadapt
#' @description Run the Cutadapt tool to remove sequencing adapters and low quality bases.
#'
#' @param mate1 List of the paths to files containing to the forward reads
#' @param mate2 List of the paths to files containing to the reverse reads
#' @param mate1.out List of paths to the files to write the trimmed forward reads
#' @param mate2.out List of paths to the files to write the trimmed reverse reads
#' @param quality The lower limit for the phred score
#' @param minimum The length at which a trimmed read will be discarded
#' @param trim.only Only keep reads that have had adapters ttrimmed
#' @param cut Remove the first 'n' bases form the 5' end of the forward read
#' @param adapter1 Sequence for the adapter for the forward read
#' @param adapter2 Sequence for the adapter for the reverse read
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
#' mate1.out <- paste(trimmed.reads.dir,(list.files(path = reads.path,
#'                                                  pattern = read1.pattern,
#'                                                  full.names = FALSE)), sep = "/")
#'
#' mate2 <- list.files(path = reads.path,
#'                     pattern = read2.pattern,
#'                     full.names = TRUE)
#' mate2.out <- paste(trimmed.reads.dir,(list.files(path = reads.path,
#'                                                  pattern = read2.pattern,
#'                                                  full.names = FALSE)), sep = "/")
#'
#' # Single end
#' run_cutadapt(mate1 = mate1,
#'              mate1.out = mate1.out,
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
#'              mate1.out = mate1.out,
#'              mate2.out = mate2.out,
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
                         mate1.out = NULL,
                         mate2.out = NULL,
                         quality = 15,
                         minimum = NULL,
                         trim.only = FALSE,
                         cut = NULL,
                         adapter1 = NULL,
                         adapter2 = NULL,
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
    args <- paste(args,"-q", quality, sep = " ")
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
  if (!is.null(cut)){
    args <- paste(args,"--cut", cut, sep = " ")
  }
  # Adapter1
  if (!is.null(adapter1)){
    args <- paste(args,"-a",adapter1, sep = " ")
  }
  # Adapter2
  if (!is.null(adapter2)){
    args <- paste(args,"-A",adapter2, sep = " ")
  }

  # Single end
  if (is.null(mate2)){
    cutadapt.run <- sprintf('%s %s -o %s %s',
                            cutadapt,args,mate1.out,mate1)
  }
  # Paired end
  if (!is.null(mate2) && !is.null(mate2.out)){
    cutadapt.run <- sprintf('%s %s -o %s -p %s %s %s',
                            cutadapt,args,mate1.out,mate2.out,mate1, mate2)
  }

  # Run Cutadapts commands
  lapply(cutadapt.run,function (cmd) system(cmd, intern = FALSE, wait = TRUE))


  return(cutadapt.run)
}

