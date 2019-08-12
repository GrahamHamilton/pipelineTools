#' Run FastP
#' @description Run the FastP tool to remove contaminating sequencing adapters and low quality bases.
#'
#' @param mate1 List of the paths to files containing to the forward reads
#' @param mate2 List of the paths to files containing to the reverse reads
#' @param mate1.out List of paths to the files to write the trimmed forward reads
#' @param mate2.out List of paths to the files to write the trimmed reverse reads
#' @param adapter1 Sequence for the adapter for the forward read
#' @param adapter2 Sequence for the adapter for the reverse read
#' @param sample.name List of the sample names
#' @param out.dir Name of the directory to write quality control results files. If NULL,
#'   which is the default, a directory named "fastP" is created in the current
#'   working directory.
#' @param phred_quality The lower limit for the phred score
#' @param min_length The length at which a trimmed read will be discarded
#' @param threads Number of threads for FastP to use, default set to 10
#' @param fastp Path to the FastP program, required
#' @param version Returns the version number
#'
#' @return A file with the FastP commands and creates a directory of adapter and quality trimmed reads
#' @examples
#'  \dontrun{
#' # Set the directory containing the raw fastq files
#' reads_path <- "raw_reads"
#' mate1 <- list.files(path = reads_path, pattern = "*_R1_001.fastq.gz$", full.names = TRUE)
#' mate2 <- list.files(path = reads_path, pattern = "*_R2_001.fastq.gz$", full.names = TRUE)
#'
#' # Set the directory for writing the trimmend reads to
#' trimmed_reads_dir <- "trimmed_reads"
#' mate1.out <- paste(trimmed_reads_dir,
#'              (list.files(path = path, pattern = "*_R1_001.fastq.gz$", full.names = FALSE)),
#'              sep = "/")
#' mate2.out <- paste(trimmed_reads_dir,
#'              (list.files(path = path, pattern = "*_R2_001.fastq.gz$", full.names = FALSE)),
#'              sep = "/")
#'
#' # Get the sample names from the first reads
#' sample_names <- unlist(lapply(strsplit
#'                 (list.files(path = path, pattern = "*_R1_001.fastq.gz$", full.names = FALSE),"_"),
#'                 `[[`, 1))
#'
#' # Set the adapter sequences, these are for Illumina
#' adapter1 <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
#' adapter2 <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
#'
#' fastp.cmds <- run_fastp(mate1 = mate1,
#'                         mate2 = mate2,
#'                         mate1.out = mate1.out,
#'                         mate2.out = mate2.out,
#'                         adapter1 = adapter1,
#'                         adapter2 = adapter2,
#'                         sample.name =  sample.names,
#'                         out.dir = fastp.results.dir,
#'                         fastp = "/software/bin/fastp")
#' }
#'
#' @export
run_fastp <- function(mate1 = NULL,
                      mate2 = NULL,
                      mate1.out = NULL,
                      mate2.out = NULL,
                      adapter1 = NULL,
                      adapter2 = NULL,
                      sample.name = NULL,
                      out.dir = NULL,
                      phred_quality = 15,
                      min_length = 54,
                      threads = 10,
                      fastp = NULL,
                      version= FALSE){
  # Check FastP program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", fastp)

  # Version
  if (isTRUE(version)){
    fastp_run <- sprintf('%s --version 2>&1',
                         fastp)
    result <- system(fastp_run, intern = TRUE)
    return(result)
    }

  # Create the directory to write the data, if it is not present
  if (is.null(out.dir)){
    out.dir <- "FastpQC"
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
    }

  # Set the additional arguments
  args <- ""
  # Phred
  if (!is.null(phred_quality)){
    args <- paste(args,"--qualified_quality_phred",phred_quality, sep = " ")
  }
  # Minimum length
  if (!is.null(min_length)){
    args <- paste(args,"--length_required",min_length, sep = " ")
  }
  # Threads
  if (!is.null(threads)){
    args <- paste(args,paste("--thread",threads,sep = " "),sep = " ")
  }
  # Forward adapter
  if (!is.null(adapter1)){
    args <- paste(args,"--adapter_sequence",adapter1, sep = " ")
  }
  # Reverse adapter
  if (!is.null(adapter2)){
    args <- paste(args,"--adapter_sequence_r2",adapter2, sep = " ")
  }
  # HTML file name
  args <- paste(args,"--html",paste(paste(out.dir,sample.name, sep = "/"),"_fastp.html", sep = ""), sep = " ")
  # JSON file name
  args <- paste(args,"--json",paste(paste(out.dir,sample.name, sep = "/"),"_fastp.json", sep = ""), sep = " ")
  # Report title
  args <- paste(args,"--report_title",sample.name, sep = " ")

  # Paired end
  if (!is.null(mate2)){
    fastp.run <- sprintf('%s %s --in1 %s --in2 %s --out1 %s --out2 %s',
                         fastp,args,mate1,mate2,mate1.out,mate2.out)
  }
  # Single end
  else if (is.null(mate2)){
    fastp.run <- sprintf('%s %s --in1 %s --out1 %s',
                         fastp,args,mate1,mate1.out)
  }

  # Run the fastp commands on the command line
  lapply(fastp.run, function (cmd) system(cmd, intern = FALSE))

  return(fastp.run)
}
