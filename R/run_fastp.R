#' Run FastP
#' @description Run the FastP tool to remove contaminating sequencing adapters and low quality bases.
#'
#' @param input1 List of the paths to files containing to the forward reads
#' @param input2 List of the paths to files containing to the reverse reads
#' @param output1 List of paths to the files to write the trimmed forward reads
#' @param output2 List of paths to the files to write the trimmed reverse reads
#' @param adapter1 Sequence for the adapter for the forward read
#' @param adapter2 Sequence for the adapter for the reverse read
#' @param sample.name List of the sample names
#' @param out.dir Name of the directory to write quality control results files. If NULL,
#'   which is the default, a directory named "fastP" is created in the current
#'   working directory.
#' @param phred.quality The lower limit for the phred score
#' @param min.length The length at which a trimmed read will be discarded
#' @param trim.front.1 Trim 'n' bases from front of read1, default is 0
#' @param trim.tail.1 Trim 'n' bases from tail of read1, default is 0
#' @param trim.front.2 Trim 'n' bases from front of read2, default is 0
#' @param trim.tail.2 Trim 'n' bases from tail of read2, default is 0
#' @param threads Number of threads for FastP to use, default set to 10
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
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
#' fastp.cmds <- run_fastp(input1 = mate1,
#'                         input2 = mate2,
#'                         output1 = mate1.out,
#'                         output2 = mate2.out,
#'                         adapter1 = adapter1,
#'                         adapter2 = adapter2,
#'                         sample.name =  sample.names,
#'                         out.dir = fastp.results.dir,
#'                         fastp = "/software/bin/fastp")
#' }
#'
#' @export
run_fastp <- function(input1 = NULL,
                      input2 = NULL,
                      output1 = NULL,
                      output2 = NULL,
                      adapter1 = NULL,
                      adapter2 = NULL,
                      sample.name = NULL,
                      out.dir = NULL,
                      phred.quality = 15,
                      min.length = NULL,
                      trim.front.1 = NULL,
                      trim.tail.1 = NULL,
                      trim.front.2 = NULL,
                      trim.tail.2 = NULL,
                      threads = 10,
                      parallel = FALSE,
                      cores = 4,
                      execute = TRUE,
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
  if (!is.null(phred.quality)){
    args <- paste(args,"--qualified_quality_phred",phred.quality, sep = " ")
  }
  # Minimum length
  if (!is.null(min.length)){
    args <- paste(args,"--length_required",min.length, sep = " ")
  }
  # Trim first n bases of read 1
  if (!is.null(trim.front.1)){
    args <- paste(args,"--trim_front1",trim.front.1, sep = " ")
  }
  # Trim last n bases of read 1
  if (!is.null(trim.tail.1)){
    args <- paste(args,"--trim_tail1",trim.tail.1, sep = " ")
  }
  # Trim first n bases of read 2
  if (!is.null(trim.front.2)){
    args <- paste(args,"--trim_front2",trim.front.2, sep = " ")
  }
  # Trim last n bases of read 2
  if (!is.null(trim.tail.2)){
    args <- paste(args,"--trim_tailt2",trim.tail.2, sep = " ")
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
  if (!is.null(input2)){
    fastp.run <- sprintf('%s %s --in1 %s --in2 %s --out1 %s --out2 %s',
                         fastp,args,input1,input2,output1,output2)
  }
  # Single end
  else if (is.null(input2)){
    fastp.run <- sprintf('%s %s --in1 %s --out1 %s',
                         fastp,args,input1,output1)
  }

  # Run the fastp commands on the command line in parallel or not.
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, fastp.run, function (cmd) system(cmd, intern = FALSE))
      stopCluster(cluster)
    }else{
      lapply(fastp.run, function (cmd) system(cmd, intern = FALSE))
    }
  }

  return(fastp.run)
}
