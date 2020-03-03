#' Run HISAT2
#'
#' @description Runs the HISAT2 tool, can be used for single end and paired end reads
#'
#' @import parallel
#'
#' @param mate1 List of the paths to files containing to the forward reads, required
#' @param mate2 List of the paths to files containing to the reverse reads, required for paired end sequence data
#' @param index Path to the reference transcriptome kallisto index, required
#' @param sample.name List of the sample names, required
#' @param strandedness Strand-specific information
#' @param no_splice Disable spliced alignment, use for Trypanosomes
#' @param known_splice File with known splice sites
#' @param assembly Reports alignments tailored for transcript assemblers
#' @param threads Number of threads for hisat2 to use, default set to 10
#' @param phred Quality score offsets, default set to illumina/sanger standard of 33
#' @param out.dir Name of the directory from the HISAT2 output. If NULL,
#'   which is the default, a directory named "hisat2_alignments" is created in the current
#'   working directory.
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param hisat2 Path to the HISAT2 program, required
#' @param version Returns the version number
#'
#' @examples
#'  \dontrun{
#'  trimmed_reads_dir <- "trimmed_reads"
#'  mate1 <- list.files(path = trimmed_reads_dir, pattern = "*_R1_001.fastq$", full.names = TRUE)
#'  mate2 <- list.files(path = trimmed_reads_dir, pattern = "*_R2_001.fastq$", full.names = TRUE)
#'
#'  sample_names <- unlist(lapply(strsplit(list.files(path = trimmed_reads_dir,
#'                         pattern = "*_R1_001.fastq$",
#'                         full.names = FALSE),"_"), `[[`, 1))
#'
#'  index <- "/export/buzz1/Genome/Homo_sapiens/Ensembl/GRCH38_p7/Sequence/Transcriptome/
#'            KallistoIndex/GRCh38_p7.kall"
#'
#'  # Paired End example
#'  strandedness <- "RF"
#'  hisat2.cmds <- run_hisat2(mate1 = mate1,
#'                            mate2 = mate2,
#'                            index = genome,
#'                            sample.name = sample.names,
#'                            strandedness = strandedness,
#'                            assembly = TRUE,
#'                            out.dir = hisat.alignments.dir,
#'                            hisat2 = "/software/hisat_v2-2.1.0/hisat2")
#'
#'  # Single End example
#'  strandedness <- "R"
#'  hisat2.cmds <- run_hisat2(mate1 = mate1,
#'                            index = genome,
#'                            sample.name = sample.names,
#'                            strandedness = strandedness,
#'                            assembly = TRUE,
#'                            out.dir = hisat.alignments.dir,
#'                            hisat2 = "/software/hisat_v2-2.1.0/hisat2")
#'  }
#'
#' @return A list with the HISAT2 commands
#' @export
#'
run_hisat2 <- function(mate1 = NULL,
                       mate2 = NULL,
                       index = NULL,
                       sample.name = NULL,
                       strandedness = NULL,
                       no_splice = NULL,
                       known_splice = NULL,
                       assembly = FALSE,
                       phred = 33,
                       threads = 10,
                       out.dir = NULL,
                       parallel = FALSE,
                       cores = 4,
                       hisat2 = NULL,
                       version = FALSE){
  # Check hisat2 program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", hisat2)

  # Version
  if (isTRUE(version)){
    hisat2_run <- sprintf('%s --version',
                          hisat2)
    result <- system(hisat2_run, intern = TRUE)
    return(result)
  }

  # Create the directory to write the data, if it is not present
  if (is.null(out.dir)){
    out.dir <- "hisat2_alignments"
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Create the sample directories for the per sample HISAT2 results
  lapply(paste(out.dir,sample.name, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))

  # Set the additional arguments
  args <- ""
  # Strandedness
  if (!is.null(strandedness)){
    args <- paste(args,"--rna-strandness",strandedness,sep = " ")
  }
  # No spliced alignment
  if (!is.null(no_splice)){
    args <- paste(args,"--no-spliced-alignment",sep = " ")
  }
  # Assembly
  if (isTRUE (assembly)){
    args <- paste(args,"--dta",sep = " ")
  }
  # Threads
  if (!is.null(threads)){
    args <- paste(args,paste("--threads",threads,sep = "="),sep = " ")
  }
  # Phred score
  if (phred != 33){
    args <- paste(args,"--phred64",sep = " ")
  }else{
    args <- paste(args,"--phred33",sep = " ")
  }
  # Known splice
  if (!is.null(known_splice)){
    args <- paste(args,"--known-splicesite-infile",known_splice,sep = " ")
  }


  # Set the names for the alignment and logfiles
  logfile <- paste(out.dir,sample.name,paste(sample.name,"log",sep = "."),sep = "/")
  samfile <- paste(out.dir,sample.name,paste(sample.name,"sam",sep = "."),sep = "/")

  # Create the HISAT2 commands
  # Paired end
  if (!is.null(mate2)){
    hisat2.run <- sprintf('%s %s -x %s -1 %s -2 %s -S %s > %s 2>&1',
                       hisat2,args,index,mate1,mate2,samfile,logfile)
  }
  # Single end
  else if (is.null(mate2)){
    hisat2.run <- sprintf('%s %s -x %s -U %s -S %s > %s 2>&1',
                         hisat2,args,index,mate1,samfile,logfile)
  }

  # Run the HISAT2 commands
  if (isTRUE(parallel)){
    cluster <- makeCluster(cores)
    parLapply(cluster, hisat2.run, function (cmd)  system(cmd))
    stopCluster(cluster)
  }else{
    lapply(hisat2.run, function (cmd)  system(cmd))
  }

  # Return the list of HISAT2 commands
  return(hisat2.run)
}

