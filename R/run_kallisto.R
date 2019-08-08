#' Run Kallisto
#'
#' @description Runs the kallisto quant tool
#'
#' @param mate1 List of the paths to files containing to the forward reads, required
#' @param mate2 List of the paths to files containing to the reverse reads, required for paired end sequence data
#' @param index Path to the reference transcriptome kallisto index, required
#' @param sample.name List of the sample names, required
#' @param fusion Search for fusions used for Pizzly
#' @param strandedness Strand spcific reads, values are "first" or "second"
#' @param fragment.length Estimated average insert fragment size, must be set for single end sequence data
#' @param std.dev Estimated standard deviation of insert fragment size, must be set for single end sequence data
#' @param out.dir Name of the directory from the Kallisto output. If NULL,
#'   which is the default, a directory named "kallisto" is created in the current working directory.
#' @param threads Number of threads for kallisto to use, default set to 10
#' @param bootstrap Number of bootstrap samples, default set to 100
#' @param kallisto Path to the kallisto program, required
#' @param version Returns the version number
#'
#' @return A list with the kallisto commands
#'
#' @examples
#'  \dontrun{
#' trimmed_reads_dir <- "trimmed_reads"
#' mate1 <- list.files(path = trimmed_reads_dir, pattern = "*_R1_001.fastq$", full.names = TRUE)
#' mate2 <- list.files(path = trimmed_reads_dir, pattern = "*_R2_001.fastq$", full.names = TRUE)
#'
#' sample_names <- unlist(lapply(strsplit(list.files(path = trimmed_reads_dir,
#'                        pattern = "*_R1_001.fastq$",
#'                        full.names = FALSE),"_"), `[[`, 1))
#'
#' index <- "/export/buzz1/Genome/Homo_sapiens/Ensembl/GRCH38_p7/Sequence/Transcriptome/
#'           KallistoIndex/GRCh38_p7.kall"
#'
#' strandedness <- "second"
#'
#' # Paired end
#' kallisto.cmds <- run_kallisto(mate1 = mate1,
#'                               mate2 = mate2,
#'                               index = transcriptome,
#'                               sample.name = sample.names,
#'                               strandedness = strandedness,
#'                               out.dir = kalisto.results.dir,
#'                               kallisto = "/software/kallisto_v0.45.1/kallisto")
#' # Single end
#' kallisto.cmds <- run_kallisto(mate1 = mate1,
#'                               index = transcriptome,
#'                               sample.name = sample.names,
#'                               strandedness = strandedness,
#'                               fragment.length = 300,
#'                               std.dev = 25,
#'                               out.dir = kalisto.results.dir,
#'                               kallisto = "/software/kallisto_v0.45.1/kallisto")
#' }
#'
#' @export
#'
run_kallisto <- function(mate1 = NULL,
                         mate2 = NULL,
                         index = NULL,
                         sample.name = NULL,
                         fusion = NULL,
                         strandedness = NULL,
                         out.dir = NULL,
                         threads = 10,
                         bootstrap = 100,
                         fragment.length = 300,
                         std.dev = 25,
                         kallisto = "",
                         version = FALSE){
  # Check kallisto program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", kallisto)

  # Version
  if (isTRUE(version)){
    kallisto.run <- sprintf('%s version',
                            kallisto)
    result <- system(kallisto.run, intern = TRUE)
    return(result)
  }

  # Create the directory to write the data, if it is not present
  if (is.null(out.dir)){
    out.dir <- "kallisto"
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
    }

  # Create the sample directories for the per sample kallisto results
  lapply(paste(out.dir,sample.name, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))

  # Logfile names
  logfile<- paste(out.dir,sample.name,paste(sample.name,"kallisto.log", sep = "_"), sep = "/")

  # Set the additional arguments
  args <- ""
  # Fusion
  if (!is.null(fusion)){
    args <- paste(args,"--fusion",sep = " ")
  }
  # Strandedness
  if (!is.null(strandedness)){
    # First strand
    if (strandedness == "first"){
      args <- paste(args,"--fr-stranded",sep = " ")
    }else if (strandedness == "second"){
      args <- paste(args,"--rf-stranded",sep = " ")
    }
  }
  # Threads
  if (!is.null(threads)){
    args <- paste(args,paste("--threads",threads,sep = "="),sep = " ")
  }
  # Bootstrap
  if (!is.null(bootstrap)){
    args <- paste(args,paste("--bootstrap-samples",bootstrap,sep = "="),sep = " ")
  }

  # Create the kallisto commands
  # Paired end
  if (!is.null(mate2)){
    kallisto.run <- sprintf('%s quant %s --index=%s --output-dir=%s %s %s > %s 2>&1',
                            kallisto,args,index,paste(out.dir,sample.name, sep = "/"),mate1,mate2,logfile)
  }
  # Single end
  else if (is.null(mate2)){
    # Single
    args <- paste(args,"--single",sep = " ")
    # Fragment length
    if (!is.null(fragment.length)){
      args <- paste(args,paste("--fragment-length",fragment.length,sep = "="),sep = " ")
    }
    # Standard deviation
    if (!is.null(fragment.length)){
      args <- paste(args,paste("--sd",std.dev,sep = "="),sep = " ")
    }
    kallisto.run <- sprintf('%s quant %s --index=%s --output-dir=%s %s > %s 2>&1',
                            kallisto,args,index,paste(out.dir,sample.name, sep = "/"),mate1,logfile)
  }


  # Run the kallisto commands
  lapply(kallisto.run, function (cmd)  system(cmd))

  # Return the list of kaliisto commands
  return(kallisto.run)
}
