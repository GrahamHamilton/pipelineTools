#' Run Macs3
#'
#'@description Macs3 is a spatial clustering approach for the identification
#'   of ChIP-enriched regions which was developed for calling narrow and broad peaks,
#'   for transcription factor binding or histone modifications from ChIP-seq data.
#'
#' @import parallel
#'
#' @param command Tool arguments, choose form callpeak,bdgpeakcall,bdgbroadcall,
#'  bdgcmp,bdgopt,cmbreps,bdgdiff,filterdup,predictd,pileup,randsample,refinepeak,callvar,hmmratac
#' @param treatment ChIP-seq treatment files list.
#' @param control Control files list.
#' @param format Format of the input files, choose from AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE
#' @param genome.size Effective genome size.
#' @param qvalue Minimum FDR (q-value) cutoff for peak detection.
#' @param broad.peaks Boolean, if set to TRUE, MACS will try to call broad peaks
#' @param broad.peaks.cutoff Minimum FDR (q-value) cutoff for broad region.
#' @param out.dir Name of output directory
#' @param experiment.name Experiment name, which will be used to generate output file names.
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param macs Path to the Macs3 program, required
#' @param version Returns the version number
#'
#' @return A list with the Macs3 commands
#'
#' @examples
#' \dontrun{
#' macs_path <- "/home/gmh5n/.localpython/bin/macs3"
#'
#' macs_version <- run_macs(macs = macs_path,
#'                          version = TRUE)
#'
#' macs_version
#'
#' macs_command <- "callpeak"
#' treatment <- c("1_sorted.bam","2_sorted.bam")
#' experiment.name <- gsub("_sorted.bam","",treatment)
#' control <- c("1-input_sorted.bam","2-input_sorted.bam")
#' size <- 1.9e9
#' macs_dir <- "macs"
#' macs_commands <- run_macs(command = command,
#'                           treatment = treatment,
#'                           control = control,
#'                           format = "BAM",
#'                           genome.size = size,
#'                           broad.peaks = TRUE,
#'                           broad.peaks.cutoff = 0.05,
#'                           out.dir = macs_dir,
#'                           parallel = TRUE,
#'                           cores = 2,
#'                           execute = FALSE,
#'                           experiment.name = experiment.name,
#'                           macs = macs_path)
#' macs_commands
#' }
#' @export
run_macs <- function(command = NULL,
                     treatment = NULL,
                     control = NULL,
                     format = NULL,
                     genome.size = NULL,
                     qvalue = NULL,
                     broad.peaks = NULL,
                     broad.peaks.cutoff = NULL,
                     out.dir = NULL,
                     experiment.name = NULL,
                     parallel = FALSE,
                     cores = 4,
                     execute = TRUE,
                     macs = NULL,
                     version = FALSE){
  # Check macs program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", macs)

  # Version
  if (isTRUE(version)){
    macs.run <- sprintf('%s  --version',
                       macs)
    version <- try(system(macs.run, intern = TRUE))
    return(version)
  }

  # Create the sample directories for the per sample macs results
  results_directories <- paste(out.dir,experiment.name, sep = "/")
  lapply(results_directories, function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))

  # Log files
  log_files <- paste(paste(results_directories,experiment.name, sep = "/"), "log", sep = ".")

  # Set the additional arguments
  args <- ""
  # Format
  if (!is.null(format)){
    args <- paste(args,"--format",format,sep = " ")
  }
  # Genome size
  if (!is.null(genome.size)){
    args <- paste(args,"--gsize",genome.size,sep = " ")
  }
  # Q value
  if(!is.null(qvalue)){
    args <- paste(args,"--qvalue",qvalue,sep = " ")
  }
  # Broad peaks
  if (!is.null(broad.peaks)){
    args <- paste(args,"--broad",sep = " ")
  }
  # Broad peaks cutoff
  if (!is.null(broad.peaks.cutoff)){
    args <- paste(args,"--broad-cutoff",broad.peaks.cutoff,sep = " ")
  }
  # Experiment name
  if (!is.null(experiment.name)){
    args <- paste(args,"--name",experiment.name,sep = " ")
  }

  macs.run <- sprintf('%s %s --treatment %s  --control %s --outdir %s %s 2> %s',
                       macs,command,treatment,control,results_directories,args,log_files)

  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, macs.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(macs.run, function (cmd)  system(cmd))
    }
  }

  return(macs.run)

}
