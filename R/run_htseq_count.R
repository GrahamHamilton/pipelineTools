#' Run HtSeq-Count
#'
#' @description Runs the htseq-count program
#'
#' @import parallel
#'
#' @param input List of aligned files in sam or bam format, required
#' @param output List of file names for output,
#' @param mode Mode to handle reads overlapping more than one feature
#' (choices: union, intersection-strict, intersection- nonempty; default: union)
#' @param type Feature type (3rd column in GTF file) to be used,
#' @param attribute GTF attribute to be used as feature ID,
#' @param format Type of <alignment_file> data, can be either sam, bam or auto
#' @param gtf Path to the GTF file
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param htseq_count Path to the htseq-count program program, required
#' @param version Returns the version number
#'
#' @return A list with the htseq_count commands
#'
#' @examples
#' \dontrun{
#' path <- "/home/gmh5n/.local/bin/htseq-count"
#'
#' # Version
#' run_htseq_count(htseq_count = path,
#'                 version = TRUE)
#'
#' hisat2.alignments.dir <- "/path/to/hisat/alignments/directory"
#' bam.files <- list.files(path = hisat2.alignments.dir,
#'                        pattern = "sorted.bam$",
#'                        full.names = TRUE,
#'                        recursive = TRUE)
#' counts.dir <- "counts_dir"
#' htseq.counts.files <- paste(counts.dir,
#'                             paste(list.files(path = hisat2.alignments.dir,recursive = FALSE),
#'                             "counts.txt",sep = "_"),
#'                             sep = "/")
#'
#' # Path to the gtf file
#' gtf.file <- "/path/to/gtffile"
#'
#' htseq.counts.cmds <- run_htseq_count(input = bam.files,
#'                                      output = htseq.counts.files,
#'                                      mode = "intersection-nonempty",
#'                                      type = "gene",
#'                                      attribute = "ID",
#'                                      format = "bam",
#'                                      gtf = gtf.file,
#'                                      htseq_count = path)
#' htseq.counts.cmds
#' }
#' @export
#'
run_htseq_count <- function(input = NULL,
                            output = NULL,
                            mode = NULL,
                            type = NULL,
                            attribute = NULL,
                            format = NULL,
                            gtf = NULL,
                            parallel = FALSE,
                            cores = 4,
                            execute = TRUE,
                            htseq_count = NULL,
                            version = FALSE){
  # Check cutadapt program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", htseq_count)

  # Version
  if (isTRUE(version)){
    htseq_count.run <- sprintf('%s --version',
                               htseq_count)
    result <- system(htseq_count.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""
  # Mode
  if (!is.null(mode)){
    args <- paste(args,paste("--mode",mode,sep = " "),sep = " ")
  }
  # Type
  if (!is.null(type)){
    args <- paste(args,paste("--type",type,sep = " "),sep = " ")
  }
  # Attribute
  if (!is.null(attribute)){
    args <- paste(args,paste("--idattr",attribute,sep = " "),sep = " ")
  }
  # Format
  if (!is.null(format)){
    args <- paste(args,paste("--format",format,sep = " "),sep = " ")
  }

  htseq_count.run <- sprintf('%s %s %s %s > %s',
                          htseq_count,args,input,gtf,output)
  # Run the htseq_count commands
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, htseq_count.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(htseq_count.run, function (cmd)  system(cmd))
    }
  }

  # Return the list of htseq-counts commands
  return(htseq_count.run)
}

