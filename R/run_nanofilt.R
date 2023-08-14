#' Run the NanoFilt program
#'
#' @description Runs the NanoFilt tool, can be used to filter nanopore reads reads
#'
#' @import parallel
#'
#' @param input List of the paths to files containing to the nanopore reads, required
#' @param sample_names List of sample names, required
#' @param out_dir Name of the directory from the BWA output
#' @param min_length Minimum length for the reads
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param nanofilt Path to the NanoFilt program, required
#' @param version Returns the version number
#'
#' @return A list with the NanoFilt commands
#'
#' @examples
#'  \dontrun{
#' path <- "/Path/To//NanoFilt"
#' input_files <- c("sample1.fastq", "sample2.fastq")
#' sample_names <- gsub(".fastq","",input_files)
#' out_directory <- "filtered_reads"
#' min_length <- 10000
#'
#' run_nanofilt(input = input_files,
#'              sample_names = sample_names,
#'              out_dir = out_directory,
#'              min_length = min_length,
#'              nanofilt = path)
#'}
#'
#'
#' @export
#'
run_nanofilt <- function(input = NULL,
                         sample_names = NULL,
                         out_dir = NULL,
                         min_length = NULL,
                         parallel = FALSE,
                         cores = 4,
                         execute = FALSE,
                         nanofilt = NULL,
                         version = FALSE){
  # Check NanoFilt program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", nanofilt)

  # Version
  if (isTRUE(version)){
    nanofilt.run <- sprintf('%s --version',
                            nanofilt)
    result <- system(nanofilt.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""
  # Minimum read length
  if (!is.null(min_length)){
    args <- paste(args,"-l",min_length, sep = " ")
  }

  # Set the output file names
  output <- paste(out_dir,paste(sample_names,"filt.fastq",sep = "_"),sep = "/")

  # Create the run commands
  nanofilt.run <- sprintf('%s %s %s > %s',
                          nanofilt,args,input,output)

  # Run the commands, if execute is true
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, nanofilt.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(nanofilt.run, function (cmd)  system(cmd))
    }
  }

  # Return the list NanoFilt commands
  return(nanofilt.run)
}
