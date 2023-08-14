#' Run the Flye assembler
#'
#' @description Run the long read assembler Flye
#'
#' @import parallel
#'
#' @param input List of the paths to files containing to the long reads, required
#' @param sample_names List of sample names, required
#' @param platform Names of the sequencing platform, choose either,
#' pacbio_raw, PacBio regular CLR reads less than 20 percent error
#' pacbio_corr, PacBio reads that were corrected with other methods less than 3 percent error
#' pacbio_hifi, PacBio HiFi reads
#' nano_raw, ONT regular reads, pre-Guppy5 less than 20 percent error
#' nano_corr, ONT reads that were corrected with other methods less than 3 percent error
#' nano_hq , ONT high-quality reads: Guppy5+ or Q20 less than 5 percent error
#' @param out_dir Name of the directory from the assembled reads
#' @param threads Number of threads for flye to use, default set to 10
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param flye Path to the flye program, required
#' @param version Returns the version number
#'
#' @return A list with the Flye commands
#'
#' @examples
#'
#' \dontrun{
#' path <- "/path/to/flye"
#' input_files <- c("sample1.fastq", "sample2.fastq")
#' sample_names <- gsub(".fastq","",input_files)
#' platform <- "nano_hq"
#' out_directory <- "Assemblies"
#'
#' run_flye(flye = path,
#'          version = TRUE)
#'
#' run_flye(input = input_files,
#'          sample_names = sample_names,
#'          platform = platform,
#'          out_dir = out_directory,
#'          threads = 10,
#'          parallel = FALSE,
#'          cores = 4,
#'          execute = FALSE,
#'          flye = path)
#' }
#'
#' @export
#'
run_flye <- function(input = NULL,
                     sample_names = NULL,
                     platform = NULL,
                     out_dir = NULL,
                     threads = 10,
                     parallel = FALSE,
                     cores = 4,
                     execute = FALSE,
                     flye = NULL,
                     version = FALSE){
  # Check flye program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", flye)

  # Version
  if (isTRUE(version)){
    flye.run <- sprintf('%s --version',
                          flye)
    result <- system(flye.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""
  # Threads
  if (!is.null(threads)){
    args <- paste(args,"--threads",threads,sep = " ")
  }
  if (!is.null(platform)){
    if (platform == "nano_raw"){
      args <- paste(args,"--nano-raw",sep = " ")
    }else if (platform == "nano_corr"){
      args <- paste(args,"--nano_corr",sep = " ")
    }else if (platform == "nano_hq"){
      args <- paste(args,"--nano-hq",sep = " ")
    }else if (platform == "pacbio_raw"){
      args <- paste(args,"--pacbio-raw",sep = " ")
    }else if (platform == "pacbio_corr"){
      args <- paste(args,"--pacbio-corr",sep = " ")
    }else if (platform == "pacbio_hifi"){
      args <- paste(args,"--pacbio-hifi",sep = " ")
    }else{
      stop("Please provide either nano_raw, nano_corr, nano_corr, pacbio_raw, pacbio_corr or pacbio_hifi for platorm variable")
    }
  }

  # Set the names of the outout directory
  output_directories <- paste(out_dir,paste(sample_names,"assembly",sep = "_"),sep = "/")
  #Set the names for the logfiles
  logfiles <- paste(output_directories,paste(sample_names,"log",sep = "."),sep = "/")

  # Create the run commands
  flye.run <- sprintf('%s %s %s --out-dir %s &> %s',
                          flye,args,input,output_directories,logfiles)

  # Run the commands, if execute is true
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, flye.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(flye.run, function (cmd)  system(cmd))
    }
  }

  # Return the list Flye commands
  return(flye.run)
}
