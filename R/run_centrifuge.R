#' Run Centrifuge
#'
#' @description Runs the Centrifuge tool
#'
#' @import parallel
#'
#' @param input1 List of the paths to files containing to the forward reads
#' @param input2 List of the paths to files containing to the reverse reads
#' @param index Path to the reference metagenome index
#' @param sample.name List of the sample names
#' @param out.dir Name of the directory from the Centrifuge output
#' @param threads Number of threads for Centrifuge to use, default set to 10
#' @param min_hit_length minimum length of partial hits, default set to the minimum allowed of 15
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param centrifuge Path to the Centrifuge program
#' @param version Returns the version number
#'
#' @return A list with the Centrifuge commands
#'
#' @examples
#' \dontrun{
#' path <- "/software/centrifuge-1.0.4/centrifuge"
#' in1 <- c("sample_1_1.fq","sample_2_1.fq")
#' in2 <- c("sample_1_2.fq","sample_2_2.fq")
#' sample_names <- c("sample_1","sample_2")
#' ref <- "/path/to/centrifuge/reference/index"
#'
#' run_centrifuge(input1 = in1,
#'                input2 = in2,
#'                index = ref,
#'                sample.name = sample_names,
#'                out.dir = "test_dir",
#'                min_hit_length = 16,
#'                execute = FALSE,
#'                centrifuge = path)
#' }
#'
#' @export
#'

run_centrifuge <- function(input1 = NULL,
                           input2 = NULL,
                           index = NULL,
                           sample.name = NULL,
                           out.dir = NULL,
                           threads = 10,
                           min_hit_length = 15,
                           parallel = FALSE,
                           cores = 4,
                           execute = TRUE,
                           centrifuge = NULL,
                           version = FALSE){
  # Check centrifuge program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", centrifuge)

  # Version
  if (isTRUE(version)){
    centrifuge.run <- sprintf('%s --version',
                              centrifuge)
    result <- system(centrifuge.run, intern = TRUE)
    return(result)
  }

  # Create the sample directories for the per sample centrifuge results
  lapply(paste(out.dir,sample.name, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))

  # Set the additional arguments
  args <- ""
  # Hit length
  if (isTRUE(min_hit_length)){
    args <- paste(args,"--min-hitlen",sep = " ")
  }
  # Threads
  if (!is.null(threads)){
    args <- paste(args,paste("--threads",threads,sep = " "),sep = " ")
  }

  # Set the names for the alignment and logfiles
  log.files <- paste(out.dir,sample.name,paste(paste(sample.name,"hitlen",min_hit_length,sep = "_"),"log",sep = "."),sep = "/")
  sam.files <- paste(out.dir,sample.name,paste(paste(sample.name,"hitlen",min_hit_length,sep = "_"),"sam",sep = "."),sep = "/")

  # Paired end
  if(!is.null(input2)){
   centrifuge.run <- sprintf('%s %s -x %s -1 %s -2 %s --report-file %s -S %s ',
                           centrifuge,args,index,input1,input2,log.files,sam.files)
  }else{
    centrifuge.run <- sprintf('%s %s -x %s -S %s -U %s > %s 2>&1',
                              centrifuge,args,index,sam.files,input1,log.files)
  }

  # Run the commands, if execute is true
    if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, centrifuge.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(centrifuge.run, function (cmd)  system(cmd))
    }
  }

  # Return the list of Centrifuge commands
  return(centrifuge.run)
}
