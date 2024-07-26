#' Run the MaSuRCA assembly program
#'
#'  @description Creates a shell script to perform the hybrid assembly an correction of illumina and either nanopore or pacbio reads
#'
#' @param input1 List of the paths to files containing to the forward reads, can be gzipped
#' @param input2 List of the paths to files containing to the reverse reads, can be gzipped
#' @param input_long List of the paths to files containing the nanopore ot pacbio read, can be gzipped
#' @param output Name of assembly script, default assemble.sh
#' @param threads Number of threads
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param masurca Path to the masurca program, required
#' @param version Returns the version number
#'
#' @return A list with the masurca commands
#'
#' @examples
#' \dontrun{
#' # Path to masurca executable file
#'  path <- "/software/MaSuRCA-4.1.1/bin/masurca"
#'
#'  # Check the program version
#'  run_masurca(masurca = path,
#'              version = TRUE)
#'
#'  # Set paths to the reads files
#'  fastq1 <- "illumina1.fq.gz"
#'  fastq2 <- "illumina2.fq.gz"
#'  nanopore <- "nanopore.fq.gz"
#'
#'  # Run the masurca commands to create the assembly script
#'  run_masurca(input1 = fastq1,
#'              input2 = fastq2,
#'              input_long = nanopore,
#'              threads = 40,
#'              masurca = path)
#' }
#'
#' @export
#'
run_masurca <- function(input1 = NULL,
                     input2 = NULL,
                     input_long = NULL,
                     output = NULL,
                     threads = 10,
                     parallel = FALSE,
                     cores = 4,
                     execute = TRUE,
                     masurca = NULL,
                     version = FALSE){
  # Check masurca program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", masurca)

  # Version
  if (isTRUE(version)){
    masurca.run <- sprintf('%s --version',
                           masurca)
    result <- system(masurca.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""
  # Threads
  if (!is.null(threads)){
    args <- paste(args,"--threads",threads,sep = " ")
  }
  # Output
  if (!is.null(output)){
    args <- paste(args,"--output",output,sep = " ")
  }

  # Create the run commands
  if (!is.null(input_long)){
    masurca.run <- sprintf('%s %s --illumina %s,%s --reads %s',
                           masurca,args,input1,input2,input_long)
  }else{
    masurca.run <- sprintf('%s %s --illumina %s,%s',
                           masurca,args,input1,input2)
  }

  # Run the commands, if execute is true
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, masurca.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(masurca.run, function (cmd)  system(cmd))
    }
  }

  # Return the list masurca commands
  return(masurca.run)
}

path <- "/software/MaSuRCA-4.1.1/bin/masurca"

run_masurca(masurca = path,
            version = TRUE)

fastq1 <- "illumina1.fq.gz"
fastq2 <- "illumina2.fq.gz"
nanopore <- "nanopore.fq.gz"

run_masurca(input1 = fastq1,
            input2 = fastq2,
            input_long = nanopore,
            threads = 40,
            masurca = path)
