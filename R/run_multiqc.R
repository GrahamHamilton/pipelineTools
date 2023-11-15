#' Run MultiQC
#'
#' @param workingDir Current working directory, parent directory for the logfiles, required
#'
#' @param multiqc Path to the MultiQC program, required
#' @param version Returns the version number
#'
#' @examples
#'  \dontrun{
#'  workingDir <- getwd()
#'  run_multiqc(workingDir,
#'              multiqc = "/usr/local/bin/multiqc")
#'  }
#' @export
#'
run_multiqc <- function(workingDir = NULL,
                        multiqc = NULL,
                        version = FALSE){
  # Check multiqc program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", multiqc)

  # Version
  if (isTRUE(version)){
    multiqc.run <- sprintf('%s --version',
                           multiqc)
    result <- system(multiqc.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""
  # Force
  args <- paste(args,"--force",sep = " ")

  # Create the multiqc command
  multiqc.run <- sprintf('%s %s %s',
                         multiqc,args,workingDir)

  # Run the multiqc command
  system(multiqc.run)
}
