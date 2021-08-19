#' Run UMITools
#'
#' @description UMITools are a set of tools for dealing with Unique Molecular Identifiers.
#' Runs the commands whitelist, extract, group, dedup, count and count_tab.
#'
#' @param command Umitools command
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param umitools Path to the Umitools program, required
#' @param version Returns the version number
#'
#' @return A file with the Umitools commands
#' @export
#'
#' @examples
#'  \dontrun{
#'  # Version number
#'  run_umitools(umitools = umitools.path,
#'                          version = TRUE)
#'}
#'
run_umitools <- function(command = NULL,

                         parallel = FALSE,
                         cores = 4,
                         umitools = NULL,
                         version = FALSE){
  # Check umitools program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", umitools)

  # Version
  if (isTRUE(version)){
    umitools.run <- sprintf('%s --version',
                            umitools)
    result <- system(umitools.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""

  if (isTRUE(parallel)){
    cluster <- makeCluster(cores)
    parLapply(cluster, umitools.run, function (cmd)  system(cmd))
    stopCluster(cluster)
  }else{
    lapply(umitools.run, function (cmd)  system(cmd))
  }

  return(umitools.run)
}


umitools.path <- "/software/anaconda3/bin/umi_tools"

run_umitools(umitools = umitools.path,
             version = TRUE)
