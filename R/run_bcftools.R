#' Run bcftools
#'
#' @param command bcftools command to run, at present can only choose isec
#' @param input vcffiles
#' @param output ouput directory or file
#' @param regions Bed formatted regions of interest file
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param bcftools Path to the bcftools program, required
#' @param version Returns the version number
#'
#' @return A list with the bcftools commands
#'
#' @examples
#' \dontrun{
#' bcftools.path <- "/software/bcftools-v1.9/bin/bcftools"
#'
#' bcftools.version <- run_bcftools(bcftool = bcftools.path,
#'                                  version = TRUE)
#'
#' bcftools.version[1]
#'
#' # bcftools intersect
#' command <- "isec"
#'
#' vcfs <- c("test1.vcf","test2.vcf")
#' vcfs <- paste(vcfs, collapse = " ")
#'
#' out.directory <- "test_directory"
#'
#'
#'
#' bcftools.cmd <- run_bcftools(command = command,
#'                              input = vcfs,
#'                              output = out.directory,
#'                              execute = FALSE,
#'                              bcftools = bcftools.path)
#' bcftools.cmd
#' }
#'
#' @export
run_bcftools <- function(command = NULL,
                         input = NULL,
                         output = NULL,
                         regions = NULL,
                         parallel = FALSE,
                         cores = 4,
                         execute = TRUE,
                         bcftools = NULL,
                         version = NULL){
  # Check bcftools program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", bcftools)

  # Version
  if (isTRUE(version)){
    bcftools.run <- sprintf('%s --version',
                            bcftools)
    result <- system(bcftools.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""

  # Regions
  if (!is.null(regions)){
    args <- paste(args,"-R",regions,sep = " ")
  }


  if(command == "isec"){
    bcftools.run <- sprintf('%s %s -p %s  %s %s',
                            bcftools,command,output,input,args)
  }

  if(command == "view"){
    bcftools.run <- sprintf('%s %s %s %s > %s',
                            bcftools,command,args,input,output)
  }


  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, bcftools.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(bcftools.run, function (cmd)  system(cmd))
    }
  }

  return(bcftools.run)
}
