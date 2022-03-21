#' Run Freebayes
#'
#' @param input List of aligned bam files
#' @param reference Reference genome sequence in fasta format
#' @param output List of file names for out put
#' @param ploidy Sets the default ploidy for the analysis, default st to 2
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param freebayes Path to the freebayes program, required
#' @param version Returns the version number
#'
#' @examples
#' \dontrun{
#'  freebayes <- "/usr/bin/freebayes"
#'
#'  freebayes.verion <- run_freebayes.1(freebayes = freebayes,
#'                                    version = TRUE)
#'
#'  freebayes.cmds <- run_freebayes(input = "input.bam",
#'                                  reference = "reference.fa",
#'                                  output = "out.vcf",
#'                                  ploidy = 2,
#'                                  execute = FALSE,
#'                                  freebayes = freebayes)
#' }
#'
#' @return A list of freebayes commands
#'
#' @export

run_freebayes <- function(input = NULL,
                          reference = NULL,
                          output = NULL,
                          ploidy = NULL,
                          parallel = FALSE,
                          cores = 4,
                          execute = TRUE,
                          freebayes = NULL,
                          version = FALSE){
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", freebayes)

  # Version
  if (isTRUE(version)){
    freebayes.run <- sprintf('%s --version',
                             freebayes)
    result <- system(freebayes.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""
  # Ploidy
  if (!is.null(ploidy)){
    args <- paste(args,"--ploidy",ploidy,sep = " ")
  }

  freebayes.run <- sprintf('%s %s --fasta-reference %s --bam %s > %s',
                           freebayes,args,reference,input,output)

  # Run the Freebayes commands
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, freebayes.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(freebayes.run, function (cmd)  system(cmd))
    }
  }

  # Return the list of freebayes commands
  return(freebayes.run)
}
