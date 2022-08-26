#' Run SnpEff
#'
#' @description Runs the SnpEff program.
#'
#' @param command SnpEff command to run
#' @param input Input file, vcf format
#' @param output Name of output file, vcf format
#' @param genome Name of the genome
#' @param config Path to the configuration file
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param snpeff Path to the SnpEff program, required
#' @param version Returns the version number
#'
#' @return A list with the snpEff commands
#'
#' @examples
#' \dontrun{
#'  snpeff.path <- "/software/snpEff-v5.0e/snpEff.jar"
#'
#'  # Version
#'  res <- run_snpeff(snpeff = snpeff.path,
#'                    version = TRUE)
#'
#'  # Download a premade genome
#'  command <- "download"
#'  genome <- "GRCh38.99"
#'
#'  res <- run_snpeff(command = command,
#'                    genome = genome,
#'                    execute = FALSE,
#'                    snpeff = snpeff.path)
#'
#'  # Run snpEff
#'  res <- run_snpeff(command = command,
#'                    genome = genome,
#'                    input = "input.vcf",
#'                    output = "annotated.vcf",
#'                    execute = FALSE,
#'                    snpeff = snpeff.path)
#' }
#'
#' @export
#'
run_snpeff <- function(command = NULL,
                       input = NULL,
                       output = NULL,
                       genome = NULL,
                       config = NULL,
                       parallel = FALSE,
                       cores = 4,
                       execute = TRUE,
                       snpeff = NULL,
                       version = FALSE){
  # Version
  if (isTRUE(version)){
    snpeff.run <- sprintf('java -jar %s -version',
                          snpeff)
    result <- system(snpeff.run, intern = TRUE)
    return(result)
  }

  # Download database
  if (command == "download"){
    snpeff.run <- sprintf('java -jar %s %s -v %s',
                          snpeff,command,genome)
  }

  # Annotate vcf files
  if (command == "ann"){
    # Name of stats file
    stats <- gsub(".vcf","_snpEff_summary.html",input)

    snpeff.run <- sprintf('java -jar %s -config %s %s %s -s %s > %s',
                          snpeff,config,genome,input,stats,output)
  }

  # Run the commands, if execute is true
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, snpeff.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(snpeff.run, function (cmd)  system(cmd))
    }
  }

  # Return the list of snpeff commands
  return(snpeff.run)
}

