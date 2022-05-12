#' Run SnpSift
#'
#' @description Runs the SnpSift program.
#'
#' @param command SnpSift command to run, currently only fiter, annotate and vartype
#' @param input Input file, vcf format
#' @param output Name of output file, vcf format
#' @param filter Filter expression
#' @param dbsnp Path to the dbsnp file
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param snpsift Path to the SnpSift program, required
#'
#' @return List of SnpSift commands
#' @examples
#' \dontrun{
#' snpsift.path <- "/software/snpEff-v5.0e/SnpSift.jar"
#'
#' run_snpsift(command = "filter",
#'             input = "file.vcf",
#'             output = "filtered.file.vcf",
#'             filter = "'( QUAL >= 30 )'",
#'             execute = FALSE,
#'             snpsift = snpsift.path)
#'
#' run_snpsift(command = "annotate",
#'             input = "file.vcf",
#'             output = "annotated.file.vcf",
#'             dbsnp = "/datastore/GATK_resource_bundle/GRCh38/dbsnp_146.GRCh38.vcf.gz",
#'             execute = FALSE,
#'             snpsift = snpsift.path)
#'
#' run_snpsift(command = "vartype",
#'             input = "file.vcf",
#'             output = "annotated.file.vcf",
#'             execute = FALSE,
#'             snpsift = snpsift.path)
#'
#' }
#' @export
#'
run_snpsift <- function(command = NULL,
                        input = NULL,
                        output = NULL,
                        filter = NULL,
                        dbsnp = NULL,
                        parallel = FALSE,
                        cores = 4,
                        execute = TRUE,
                        snpsift = NULL){

  # Filter
  if (command == "filter"){
    snpsift.run <- sprintf('java -jar %s %s %s %s > %s',
                          snpsift,command,filter,input,output)
  }
  # Annotate
  if (command == "annotate"){
    snpsift.run <- sprintf('java -jar %s %s %s %s > %s',
                           snpsift,command,dbsnp,input,output)
  }
  # VarType
  if (command == "vartype"){
    snpsift.run <- sprintf('java -jar %s %s %s > %s',
                           snpsift,command,input,output)
  }

  # Run the commands, if execute is true
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, snpsift.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(snpsift.run, function (cmd)  system(cmd))
    }
  }

  # Return the list of snpeff commands
  return(snpsift.run)
}
