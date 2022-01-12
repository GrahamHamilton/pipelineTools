#' Run GATK
#'
#' @description Runs the GATK suite of programs.
#'
#' @param command GATK command to run, required
#' @param input List of sorted bam files, required
#' @param output List of output files
#' @param reference Path to the fasta formatted reference
#' @param intervals Path to the intrvals list file, usually the exome coordiantes file
#' @param known.sites List of paths to the files containing known polymorphic sites
#' @param bsqr List of base quality recalibration files
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param gatk Path to the GATK suit of programs, required
#'
#' @examples
#'
#' \dontrun{
#' known.site.list <- c("Mills_and_1000G_gold_standard.indels.GRCh38.vcf.gz",
#'                       "1000G_phase1.snps.high_confidence.GRCh38.vcf.gz")
#'
#' recalibration.files <- gsub(".bam","_recal_data.table",rg.bam.files)
#' command <-  "BaseRecalibrator"
#' fasta <- Path the genome fats
#' # BaseRecalibration
#' base.recalibration.cmds <- run_gatk(command = command,
#'                                     input = rg.bam.files,
#'                                     output = recalibration.files,
#'                                     reference = fasta,
#'                                     intervals = intervals.file,
#'                                     known.sites = known.site.list,
#'                                     execute = TRUE,
#'                                     gatk = gatk)
#' }
#'
#' @return List of GATK commands
#'
#' @export
#'

run_gatk <- function(command = NULL,
                     input = NULL,
                     output = NULL,
                     reference = NULL,
                     intervals = NULL,
                     known.sites = NULL,
                     bsqr = NULL,
                     parallel = FALSE,
                     cores = 4,
                     execute = TRUE,
                     gatk = NULL){
  # Check gatk program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", gatk)

  # Set the additional arguments
  args <- ""

  if (!is.null(known.sites)){
    for (site in known.sites){
      args <- paste(args,"--known-sites",site, sep = " ")
    }
  }

  # BaseRecalibrator
  if (command == "BaseRecalibrator"){
    if (is.null(intervals)){
      gatk.run <- sprintf('%s %s -I %s -O %s -R %s %s',
                          gatk,command,input,output,reference,args)
    }else{
      gatk.run <- sprintf('%s %s -I %s -O %s -R %s -L %s  %s',
                          gatk,command,input,output,reference,intervals,args)
      }
  }

  # ApplyBSQR
  if (command == "ApplyBQSR"){
    if (is.null(intervals)){
      gatk.run <- sprintf('%s %s -I %s -O %s -R %s --bsqr-recal-file %s %s',
                          gatk,command,input,output,reference,bsqr,args)
    }else{
      gatk.run <- sprintf('%s %s -I %s -O %s -R %s --bsqr-recal-file %s -L %s  %s',
                          gatk,command,input,output,reference,bsqr,intervals,args)
    }
  }

  # HaplotypeCaller
  if (command == "HaplotypeCaller"){
    if (is.null(intervals)){
      gatk.run <- sprintf('%s %s -I %s -O %s -R %s  %s',
                          gatk,command,input,output,reference,args)
    }else{
      gatk.run <- sprintf('%s %s -I %s -O %s -R %s -L %s -ERC GVCF  %s',
                          gatk,command,input,output,reference,intervals,args)
    }
  }

  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, gatk.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(gatk.run, function (cmd)  system(cmd))
    }
  }

  # Return the list of Picard commands
  return(gatk.run)
}
