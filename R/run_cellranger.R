#' Run Cellranger
#'
#' @description Runs the Cellranger program, currently only runs cellranger count
#'
#' @param command Cellranger command
#' @param input Path to input FASTQ directory
#' @param run.id A unique run id and output folder name
#' @param description Sample description to embed in output files
#' @param reference Path of folder containing 10x-compatible reference
#' @param out.dir Output the results to this directory
#' @param sample Prefix of the filenames of FASTQs to select
#' @param create.bam Enable or disable BAM file generation,"true or false", default set to "false"
#' @param local.cores Set max cores the pipeline may request at one time, default set to 20
#' @param local.memory Set max GB the pipeline may request at one time, defaults set to 96Gb
#' @param parallel Run in parallel, default set to "false", quoted in lowercase
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param cellranger Path to the cellranger program, required
#' @param version Returns the version number
#'
#' @returns A list with the cellranger commands
#'
#' @examples
#' \dontrun{
#' command <- "count"
#' fastq_paths <- c("/export/sequencers/nextseq02/processed_data/2026_04/AurelieNajm/3421_VH00574_364_AurelieNajm_20260413/fastq/",
#'                  "/export/sequencers/nextseq02/processed_data/2026_04/AurelieNajm/3421_VH00574_366_AurelieNajm_20260427/fastq/")
#' ids <- c("GSYNC-067L","GSYNC-067L","GSYNC-067R","GSYNC-067R","GSYNC-100L","GSYNC-100L","GSYNC-65L","GSYNC-65L","GSYNC-65R",
#'          "GSYNC-65R","GSYNC-66R","GSYNC-66R","GSYNC-66L","GSYNC-66L","GSYNC-66R","GSYNC-66R","GSYNC-88L","GSYNC-88L","GSYNC-88R","GSYNC-88R")
#'
#' genome <- "/datastore/Genomes/Homo_sapiens/Ensembl/release-112/indexes/cellranger"
#'
#' path <- "/software/cellranger-10.0.0/cellranger"
#'
#' run_cellranger(command = command,
#'                input = fastq_paths,
#'                run.id = ids,
#'                reference = genome,
#'                sample = ids,
#'                cellranger = path,
#'                version = FALSE)
#' }
#' @export
#'

run_cellranger <- function(command = NULL,
                           input = NULL,
                           run.id = NULL,
                           description = NULL,
                           reference = NULL,
                           out.dir = NULL,
                           sample = NULL,
                           create.bam = "false",
                           local.cores = 20,
                           local.memory = 96,
                           parallel = FALSE,
                           cores = 4,
                           execute = FALSE,
                           cellranger = NULL,
                           version = FALSE){
  # Check cellranger program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", cellranger)

  # Version
  if (isTRUE(version)){
    cellranger.run <- sprintf('%s --version',
                            cellranger)
    result <- system(cellranger.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""
  # Sample
  if (!is.null(sample)){
    args <- paste(args,paste("--sample",sample,sep = "="), sep = " ")
  }
  # Run ID
  if (!is.null(run.id)){
    args <- paste(args,paste("--id",run.id,sep = "="), sep = " ")
  }
  #  Create bam files
  if (!is.null(create.bam)){
    args <- paste(args,"--create-bam",create.bam, sep = " ")
  }
  #  Local cores
  if (!is.null(local.cores)){
    args <- paste(args,paste("--localcore",local.cores,sep = "="), sep = " ")
  }
  # Local memory
  if (!is.null(local.memory)){
    args <- paste(args,paste("--localmem",local.memory,sep = "="), sep = " ")
  }

  if (command == "count"){
    cellranger.run <- sprintf('%s %s %s --transcriptome=%s --fastqs=%s',
                            cellranger,command,args,reference,input)
  }

  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, cellranger.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(cellranger.run, function (cmd)  system(cmd))
    }
  }

  return(cellranger.run)
}
