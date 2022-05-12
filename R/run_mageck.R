#' Run Mageck
#'
#' @param command Mageck command to run, currently only count is implemented
#' @param umi.list file contating the list of UMIs
#' @param fastq Trimmed sample files in fastq format, should only contain the umi sequence after trimming with cutadapt
#' @param sample.label Sample name
#' @param prefix Name to add to results files
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param version Returns the version number
#' @param mageck Pathe to the mageck program, required
#'
#' @return List of Mageck commands
#'
#' @examples
#' \dontrun{
#' mageck.version <- run_mageck(version = TRUE,
#'                             mageck = path)
#' mageck.version
#'
#' command <- "count"
#' list <- "umi-list.tsv"
#' reads <- "MPRA_R1.fastq.gz"
#' label <- "MPRA"
#' prefix <- "MPRA"
#'
#' mageck.cmds <- run_mageck(command = command,
#'                           umi.list = list,
#'                           fastq = reads,
#'                           sample.label = label,
#'                           prefix = prefix,
#'                           execute = FALSE,
#'                           mageck = path)
#'
#' mageck.cmds
#'
#' }
#' @export
#'

run_mageck <- function(command = NULL,
                       umi.list = NULL,
                       fastq = NULL,
                       sample.label = NULL,
                       prefix = NULL,
                       parallel = FALSE,
                       cores = 4,
                       execute = TRUE,
                       version = NULL,
                       mageck = NULL){
  # Version
  if (isTRUE(version)){
    mageck.run <- sprintf('%s --version',
                            mageck)
    result <- system(mageck.run, intern = TRUE)
    return(result)
  }

  # Count
  if (command == "count"){
    mageck.run <- sprintf(' %s %s --list-seq %s --sample-label %s --fastq %s --output-prefix %s',
                          mageck,command,umi.list,sample.label,fastq,prefix)
  }
}

