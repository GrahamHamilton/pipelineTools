#' Run BLAST programs
#'
#' @description Runs the BLAST progam and creates BLAST databases
#'
#' @import parallel
#'
#' @param input Fasta formatted file for input to makeblastdb or query input fasta for blast
#' @param output Name of output file for blast
#' @param title Title for the blast database
#' @param database Database type e.g. "nucl" 0r "prot", for makeblastdb or database name for blast
#' @param format Output format for blast, numeric e.g 6 is a tabular format
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param blast Path to the BLAST program
#' @param version Returns the version number
#'
#' @return A list with the BLAST commands

#'
#' @examples
#' \dontrun{
#' path_db <- "/software/blast-v2.13.0/bin/makeblastdb"
#' path_blast <- "/software/blast-v2.13.0/bin/blastn" or "/software/blast-v2.13.0/bin/blastp"
#'
#' version <- run_blast(blast = path_blast,
#'                      version =TRUE)
#' version[1]
#'
#' run_blast(input = "sample1.fasta",
#'           title = "sample_name",
#'           database = "nucl",
#'           blast = path_db)
#'
#' run_blast(input = "cbp1.fasta",
#'           output = "sample1.blast",
#'           database = "blastDB",
#'           format = 6,
#'           blast = path_blast)
#' }
#'
#' @export

run_blast <- function(input = NULL,
                      output = NULL,
                      title = NULL,
                      database = NULL,
                      format = NULL,
                      parallel = FALSE,
                      cores = 4,
                      execute = TRUE,
                      blast = NULL,
                      version = FALSE){
  # Check blast program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", blast)

  # Version
  if (isTRUE(version)){
    blast.run <- sprintf('%s -version',
                            blast)
    result <- system(blast.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments, not in usee at present
  args <- ""

  #
  if (grepl("blastn | blastp",blast)){
    blast.run <- sprintf('%s -db %s -out %s -query %s -outfmt %s',
                         blast,database,output,input,format)
  }else if(grepl("makeblastdb",blast)){
    blast.run <- sprintf('%s -in %s -parse_seqids -title "%s" -dbtype %s',
                         blast,input,title,database)
  }

  # Run the commands, if execute is true
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, blast.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(blast.run, function (cmd)  system(cmd))
    }
  }

  return(blast.run)
}
