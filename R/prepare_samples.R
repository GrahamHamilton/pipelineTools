#' Prepare samples
#'
#' @description Reads in the sample file names for the directory or directories. Ensures the files are in order, returns a dataframe with the paths
#'              to the raw reads files, the trimmed reads files and the sample names. Can take paired end or single end data.
#'
#' @param path List of full paths to directory or directories containing the raw reads data
#' @param patt List of suffix patterns for the raw reads data for forward and (optionally) reverse reads
#' @param trimmed.reads Name of the directory for the quality and adapter trimmed reads
#'
#' @examples
#' \dontrun{
#' # Paths to reads
#' reads.path <- c("/Path to reads 1",
#'                 "/Path to reads 2")
#'
#' # Standard suffixes for reads files for forward and reverse files
#' reads.patt.1 <- "*_R1_001.fastq.gz$"
#' reads.patt.2 <- "*_R2_001.fastq.gz$"
#'
#' # Trimmed reads directory
#' trimmed.reads <- "trimmed_reads"
#'
#' # Paired end data
#' sample.dataframe <- prepare_samples(reads.path, c(reads.patt.1,reads.patt.2),trimmed.reads)
#' mate1 <- as.character(sample.dataframe$reads.path.1)
#' mate1.trim <- as.character(sample.dataframe$trimmed.reads.path.1)
#' mate2 <- as.character(sample.dataframe$reads.path.2)
#' mate2.trim <- as.character(sample.dataframe$trimmed.reads.path.2)
#'
#' # Single end data
#' sample.dataframe <- prepare_samples(reads.path, reads.patt.1,trimmed.reads)
#' mate1 <- as.character(sample.dataframe$reads.path.1)
#' mate1.trim <- as.character(sample.dataframe$trimmed.reads.path.1)
#'
#' # Sample names
#' sample.names <- as.character(sample.dataframe$sample.names)
#' }
#'
#' @return Dataframe with sorted file paths and sample names
#'
#' @export
#'

prepare_samples<- function(path,patt,trimmed.reads){
  reads1 <- as.data.frame(matrix(list.files(path = path, pattern = patt[1], full.names = TRUE, recursive = TRUE)))
  colnames(reads1) <- c("reads.path.1")
  reads1$trimmed.reads.path.1 <- paste(trimmed.reads,(basename(as.character(reads1$reads.path.1))),sep = "/")
  reads1$sample.names <- unlist(lapply(strsplit(basename(as.character(reads1$reads.path.1)),"_"), `[[`, 1))

  if (length(patt) > 1){
    reads2 <- as.data.frame(matrix(list.files(path = path, pattern = patt[2], full.names = TRUE, recursive = TRUE)))
    colnames(reads2) <- c("reads.path.2")
    reads2$trimmed.reads.path.2 <- paste(trimmed.reads,(basename(as.character(reads2$reads.path.2))),sep = "/")
    reads2$sample.names <- unlist(lapply(strsplit(basename(as.character(reads2$reads.path.2)),"_"), `[[`, 1))
    df <- merge(reads1,reads2, by = "sample.names")
    return(df)
  }else{
    return(reads1)
  }
}
