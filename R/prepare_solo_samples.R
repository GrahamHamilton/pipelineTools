#' Title
#'
#' @description Reads in the sample file names for the directory or directories, for scRNA samples. Requires paired end reads.
#' @param path List of full paths to directory or directories containing the raw reads data
#' @param patt List of suffix patterns for the raw reads data for forward and reverse reads
#'
#' @examples
#' reads.path <- c("/Path to reads 1",
#'                 "/Path to reads 2")
#'
#' # Standard suffixes for reads files for forward and reverse files
#' reads.patt.1 <- "*_R1_001.fastq.gz$"
#' reads.patt.2 <- "*_R2_001.fastq.gz$"
#'
#' sample.dataframe <- prepare_solo_samples(reads.path, c(reads.patt.1,reads.patt.2))
#'
#' mate1 <- as.character(sample.dataframe$reads.path.1)
#' mate2 <- as.character(sample.dataframe$reads.path.2)
#'
#' sample.names <- as.character(sample.dataframe$sample.names)
#'
#' @return Dataframe with sorted file paths and sample names
#'
#' @export

prepare_solo_samples <- function(path = path, patt) {
  reads1 <-
    as.data.frame(matrix(
      list.files(
        path = path,
        pattern = patt[1],
        full.names = TRUE,
        recursive = TRUE
      )
    ))
  colnames(reads1) <- c("reads.path.1")
  reads1$sample.names <-
    unlist(lapply(strsplit(basename(
      as.character(reads1$reads.path.1)
    ), "_"), `[[`, 1))
  samples <- unique(reads1$sample.names)

  reads2 <-
    as.data.frame(matrix(
      list.files(
        path = path,
        pattern = patt[2],
        full.names = TRUE,
        recursive = TRUE
      )
    ))
  colnames(reads2) <- c("reads.path.2")
  reads2$sample.names <-
    unlist(lapply(strsplit(basename(
      as.character(reads2$reads.path.2)
    ), "_"), `[[`, 1))

  df <- data.frame()
  for (samp in samples) {
    d <-
      data.frame(
        "sample.names" = samp,
        "reads.path.1" = paste(as.character(reads1[reads1$sample.names == samp, ]$reads.path.1), collapse =
                                 ","),
        "reads.path.2" = paste(as.character(reads2[reads2$sample.names == samp, ]$reads.path.2), collapse =
                                 ",")
      )
    df <- rbind(df, d)
  }
  return(df)
}
