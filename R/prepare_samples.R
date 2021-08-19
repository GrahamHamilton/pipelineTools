#' Prepare samples
#'
#' @description Reads in the sample file names for the directory or directories. Ensures the files are in order, returns a dataframe with the paths
#'              to the raw reads files, the trimmed reads files and the sample names. Can take paired end or single end data.
#'              If there are duplicated sample names the files will be merged into a temporary directory.
#'              For scRNA samples, require paired end reads and do not require the trimmed reads directory. Call the function prepare_solo_samples.
#'
#' @param path List of full paths to directory or directories containing the raw reads data
#' @param patt List of suffix patterns for the raw reads data for forward and (optionally) reverse reads
#' @param trimmed.reads Name of the directory for the quality and adapter trimmed reads, optional
#' @param merge Boolean for whether to merge samples, default set to FALSE
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

prepare_samples <- function(path = NULL,
                            patt = NULL,
                            trimmed.reads = NULL,
                            merge = FALSE) {

  sample.names<- NULL

  df <-
    as.data.frame(matrix(
      list.files(
        path = path,
        pattern = patt[1],
        full.names = TRUE,
        recursive = FALSE
      )
    ))

  colnames(df) <- c("reads.path.1")
  df$trimmed.reads.path.1 <-
    paste(trimmed.reads, (basename(as.character(
      df$reads.path.1
    ))), sep = "/")
  df$sample.names <-
    unlist(lapply(strsplit(basename(
      as.character(df$reads.path.1)
    ), "_"), `[[`, 1))

  if (length(patt) > 1) {
    reads2 <-
      as.data.frame(matrix(
        list.files(
          path = path,
          pattern = patt[2],
          full.names = TRUE,
          recursive = FALSE
        )
      ))
    colnames(reads2) <- c("reads.path.2")
    reads2$trimmed.reads.path.2 <-
      paste(trimmed.reads, (basename(as.character(
        reads2$reads.path.2
      ))), sep = "/")
    reads2$sample.names <-
      unlist(lapply(strsplit(basename(
        as.character(reads2$reads.path.2)
      ), "_"), `[[`, 1))

    if (anyDuplicated(df$sample.names)) {
      df <- cbind(df, reads2)
    }else{
      df <- merge(df, reads2, by = "sample.names")
    }
  }

  if ((anyDuplicated(df$sample.names)) && merge){
    # Get duplicate names
    dup.names <- df[duplicated(df$sample.names),]$sample.names
    # Create a temporary directory for combined files
    temp.out.dir <- "tmp"
    dir.create(temp.out.dir, showWarnings = FALSE)
    for(name in dup.names){
      sub <- subset(df, sample.names == name)
      join <- paste(as.character(sub$reads.path.1), collapse=' ' )
      combine <- file.path(temp.out.dir,
                           unlist(strsplit(as.character(sub$reads.path.1[1]),"/"))[length(unlist(strsplit(as.character(sub$reads.path.1[1]),"/")))],
                           fsep = .Platform$file.sep)
      cat.files.run <- sprintf('%s %s > %s',
                               "cat",join,combine )
      lapply(cat.files.run, function (cmd)  system(cmd))

      df$reads.path.1 <- as.character(df$reads.path.1)
      df[df$sample.names == name,"reads.path.1"] <- combine

      # If the reads are paired end
      if (length(sub$reads.path.2) > 0){
        join <- paste(as.character(sub$reads.path.2), collapse=' ' )
        combine <- file.path(temp.out.dir,
                             unlist(strsplit(as.character(sub$reads.path.2[1]),"/"))[length(unlist(strsplit(as.character(sub$reads.path.2[1]),"/")))],
                             fsep = .Platform$file.sep)
        cat.files.run <- sprintf('%s %s > %s',
                                 "cat",join,combine )
        lapply(cat.files.run, function (cmd)  system(cmd))

        df$reads.path.2 <- as.character(df$reads.path.2)
        df[df$sample.names == name, "reads.path.2"] <- combine
      }
    }
    df <- df[!duplicated(df$sample.names),]
  }

  return(df)
}
