#' Run Picard
#'
#' @description Runs the Picard program. Curently only works for CollectRnaSeqMetrics & MarkDuplicates command
#'
#' @param command Picard command to run, required
#' @param input List of sorted bam files, required
#' @param output List of output bam files, required for MarkDuplicates command
#' @param out.dir Name of the output directory, required
#' @param refFlat Path to the refFlat file, required for CollectRnaSeqMetrics commmand
#' @param rRNA.intervals Path to the rRNAintrvals list file
#' @param strand Strand-specific information, "FR" for first strand or "RF" for reverse strand,
#'               required for CollectRnaSeqMetrics commmand
#' @param remove.duplicates Set for removing marked duplicates, boolean default set to FALSE
#' @param sample.name List of the sample names, required
#' @param picard Path to the Picard program, required
#'
#' @examples
#' \dontrun{
#' sorted.bam.files <- list.files(path = hisat2.alignments.dir, pattern = "sorted.bam$",
#'                                full.names = TRUE, recursive = TRUE)
#' sample_names <- unlist(lapply(strsplit(list.files(path = trimmed_reads_dir,
#'                        pattern = "*_R1_001.fastq$",
#'                        full.names = FALSE),"_"), `[[`, 1))
#' picard.dir <- "picard"
#' refFlat.file <- "refFlatBDGP6.txt"
#'
#' run_picard(command = "CollectRnaSeqMetrics",
#'            input = sorted.bam.files,
#'            out.dir = picard.dir,
#'            refFlat = refFlat.file,
#'            strand = strandedness,
#'            sample.name = sample.names,
#'            picard = picard.path)
#' }
#'
#' @return A list with the Picard commands
#'
#' @export
#'

run_picard <- function(command = NULL,
                       input = NULL,
                       output = NULL,
                       out.dir = NULL,
                       refFlat = NULL,
                       rRNA.intervals = NULL,
                       strand = NULL,
                       remove.duplicates = FALSE,
                       sample.name = NULL,
                       picard = NULL
                       ){
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", picard)

  # Set the additional arguments
  args <- ""

  if (command == "CollectRnaSeqMetrics"){
    # Strand
    if (!is.null(strand)){
      if (strand == "FR" || strand == "F"){
        args <- paste(args,"STRAND=FIRST_READ_TRANSCRIPTION_STRAND",sep = " ")
      }
      if (strand == "RF" || strand == "R"){
        args <- paste(args,"STRAND=SECOND_READ_TRANSCRIPTION_STRAND",sep = " ")
      }
    }else{
      args <- paste(args,"STRAND=NONE",sep = " ")
    }
    # RNA intervals
    if (!is.null(rRNA.intervals)){
      args <- paste(args,paste("RIBOSOMAL_INTERVALS=",rRNA.intervals,sep = ""),sep = " ")
    }

    # Create the output files list
    picard.file <- paste(out.dir,paste(sample.name,"metrics.txt",sep = "."),sep = "/")

    # Create the Picard commands
    picard.run <- sprintf('java -jar %s %s I=%s O=%s REF_FLAT=%s %s',
                          picard,command,input,picard.file,refFlat,args)
  }

  if (command == "MarkDuplicates"){
   # Create the output metric files list
    metric.files <- paste(out.dir,paste(sample.name,"dup","metrics.txt",sep = "."),sep = "/")

    # Create the Picard commands
    picard.run <- sprintf('java -jar %s %s I=%s O=%s M=%s REMOVE_DUPLICATES=%s',
                          picard,command,input,output,metric.files,remove.duplicates)
  }

  # Run the Picard commands
  lapply(picard.run, function (cmd)  system(cmd))

  # Return the list of Picard commands
  return(picard.run)
}
