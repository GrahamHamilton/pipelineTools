#' Run Stringtie
#'
#' @description Runs the stringtie tool
#'
#' @param input List of aligned files in sorted bam format, required
#' @param threads Number of threads for stringtie to use, default set to 10
#' @param trimming Disable trimming of predicted transcripts based on coverage, by default coverage trimming is enabled
#' @param strandedness Strand spcific reads, values are "first" or "second"
#' @param reference.gtf Path to the gtf file for guiding the assembly process
#' @param estimate Only estimate the abundance of given reference transcripts, requires the reference GTF
#' @param merge Assemble transcripts from multiple input files generating a unified non-redundant set of isoforms
#' @param out Name of the directory from the Kallisto output. If NULL,
#'            which is the default, a directory named "stringtie" is created in the current working directory.
#' @param sample.name List of the sample names, required
#' @param ballgown Enable output of Ballgown table files which will be created in the
#'                 same directory as the output GTF (requires -G, -o recommended)
#' @param stringtie Path to the stringtie program, required
#' @param version Returns the version number
#'
#' @return A list with the stringtie commands
#'
#' @examples
#'  \dontrun{
#' # Version
#' stringtie.commands <- run_stringtie(stringtie = "/software/stringtie-1.3.6/stringtie",
#'                                     version = TRUE)
#' stringtie.commands
#'
#' # Test command
#' sample_names <- unlist(lapply(strsplit(list.files(path = "trimmed_reads",
#'                        pattern = "*_R1_001.fastq$", full.names = FALSE),"_"), `[[`, 1))
#' gtf <- "/export/buzz1/Genome/Homo_sapiens/Ensembl/GRCh38_release_79/Sequence/Transcriptome/
#'          WholeTranscriptome/Homo_sapiens_GRCh38_79_primary_assembly_whole_transcriptome.gff"
#' bam.files <- list.files(path = hisat.alignments.dir, pattern = "sorted.bam$",
#'                         full.names = TRUE, recursive = TRUE)
#' stringtie.commands <- run_stringtie(input = bam.files,
#'                                     trimming = FALSE,
#'                                     strandedness = "first",
#'                                     reference.gtf = gtf,
#'                                     estimate = TRUE,
#'                                     out = "stringtie",
#'                                     sample.name = sample_names,
#'                                     stringtie = "/software/stringtie-1.3.6/stringtie")
#' stringtie.commands
#' }
#'
#' @export
#'
run_stringtie <- function(input = NULL,
                          threads = 10,
                          trimming = TRUE,
                          strandedness = NULL,
                          reference.gtf = NULL,
                          estimate = FALSE,
                          merge = FALSE,
                          out = NULL,
                          sample.name = NULL,
                          ballgown = FALSE,
                          stringtie = NULL,
                          version = FALSE){
  # Check stringtie program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", stringtie)

  # Version
  if (isTRUE(version)){
    stringtie.run <- sprintf('%s --version',
                            stringtie)
    result <- system(stringtie.run, intern = TRUE)
    return(result)
  }

  # Create the directory to write the data, if it is not present
  if (is.null(out) && !isTRUE(merge)){
    out <- "stringtie"
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
  }

  # Create the sample directories for the per sample stringtie results
  if (!isTRUE(merge)){
    lapply(paste(out,sample.name, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))
    }

  # Set the additional arguments
  args <- ""
  # Merge
  if (isTRUE(merge)){
    args <- paste(args,"--merge",sep = " ")
    }
  # Threads
  if (!is.null(threads)){
    args <- paste(args,"-p",threads,sep = " ")
    }
  # Trimming transcripts based on coverage
  if (!isTRUE(trimming)){
    args <- paste(args,"-t",sep = " ")
    }
  # Strandedness
  if (!is.null(strandedness)){
    # First strand
    if (strandedness == "first"){
      args <- paste(args,"--rf",sep = " ")
      }else if (strandedness == "second"){
      args <- paste(args,"--fr",sep = " ")
      }
    }
  # Reference GTF as guide
  if (!is.null(reference.gtf)){
    args <- paste(args,"-G",reference.gtf,sep = " ")
    }
  # Estimate abundances of known transcripts
  if (isTRUE(estimate) && !is.null(reference.gtf)){
    args <- paste(args,"-e",sep = " ")
  }
  # Ballgown
  if (isTRUE(ballgown)){
    args <- paste(args,"-B",sep = " ")
  }

  finalGTF <- paste(out,sample.name,paste(sample.name,"transcripts.gtf",sep = "_"),sep = "/")

  # Assemble commands
  if (!isTRUE(merge)){
    stringtie.run <- sprintf('%s %s -o %s %s',
                             stringtie,args,finalGTF,input)
    }
  # Merge commands
  else if (isTRUE(merge)){
    stringtie.run <- sprintf('%s %s -o %s %s',
                             stringtie,args,out,input)
    }


  lapply(stringtie.run, function (cmd)  system(cmd))

  return(stringtie.run)
}
