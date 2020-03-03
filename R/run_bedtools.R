#' Run the Bedtools program
#'
#' @description Runs the bedtools program, currently only supports coverage
#'
#' @param command Bedtools command to run, at present can only choose from 'coverage' or 'bamtobed', required
#' @param input List of aligned files in bam format, required
#' @param reference GTF file fro calculating the depth of coverage intervals
#' @param out.dir Name of the directory from the Bedtools output
#' @param sample.names List of sample names, required
#' @param counts Only report the count of overlaps, default set to FALSE
#' @param bedtools Path to the bedtools program, required
#' @param version Returns the version number
#'
#' @examples
#' \dontrun{
#' # Version
#' bedtools.path <- "/software/bedtools2/bin/bedtools"
#' bedtools.version <- run_bedtools(bedtools = bedtools.path,
#'                                  version = TRUE)
#'
#' # Bedtools coverage
#' command = "coverage"
#' outputDirectory <- "coverage"
#' sam.files <- list.files(path = alignments.path, pattern = "sam$",
#'                         full.names = TRUE,
#'                         recursive = TRUE)
#' sorted.bam.files <- gsub(sam.files,
#'                          pattern = ".sam",
#'                          replacement = "_sorted.bam")
#' sample_names <- unlist(lapply(strsplit(list.files(path = reads_path,
#'                        pattern = "*_R1_001.fastq$",
#'                        full.names = FALSE),"_"), `[[`, 1))
#' gtfFile <- "mirBase/hsa.ensembl.gff3"
#'
#' bedtools.cmds <- run_bedtools(command = command,
#'                               input = sorted.bam.files,
#'                               reference = gtfFile,
#'                               out.dir = outputDirectory,
#'                               sample.names = sample.names,
#'                               counts = TRUE,
#'                               bedtools = bedtools.path)
#' bedtools.cmds
#' }
#'
#' @return A list with the bedtools commands
#'
#' @export
#'

run_bedtools <- function(command = NULL,
                         input = NULL,
                         reference = NULL,
                         out.dir = NULL,
                         sample.names = NULL,
                         counts = FALSE,
                         bedtools = NULL,
                         version = FALSE){
  # Check bedtools program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", bedtools)

  # Version
  if (isTRUE(version)){
    bedtools.run <- sprintf('%s --version',
                           bedtools)
    result <- system(bedtools.run, intern = TRUE)
    return(result)
  }

  # Create the sample directories for the per sample bedtools results
  lapply(paste(out.dir,sample.names, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))

  # Set the additional arguments
  args <- ""
  # Counts
  if (isTRUE(counts)){
    args <- paste(args,"-counts", sep = " ")
  }

 if (command == "coverage"){
   # Set the names and paths for the bedtools counts files
   out.files <- paste(out.dir,sample.names,paste(sample.names,"txt",sep = "."),sep = "/")

   bedtools.run <- sprintf('%s %s %s -a %s -b %s > %s',
                          bedtools,command,args,reference,input,out.files)
   }

  if (command == "bamtobed"){
    # Set the names and paths for the bedtools bamtobed files
    out.files <- paste(out.dir,sample.names,paste(sample.names,"bed",sep = "."),sep = "/")

    bedtools.run <- sprintf('%s %s %s -i %s > %s',
                            bedtools,command,args,input,out.files)
  }

  lapply(bedtools.run, function (cmd)  system(cmd))
  return(bedtools.run)

}

