#' Run the Star program
#'
#' @description Runs the Star alignment program, can be run in parallel on multiple cores.
#'
#' @import parallel
#'
#' @param mate1 List of the paths to files containing to the forward reads
#' @param mate2 List of the paths to files containing to the reverse reads
#' @param genome Path to the reference genome index directory
#' @param sample.name List of the sample names
#' @param out.dir Name of the directory from the Star output
#' @param out.format Format of output file, default set to BAM
#' @param unmapped How to deal with unmapped reads, default set to "Within" which write unmapped reads to SAM/BAM file
#' @param sam.attributes Alignment attributes for the SAM/BAM file, default set to "Standard"
#' @param quant.mode Type of quantification required, default set to "GeneCounts"
#' @param compressed Compression mode for input reads files, default set to "zcat" for gzipped files, can use "bzcat" for bz2 files
#' @param threads Number of threads
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param star Path to the Star program
#' @param version Returns the version number
#'
#' @return A list with the Star commands
#'
#' @examples
#' \dontrun{
#' path <- "/full/path/to/program"
#' genome <- "/full/path/to/genome"
#'
#' mate1.trim <- List of paths to trimmed forward reads for alignment
#' mate2.trim <- List of paths to trimmed reverse reads for alignment
#' sample.names <- List os sample names
#'
#' cmds <- run_star(mate1 = mate1.trim,
#'                  mate2 = mate2.trim,
#'                  genome = genome,
#'                  sample.name = sample.names,
#'                  out.dir = results.dir,
#'                  unmapped = "Within",
#'                  sam.attributes = "Standard",
#'                  quant.mode = "GeneCounts",
#'                  parallel = TRUE,
#'                  cores = 4,
#'                  star = path)
#'
#' # Version number
#' version <- run_star(star = path,
#'                     version = TRUE)
#' }
#'
#' @export
#'

run_star <- function(mate1 = NULL,
                     mate2 = NULL,
                     genome = NULL,
                     sample.name = NULL,
                     out.dir = NULL,
                     out.format = "BAM",
                     unmapped = "Within",
                     sam.attributes = "Standard",
                     quant.mode = "GeneCounts",
                     compressed = "zcat",
                     threads = 10,
                     parallel = FALSE,
                     cores = 4,
                     star = NULL,
                     version = FALSE
                     ){
  # Check star program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", star)

  # Version
  if (isTRUE(version)){
    star.run <- sprintf('%s --version',
                           star)
    result <- system(star.run, intern = TRUE)
    return(result)
  }

  # Create the sample directories for the per sample bowtie2 results
  lapply(paste(out.dir,sample.name, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))

  # Set the additional arguments
  args <- ""
  # Aligment output format
  if (!is.null(out.format)){
    args <- paste(args, "--outSAMtype", out.format, "SortedByCoordinate", sep = " ")
  }
  # Unmapped reads
  if (!is.null(unmapped)){
    args <- paste(args, "--outSAMunmapped", unmapped, sep = " ")
  }
  # SAM attributes
  if (!is.null(sam.attributes)){
    args <- paste(args, "--outSAMattributes", sam.attributes, sep = " ")
  }
  # Gene quantification
  if (!is.null(quant.mode)){
    args <- paste(args, "--quantMode", quant.mode, sep = " ")
  }
  # Compressed files
  if (!is.null(compressed)){
    args <- paste(args, "--readFilesCommand zcat", sep = " ")
  }

  # Results prefix
  results.prefix <- paste(out.dir,sample.name,sample.name, sep = "/")

  # Paired end reads
  if (!is.null(mate2)){
    star.run <- sprintf('%s --genomeDir %s --runThreadN %d --readFilesIn %s %s --outFileNamePrefix %s %s',
                        star,genome,threads,mate1,mate2,results.prefix,args)
  }else{
    star.run <- sprintf('%s --genomeDir %s --runThreadN %d --readFilesIn %s --outFileNamePrefix %s %s',
                        star,genome,threads,mate1,results.prefix,args)
  }

  # Run the Star commands
  if (isTRUE(parallel)){
    cluster <- makeCluster(cores)
    parLapply(cluster, star.run, function (cmd)  system(cmd))
    stopCluster(cluster)
  }else{
    lapply(star.run, function (cmd)  system(cmd))
  }

  return(star.run)
}
