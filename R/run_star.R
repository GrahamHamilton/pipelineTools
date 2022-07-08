#' Run the STAR/STARsolo program
#'
#' @description Runs the STAR alignment program, can be run in parallel on
#'   multiple cores. STARsolo also implmented.
#'
#' @import parallel
#'
#' @param mate1 List of the paths to files containing to the forward reads
#' @param mate2 List of the paths to files containing to the reverse reads
#' @param genome.dir Path to the directory where genome files are stored
#' @param sample.name List of the sample names
#' @param out.dir Name of the directory from the Star output
#' @param out.format Format of output file. Can select "BAM SortedByCoordinate",
#'   "BAM Unsorted" or "BAM Unsorted SortedByCoordinate"
#' @param unmapped Fastx will output unmapped and partially mapped (i.e. mapped
#'   only one mate of a paired end read) reads into separate file(s)
#'   Unmapped.out.mate1(2), formatted the same way as input read files.
#' @param sam.attributes Alignment attributes for the SAM/BAM file, default set
#'   to "Standard"
#' @param quant.mode Type of quantification required, recommend set to
#'   "GeneCounts"
#' @param compressed Compression mode for input reads files, recommend set to
#'   "zcat" for gzipped files, can use "bzcat" for bz2 files
#' @param filter.type Filtering to reduce the number of spurious junctions,
#'   default is Normal BySJout for filtering
#' @param filter.multi Set maximum number of multiple alignments for a read, if
#'   exceeded read considered unmapped
#' @param filter.mismatch Maximum number of mismatches per pair. Default 10,
#'   large number switches off this filter e.g. 999
#' @param filter.mismatch.pair Max number of mismatches per pair relative to
#'   read length
#' @param intron.min Minimum intron length, default 21
#' @param intron.max Maximum intron lenght, default 0
#' @param mate.gap Maximum gap between read pairs, default 0
#' @param min.overhang.annotated minimum overhang for annotated junctions,
#'   default 3
#' @param min.overhang.unannotated minimum overhang for unannotated junctions,
#'   default 5
#' @param solo.type Type of single-cell RNASeq, for 10x Chromium or DropSeq use
#'   "CB_UMI_Simple"
#' @param solo.cell.filtering Cell filtering type and parameters
#' @param white.list Path to the file with the whitelist of cell barcodes
#' @param solo.cb.start Cell barcode start base
#' @param solo.cb.len Cell barcode length
#' @param solo.umi.start UMI start base
#' @param solo.umi.len UMI length
#' @param solo.barcode.read.length Length of the barcode read. Set to 1 equal to
#'   sum of soloCBlen+soloUMIlen, set to 0 for do not check
#' @param solo.strand Strandedness of the scRNA libraries
#' @param solo.features Genomic features for which the UMI counts per Cell
#'   Barcode are collected
#' @param solo.multi.mappers Counting method for reads mapping to multiple genes.
#'   Set to Unique, Uniform, Rescue, PorpUnique or EM
#' @param solo.umi.dedup Type of UMI deduplication (collapsing) algorithm.
#'   1MM_All - all UMIs with 1 mismatch distance to each other are collapsed.
#'   1MM_Directional - follows the "directional" method from the UMI-tools by
#'   Smith, Heger and Sudbery (Genome Research 2017). 1MM_NotCollapsed - UMIs
#'   with 1 mismatch distance to others are not collapsed (i.e. all counted)
#' @param solo.umi.filter Type of UMI filtering
#' @param solo.cb.wl.match Matching the Cell Barcodes to the WhiteList
#' @param solo.out.filenames File names for STARsolo output
#' @param threads Number of threads
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default
#'   set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
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
                     genome.dir = NULL,
                     sample.name = NULL,
                     out.dir = NULL,
                     out.format = NULL,
                     unmapped = NULL,
                     sam.attributes = NULL,
                     quant.mode = NULL,
                     compressed = NULL,
                     filter.type = NULL,
                     filter.multi = NULL,
                     filter.mismatch = NULL,
                     filter.mismatch.pair = NULL,
                     intron.min = NULL,
                     intron.max = NULL,
                     mate.gap = NULL,
                     min.overhang.annotated = NULL,
                     min.overhang.unannotated = NULL,
                     solo.type = NULL,
                     solo.cell.filtering = NULL,
                     white.list = NULL,
                     solo.cb.start = NULL,
                     solo.cb.len = NULL,
                     solo.umi.start = NULL,
                     solo.umi.len = NULL,
                     solo.barcode.read.length = NULL,
                     solo.strand = NULL,
                     solo.features = NULL,
                     solo.multi.mappers = NULL,
                     solo.umi.dedup = NULL,
                     solo.umi.filter = NULL,
                     solo.cb.wl.match = NULL,
                     solo.out.filenames = NULL,
                     threads = 10,
                     parallel = FALSE,
                     cores = 4,
                     execute = TRUE,
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
    args <- paste(args, "--outSAMtype", out.format, sep = " ")
  }
  # Unmapped reads
  if (!is.null(unmapped)){
    args <- paste(args, "--outReadsUnmapped", unmapped, sep = " ")
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
    args <- paste(args, "--readFilesCommand", compressed, sep = " ")
  }
  # Filter type
  if (!is.null(filter.type)){
    args <- paste(args, "--outFilterType", filter.type, sep = " ")
  }else{
    args <- paste(args, "--outFilterType", "Normal", sep = " ")
  }
  # Max multi mappers
  if (!is.null(filter.multi)){
    args <- paste(args, "--outFilterMultimapNmax", filter.multi, sep = " ")
  }
  # Max mismatches
  if (!is.null(filter.mismatch)){
    args <- paste(args, "--outFilterMismatchNmax", filter.mismatch, sep = " ")
  }
  # Max mismatches per paired read length
  if (!is.null(filter.mismatch.pair)){
    args <- paste(args, "--outFilterMismatchNoverReadLmax", filter.mismatch.pair, sep = " ")
  }
  # Minimum intron size
  if (!is.null(intron.min)){
    args <- paste(args, "--alignIntronMin", intron.min, sep = " ")
  }
  # Maximum intron size
  if (!is.null(intron.max)){
    args <- paste(args, "--alignIntronMax", intron.max, sep = " ")
  }
  # Maximum gap between read pairs
  if (!is.null(mate.gap)){
    args <- paste(args, "--alignMatesGapMax", mate.gap, sep = " ")
  }
  # Minimum overhang of annotated junctions
  if (!is.null(min.overhang.annotated)){
    args <- paste(args, "--alignSJDBoverhangMin", min.overhang.annotated, sep = " ")
  }
  # Minimum overhang of unannotated junctions
  if (!is.null(min.overhang.unannotated)){
    args <- paste(args, "--alignSJoverhangMin", min.overhang.unannotated, sep = " ")
  }
  #####################################################################
  # Solo
  #####################################################################
  # scRNA type
  if (!is.null(solo.type)){
    args <- paste(args, "--soloType", solo.type, sep = " ")
  }
  # Cell filtering
  if (!is.null(solo.cell.filtering)){
    args <- paste(args, "--soloCellFilter", solo.cell.filtering, sep = " ")
  }
  # Whitelist
  if (!is.null(white.list)){
    args <- paste(args, "--soloCBwhitelist", white.list, sep = " ")
  }
  # Cell barcode start base
  if (!is.null(solo.cb.start)){
    args <- paste(args,"--soloCBstart",solo.cb.start, sep = " ")
  }
  # Cell barcode length
  if (!is.null(solo.cb.len)){
    args <- paste(args,"--soloCBlen",solo.cb.len, sep = " ")
  }
  # UMI start base
  if (!is.null(solo.umi.start)){
    args <- paste(args,"--soloUMIstart",solo.umi.start, sep = " ")
  }
  # UMI length
  if (!is.null(solo.umi.len)){
    args <- paste(args,"--soloUMIlen",solo.umi.len, sep = " ")
  }
  # Solo barcode read length
  if (!is.null(solo.barcode.read.length)){
    args <- paste(args,"--soloBarcodeReadLength",solo.barcode.read.length, sep = " ")
  }
  # Strandedness
  if (!is.null(solo.strand)){
    args <- paste(args,"--soloStrand",solo.strand, sep = " ")
  }
  # Genomic features
  if (!is.null(solo.features)){
    args <- paste(args,"--soloFeatures",solo.features, sep = " ")
  }
  # Multiple gene reads
  if (!is.null(solo.multi.mappers)){
    args <- paste(args,"--soloMultiMappers",solo.multi.mappers, sep = " ")
  }
  # UMI deduplication
  if (!is.null(solo.umi.dedup)){
    args <- paste(args,"--soloUMIdedup",solo.umi.dedup, sep = " ")
  }
  # UMI filtering
  if (!is.null(solo.umi.filter)){
    args <- paste(args,"--soloUMIfiltering",solo.umi.filter, sep = " ")
  }
  # Matching cell barcodes to white list
  if (!is.null(solo.cb.wl.match)){
    args <- paste(args,"--soloCBmatchWLtype",solo.cb.wl.match, sep = " ")
  }

  # Output file names
  if (!is.null(solo.out.filenames)){
    args <- paste(args,"--soloOutFileNames",solo.out.filenames, sep = " ")
  }

  # Results prefix
  results.prefix <- paste(out.dir,sample.name,sample.name, sep = "/")

  # Paired end reads
  if (!is.null(mate2)){
    star.run <- sprintf('%s --genomeDir %s --runThreadN %d --readFilesIn %s %s --outFileNamePrefix %s %s',
                        star,genome.dir,threads,mate1,mate2,results.prefix,args)
  }else{
    star.run <- sprintf('%s --genomeDir %s --runThreadN %d --readFilesIn %s --outFileNamePrefix %s %s',
                        star,genome.dir,threads,mate1,results.prefix,args)
  }

  # Run the Star commands
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, star.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(star.run, function (cmd)  system(cmd))
    }
  }

  return(star.run)
}
