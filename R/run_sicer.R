#' Run Sicer2
#'
#' @description Sicer2 is a spatial clustering approach for the identification
#'   of ChIP-enriched regions which was developed for calling broad peaks of
#'   histone modifications from ChIP-seq data.
#'
#' @param treatment List of the paths to files containing the treatment files.
#'   This can either be the relative or the absolute path of the file. Must be
#'   in BED or BAM format.
#' @param control List of the paths to files containing the control files. This
#'   can either be the relative or the absolute path of the file. Must be in BED
#'   or BAM format. OPTIONAL.
#' @param comparison.names List of the comparisons to be made. OPTIONAL.
#' @param species The species/genome used (ex: hg38).
#' @param redundancy_threshold The number of copies of indentical reads allowed
#'   in a library. Default is 1.
#' @param window_size Resolution of SICER. Default value is 200 (bp)
#' @param fragment_size The amount of shift from the beginning of a read to the
#'   center of the DNA fragment represented by the read. Default is 150 (bp).
#' @param effective_genome_fraction Effective genome as fraction of the genome
#'   size. Default is 0.74.
#' @param false_discovery_rate Remove all islands with an false_discovery_rate
#'   below cutoff. Default is 0.01.
#' @param gap_size The minimum length of a "gap" such that neighboring window is
#'   an "island." This value must be a multiple of the window size. Default is
#'   600 (bp).
#' @param e_value Requires user input when no control library is provided.
#'   Default is 1000.
#' @param step_score Step Score: The minimum number of positive elements in the
#'   graining unit to call the unit positive. Used for RECOGNICER algorithm.
#' @param cpu The number of CPU cores SICER program will use when executing
#'   multi-processing tasks. Optimal core count is the species' number of
#'   chromosomes. Default value is the maximum number of cores avaiable in the
#'   system.
#' @param significant_reads SICER produces a BED file of treatment reads
#'   filtered by significant islands and WIG file of filtered reads binned into
#'   windows
#' @param sicer Path to the Sicer2 program
#'
#' @return A list with the Sicer2 commands
#'
#' @examples
#' \dontrun{
#'  sicer.cmds <- run_sicer(treatment = treatment.list,
#'                          control = control.list,
#'                          comparison.names = comp.list,
#'                          significant_reads = TRUE,
#'                          sicer = sicer.path
#'                         )
#' }
#' @export
#'

run_sicer <- function(treatment = NULL,
                      control = NULL,
                      comparison.names = NULL,
                      species = NULL,
                      redundancy_threshold = 1,
                      window_size = 200,
                      fragment_size = 150,
                      effective_genome_fraction = 0.74,
                      false_discovery_rate = 0.01,
                      gap_size = 600,
                      e_value = 1000,
                      step_score = NULL,
                      cpu = NULL,
                      significant_reads = FALSE,
                      sicer = NULL){
  # Check sicer program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", sicer)

  # Create the directories to write the data, if it is not present
  if(!is.null(comparison.names)){
    out.dir <- "sicer_data"
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
    # Create the sample directories for the per sample Sicer2 results
    lapply(paste(out.dir,comparison.names, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))
    output.directories <- paste(out.dir,comparison.names, sep = "/")
  }

  # Set the additional arguments
  args <- ""
  # Species
  if(!is.null(species)){
    args <- paste(args,"--species",species,sep = " ")
  }
  # Redundancy threshold
  if(!is.null(redundancy_threshold)){
    args <- paste(args,"--redundancy_threshold",redundancy_threshold,sep = " ")
  }
  # Window size
  if(!is.null(window_size)){
    args <- paste(args,"--window_size",window_size,sep = " ")
  }
  # Fragment size
  if(!is.null(fragment_size)){
    args <- paste(args,"--fragment_size",fragment_size,sep = " ")
  }
  # Effective genome fraction
  if(!is.null(effective_genome_fraction)){
    args <- paste(args,"--effective_genome_fraction",effective_genome_fraction,sep = " ")
  }
  # False discovery rate
  if(!is.null(false_discovery_rate)){
    args <- paste(args,"--false_discovery_rate",false_discovery_rate,sep = " ")
  }
  # Gap size
  if(!is.null(gap_size)){
    args <- paste(args,"--gap_size",gap_size,sep = " ")
  }
  # E value
  if(!is.null(e_value)){
    args <- paste(args,"--e_value",e_value,sep = " ")
  }
  # Step score
  if(!is.null(step_score)){
    args <- paste(args,"--step_score",step_score,sep = " ")
  }
  # CPU
  if(!is.null(cpu)){
    args <- paste(args,"--cpu",cpu,sep = " ")
  }
  # Significant reads
  if(isTRUE(significant_reads)){
    args <- paste(args,"--significant_reads",sep = " ")
  }
  #Output directories
  if(!is.null(output.directories)){
    args <- paste(args,"--output_directory",output.directories,sep = " ")
  }

  if(!is.null(control)){
    sicer.run <- sprintf('%s --treatment_file %s --control_file %s %s',
                            sicer,treatment,control,args,output.directories)
  }else{
    sicer.run <- sprintf('%s --treatment_file %s  %s',
                         sicer,treatment,args)
  }

  lapply(sicer.run, function (cmd)  system(cmd))

  return(sicer.run)
}

# comparisons <- read.csv("comparisons.tsv", sep = "\t", header = TRUE)
# class(comparisons)
# for (index in 1:nrow(comparisons)){
#   test.cond <- as.character(comparisons[index,]$test.condition)
#   print(test.cond)
#
#   base.cond <- as.character(comparisons[index,]$base.condition)
#   print(base.cond)
# }
#
# treatment1 <- paste("treatment1","treatment2","treatment3",sep = " ")
# treatment2 <- paste("treatment4","treatment5","treatment6",sep = " ")
# treatment.list <- list(treatment1,treatment2)
#
# control1 <- paste("control1","control2","control3",sep = " ")
# control2 <- paste("control4","control5","control6",sep = " ")
# control.list <- list(control1,control2)
#
# sicer.path <- "/usr/local/bin/recognicer"
#
# comp.list <- list()
# for (index in 1:length(treatment.list)){
#   comp.list[index] <- paste(treatment.list[index],control.list[index], sep = "_vs_")
# }
#
#
# sicer.cmds <- run_sicer(treatment = treatment.list,
#                         #control = control.list,
#                         comparison.names = comp.list,
#                         significant_reads = TRUE,
#                         sicer = sicer.path
#                         )
# print(sicer.cmds)
