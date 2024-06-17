#' Run Deeptools
#'
#' @description Runs the deeptools suite of programs
#'
#' @import parallel
#'
#' @param command deeptools command to run, at present can choose from
#'   bamCoverage, computeMatrix and plotHeatmap, required
#' @param input List of files to be processed, required
#' @param control List of input control files, for bamCompare only
#' @param output Name of output directory
#' @param output.format Output format, bigwig or bedgraph, default bigwig
#' @param binsize Size of the bins, in bases, default 50
#' @param sample.names list of the sample name
#' @param gtf Path to the gtf file
#' @param region Region of the genome to limit the operation to,format is chr:start:end
#' @param before.region Distance upstream of the reference-point selected
#' @param region.bases Distance in bases to which all regions will be fit
#' @param after.region Distance downstream of the reference-point selected
#' @param plotType Type of plot can select lines, fill, se, std,
#'   overlapped_lines or heatmap
#' @param kmeans Number of kmeans clusters to compute.
#' @param effective.genome.size The effective genome size is the portion of the
#'   genome that is mappable. Large fractions of the genome are stretches of
#'   NNNN that should be discarded. A table of values is available here:
#'   http://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
#'
#' @param normalization Possible choices: RPKM (Reads Per Kilobase per Million
#'   mapped reads), CPM (Counts Per Million mapped reads), BPM (Bins Per Million
#'   mapped reads, similar to TPM), RPGC (reads per genomic content)

#' @param scale.factors Method to use to scale the samples. Possible choices: readCount, SES and None
#' @param threads Number of threads for each instance of deeptools to use,
#'   default set to 10
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default
#'   set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param deeptools Path to the where the deeptools programs are sorted (usually
#'   /usr/local/bin), required
#' @param version Returns the version number
#'
#' @return A list with the deeptools commands
#'
#' @examples
#' \dontrun{
#'
#' # Version
#' command <- "deeptools"
#'
#' deeptools.cmd <- run_deeptools(command = command,
#'                                deeptools.path = path,
#'                                version = TRUE)
#' deeptools.cmd
#'
#' # Bam to bigwig
#' command <- "bamCoverage"
#' input.names <- "list of input file paths"
#' output.dir <- "directoryname"
#' output.names <- "list of names for samples"
#'
#' deeptools.cmd <- run_deeptools(command = command,
#'                                input = input,
#'                                output = output,
#'                                sample.names = sample.names,
#'                                deeptools = path)
#' deeptools.cmd
#'
#' # Bigwig to matrix
#' command <- "computeMatrix"
#' gtf <- "path/to/gtf"
#' input.names <- "list of input file paths"
#' output.dir <- "directoryname"
#' output.names <- "list of names for samples"
#'
#' deeptools.cmd <- run_deeptools(command = command,
#'                                input = matrix.lists,
#'                                output = output,
#'                                sample.names = output.names,
#'                                gtf = gtf,
#'                                # parallel = TRUE,
#'                                # cores = 6,
#'                                deeptools = path)
#' deeptools.cmd
#'
#' # plot heatmap
#' command <- "plotHeatmap"
#' input.names <- "list of input file paths"
#' output.dir <- "directoryname"
#' output.names <- "list of names for samples"
#'
#' deeptools.cmd <- run_deeptools(command = command,
#'                                input = input.names,
#'                                output = output.dir
#'                                sample.names = output.names,
#'                                # parallel = TRUE,
#'                                # cores = 6,
#'                                deeptools = path)
#' deeptools.cmd
#' }
#' @export
#'
run_deeptools <- function(command = NULL,
                          input = NULL,
                          control = NULL,
                          output = NULL,
                          output.format = NULL,
                          binsize = NULL,
                          sample.names = NULL,
                          gtf = NULL,
                          region = NULL,
                          before.region = NULL,
                          region.bases = NULL,
                          after.region = NULL,
                          plotType = NULL,
                          kmeans = NULL,
                          effective.genome.size = NULL,
                          normalization = NULL,
                          scale.factors = NULL,
                          threads = 10,
                          parallel = FALSE,
                          cores = 4,
                          execute = TRUE,
                          deeptools = NULL,
                          version = FALSE){
  # Version
  if (isTRUE(version)){
    deeptools.run <- sprintf('%s%s --version',
                             deeptools,command)
    result <- system(deeptools.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""
  # Format
  if (!is.null(output.format)){
    args <- paste(args,"--outFileFormat",output.format,sep = " ")
  }
  # Binsize
  if (!is.null(binsize)){
    args <- paste(args,"--binSize",binsize,sep = " ")
  }
  # Region
  if (!is.null(region)){
    args <- paste(args,"--region",region,sep = " ")
  }
  # Before region
  if (!is.null(before.region)){
    args <- paste(args,"--beforeRegionStartLength",before.region,sep = " ")
  }
  # Region bases
  if (!is.null(region.bases)){
    args <- paste(args,"--regionBodyLength",region.bases,sep = " ")
  }
  # After region
  if (!is.null(after.region)){
    args <- paste(args,"--afterRegionStartLength",after.region,sep = " ")
  }
  # Plot type
  if (!is.null(plotType)){
    args <- paste(args,"--plotType",plotType,sep = " ")
  }
  # K means
  if (!is.null(kmeans)){
    args <- paste(args,"--kmeans",kmeans,sep = " ")
  }
  # Effective genome size
  if (!is.null(effective.genome.size)){
    args <- paste(args,"--effectiveGenomeSize",effective.genome.size,sep = " ")
  }
  # Normalisation
  if (!is.null(normalization)){
    args <- paste(args,"--normalizeUsing",normalization,sep = " ")
  }
  # Scale factors
  if (!is.null(scale.factors)){
    args <- paste(args,"--scaleFactorsMethod",scale.factors,sep = " ")
  }

  # Bam to BigWig
  if (command == "bamCoverage"){
    deeptools.run <- sprintf('%s%s --numberOfProcessors %s --bam %s --outFileName %s %s',
                             deeptools,command,threads,input,file.path(output,(paste(sample.names,"bw",sep = ".")),fsep = .Platform$file.sep),args)
  }
  if (command == "bamCompare"){
    deeptools.run <- sprintf('%s%s --numberOfProcessors %s --bamfile1 %s --bamfile2 %s --outFileName %s %s',
                             deeptools,command,threads,input,control,file.path(output,(paste(sample.names,"bw",sep = ".")),fsep = .Platform$file.sep),args)
  }

  # Compute Matrix
  if (command == "computeMatrix"){
      deeptools.run <- sprintf('%s%s scale-regions --numberOfProcessors %s --scoreFileName %s --regionsFileName %s --outFileName %s %s',
                               deeptools,command,threads,input,gtf,file.path(output,(paste(sample.names,"mat","gz",sep = ".")),fsep = .Platform$file.sep),args)
  }

  # plotHeatmap
  if (command == "plotHeatmap"){
    deeptools.run <- sprintf('%s%s --matrixFile %s --outFileName %s',
                             deeptools,command,input,file.path(output,(paste(sample.names,"heatmap","png",sep = ".")),fsep = .Platform$file.sep))
  }

  # profile plots
  if (command == "plotProfile"){
    if (plotType == "heatmap"){
      deeptools.run <- sprintf('%s%s --matrixFile %s --outFileName %s --perGroup %s',
                               deeptools,command,input,file.path(output,(paste(sample.names,"profile","heatmap","png",sep = ".")),fsep = .Platform$file.sep),args)
    }else{
      deeptools.run <- sprintf('%s%s --matrixFile %s --outFileName %s --perGroup %s',
                               deeptools,command,input,file.path(output,(paste(sample.names,"profile","png",sep = ".")),fsep = .Platform$file.sep),args)
    }
  }

  # Run the deeptools commands
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, deeptools.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(deeptools.run, function (cmd)  system(cmd))
    }
  }

  return(deeptools.run)
}
