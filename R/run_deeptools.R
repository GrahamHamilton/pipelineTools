#' Run Deeptools
#'
#' @description Runs the deeptools suite of programs
#'
#' @import parallel
#'
#' @param command deeptools command to run, at present can choose from
#'   bamCoverage, computeMatrix and plotHeatmap, required
#' @param input List of files to be processed, required
#' @param output Name of output direcotry
#' @param sample.names list of the sample name
#' @param gtf Path to the gtf file
#' @param before.region Distance upstream of the reference-point selected
#' @param region.bases Distance in bases to which all regions will be fit
#' @param after.region Distance downstream of the reference-point
#'   selected
#' @param plotType Type of plot can select lines, fill, se, std, overlapped_lines or heatmap
#' @param kmeans Number of kmeans clusters to compute.
#' @param threads Number of threads for each instance of deeptools to use,
#'   default set to 10
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default
#'   set to 4
#' @param deeptools Path to the where the deeptools programs are sorted
#'   (usually /usr/local/bin), required
#' @param version Returns the version number
#'
#' @return A list with the deeptools commands
#'
#' @examples
#' \dontrun{
#' path <- "/usr/local/bin/"
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
                          output = NULL,
                          sample.names = NULL,
                          gtf = NULL,
                          before.region = NULL,
                          region.bases = NULL,
                          after.region = NULL,
                          plotType = NULL,
                          kmeans = NULL,
                          threads = 10,
                          parallel = FALSE,
                          cores = 4,
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
  # Bam to BigWig
  if (command == "bamCoverage"){
    deeptools.run <- sprintf('%s%s --numberOfProcessors %s --bam %s --outFileName %s',
                             deeptools,command,threads,input,file.path(output,(paste(sample.names,"bw",sep = ".")),fsep = .Platform$file.sep))
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
  if (isTRUE(parallel)){
    cluster <- makeCluster(cores)
    parLapply(cluster, deeptools.run, function (cmd)  system(cmd))
    stopCluster(cluster)
  }else{
    lapply(deeptools.run, function (cmd)  system(cmd))
  }

  return(deeptools.run)
}
