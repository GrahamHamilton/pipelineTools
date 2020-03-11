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
#' @param before.region Distance upstream of the reference-point selected,
#'   default 500
#' @param region.bases Distance in bases to which all regions will be fit,
#'   default 1000
#' @param after.region Distance downstream of the reference-point
#'   selected,default 1500
#' @param threads Number of threads for each instance of deeptools to use,
#'   default set to 10
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default
#'   set to 4
#' @param deeptools.path Path to the where the deeptools programs are sorted
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
#'                                deeptools.path = path)
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
#'                                deeptools.path = path)
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
#'                                deeptools.path = path)
#' deeptools.cmd
#' }
#' @export
#'
run_deeptools <- function(command = NULL,
                          input = NULL,
                          output = NULL,
                          sample.names = NULL,
                          gtf = NULL,
                          before.region = 500,
                          region.bases = 5000,
                          after.region = 1500,
                          threads = 10,
                          parallel = FALSE,
                          cores = 4,
                          deeptools.path = NULL,
                          version = FALSE){
  # Version
  if (isTRUE(version)){
    deeptools.run <- sprintf('%s%s --version',
                             deeptools.path,command)
    result <- system(deeptools.run, intern = TRUE)
    return(result)
  }

  # Bam to BigWig
  if (command == "bamCoverage"){
    deeptools.run <- sprintf('%s%s --numberOfProcessors %s --bam %s --outFileName %s',
                             deeptools.path,command,threads,input,file.path(output,(paste(sample.names,"bw",sep = ".")),fsep = .Platform$file.sep))
  }

  # Compute Matrix
  if (command == "computeMatrix"){
      deeptools.run <- sprintf('%s%s scale-regions --numberOfProcessors %s --scoreFileName %s --regionsFileName %s --outFileName %s --beforeRegionStartLength %s --regionBodyLength %s --afterRegionStartLength %s',
                               deeptools.path,command,threads,input,gtf,file.path(output,(paste(sample.names,"mat","gz",sep = ".")),fsep = .Platform$file.sep),before.region,region.bases,after.region)
  }

  # plotHeatmap
  if (command == "plotHeatmap"){
    deeptools.run <- sprintf('%s%s --matrixFile %s --outFileName %s',
                             deeptools.path,command,input,file.path(output,(paste(sample.names,"heatmap","png",sep = ".")),fsep = .Platform$file.sep))
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
