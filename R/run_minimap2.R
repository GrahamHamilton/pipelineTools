#' Run Minimap2
#'
#' @description Runs the Minimap2 tool, can be used for single end and paired end reads from Illumina,
#' Oxford Nanopore and PacBio platforms
#'
#' @import parallel
#'
#' @param input1 input1 List of the paths to files containing to the forward reads, required
#' @param input2 input2 List of the paths to files containing to the reverse reads, required for paired end sequence data
#' @param index Path to the reference genome minimap index
#' @param genome Path to the reference genome fasta
#' @param sample.name List of the sample names, required
#' @param out.dir Name of the directory from the Minimap2 alignments.
#' @param platform Names of the sequencing platform used for the reads, choose either "ONT", "PacBio" or "Illumina"
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param minimap2 Path to the Minimap2 program, required
#' @param version Returns the version number
#'
#' @examples
#' \dontrun{
#'   # Set the variables
#'   minimap2 <- "/software/minimap2-v2.21/minimap2"
#'   out.dir <- "minmap2_alignments"
#'   input1 <- c("sample1_R1.fq", "sample2_R1.fq")
#'   input2 <- c("sample1_R2.fq", "sample2_R2.fq")
#'   sample.names <- c("sample1", "sample2")
#'   genome <- "reference.fa"
#'   index <- "reference.mmi"
#'
#'   # Get the version number
#'   run_minimap2(minimap2 = minimap2,
#'                version = TRUE)
#'
#'   # Alignment commands
#'   minimap_cmds <- run_minimap2(input1 = input,
#'                                sample.name = sample.names,
#'                                index = index,
#'                                out.dir = out.dir,
#'                                platform = "ONT",
#'                                parallel = TRUE,
#'                                cores = 2,
#'                                execute = FALSE,
#'                                minimap2 = minimap2)
#'
#'   minimap_cmds
#' }
#'
#' @return A list with the Minimap2 commands
#' @export
#'
run_minimap2 <- function(input1 = NULL,
                         input2 = NULL,
                         index = NULL,
                         genome = NULL,
                         sample.name = NULL,
                         out.dir = NULL,
                         platform = NULL,
                         parallel = FALSE,
                         cores = 4,
                         execute = TRUE,
                         minimap2 = NULL,
                         version = FALSE){
  # Check minimap2 program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", minimap2)

  # Version
  if (isTRUE(version)){
    minimap2.run <- sprintf('%s --version',
                            minimap2)
    result <- system(minimap2.run, intern = TRUE)
    return(result)
  }

  # Create the sample directories for the per sample bwa results
  lapply(paste(out.dir,sample.name, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))

  # Set the additional arguments
  args <- ""
  # Platform
  if (!is.null(platform)){
    if (platform == "ONT"){
      args <- paste(args,"--ax map-ont",sep = " ")
    }else if (platform == "PacBio"){
      args <- paste(args,"--ax map-opb",sep = " ")
    }else if (platform =="Illumina"){
      args <- paste(args,"--ax map-iclr",sep = " ")
    }else{
      stop("Please provide either ONT, PacBio or Illumina for platorm variable")
    }
  }

  # Set the names for the alignment
  sam.files <- paste(out.dir,sample.name,paste(sample.name,"sam",sep = "."),sep = "/")

  # Create the sample directories for the per sample bowtie2 results
  lapply(paste(out.dir,sample.name, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))

  # Create the Minimap2 commands
  if (is.null(input2)){
    # Single end fasta reference genome
    if (!is.null(genome)){
      minimap2.run <- sprintf('%s %s %s %s > %s',
                            minimap2,args,genome,input1,sam.files)
    }
    # Single end indexed reference geneome
    else if (!is.null(index)){
      minimap2.run <- sprintf('%s %s %s %s > %s',
                            minimap2,args,index,input1,sam.files)
    }
  }else if (!is.null(input2)){
    # Paired end fasta reference genome
    if (!is.null(genome)){
      minimap2.run <- sprintf('%s %s %s %s %s > %s',
                              minimap2,args,genome,input1,input2,sam.files)
    }
    # Paired end indexed reference geneome
    else if (!is.null(index)){
      minimap2.run <- sprintf('%s %s %s %s %s > %s',
                              minimap2,args,index,input1,input2,sam.files)
    }
  }

  # Run the minimap2 commands
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, minimap2.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(minimap2.run, function (cmd)  system(cmd))
    }
  }

  # Return the list of Minimap2 commands
  return(minimap2.run)

}
