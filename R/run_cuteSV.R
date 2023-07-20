#' Run cuteSV
#'
#' @description cuteSV, a long-read-based structural variant detection program
#'
#' @import parallel
#'
#' @param input List of sorted bam files, required
#' @param reference Path to the fasta formatted reference
#' @param sample.names List of the sample names
#' @param bed Only detect SVs in regions in the BED file
#' @param out.dir Name of the directory from the cuteSV output
#' @param platform Names of the sequencing platform, choose either "ONT", "PacBio_CLR",
#' "PacBio_CCS" or "FORCE"
#' @param max_cluster_bias_INS Maximum distance to cluster reads together for insertion.
#' @param diff_ratio_merging_INS Do not merge breakpoints with basepair identity more
#' than [0.3] for insertion.
#' @param max_cluster_bias_DEL Maximum distance to cluster read together for
#' deletion.
#' @param diff_ratio_merging_DEL Do not merge breakpoints with basepair identity more
#' than [0.5] for deletion.
#' @param min_mapq Minimum mapping quality value of alignment to be taken
#' into account.
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param cuteSV Path to the cuteSV program, required
#' @param version Returns the version number
#'
#' @examples
#' \dontrun{
#' cuteSV.path <- "path/to/cuteSV"
#'
#' run_cuteSV(cuteSV = cuteSV.path,
#'            version = TRUE)
#'
#' out.dir <- "VCF"
#' sample.names <- c("sample_1","sample_2","sample_3")
#' sorted.bam.files <- c("sample_1.bam","sample_2.bam","sample_3.bam")
#' genome = "/datastore/Gneomes/Reference/genome.fa"
#'
#' run_cuteSV(input = sorted.bam.files,
#'            reference = genome,
#'            sample.names = sample.names,
#'            out.dir =out.dir,
#'            platform = "ONT",
#'            cuteSV = cuteSV.path)
#' }
#'
#' @return A list with the Minimap2 commands
#' @export
#'
run_cuteSV <- function(input = NULL,
                       reference = NULL,
                       sample.names = NULL,
                       bed = NULL,
                       out.dir = NULL,
                       platform = NULL,
                       max_cluster_bias_INS = NULL,
                       diff_ratio_merging_INS = NULL,
                       max_cluster_bias_DEL = NULL,
                       diff_ratio_merging_DEL = NULL,
                       min_mapq = NULL,
                       parallel = FALSE,
                       cores = 4,
                       execute = FALSE,
                       cuteSV = NULL,
                       version = FALSE){
  # Check cuteSV program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", cuteSV)

  # Version
  if (isTRUE(version)){
    cuteSV.run <- sprintf('%s --version',
                              cuteSV)
    result <- system(cuteSV.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""
  # Set cuteSV suggested parameters
  if (!is.null(platform)){
    if (platform == "ONT"){
      args <- paste(args,"--max_cluster_bias_INS 100",sep = " ")
      args <- paste(args,"--diff_ratio_merging_INS 0.3",sep = " ")
      args <- paste(args,"--max_cluster_bias_DEL 100",sep = " ")
      args <- paste(args,"--diff_ratio_merging_DEL 0.3",sep = " ")
    }else if (platform == "PacBio_CLR"){
      args <- paste(args,"--max_cluster_bias_INS 100",sep = " ")
      args <- paste(args,"--diff_ratio_merging_INS 0.3",sep = " ")
      args <- paste(args,"--max_cluster_bias_DEL 200",sep = " ")
      args <- paste(args,"--diff_ratio_merging_DEL 0.5",sep = " ")
    }else if (platform == "PacBio_CCS"){
      args <- paste(args,"--max_cluster_bias_INS 1000",sep = " ")
      args <- paste(args,"--diff_ratio_merging_INS 0.3",sep = " ")
      args <- paste(args,"--max_cluster_bias_DEL 100",sep = " ")
      args <- paste(args,"--diff_ratio_merging_DEL 0.3",sep = " ")
    }else if (platform == "FORCE"){
      args <- paste(args,"--min_mapq 10",sep = " ")
    }else{
      stop("Please provide either ONT, PacBio_CLR, PacBio_CCS of FORCE for platorm variable")
    }
  }
  # Maximum disatnce for insertion clustering
  if (!is.null(max_cluster_bias_INS)){
    args <- paste(args,"--max_cluster_bias_INS",max_cluster_bias_INS, sep = " ")
  }
  # Insertion breakpoint merging
  if (!is.null(diff_ratio_merging_INS)){
    args <- paste(args,"--diff_ratio_merging_INS",diff_ratio_merging_INS, sep = " ")
  }
  # Maximum disatnce for deletion clustering
  if (!is.null(max_cluster_bias_DEL)){
    args <- paste(args,"--max_cluster_bias_DEL",max_cluster_bias_DEL, sep = " ")
  }
  # Deletion breakpoint merging
  if (!is.null(diff_ratio_merging_DEL)){
    args <- paste(args,"--diff_ratio_merging_DEL",diff_ratio_merging_DEL, sep = " ")
  }
  # Minimum mapping quality
  if (!is.null(min_mapq)){
    args <- paste(args,"--min_mapq",min_mapq, sep = " ")
  }
  # bed file
  if (!is.null(bed)){
    args <- paste(args,"-include_bed",bed, sep = " ")
  }

  # Set the vcf output file names
  vcf.files <- paste(sample.names,"vcf",sep = ".")

  # Create the run commands
  cuteSV.run <- sprintf('%s %s %s  %s %s %s',
                            cuteSV,args,input,reference,vcf.files,out.dir)

  # Run the commands, if execute is true
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, cuteSV.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(cuteSV.run, function (cmd)  system(cmd))
    }
  }

  # Return the list of cuteSV commands
  return(cuteSV.run)
}
