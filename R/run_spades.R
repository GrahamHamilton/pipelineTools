#' Run SPAdes tool
#'
#' @description Runs the SPAdes assembler
#'
#' @import parallel
#'
#' @param input1 List of the paths to files containing to the forward reads
#' @param input2 List of the paths to files containing to the reverse reads
#' @param out.dir Name of the directory from the SPAdes output
#' @param sample.name List of the sample names
#' @param isolate This flag is highly recommended for high-coverage isolate and multi-cell data
#' @param sc This flag is required for MDA (single-cell) data
#' @param meta This flag is required for metagenomic data
#' @param bio This flag is required for biosyntheticSPAdes mode
#' @param sewage This flag is required for sewage mode
#' @param corona This flag is required for coronaSPAdes mode
#' @param rna This flag is required for RNA-Seq data
#' @param plasmid Runs plasmidSPAdes pipeline for plasmid detection
#' @param metaviral Runs metaviralSPAdes pipeline for virus detection
#' @param metaplasmid Runs metaplasmidSPAdes pipeline for plasmid detection in metagenomic datasets
#'                    (equivalent for meta & plasmid)
#' @param rnaviral This flag enables virus assembly module from RNA-Seq data
#' @param iontorrent This flag is required for IonTorrent data
#' @param sanger Sanger sequence files
#' @param pacbio PacBIo sequence files
#' @param nanopore Nanopore sequence files
#' @param contigs Trusted cotig files
#' @param error_correction_only Run read error correction only
#' @param assembly_only Run assembly only
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param spades Path to the SPAdes program
#' @param version Returns the version number
#'
#' @return A list of SPAdes commands
#'
#' @examples
#'  \dontrun{
#'  path <- "/software/SPAdes-4.0.0-Linux/bin/spades.py"
#'
#'  run_spades(spades = path,
#'             version = TRUE)
#'  fastq1 <- c("sample1_R1.fq", "sample2_R1.fq")
#'  fastq2 <- c("sample1_R2.fq", "sample2_R2.fq")
#'  sample_names <- c("sample1", "sample2")
#'  run_spades(input1 = fastq1,
#'             input2 = fastq2,
#'             out.dir = "Out",
#'             isolate = TRUE,
#'             sample.name = sample_names,
#'             execute = FALSE,
#'             spades = path)
#' }
#' @export
#'

run_spades <- function(input1 = NULL,
                       input2 = NULL,
                       out.dir = NULL,
                       sample.name = NULL,
                       isolate	= FALSE,
                       sc	= FALSE,
                       meta	= FALSE,
                       bio	= FALSE,
                       sewage	= FALSE,
                       corona	= FALSE,
                       rna	= FALSE,
                       plasmid	= FALSE,
                       metaviral	= FALSE,
                       metaplasmid	= FALSE,
                       rnaviral	= FALSE,
                       sanger = NULL,
                       pacbio = NULL,
                       nanopore = NULL,
                       contigs = NULL,
                       error_correction_only = NULL,
                       assembly_only = NULL,
                       parallel = FALSE,
                       cores = 4,
                       execute = TRUE,
                       spades = NULL,
                       version = FALSE){
  # Checkspades program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", spades)

  # Version
  if (isTRUE(version)){
    spades.run <- sprintf('%s --version',
                           spades)
    result <- system(spades.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""
  # Isolate
  if (isTRUE(isolate)){
    args <- paste(args,"--isolate",sep = " ")
  }
  # Single cell
  if (isTRUE(sc)){
    args <- paste(args,"--sc",sep = " ")
  }
  # Metagenomics
  if (isTRUE(meta)){
    args <- paste(args,"--meta",sep = " ")
  }
  # Biosynthetic
  if (isTRUE(bio)){
    args <- paste(args,"--bio",sep = " ")
  }
  # Sewage
  if (isTRUE(sewage)){
    args <- paste(args,"--sewage",sep = " ")
  }
  # Corona
  if (isTRUE(corona)){
    args <- paste(args,"--corona",sep = " ")
  }
  # RNA
  if (isTRUE(rna)){
    args <- paste(args,"--rna",sep = " ")
  }
  # Plasmid
  if (isTRUE(plasmid)){
    args <- paste(args,"--plasmid",sep = " ")
  }
  # Metaviral
  if (isTRUE(metaviral)){
    args <- paste(args,"--metaviral",sep = " ")
  }
  # IMetaplasmid
  if (isTRUE(metaplasmid)){
    args <- paste(args,"--metaplasmid",sep = " ")
  }
  # RNA viral
  if (isTRUE(rnaviral)){
    args <- paste(args,"--rnaviral",sep = " ")
  }
  # Sanger
  if (isTRUE(sanger)){
    args <- paste(args,"--sanger", sanger,sep = " ")
  }
  # Pacbio
  if (isTRUE(pacbio)){
    args <- paste(args,"-pacbio", pacbio,sep = " ")
  }
  # Nanopore
  if (isTRUE(nanopore)){
    args <- paste(args,"-nanopore", nanopore,sep = " ")
  }

  # Create the sample directories for the per sample spades results
  lapply(paste(out.dir,sample.name, sep = "/"), function(cmd) dir.create(cmd, showWarnings = FALSE, recursive = TRUE))

  out_dirs <- paste(out.dir,sample.name, sep = "/")

  spades.run <- sprintf('%s %s -1 %s -2 %s -o %s',
                        spades,args,input1,input2,out_dirs)

  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, spades.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(spades.run, function (cmd)  system(cmd))
    }
  }

  return(spades.run)
}
