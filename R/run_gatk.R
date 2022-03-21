#' Run GATK
#'
#' @description Runs the GATK suite of programs.
#'
#' @param command GATK command to run, required
#' @param input List of sorted bam files, required
#' @param output List of output files or empty/non-existant directory
#' @param reference Path to the fasta formatted reference
#' @param intervals Path to the intrvals list file, usually the exome coordiantes file
#' @param known.sites List of paths to the files containing known polymorphic sites
#' @param sample.map Path to tab seperated file mapping sample names to the gvcf file
#' @param erc Mode for emitting reference confidence scores, can be "NONE", "BP_RESOLUTION" and "GVCF"
#' @param temp Path to temporary directory
#' @param batch Batch size fior number of readers open at once,GenomicsDBImport only
#' @param threads Number of threads for opening VCFs in batches
#' @param database Name of the database directory
#' @param bqsr List of base quality recalibration files
#' @param vqsr Varian recalibration file
#' @param tranches List of levels of truth sensitivity at which to slice the data
#' @param resources Pre defined list of sites for which to apply a prior probability of being correct
#' @param annotations List ofnames of the annotations to be used for calculations
#' @param mode Recalibration mode to employ, SNP or INDEL
#' @param tranches.file The input tranches file describing where to cut the data, from VariantRecalibrator
#' @param sensitivity.filter The truth sensitivity level at which to start filtering
#' @param variant.index Create a VCF index when writing a coordinate-sorted VCF file, boolean, default set to TRUE
#' @param variant.type Variant type to include in output, SNP or INDEL.
#' @param select.method Method to select filtered vaiants, to select only variant that pass all filters use 'vc.isNotFiltered()'
#' @param filter.expression String of filters and values
#' @param filter.name Name to identify the filtered variants
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param gatk Path to the GATK suit of programs, required
#'
#' @examples
#'
#' \dontrun{
#' known.site.list <- c("Mills_and_1000G_gold_standard.indels.GRCh38.vcf.gz",
#'                       "1000G_phase1.snps.high_confidence.GRCh38.vcf.gz")
#'
#' recalibration.files <- gsub(".bam","_recal_data.table",rg.bam.files)
#' command <-  "BaseRecalibrator"
#' fasta <- Path the genome fats
#' # BaseRecalibration
#' base.recalibration.cmds <- run_gatk(command = command,
#'                                     input = rg.bam.files,
#'                                     output = recalibration.files,
#'                                     reference = fasta,
#'                                     intervals = intervals.file,
#'                                     known.sites = known.site.list,
#'                                     execute = TRUE,
#'                                     gatk = gatk)
#' }
#'
#' @return List of GATK commands
#'
#' @export
#'

run_gatk <- function(command = NULL,
                     input = NULL,
                     output = NULL,
                     reference = NULL,
                     intervals = NULL,
                     known.sites = NULL,
                     sample.map = NULL,
                     erc = NULL,
                     temp = NULL,
                     batch = NULL,
                     threads = NULL,
                     database = NULL,
                     bqsr = NULL,
                     vqsr = NULL,
                     tranches = NULL,
                     resources = NULL,
                     annotations = NULL,
                     mode = NULL,
                     tranches.file = NULL,
                     sensitivity.filter = NULL,
                     variant.index = TRUE,
                     variant.type = NULL,
                     select.method = NULL,
                     filter.expression = NULL,
                     filter.name = NULL,
                     parallel = FALSE,
                     cores = 4,
                     execute = TRUE,
                     gatk = NULL){
  # Check gatk program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", gatk)

  # Set the additional arguments
  args <- ""

  # Known sites
  if (!is.null(known.sites)){
    for (site in known.sites){
      args <- paste(args,"--known-sites",site, sep = " ")
    }
  }
  # Temporary directory
  if (!is.null(temp)){
    args <- paste(args,paste("--tmp-dir",temp, sep = "="), sep = " ")
  }
  # Batch size
  if (!is.null(batch)){
    args <- paste(args,"--batch-size",batch, sep = " ")
  }
  # Reader threads
  if (!is.null(threads)){
    args <- paste(args,"--reader-threads",threads, sep = " ")
  }
  # GenomicsDB
  if (!is.null(database)){
    args <- paste(args,"-V",paste("gendb",database, sep = "://"), sep = " ")
  }
  # Tranches
  if (!is.null(tranches)){
    list <- sapply(tranches, function(l) paste("-tranche",l,sep = " "))
    args <- paste(args,paste(list, collapse = " "), sep = " ")
  }
  # Annotations
  if (!is.null(annotations)){
    list <- sapply(annotations, function(l) paste("-an",l,sep = " "))
    args <- paste(args,paste(list, collapse = " "), sep = " ")
  }
  # EmitReferenceConfidence
  if (!is.null(erc)){
    args <- paste(args,"-ERC", erc, sep = " ")
  }
  # Select type to include
  if (!is.null(variant.type)){
    args <- paste(args, "--select-type-to-include", variant.type, sep = " ")
  }
  # Selection method
  if (!is.null(select.method)){
    args <- paste(args, "-select", select.method, sep = " ")
  }
  # Filter expression
  if (!is.null(filter.expression)){
    args <- paste(args, "--filter-expression",filter.expression, sep = " ")
  }
  # Filter name
  if (!is.null(filter.name)){
    args <- paste(args, "--filter-name", filter.name, sep = " ")
  }

  ## Set up commands
  # BaseRecalibrator
  if (command == "BaseRecalibrator"){
    if (is.null(intervals)){
      gatk.run <- sprintf('%s %s -I %s -O %s -R %s %s',
                          gatk,command,input,output,reference,args)
    }else{
      gatk.run <- sprintf('%s %s -I %s -O %s -R %s -L %s  %s',
                          gatk,command,input,output,reference,intervals,args)
      }
  }

  # ApplyBQSR
  if (command == "ApplyBQSR"){
    if (is.null(intervals)){
      gatk.run <- sprintf('%s %s -I %s -O %s -R %s --bqsr-recal-file %s %s',
                          gatk,command,input,output,reference,bqsr,args)
    }else{
      gatk.run <- sprintf('%s %s -I %s -O %s -R %s --bqsr-recal-file %s -L %s  %s',
                          gatk,command,input,output,reference,bqsr,intervals,args)
    }
  }

  # HaplotypeCaller
  if (command == "HaplotypeCaller"){
    if (is.null(intervals)){
      gatk.run <- sprintf('%s %s -I %s -O %s -R %s %s',
                          gatk,command,input,output,reference,args)
    }else{
      gatk.run <- sprintf('%s %s -I %s -O %s -R %s -L %s %s',
                          gatk,command,input,output,reference,intervals,args)
    }
  }

  # GenomicsDBImport
  if (command == "GenomicsDBImport"){
    if (is.null(intervals)){
      gatk.run <- sprintf('%s %s --genomicsdb-workspace-path %s --sample-name-map $s %s',
                          gatk,command,output,sample.map,args)
    }else{
      gatk.run <- sprintf('%s %s --genomicsdb-workspace-path %s --sample-name-map %s -L %s %s',
                          gatk,command,output,sample.map,intervals,args)
    }
  }

  # GenotypeGVCFs
  if (command == "GenotypeGVCFs"){
    gatk.run <- sprintf('%s %s -R %s -O %s %s',
                          gatk,command,reference,output,args)
  }

  # VariantRecalibrator
  if (command == "VariantRecalibrator"){
    tranches.file <- gsub(".vcf",paste("",mode,"tranches",sep = "."),input)
    rscript.file <- gsub(".vcf",paste("",mode,"plots.R",sep = "."),input)
    gatk.run <- sprintf('%s %s -R %s -V %s %s -mode %s -O %s --tranches-file %s  --rscript-file %s %s',
                        gatk,command,reference,input,resources,mode,output,tranches.file,rscript.file,args)
  }

  # ApplyVSQR
  if (command == "ApplyVQSR"){
    gatk.run <- sprintf('%s %s -V %s --recal-file %s -mode %s --tranches-file %s --truth-sensitivity-filter-level %s --create-output-variant-index %s -O %s',
                       gatk,command,input,vqsr,mode,tranches.file,sensitivity.filter,variant.index,output)
  }

  # SelectVariants
  if (command == "SelectVariants"){
    gatk.run <- sprintf('%s %s -R %s  -V %s -O %s %s',
                        gatk,command,reference,input,output,args)
  }

  # VariantFiltration
  if (command == "VariantFiltration"){
    gatk.run <- sprintf('%s %s  -V %s -O %s %s',
                        gatk,command,input,output,args)
  }

  # Run the commands, if execute is true
  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, gatk.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(gatk.run, function (cmd)  system(cmd))
    }
  }

  # Return the list of Picard commands
  return(gatk.run)
}
