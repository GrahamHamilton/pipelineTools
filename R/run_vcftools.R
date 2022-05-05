#' Run vcftools
#'
#' @description Runs the vcftool, currently only set for creating PLINK Map and Ped files
#'
#' @param input Path to the input vcf file
#' @param output Path to the output files
#' @param plink Create plink files, boolean, default st to FALSE
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param vcftools Path to the vcftools program, required
#' @param version Returns the version number
#'
#' @return Alist of vcftools commnads
#' @export
#'
#' @examples
#' \dontrun{
#' vcftools.path <- "/software/vcftools-v0.1.16/bin/vcftools"
#'
#' run_vcftools(vcftools = vcftools.path,
#'              version = TRUE)
#'
#' input.file <- "variants.vcf"
#' output.file <- gsub("vcf","plink",input.file)
#'
#' run_vcftools(input = input.file,
#'              output = output.file,
#'              plink = TRUE,
#'              execute = FALSE,
#'              vcftools = vcftools.path)
#' }
#'
#' @export
#'
run_vcftools <- function(input = NULL,
                         output = NULL,
                         plink = FALSE,
                         execute = TRUE,
                         vcftools = NULL,
                         version = NULL){
  # Version
  if (isTRUE(version)){
    vcftools.run <- sprintf('%s --version',
                            vcftools)
    result <- system(vcftools.run, intern = TRUE)
    return(result)
  }

  if(isTRUE(plink)){
    vcftools.run <- sprintf('%s --vcf %s --plink --out %s',
                            vcftools,input,output)
  }

  if (isTRUE(execute)){
    lapply(vcftools.run, function (cmd)  system(cmd))
  }

  return(vcftools.run)
}
