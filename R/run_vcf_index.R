#' Compress and index vcf files
#'
#'
#' @import parallel
#'
#' @param input List of vcf files to index
#' @param type File type for indexing, e.g. "vcf"
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param execute Whether to execute the commands or not, default set to TRUE
#' @param bgzip Path to the bgzip program, required
#' @param tabix Path to the tabix program, required
#'
#' @examples
#' \dontrun{
#' bgzip.path <- "/software/htslib-1.9/bin/bgzip"
#' tabix.path <- "/software/htslib-1.9/bin/tabix"
#'
#' vcf.files <- c("file1.vcf","file2.vcf")
#'
#' file.type <- "vcf"
#'
#' index.vcf.cmds <- run_vcf_index(input = vcf.files,
#'                                 type = file.type,
#'                                 execute = FALSE,
#'                                 bgzip = bgzip.path,
#'                                 tabix = tabix.path)
#' }
#'
#' @return A list with the indexing commands
#'
#' @export

run_vcf_index <- function(input = NULL,
                          type = NULL,
                          parallel = FALSE,
                          cores = 4,
                          execute = TRUE,
                          bgzip = NULL,
                          tabix = NULL){
  output <- paste(input,".gz",sep = "")
  vcf.index.run <- sprintf('%s -c %s > %s%s%s -p %s %s',
                          bgzip,input,output,"\n",tabix,type,output)

  if (isTRUE(execute)){
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, vcf.index.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(vcf.index.run, function (cmd)  system(cmd))
    }
  }

  return(vcf.index.run)
}
