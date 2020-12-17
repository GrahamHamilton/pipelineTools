#' Run Samtools
#'
#' @description Runs the samtools program
#'
#' @import parallel
#'
#' @param command Samtools command to run, at present can choose from 'view', 'sort', 'index' and 'depth', required
#' @param input List of aligned files in sam or bam format, required
#' @param output List of file names for output,
#' @param sample.name List of sample names, required
#' @param threads Number of threads for samtools to use, default set to 10
#' @param memory Set maximum memory per thread, default set to 5Gb
#' @param mapq Set to minimum mapping quality, for filtering bam files
#' @param parallel Run in parallel, default set to FALSE
#' @param cores Number of cores/threads to use for parallel processing, default set to 4
#' @param samtools Path to the samtools program, required
#' @param version Returns the version number
#'
#' @return A list with the samtools commands
#'
#' @examples
#' \dontrun{
#' alignments_path <- "hisat2_alignments"
#' sam.files <- list.files(path = alignments_path, pattern = "sam$",
#'                        full.names = TRUE, recursive = TRUE)
#' bam.files <- gsub(sam.files, pattern = ".sam", replacement = ".bam")
#' sorted.bam.files <- gsub(sam.files, pattern = ".sam", replacement = "_sorted.bam")
#' sample_names <- unlist(lapply(strsplit(list.files(path = reads_path,
#'                        pattern = "*_R1_001.fastq$",
#'                        full.names = FALSE),"_"), `[[`, 1))
#'
#' command <- "view"
#' samtools.cmds <- run_samtools(command = command,
#'                               input = sam.files,
#'                               output = bam.files,
#'                               sample.name = sample_names
#'                               samtools = samtools.path))
#' samtools.cmds
#'
#' command <- "sort"
#' samtools.cmds <- run_samtools(command = command,
#'                               input = bam.files,
#'                               output = sorted.bam.files
#'                               samtools = samtools.path)
#' samtools.cmds
#'
#' command <- "index"
#' samtools.cmds <- run_samtools(command = command,
#'                               input = sorted.bam.files,
#'                               samtools = samtools.path)
#' samtools.cmds
#' command <- "depth"
#' samtools.cmds <- run_samtools(command = command,
#'                               input = sorted.bam.files,
#'                               output = outfile.names,
#'                               samtools = samtools.path)
#' samtools.cmds
#'}
#'
#' @export
#'
run_samtools <- function(command = NULL,
                         input = NULL,
                         output = NULL,
                         sample.name = NULL,
                         threads = 10,
                         mapq = NULL,
                         memory = "5G",
                         parallel = FALSE,
                         cores = 4,
                         samtools = NULL,
                         version = FALSE){
  # Check samtools program can be found
  sprintf("type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'", samtools)

  # Version
  if (isTRUE(version)){
    samtools.run <- sprintf('%s --version',
                            samtools)
    result <- system(samtools.run, intern = TRUE)
    return(result)
  }

  # Set the additional arguments
  args <- ""

  # SAMTOOLS VIEW
  if (command == "view"){
    args <- paste(args,"-bS", sep = " ")
    # Threads
    if (!is.null(threads)){
      args <- paste(args,"--threads",threads,sep = " ")
    }
    if (!is.null(mapq)){
      args <- paste(args,"-q",mapq,sep = " ")
    }

    samtools.run <- sprintf('%s %s %s %s -o %s',
                            samtools,command,args,input,output)

    # Run the samtools commands
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, samtools.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(samtools.run, function (cmd)  system(cmd))
    }

    return(samtools.run)
  }

  # SAMTOOLS SORT
  else if (command == "sort"){
    # Threads
    if (!is.null(threads)){
      args <- paste(args,"--threads",threads,sep = " ")
    }
    # Memory
    if (!is.null(memory)){
      args <- paste(args,"-m",memory, sep = " ")
    }

    samtools.run <- sprintf('%s %s %s %s -o %s',
                            samtools,command,args,input,output)

    # Run the samtools commands
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, samtools.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(samtools.run, function (cmd)  system(cmd))
    }

    return(samtools.run)
  }

  # SAMTOOLS INDEX
  else if (command == "index"){
    # Threads
    if (!is.null(threads)){
      args <- paste(args,"-@",threads, sep = " ")
    }

    samtools.run <- sprintf('%s %s %s %s',
                            samtools,command,args,input)

    # Run the samtools commands
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, samtools.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(samtools.run, function (cmd)  system(cmd))
    }

    return(samtools.run)
  }

  # SAMTOOLS DEPTH
  else if (command == "depth"){

    samtools.run <- sprintf('%s %s %s > %s',
                            samtools,command,input)

    # Run the samtools commands
    if (isTRUE(parallel)){
      cluster <- makeCluster(cores)
      parLapply(cluster, samtools.run, function (cmd)  system(cmd))
      stopCluster(cluster)
    }else{
      lapply(samtools.run, function (cmd)  system(cmd))
    }

    return(samtools.run)
  }
  else{
    result <- sprintf('%s %s %s',"WARNING: Samtools command",command,"not found or misspelled")
    return(result)
  }
}

