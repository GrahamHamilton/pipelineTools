% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_htseq_count.R
\name{run_htseq_count}
\alias{run_htseq_count}
\title{Run HtSeq-Count}
\usage{
run_htseq_count(
  input = NULL,
  output = NULL,
  mode = NULL,
  type = NULL,
  attribute = NULL,
  format = NULL,
  gtf = NULL,
  parallel = FALSE,
  cores = 4,
  execute = TRUE,
  htseq_count = NULL,
  version = FALSE
)
}
\arguments{
\item{input}{List of aligned files in sam or bam format, required}

\item{output}{List of file names for output,}

\item{mode}{Mode to handle reads overlapping more than one feature
(choices: union, intersection-strict, intersection- nonempty; default: union)}

\item{type}{Feature type (3rd column in GTF file) to be used,}

\item{attribute}{GTF attribute to be used as feature ID,}

\item{format}{Type of <alignment_file> data, can be either sam, bam or auto}

\item{gtf}{Path to the GTF file}

\item{parallel}{Run in parallel, default set to FALSE}

\item{cores}{Number of cores/threads to use for parallel processing, default set to 4}

\item{execute}{Whether to execute the commands or not, default set to TRUE}

\item{htseq_count}{Path to the htseq-count program program, required}

\item{version}{Returns the version number}
}
\value{
A list with the htseq_count commands
}
\description{
Runs the htseq-count program
}
\examples{
\dontrun{
path <- "/home/gmh5n/.local/bin/htseq-count"

# Version
run_htseq_count(htseq_count = path,
                version = TRUE)

hisat2.alignments.dir <- "/path/to/hisat/alignments/directory"
bam.files <- list.files(path = hisat2.alignments.dir,
                       pattern = "sorted.bam$",
                       full.names = TRUE,
                       recursive = TRUE)
counts.dir <- "counts_dir"
htseq.counts.files <- paste(counts.dir,
                            paste(list.files(path = hisat2.alignments.dir,recursive = FALSE),
                            "counts.txt",sep = "_"),
                            sep = "/")

# Path to the gtf file
gtf.file <- "/path/to/gtffile"

htseq.counts.cmds <- run_htseq_count(input = bam.files,
                                     output = htseq.counts.files,
                                     mode = "intersection-nonempty",
                                     type = "gene",
                                     attribute = "ID",
                                     format = "bam",
                                     gtf = gtf.file,
                                     htseq_count = path)
htseq.counts.cmds
}
}
