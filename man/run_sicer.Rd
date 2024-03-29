% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_sicer.R
\name{run_sicer}
\alias{run_sicer}
\title{Run Sicer2}
\usage{
run_sicer(
  treatment = NULL,
  control = NULL,
  comparison.names = NULL,
  species = NULL,
  redundancy_threshold = 1,
  window_size = 200,
  fragment_size = 150,
  effective_genome_fraction = 0.74,
  false_discovery_rate = 0.01,
  gap_size = 600,
  e_value = 1000,
  step_score = NULL,
  cpu = NULL,
  significant_reads = FALSE,
  parallel = FALSE,
  cores = 4,
  execute = TRUE,
  sicer = NULL
)
}
\arguments{
\item{treatment}{List of the paths to files containing the treatment files.
This can either be the relative or the absolute path of the file. Must be
in BED or BAM format.}

\item{control}{List of the paths to files containing the control files. This
can either be the relative or the absolute path of the file. Must be in BED
or BAM format. OPTIONAL.}

\item{comparison.names}{List of the comparisons to be made. OPTIONAL.}

\item{species}{The species/genome used (ex: hg38).}

\item{redundancy_threshold}{The number of copies of indentical reads allowed
in a library. Default is 1.}

\item{window_size}{Resolution of SICER. Default value is 200 (bp)}

\item{fragment_size}{The amount of shift from the beginning of a read to the
center of the DNA fragment represented by the read. Default is 150 (bp).}

\item{effective_genome_fraction}{Effective genome as fraction of the genome
size. Default is 0.74.}

\item{false_discovery_rate}{Remove all islands with an false_discovery_rate
below cutoff. Default is 0.01.}

\item{gap_size}{The minimum length of a "gap" such that neighboring window is
an "island." This value must be a multiple of the window size. Default is
600 (bp).}

\item{e_value}{Requires user input when no control library is provided.
Default is 1000.}

\item{step_score}{Step Score: The minimum number of positive elements in the
graining unit to call the unit positive. Used for RECOGNICER algorithm.}

\item{cpu}{The number of CPU cores SICER program will use when executing
multi-processing tasks. Optimal core count is the species' number of
chromosomes. Default value is the maximum number of cores avaiable in the
system.}

\item{significant_reads}{SICER produces a BED file of treatment reads
filtered by significant islands and WIG file of filtered reads binned into
windows}

\item{parallel}{Run in parallel, default set to FALSE}

\item{cores}{Number of cores/threads to use for parallel processing, default set to 4}

\item{execute}{Whether to execute the commands or not, default set to TRUE}

\item{sicer}{Path to the Sicer2 or recognicer programs}
}
\value{
A list with the Sicer2 commands
}
\description{
Sicer2 is a spatial clustering approach for the identification
  of ChIP-enriched regions which was developed for calling broad peaks of
  histone modifications from ChIP-seq data.
}
\examples{
\dontrun{
 sicer.cmds <- run_sicer(treatment = treatment.list,
                         control = control.list,
                         comparison.names = comp.list,
                         significant_reads = TRUE,
                         sicer = sicer.path
                        )
}
}
