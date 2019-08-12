# pipelineTools

<!-- badges: start -->
<!-- badges: end -->

The goal of pipelineTools is to streamline the NGS analysis pipelines and result reporting within RStudio. PipelineTools
provides packages to run standard open source NGS tools

## Installation
Can be installed using devtools directly from GitHub using the following commands.
```{r}
# Install devtools
install.packages("devtools")
#Install package from GitHub
install_github("GrahamHamilton/pipelineTools")
```

## Example
Load the required libraries
```{r load libraries}
library("pipelineTools")
library("biomaRt")
library("tximport")
```
Set the paths to the software installed on the system
```{r software paths, echo = FALSE}
fastq.screen.path <- "/software/fastq_screen_v0.13.0/fastq_screen"
fastp.path <- "/software/bin/fastp"
kallisto.path <- "/software/kallisto_v0.45.1/kallisto"
hisat.path <- "/software/hisat_v2-2.1.0/hisat2"
multiqc.path <- "/usr/local/bin/multiqc"
samtools.path <- "/software/samtools_v1.9/samtools"
stringtie.path <- "/software/stringtie-1.3.6/stringtie"
```

Version numbers for the software used
```{r versions, warning=FALSE, echo = FALSE}
fastqscreen.version <- run_fastq_screen(fastq_screen = fastq.screen.path, version = TRUE)
fastp.version <- run_fastp(fastp = fastp.path, version = TRUE)
kallisto.version <- run_kallisto(kallisto = kallisto.path, version = TRUE)
hisat.version <- run_hisat2(hisat2 = hisat.path, version = TRUE)
multiqc.version <- run_multiqc(multiqc = multiqc.path, version = TRUE)
samtools.version <- run_samtools(samtools = samtools.path, version = TRUE)
stringtie.version <- run_stringtie(stringtie = stringtie.path, version = TRUE)
r.version <- getRversion()
pipelineTools.version <- packageDescription("pipelineTools")$Version
deseq2.version <- packageDescription("DESeq2")$Version
```

Create the results directories
```{r results directories}
# Trimmed reads directory
trimmed.reads.dir <- "trimmed_reads"
#Create the directory for the trimmed reads
dir.create(trimmed.reads.dir, showWarnings = FALSE)

# FastQScreen results directory
fastq.screen.dir <- "Screen"
# Create the directory for the fastq screen results
dir.create(fastq.screen.dir, showWarnings = FALSE)

# FastP results directory
fastp.results.dir <- "FastpQC"
# Create the directory for the FastP results
dir.create(fastp.results.dir, showWarnings = FALSE)

# Kallisto results directory
kalisto.results.dir <- "kallisto"
#Create the directory for the Kallisto results
dir.create(kalisto.results.dir, showWarnings = FALSE)

# Hisat2 alignment results directory
hisat2.alignments.dir <- "hisat2_alignments"
#Create the directory for the HiSat2 results
dir.create(hisat2.alignments.dir, showWarnings = FALSE)

# Stringtie results directory
stringtie.dir <- "stringtie"
#Create the directory for the HiSat2 results
dir.create(stringtie.dir, showWarnings = FALSE)
```

Read in the fastq file paths to lists and then create sample names lists and trimmed read file paths
```{r setup files}
reads.path <- "raw_reads"

mate1 <- list.files(path = reads.path, pattern = "*_R1_001.fastq.gz$", full.names = TRUE)
mate2 <- list.files(path = reads.path, pattern = "*_R2_001.fastq.gz$", full.names = TRUE)

mate1.out <- paste(trimmed.reads.dir,(list.files(path = reads.path, pattern = "*_R1_001.fastq.gz$", full.names = FALSE)), sep = "/")
mate2.out <- paste(trimmed.reads.dir,(list.files(path = reads.path, pattern = "*_R2_001.fastq.gz$", full.names = FALSE)), sep = "/")

sample.names <- unlist(lapply(strsplit(list.files(path = reads.path, pattern = "*_R1_001.fastq.gz$", full.names = FALSE),"_"), `[[`, 1))
```

Sequence adapters
```{r sequence adapters}
adapter1 <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter2 <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
```

The full paths to the reference genome/transcriptome indexes and annotation files
```{r references}
# Path to the reference transcriptome
transcriptome <- "/path/to/kallisto/index"

# Path to the reference genome
genome <- "/path/to/hisat2/index"

# Path to the gtf file
gtf <- "/path/to/gtf/file"
```
### Fastq Screen
Results are stored in the Screen folder
```{r fastqscreen, echo = FALSE, eval = eval}
fastq.screen.conf <- "/export/jessie3/gmh5n/PipelineTest/fastq_screen.conf"
fastq.screen.cmds <- run_fastq_screen(fq.files = mate1,
                                      out.dir = fastq.screen.dir,
                                      conf = fastq.screen.conf,
                                      fastq_screen = fastq.screen.path)
write.table(fastq.screen.cmds,"FastqScreen_commands.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

### FastP
Adapter and quality trimmed reads are stored in the `r trimmed.reads.dir` directory, QC files are stored in the FastpQC folder
```{r fastp, echo = FALSE, eval = eval}
fastp.cmds <- run_fastp(mate1 = mate1,
                        mate2 = mate2,
                        mate1.out = mate1.out,
                        mate2.out = mate2.out,
                        adapter1 = adapter1,
                        adapter2 = adapter2,
                        sample.name =  sample.names,
                        out.dir = fastp.results.dir,
                        fastp = fastp.path)

write.table(fastp.cmds,"FastP_commands.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

### Kallisto
Pseudo align the reads to the reference transcriptome with Kallisto
```{r kallisto, echo = FALSE, eval = eval}
strandedness <- "second"
# Paired end
kallisto.cmds <- run_kallisto(mate1 = mate1,
                              mate2 = mate2,
                              index = transcriptome,
                              sample.name = sample.names,
                              strandedness = strandedness,
                              out.dir = kalisto.results.dir,
                              kallisto = kallisto.path)

write.table(kallisto.cmds,"Kallisto_commands.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

### TXImport for Kallisto data
```{r include=FALSE}
mart<-"ensembl"
db<-"mmusculus_gene_ensembl" # Change to organism in study
filt<-"ensembl_gene_id"

# Create the biomaRt object
ensembl = useEnsembl(biomart=mart, dataset=db)

# Get all the transcript ids and corresponding gene ids from BiomaRt
att<-c("ensembl_transcript_id","ensembl_gene_id")
txTable<-getBM(attributes=att,mart=ensembl)

# Set the column names for the transcript to gene table
colnames(txTable)<-c("tx_id","gene_id")

# Read in the experimental design file, tab seperated
sampleTable <- read.csv("SampleDescription.txt", sep="\t", row.names=1)

# Read in the file names form the kallisto results directory
dir <- getwd()
files <- file.path(dir,kalisto.results.dir,row.names(sampleTable),"abundance.h5", fsep = .Platform$file.sep)
names(files)<-row.names(sampleTable)
 
txi <- tximport(files, type = "kallisto", tx2gene = txTable,ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
```
