# bioconductor RNA-seq lab
# http://www.bioconductor.org/help/workflows/rnaseqGene/
library(airway)
dir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(dir)
csvfile = file.path(dir, "sample_table.csv")
# csvfile
(sampleTable = read.csv(csvfile, row.names = 1))
filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))
# filenames
library("Rsamtools")
#
# specifies only counting 2M reads at a time
bamfiles <- BamFileList(filenames, yieldSize=2000000)
# bamfiles
# 
# making sure the chromosome names are the same in the bam files and the gff/gtf
seqinfo(bamfiles[1])
# 
# Next, we need to read in the gene model which will be used for counting reads. We will read the gene model from a GTF file, using makeTranscriptDbFromGFF from the GenomicFeatures package. GTF files can be downloaded from Ensembl's FTP site or other gene model repositories. A TranscriptDb object is a database which can be used to generate a variety of range-based objects, such as exons, transcripts, and genes. We will want to make a list of exons grouped by gene.
# There are other options for constructing a TranscriptDB. For the known genes track from the UCSC Genome Browser, one can use the pre-built Transcript DataBase: TxDb.Hsapiens.UCSC.hg19.knownGene. The makeTranscriptDbFromBiomart function can be used to automatically pull a gene model from Biomart.
# 
library("GenomicFeatures")
gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
(txdb <- makeTranscriptDbFromGFF(gtffile, format="gtf"))
(genes <- exonsBy(txdb, by="gene"))
#
# the actual counting
library("GenomicAlignments")
se <- summarizeOverlaps(features=genes, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
se
# the actual counts are accessed by the assay function
head(assay(se))
colSums(assay(se))
# colData is sample information, but this has none
colData(se)
# rowData is gene information - GRanges object
rowData(se)
# 
# Note that the rowData slot is a GRangesList, which contains all the information about the exons for each gene, i.e., for each row of the count matrix. It also contains metadata about the construction of the gene model in the metadata slot.
str(metadata(rowData(se)))
# The colData slot, so far empty, should contain all the metadata. We hence assign our sample table to it:
(colData(se) <- DataFrame(sampleTable))
#
# Should come back to the rest of this - has some good ideas for eda and clustering etc
# I never did a straight replicate v replicate scatter plot of the size-normalized and/or rlogged counts

#############
# Exploring the alternative: featureCounts() from package Rsubread
# 
# Need some mapped reads first, from the vignette
source("http://bioconductor.org/biocLite.R")
biocLite("Rsubread")
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
ref
buildindex(basename="reference_index",reference=ref)
reads <- system.file("extdata","reads.txt.gz",package="Rsubread")
align(index="reference_index",readfile1=reads,output_file="alignResults.BAM")
ann <- data.frame(GeneID=c("gene1","gene1","gene2","gene2"),
                  Chr="chr_dummy",
                  Start=c(100,1000,3000,5000),
                  End=c(500,1800,4000,5500),
                  Strand=c("+","+","-","-"),
                  stringsAsFactors=FALSE)
ann
fc_SE <- featureCounts("alignResults.BAM",annot.ext=ann)
fc_SE
# Yeah, this is cool because the annotation is so flexible, but making sure only the right alignments 
# are counted seems more straightforward with summerizeOverlaps




