# DESEQ2 VIGNETTE
# 201500113
# 
# THESE COMMANDS DIDN'T WORK
# install.packages("airway")
# source("http://bioconductor.org/workflows.R")
# workflowInstall("rnaseqGene")
# 
source("http://bioconductor.org/biocLite.R")
biocLite("airway")
# 
# READING IN DATA FROM A SummarizedExperiment OBJECT
library("airway")
data("airway")
se <- airway
se
#
# ADDING DESIGN INFORMATION
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
ddsSE <- DESeqDataSet(se, design = ~ cell + dex)
ddsSE
# 
# COUNT MATRIX INPUT
source("http://bioconductor.org/biocLite.R")
biocLite("pasilla")
library("pasilla")
library("Biobase")
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq")
data("pasillaGenes")
countData <- counts(pasillaGenes)
head(countData)
class(countData)
colData <- pData(pasillaGenes)[,c("condition","type")]
colData
class(colData)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
dds
# 
# ADDING REDUNDANT METADATA JUST TO SHOW THAT YOU CAN
featureData <- data.frame(gene=rownames(pasillaGenes))
head(featureData)
(mcols(dds) <- DataFrame(mcols(dds), featureData))
# ?DataFrame
# mcols(dds)
# THERE'S ALSO A WORKFLOW TO CONSTRUCT THIS FROM THE HTSeq PYTHON PACKAGE
# 
# REMEMBER TO SET THE CONTROL AS THE BASE FACTOR, AND DROP LEVELS WHEN YOU REMOVE STUFF
# dds$condition <- relevel(dds$condition, "untreated")
# dds$condition <- droplevels(dds$condition)
# 
# CONSIDER collapseReplicates...
# 
# Results tables are generated using the function results, which extracts a results table with log2 fold
# changes, p values and adjusted p values. With no arguments to results, the results will be for the
# last variable in the design formula, and if this is a factor, the comparison will be the last level of this
# variable over the first level.
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
head(resOrdered)
# head(as.data.frame(res))
# 
# library(dplyr)
# res %>%
#   arrange(padj)
# # THIS DIDN'T WORK SOMEHOW, WRONG CLASS
# 
summary(res)
# The results function contains a number of arguments to customize the results table which is generated.
# Note that the results function automatically performs independent filtering based on the mean of
# counts for each gene, optimizing the number of genes which will have an adjusted p value below a given
# threshold. This will be discussed further in Section 3.8.
# 
# In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the
# mean of normalized counts. Points will be colored red if the adjusted p value is less than 0.1. Points
# which fall out of the window are plotted as open triangles pointing either up or down
plotMA(res, main="DESeq2", ylim=c(-2,2))
# The MA-plot of log2 fold changes returned by DESeq2 allows us to see how the shrinkage of fold changes
# works for genes with low counts. You can still obtain results tables which include the “unshrunken”
# log2 fold changes (for a simple comparison, the ratio of the mean normalized counts in the two groups).
# A column lfcMLE with the unshrunken maximum likelihood estimate (MLE) for the log2 fold change
# will be added with an additional argument to results:
resMLE <- results(dds, addMLE=TRUE)
head(resMLE, 4)
## log2 fold change (MAP): condition treated vs untreated
plotMA(resMLE, main = "Unshrunken", ylim=c(-2,2))
# 
# 1.4.2 Plot counts
# It can also be useful to examine the counts of reads for a single gene across the groups. A simple
# function for making this plot is plotCounts, which normalizes counts by sequencing depth and adds a
# pseudocount of 1/2 to allow for log scale plotting. The counts are grouped by the variables in intgroup,
# where more than one variable can be specified. Here we specify the gene which had the smallest p value
# from the results table created above. You can select the gene to plot by rowname or by numeric index.

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# For customized plotting, an argument returnData specifies that the function should only return a
# data.frame for plotting with ggplot.

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition",
                returnData=TRUE)
d$condition = relevel(d$condition, "untreated")
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))
# 
# VERY IMPORTANT ONE
mcols(res)$description
# 
# ADDING THE UNSHRUNKEN LFCS BACK IN
head(results(dds, addMLE=TRUE),4)
# 
#  LOOK INTO see the “RNA-seq differential expression” vignette at the ReportingTools page, or the
# manual page for the publish method for the DESeqDataSet class.
# 
write.csv(as.data.frame(resOrdered),
          file="condition_treated_results.csv")
resSig <- subset(resOrdered, padj < 0.1)
resSig
# 
# SKIPPED THE SECTION ON MULT-FACTOR DESIGN
# 
# DATA TRANSFORMATIONS
# 2.1.2 Extracting transformed values
# The two functions return SummarizedExperiment objects, as the data are no longer counts. The assay
# function is used to extract the matrix of normalized values:
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)
# rld
# vsd
# head(rlogMat)
# head(vstMat)
# The regularized log transformation is
# preferable to the variance stabilizing transformation if the size factors vary widely
# 
# GRAPHING FOR QC
library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
# 
# 2.2.2 Heatmap of the sample-to-sample distances
# Another use of the transformed data is sample clustering. Here, we apply the dist function to the
# transpose of the transformed count matrix to get sample-to-sample distances. We could alternatively
# use the variance stabilized transformation here.
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition, type, sep=" : "))
hc <- hclust(distsRL)
hc
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(13, 13))
# 
# 2.2.3 Principal component plot of the samples
# Related to the distance matrix of Section 2.2.2 is the PCA plot of the samples, which we obtain as
# follows (Figure 7).

plotPCA(rld, intgroup=c("condition", "type"))

# It is also possible to customize the PCA plot using the ggplot function.
data <- plotPCA(rld, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
# 
# QC
plotDispEsts(dds)
# 
# 3.7.1 Local or mean dispersion fit
# A local smoothed dispersion fit is automatically substitited in the case that the parametric curve doesn’t
# fit the observed dispersion mean relationship. This can be prespecified by providing the argument
# fitType="local" to either DESeq or estimateDispersions. Additionally, using the mean of genewise
# disperion estimates as the fitted value can be specified by providing the argument fitType="mean".
# 3.7.2 Supply a custom dispersion fit
# Any fitted values can be provided during dispersion estimation, using the lower-level functions described
# in the manual page for estimateDispersionsGeneEst. In the code chunk below, we store the genewise
# estimates which were already calculated and saved in the metadata column dispGeneEst. Then
# we calculate the median value of the dispersion estimates above a threshold, and save these values
# as the fitted dispersions, using the replacement function for dispersionFunction. In the last line,
# the function estimateDispersionsMAP, uses the fitted dispersions to generate maximum a posteriori
# (MAP) estimates of dispersion.
ddsCustom <- dds
useForMedian <- mcols(ddsCustom)$dispGeneEst > 1e-7
medianDisp <- median(mcols(ddsCustom)$dispGeneEst[useForMedian],na.rm=TRUE)
dispersionFunction(ddsCustom) <- function(mu) medianDisp
ddsCustom <- estimateDispersionsMAP(ddsCustom)
plotDispEsts(ddsCustom)
# 
# SKIPPED AHEAD TO SECTION ON INDEPENDANT FILTERING:
plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
# 
# I DON'T GET SECTION 4.6





