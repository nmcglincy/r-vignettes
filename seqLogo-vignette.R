source("http://bioconductor.org/biocLite.R")
biocLite("seqLogo")
library(seqLogo)
#
mFile <- system.file("Exfiles/pwm1", package="seqLogo")
m <- read.table(mFile)
m
str(m)
# 
p <- makePWM(m)
# makePWM() checks that all column probabilities add up to 1.0 and also obtains the information
# content profile and consensus sequence for the position weight matrix. These can then be
# accessed through the corresponding slots of the created object:
slotNames(p)
# [1] "pwm" "consensus" "ic" "width" "alphabet"
p@pwm
# 1 2 3 4 5 6 7 8
# A 0.0 0.0 0.0 0.3 0.2 0.0 0.0 0.0
# C 0.8 0.2 0.8 0.3 0.4 0.2 0.8 0.2
# G 0.2 0.8 0.2 0.4 0.3 0.8 0.2 0.8
# T 0.0 0.0 0.0 0.0 0.1 0.0 0.0 0.0
p@ic
# [1] 1.2780719 1.2780719 1.2780719 0.4290494 0.1535607 1.2780719 1.2780719
# [8] 1.2780719
p@consensus
# [1] "CGCGCGCG"
p@width
p@alphabet
# 
seqLogo(p)
seqLogo(p, ic.scale = FALSE)
# 
# 
# AN ALTERNATIVE IS THE WEBLOGO PACKAGE ON CRAN
install.packages("RWebLogo")
# 
# Make a sequence logo using an external alignment file format
# In this example we'll use the EMBOSS alignment format or msf
# However, you can use any format supported by WebLogo e.g. fasta
library(RWebLogo)
fpath = system.file("extdata", "example_data.msf", package="RWebLogo")
weblogo(file.in=fpath)
# Now for an example using an alignment as an R character vector
aln <- c('CCAACCCAA', 'CCAACCCTA', 'AAAGCCTGA', 'TGAACCGGA')
aln
# Simple WebLogo
weblogo(seqs=aln)
# Lets get rid of those ugly error bars and add some text!
weblogo(seqs=aln, errorbars=FALSE, title='Yay, No error bars!',
        fineprint='RWebLogo 1.0', label='1a')
# We can also change the format of the output like this
weblogo(seqs=aln, format='png', resolution=500)
# You can change the axis labels like this
weblogo(seqs=aln, xlabel='My x-axis', ylabel='Awesome bits')
# You get the idea! See ?weblogo for more awesome options!
?weblogo
