# bioconductor RNA-seq lab
# 
library(airway)
dir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(dir)
csvfile = file.path(dir, "sample_table.csv")
# csvfile
(sampleTable = read.csv(csvfile, row.names = 1))
