# VINGETTE FOR MGSA
library(mgsa)
data("example")
# 
# GO ANNOTATIONS
example_go
# 
# THE IDENTIFIED GENES
example_o

set.seed(0)
# 
# FIT THE MODEL
fit = mgsa(example_o, example_go)
fit
plot(fit)
# 
# EXTRACTING INTERESTING RESULTS AS A DATAFRAME
res = setsResults(fit)
subset(res, estimate > 0.5)
# 
# READING FULL GAF FILE
go.full = readGAF(filename = "gene_association.sgd")
