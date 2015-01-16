# VINGETTE FOR MGSA
library(mgsa)
data("example")
example_go
example_o
# 
# FIT THE MODEL
set.seed(0)
fit = mgsa(example_o, example_go)
fit
plot(fit)
# 
# READING FULL GAF FILE
# DOWNLOADED FROM http://geneontology.org/page/download-annotations ON 20150115
go.full = readGAF(filename = "gene_association.sgd")
go.full
fit.full = mgsa(o = example_o, 
                sets = go.full)
# 
# Error in matrix(raw$marg, ncol = restarts) : 
#   'data' must be of a vector type, was 'NULL'
# In addition: Warning message:
#   In mgsa(o = example_o, sets = go.full) :
#   Specified observation contains 43 unmapple items. Excluded them from calculation.
# 
traceback()
# 4: matrix(raw$marg, ncol = restarts)
# 3: mgsa.wrapper(items, sets@sets, sets@numberOfItems, alpha, beta, 
#                 p, steps, restarts, threads)
# 2: mgsa(o = example_o, sets = go.full)
# 1: mgsa(o = example_o, sets = go.full)
# 
library(devtools)
devtools::session_info()
# 
# Session info-----------------------------------------------------------------------------------------------
#   setting  value                       
# version  R version 3.1.2 (2014-10-31)
# system   x86_64, darwin13.4.0        
# ui       RStudio (0.98.1091)         
# language (EN)                        
# collate  en_US.UTF-8                 
# tz       America/Los_Angeles         
# 
# Packages---------------------------------------------------------------------------------------------------
#   package       * version   date       source        
# AnnotationDbi * 1.28.1    2014-10-28 Bioconductor  
# Biobase       * 2.26.0    2014-10-14 Bioconductor  
# BiocGenerics  * 0.12.1    2014-11-14 Bioconductor  
# bitops          1.0.6     2013-08-17 CRAN (R 3.1.0)
# caTools         1.17.1    2014-09-10 CRAN (R 3.1.1)
# DBI           * 0.3.1     2014-09-24 CRAN (R 3.1.1)
# devtools      * 1.6.1     2014-10-07 CRAN (R 3.1.1)
# gdata           2.13.3    2014-04-06 CRAN (R 3.1.0)
# GenomeInfoDb  * 1.2.4     2014-12-19 Bioconductor  
# GO.db         * 3.0.0     2014-09-26 Bioconductor  
# gplots        * 2.16.0    2015-01-07 CRAN (R 3.1.2)
# gtools          3.4.1     2014-05-28 CRAN (R 3.1.0)
# IRanges       * 2.0.1     2014-12-12 Bioconductor  
# KernSmooth      2.23.13   2014-09-14 CRAN (R 3.1.2)
# mgsa          * 1.14.2    2014-11-08 Bioconductor  
# RSQLite       * 1.0.0     2014-10-25 CRAN (R 3.1.2)
# rstudio         0.98.1091 2015-01-06 local         
# rstudioapi      0.2       2014-12-31 CRAN (R 3.1.2)
# S4Vectors     * 0.4.0     2014-10-14 Bioconductor  

# POSTED ONTO BIOCONDUCTOR SUPPORT PAGE
