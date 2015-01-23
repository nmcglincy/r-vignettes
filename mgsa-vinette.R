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
# 
# 
vignette("mgsa")
mgsa( c("A", "B"), list(set1=LETTERS[1:3], set2=LETTERS[2:5]) )

gaf.data = read.table("gene_association.sgd",
                      header = FALSE,
                      skip = 26,
                      stringsAsFactors = FALSE,
                      sep = "\t",
                      fill = TRUE,
                      row.names = NULL,
                      quote = "")
# Error in scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  : 
#                 line 16253 did not have 17 elements
# 
# NO IDEA WHAT THIS IS REFERING TO, LOOKS FINE TO ME - try adding fill = TRUE
# 
# Warning message:
#   In scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  :
#             EOF within quoted string
# http://stackoverflow.com/questions/17414776/read-csv-warning-eof-within-quoted-string-prevents-complete-reading-of-file
# dim(gaf.data)
# 
#### 
# Columns are:
#   
#   1: DB, database contributing the file (always "SGD" for this file).
# 2: DB_Object_ID, SGDID (SGD's unique identifier for genes and
#     features).
#  3: DB_Object_Symbol, see below
#  4: Qualifier (optional), one or more of 'NOT', 'contributes_to',
#     'colocalizes_with' as qualifier(s) for a GO annotation, when needed,
#     multiples separated by pipe (|)
#  5: GO ID, unique numeric identifier for the GO term
#  6: DB:Reference(|DB:Reference), the reference associated with the GO
#     annotation
#  7: Evidence, the evidence code for the GO annotation
#  8: With (or) From (optional), any With or From qualifier for the GO
#     annotation
#  9: Aspect, which ontology the GO term belongs (Function, Process or
#     Component)
# 10: DB_Object_Name(|Name) (optional), a name for the gene product in
#     words, e.g. 'acid phosphatase'
# 11: DB_Object_Synonym(|Synonym) (optional), see below
# 12: DB_Object_Type, type of object annotated, e.g. gene, protein, etc.
# 13: taxon(|taxon), taxonomic identifier of species encoding gene
#     product
# 14: Date, date GO annotation was defined in the format YYYYMMDD
# 15: Assigned_by, source of the annotation (always "SGD" for this file)
# 16: Annotation Extension, optional, Contains cross references to other
#     ontologies that can be used to qualify or enhance the
#     annotation. The cross-reference is prefaced by an appropriate GO
#     relationship; references to multiple ontologies can be
#     entered. For example, if a gene product is localized to the
#     mitochondria of lymphocytes, the GO ID (column 5) would be
#     mitochondrion ; GO:0005439, and the annotation extension column
#     would contain a cross-reference to the term lymphocyte from the
#     Cell Type Ontology.
# 17: Gene Product Form ID, optional, this field allows the annotation
#     of specific variants of that gene or gene
#     product. UniProtKB:P12345-2
# 
# Note on SGD nomenclature, pertains to columns 3 and 11:
# 
# Column 3 - When a Standard Gene Name (e.g. CDC28, COX2) has been
# conferred, it will be present in Column 3. When no Gene Name has been
# conferred, the Systematic Name (e.g. YAL001C, YGR116W, YAL034W-A) will
# be present in column 3.
# 
# Column 11 - The Systematic Name (e.g. YAL001C, YGR116W, YAL034W-A,
# Q0010) will be the first name present in Column 11. Any other names
# (except the Standard Name, which will be in Column 3 if one exists),
# including Aliases used for the gene will also be present in this
# column.
# ###
names(gaf.data) = c("DB",
					"SGDID",
					"DB_Symbol",
					"Qualifier",
					"GO_ID",
					"DB_Reference",
					"Evidence",
					"With",
					"Aspect", 
					"DB_Object_Name",
					"DB_Object_Synonym",
					"DB_Object_Type",
					"taxon",
					"Date",
					"Assigned_by",
					"Annotation_Extension",
					"Gene_Product_Form_ID")
# head(gaf.data)
# 
# Split by Aspect to make separate analyses of Function, Process or Component
aspect.l = list(Function  = subset(gaf.data, Aspect == "F"),
                Process   = subset(gaf.data, Aspect == "P"),
                Component = subset(gaf.data, Aspect == "C"))
# 
# From each subset, select the SGDID & GO_ID
selectSGDIDandGOID = function(x) {
  require(dplyr)
  x %>%
    select(SGDID, GO_ID)
}
aspect.l.slm = lapply(aspect.l, selectSGDIDandGOID)
str(aspect.l.slm)
library(plyr)
aspect.l.slm.l = lapply(aspect.l.slm, dlply, .(GO_ID))
# 
# Within one GO term, SGDIDs are represented multiple times because there is one entry
# per evidence code, also were there are different sources of information for the same
# evidence code. I can't see the sense in maintaining the repitition if I'm removing 
# the rest of the information, if I want to sort by evidence code, then I'll have to 
# do it at the original subset command
# 
# Make each bit into a list of SGDID by GO_ID
aspect.l.final = lapply(aspect.l.slm.l, lapply, function(x) {x = unique(as.vector(x$SGDID))})
str(aspect.l.final)
# Ok, that looks like how I wanted
# 
# Converting from symbol or systematic name to SGDID
example_o
library(org.Sc.sgd.db)
ls("package:org.Sc.sgd.db")
# these tables look interesting:
summary(org.Sc.sgdCOMMON2ORF)
# symbol (accepted?) to systematic name
summary(org.Sc.sgdDESCRIPTION)
# systematic name to description
summary(org.Sc.sgdSGD)
# systematic name to SGDID
summary(org.Sc.sgdALIAS)
# systematic name to alias (like an alternative symbol?)
summary(org.Sc.sgdALIAS2ORF)
# the same as above?
summary(org.Sc.sgdGENENAME)
# don't understand this one - the first couple entries are the same in both tables
# 
# attempt to map example_o > systematic_name > SGDID
example_o

# Convert to a list
xx <- as.list(org.Sc.sgdCOMMON2ORF)
str(xx)
length(xx)
# Remove probes that do not map in COMMON2ORF
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
  # Gets the ORF identifiers for the first five gene names/alias
  xx[1:5]
  # Get the first one
  xx[[1]]
}
foo = xx[example_o]
# PROBLEM 
#   one has an NA
#   two have more than one association
#     looked at both manually in SGD; both times the first hit was the right one, might
#     be a general strategy
lapply(foo, length)
bar = ldply(lapply(foo, function(x) {x = x[1]}))
names(bar) = c("symbol", "sys.id")
bar
yy = as.list(org.Sc.sgdSGD)
head(yy)

summary(org.Sc.sgdGO)
tail(yy)
length(yy)
yy = yy[!is.na(yy)]
yy[bar$sys.id]

foo2 = yy[bar$sys.id]
bar2 = ldply(lapply(foo2, function(x) {x = x[1]}))
names(bar2) = c("sys.id", "sgdid")

banana = merge(bar, bar2,
               by.x = "sys.id",
               by.y = "sys.id")
banana
length(example_o)
# 
# The only problem with all of this is that I'm not sure whether this preserves all the 
# parent-child relationships in the terms, and I have this nagging feeling the lists are
# too short
?mgsa
unfilt.gaf.func = new("MgsaSets", sets = aspect.l.final$Function)
fit.func = mgsa(banana$sgdid, unfilt.gaf.func)
fit.func
res.func = setsResults(fit.func)
subset(res.func, estimate > 0.5)


unfilt.gaf.proc = new("MgsaSets", sets = aspect.l.final$Process)
fit.proc = mgsa(banana$sgdid, unfilt.gaf.proc)
plot(fit.proc)

unfilt.gaf.comp = new("MgsaSets", sets = aspect.l.final$Component)
fit.comp = mgsa(banana$sgdid, unfilt.gaf.comp)
plot(fit.comp)
fit.comp

str(unfilt.gaf.comp)
names(unfilt.gaf.comp@sets)
# 
# Would be cool to incorporate the item & set annotations, should be able to port
# them over from the GAF, or through the tables
summary(org.Sc.sgdGO)
summary(org.Sc.sgdGO2ORF)
summary(org.Sc.sgdGO2ALLORFS)

library(GO.db)
ls("package:GO.db")
xx <- as.list(GOTERM)
head(xx)
str(xx[1])
res.func.sig = subset(res.func, estimate > 0.5)
res.func.sig
sig.func.go = rownames(res.func.sig)
sig.func.go

foo = xx[sig.func.go]
ldply(foo)
str(foo)
str(lapply(foo, unlist))
lapply(foo, Term)

