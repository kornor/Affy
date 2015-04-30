setwd("~/Bioinformatics work/Affy stuff/Affy")

#source("http://bioconductor.org/biocLite.R")
#biocLite("oligo")
#biocLite("affy")
#biocLite("annotate")
#biocLite("limma")

library(affy)
library(oligo)
library(annotate)
library(limma)

##automatically reads in all .CEL files
AffyData <- ReadAffy()

#I made a pheno table from clinical info on the website
pheno <- read.table("Phenotype_MA.txt", sep = "\t", header = TRUE, row.names = 1)

status <- as.factor(pheno$Status)


## Make exprs set of data

