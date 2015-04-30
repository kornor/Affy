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



## Read in and make expression set of data 
# function automatically finds all .CEL files in directory
eset <- justRMA()
exp <- exprs(eset)

#you can save out the exp file as you like (big tho)

#I made a pheno table from clinical info on the website
pheno <- read.table("Phenotype_MA.txt", sep = "\t", header = TRUE, row.names = 1)

status <- as.factor(pheno$Status)



