setwd("~/Bioinformatics work/Affy stuff/Affy")

#source("http://bioconductor.org/biocLite.R")
#biocLite("oligo")
#biocLite("affy")
#biocLite("annotate")
#biocLite("limma")
#biocLite("affycoretools")



library(affy)
#library(oligo)
library(annotate)
library(limma)
library(affycoretools)

#I made a pheno table from clinical info on the website
pheno <- read.AnnotatedDataFrame("Phenotype_MA.txt", sep = "\t", header = TRUE, row.names = 1)
##makes this file a "pheno" class file for affy package
pData(pheno)

## Read in and make expression set of data 
# function automatically finds all .CEL files in directory
##RMA is normalising function also
eset <- justRMA(phenoData = pheno)
#with no pehno data to include
#eset <- justRMA()
exp <- exprs(eset)

#you can save out the exp file as you like (big tho)

write.table(exp, "Expression_BRCA_set.txt", sep = "\t")

### set some subgroups for maybe later use
group <- (factor(eset$Status))

Non <- which(eset$Status == "Non-carrier")
BRCA1 <- which(eset$Status == "BRCA1 mutation carrier")
BRCA2 <- which(eset$Status == "BRCA2 mutation carrier")


### Check the effect of rma (normalization)
boxplot(exp)

##looks ok

#Do a PCA
plotPCA(eset, groups = group, groupnames = levels(group))

## processing
### This makes a model fit among groups and tells you top most diff expr

design <- model.matrix(~group)
fit <- lmFit(eset, design)
ebayes <- eBayes(fit)

##Get the top most differentially expressed genes
## change "n = " to whatever number you want; top10, top50 for eg
tab <- topTable(ebayes, coef = 2, adjust = "fdr", n = 50)

##only look at the top ones with a p.value of <0.05
select <- p.adjust(ebayes$p.value[,2])<0.05
esetSel <- eset[select,]

esetSel

heatmap(exprs(esetSel), labCol = group)

###Annotation time!   This is very variable depending on what you are wanting to do

##find the annoation level of this array
eset@annotation

##Now download the annoation package for that (to go with "annotate" library from above)

biocLite("hgu133a.db")
biocLite("R2HTML")
library(hgu133a.db)
library(R2HTML)


### Make a subset of only the top
list <- match(rownames(tab), rownames(eset))
eset2 <- eset [list,] 

ID <- featureNames(eset2)
Symbol <- getSYMBOL(ID, "hgu133a.db")
Name <- as.character(lookUp(ID, "hgu133a.db", "GENENAME"))
Ensembl <- as.character(lookUp(ID, "hgu133a.db", "ENSEMBL"))

Ensembl <- ifelse(Ensembl=="NA", NA,
                  paste("<a href='http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=",
                        Ensembl, "'>", Ensembl, "</a>", sep=""))
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, Ensembl=Ensembl, stringsAsFactors=F)
tmp[tmp=="NA"] <- NA # The stringsAsFactors makes "NA" characters. This fixes that problem.

HTML(tmp, "out.html", append=F)
write.table(tmp ,file="target.txt",row.names=F, sep="\t")

##add the gene names into the heatmap
##make a genenames column in tab, from the tmp symbols
tab$genenames <- tmp$Symbol

sizeGrWindow(12,9)
par(mar=c(5,4,4,4))
heatmap(exprs(esetSel), labCol = group, labRow = tab$genenames)


