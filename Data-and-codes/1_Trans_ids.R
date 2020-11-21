# trans_ids
setwd("./reference/gff")
rm(list = ls())
a <- read.table("KM_Y179_gff_final.GFF3",header = F,sep = "\t",quote ="")
library(stringr)
k <- which(a[,3]=="gene", arr.ind = TRUE)
k2 <- which(a[,3]=="exon", arr.ind = TRUE)
k3 <- k+2
d <- intersect(k2,k3)
gene <- a[k,]
exon <- a[k3,]
gene <- str_split(gene[,9],";",simplify = T)
gene <- gene[,2]
gene <- str_split(gene,"=",simplify = T)[,2]
gene <- as.data.frame(gene)
id <- str_split(exon[,9],";",simplify = T)
id <- id[,1]
id <- str_split(id,"=",simplify = T)[,2]
id <- as.data.frame(id)
ids <- cbind(id,gene)
load("../../countdata2.Rdata")
table(rownames(countdata2) %in% ids[,1])
exprSet=countdata2[rownames(countdata2) %in% ids[,1],]
tmp = by(exprSet,
         ids[,2],
         function(x) rownames(x)[which.max(rowMeans(x))])
probe <- as.character(tmp)
exprSet=exprSet[rownames(exprSet) %in% probe ,]
ids = ids[match(rownames(exprSet),ids[,1]),]
rownames(exprSet) <- ids[match(rownames(exprSet),ids[,1]),2]
countdata_2 <- exprSet
save(countdata_2,file = "countdata_2.Rdata")
