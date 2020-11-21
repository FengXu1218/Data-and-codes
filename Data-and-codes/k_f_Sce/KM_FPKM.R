# 将KM酵母RNAseq数据经过FPKM处理；
rm(list = ls())
setwd("/Users/xufeng/Desktop/others/AS_Kluyveromyces/countdata2/")
b <- read.table(file = "all.idKIM1.txt",sep = "\t",header = T)
rownames(b)<-b[,1]
b<-b[,-1]
kb <- b$Length / 1000
kb
countdata <- b[,6:ncol(b)]
rpk <- countdata / kb
rpk
fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
k_f_fpkm <- fpkm[,c(1,2,15,12,13,14)]
all_fpkm <- fpkm[,c(1,2,15,16,17,3,4,5,18,6,7,8,9,10,11,19,20,21,12,13,14)]
colnames(k_f_fpkm) <- c("k_1","k_2","k_3","f_1","f_2","f_3")
colnames(all_fpkm) <- c("k_1","k_2","k_3","aa_1","aa_2","aa_3","bb_2","bb_3","bb_1",
                        "ab_1","ab_2","ab_3","bba_1","bba_2","bba_3","aab_1","aab_2",
                        "aab_3","f_1","f_2","f_3")
write.csv(k_f_fpkm,file="k_f_fpkm.csv")
write.csv(all_fpkm,file="all_fpkm.csv")

# ID转换：
setwd("/Users/xufeng/Desktop/others/reference/gff")
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
table(rownames(k_f_fpkm) %in% ids[,1])
exprSet=k_f_fpkm[rownames(k_f_fpkm) %in% ids[,1],]
tmp = by(exprSet,
         ids[,2],
         function(x) rownames(x)[which.max(rowMeans(x))])
probe <- as.character(tmp)
exprSet=exprSet[rownames(exprSet) %in% probe ,]
ids = ids[match(rownames(exprSet),ids[,1]),]
rownames(exprSet) <- ids[match(rownames(exprSet),ids[,1]),2]
k_f_fpkm <- exprSet
save(k_f_fpkm,file = "k_f_fpkm.Rdata")
number2 <- matrix(c(drugA,drugB),2,5)
write.csv(k_f_fpkm,file = "k_f_fpkm.csv")

# 合并数据集
setwd("/Users/xufeng/Desktop/k_f_Sce/")
sce <- read.table(file = "GSE53720_FPKM_NR_CR.txt",sep = "\t",header = T)
rownames(sce)<-sce[,1]
sce<-sce[,-1]
k_f_fpkm_jec <- k_f_fpkm[rownames(k_f_fpkm) %in% rownames(sce),]
k_f_fpkm_jec <- k_f_fpkm_jec[order(rownames(k_f_fpkm_jec)),]
sce_fpkm_jec <- sce[rownames(sce) %in% rownames(k_f_fpkm),]
sce_fpkm_jec <- sce_fpkm_jec[,3:6]
k_f_s_fpkm <- cbind(k_f_fpkm_jec,sce_fpkm_jec)
write.csv(k_f_s_fpkm,file = "k_f_s_fpkm.csv")
plot(hclust(dist(t(k_f_s_fpkm))))
dist <- as.matrix(dist(t(k_f_s_fpkm)))
write.csv(dist,file = "dist.csv")
