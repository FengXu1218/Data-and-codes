setwd("/Users/xufeng/Desktop/others/reference/gff")
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
ids[,2] <- toupper(ids[,2])
all_fpkm <- read.csv(file = "/Users/xufeng/Desktop/others/AS_Kluyveromyces/countdata2/all_fpkm.csv")
rownames(all_fpkm) <- all_fpkm[,1]
all_fpkm <- all_fpkm[,-1]
table(rownames(all_fpkm) %in% ids[,1])
exprSet=all_fpkm[rownames(all_fpkm) %in% ids[,1],]
tmp = by(exprSet,
         ids[,2],
         function(x) rownames(x)[which.max(rowMeans(x))])
probe <- as.character(tmp)
exprSet=exprSet[rownames(exprSet) %in% probe ,]
ids = ids[match(rownames(exprSet),ids[,1]),]
rownames(exprSet) <- ids[match(rownames(exprSet),ids[,1]),2]
all_fpkm <- exprSet
save(all_fpkm,file = "all_fpkm_nosmall.Rdata")
write.csv(all_fpkm,file = "all_fpkm_nosmall.csv")

# 合并数据集
setwd("/Users/xufeng/Desktop/k_f_Sce/")
sce <- read.table(file = "GSE53720_FPKM_NR_CR.txt",sep = "\t",header = T)
rownames(sce)<-sce[,1]
sce<-sce[,-1]
all_fpkm_jec <- all_fpkm[rownames(all_fpkm) %in% rownames(sce),]
all_fpkm_jec <- all_fpkm_jec[order(rownames(all_fpkm_jec)),]
sce_fpkm_jec <- sce[rownames(sce) %in% rownames(all_fpkm),]
sce_fpkm_jec <- sce_fpkm_jec[,3:6]
all_s_fpkm <- cbind(all_fpkm_jec,sce_fpkm_jec)
write.csv(all_s_fpkm,file = "all_s_fpkm.csv")

# PCA和htclust
library(ggfortify)
# all_s的PCA
df=as.data.frame(t(all_s_fpkm))
group_list=c(rep("K",3),rep("aa",3),rep("bb",3),rep("ab",3),rep("bba",3),
             rep("aab",3),rep("f",3),rep("S12",2),rep("S34",2))
df$group=group_list 
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')
# k_f_s的PCA
df=as.data.frame(t(k_f_s_fpkm))
group_list=c(rep("K",3),rep("f",3),rep("S12",2),rep("S34",2))
df$group=group_list 
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')
#htclust
plot(hclust(dist(t(all_s_fpkm))))
dist <- as.matrix(dist(t(all_s_fpkm)))
write.csv(dist,file = "all_dist.csv")
