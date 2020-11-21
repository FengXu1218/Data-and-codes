# KM_Sc比较（2）
rm(list = ls())
library(tidyr)
# 载入注释文件，进行id转换——获得id和gene对应关系；
setwd("/Users/xufeng/Desktop/k_f_Sce/")
a <- read.table("GCa_000218975.1_aSM21897v1_genomic.gff",header = F,sep = "\t",quote ="")
library(stringr)
k <- which(a[,3]=="gene", arr.ind = TRUE)
b <- a[k,]
d <- read.table("EC1118_PRJEa37863.gff",header = F,sep = "\t",quote ="",nrows = 15060)
k <- which(d[,3]=="gene", arr.ind = TRUE)
e <- d[k,]
b <- tidyr::unite(b, "V4_V5", V4, V5)
e <- tidyr::unite(e, "V4_V5", V4, V5)
b1 <- b[b[,4] %in% e[,4],]
b1 <- b1[order(b1$V4),]
e1 <- e[e[,4] %in% b[,4],]
e1 <- e1[order(e1$V4),]

b1 <- b1[,c(4,8)]
e1 <- e1[,c(4,8)]
id <- str_split(b1[,2],";",simplify = T)
id <- as.matrix(str_split(id[,2],"=",simplify = T)[,2])
id_num <- cbind(b1[,1],id)

gene <- str_split(e1[,2],";",simplify = T)
gene <- as.matrix(str_split(gene[,4],"=",simplify = T)[,2])
gene <- cbind(e1[,1],gene)
gene <- gene[-grep("SGD",gene[,2]),]
colnames(gene) <- c("start_num","gene")
colnames(id_num) <- c("start_num","id")
id_gene <- merge(gene,id_num,by="start_num")
id_gene <- id_gene[,-1]

# 读入数据，获得Sc各个样本生物学重复的表达量；
files <- (Sys.glob("./GSE129483_RaW/*.tab"))
f<-length(files)
for (i in 1:32) {
  a <- read.table(file = files[i],header = T)
  a <- a[,c(1,10)]
  a[,1] <- str_split(a[,1],":",simplify = T)[,2]
  rownames(a) <- a[,1]
  colnames(a)[2] <- str_split(str_split(files[i],"W/",simplify = T)[,2],".t",simplify = T)[,1]
  if (i > 1) {
    KM <- merge(KM,a,by="tracking_id")
  }
  if (i < 2) {
    KM <- a
  }
}
for (i in 33:f) {
  a <- read.table(file = files[i],header = T)
  a <- a[,c(1,10)]
  a[,1] <- str_split(a[,1],":",simplify = T)[,2]
  rownames(a) <- a[,1]
  colnames(a)[2] <- str_split(str_split(files[i],"W/",simplify = T)[,2],".t",simplify = T)[,1]
  if (i > 33) {
    Sc <- merge(Sc,a,by="tracking_id")
  }
  if (i < 34) {
    Sc <- a
  }
}
colnames(Sc)[1] <- "id"
Sc_gene <- merge(id_gene,Sc,by = "id")
Sc_gene_nodup <- Sc_gene[!duplicated(Sc_gene[,1], fromLast=TRUE), ] 
Sc_gene_nodup2 <- Sc_gene[!duplicated(Sc_gene[,2], fromLast=TRUE), ] 
rownames(Sc_gene_nodup2) <- Sc_gene_nodup2[,2]
Sc_gene_nodup2 <- Sc_gene_nodup2[,c(-1,-2)]
write.csv(Sc_gene_nodup2,file ="Sc_exprSet.csv")

# KM对应的id转换；
setwd("/Users/xufeng/Desktop/others/reference/gff")
a <- read.table("GCF_001417885.1_Kmar_1.0_genomic.gff",header = F,sep = "\t",quote ="")
k <- which(a[,3]=="gene", arr.ind = TRUE)
k2 <- which(a[,3]=="exon", arr.ind = TRUE)
k3 <- k+2
d <- intersect(k2,k3)
gene <- a[k,]
gene <- str_split(str_split(gene[,9],"Name=",simplify = T)[,2],"locus_tag=",simplify = T)
gene[,1] <- str_split(gene[,1],";",simplify = T)[,1]
gene[,2] <- str_split(gene[,2],";",simplify = T)[,1]
gene <- gene[,-3]
ids <- gene

# KM数据合并id转换
rownames(KM) <- KM[,1]
ids[,c(1,2)] <- ids[,c(2,1)]
table(rownames(KM) %in% ids[,1])
exprSet=KM[rownames(KM) %in% ids[,1],]
ids <- ids[ids[,1] %in% rownames(KM),]
exprSet <- exprSet[,-1]
tmp = by(exprSet,
         ids[,2],
         function(x) rownames(x)[which.max(rowMeans(x))])
probe <- as.character(tmp)
exprSet=exprSet[rownames(exprSet) %in% probe ,]
ids = ids[match(rownames(exprSet),ids[,1]),]
rownames(exprSet) <- ids[match(rownames(exprSet),ids[,1]),2]
KM_nodup <- exprSet
setwd("/Users/xufeng/Desktop/k_f_Sce/")
write.csv(KM_nodup,file = "KM_exprSet.csv")


###########################################################################################################3
# Sc数据与KM倍性数据合并：
load("combat_all_s_fpkm.Rdata")
all_fpkm <- b
sce <- Sc_gene_nodup2
all_fpkm_jec <- all_fpkm[rownames(all_fpkm) %in% rownames(sce),]
all_fpkm_jec <- all_fpkm_jec[order(rownames(all_fpkm_jec)),]
sce_fpkm_jec <- sce[rownames(sce) %in% rownames(all_fpkm),]
sce_fpkm_jec <- sce_fpkm_jec[order(rownames(sce_fpkm_jec)),]
all_s_fpkm <- cbind(all_fpkm_jec,sce_fpkm_jec)
write.csv(all_s_fpkm,file = "all_Sc_combat_fpkm.csv")

# PCA和htclust
library(ggfortify)
# all_s的PCA
df=as.data.frame(t(all_s_fpkm))
group_list=c(rep("K",3),rep("aa",3),rep("bb",3),rep("ab",3),rep("bba",3),
             rep("aab",3),rep("f",3),rep("S12",2),rep("S34",2),rep("SC_6h",22),rep("SC_4h",13))
df$group=group_list 
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')

#htclust
plot(hclust(dist(t(all_s_fpkm))))
dist <- as.matrix(dist(t(all_s_fpkm)))
write.csv(dist,file = "all_dist.csv")

##########################################################################################################
# KM第二次数据与之前的数据合并：
KM <- KM[,-1]
all_fpkm_jec <- all_s_fpkm[rownames(all_s_fpkm) %in% rownames(KM_nodup),]
all_fpkm_jec <- all_fpkm_jec[order(rownames(all_fpkm_jec)),]
KM_fpkm_jec <- KM_nodup[rownames(KM_nodup) %in% rownames(all_s_fpkm),]
KM_fpkm_jec <- KM_fpkm_jec[order(rownames(KM_fpkm_jec)),]
all_s_KM_fpkm <- cbind(all_fpkm_jec,KM_fpkm_jec)
write.csv(all_s_KM_fpkm,file = "all_Sc_KM_fpkm.csv")

# PCA和htclust
library(ggfortify)
# all_s的PCA
df=as.data.frame(t(all_s_KM_fpkm))
group_list=c(rep("K",3),rep("aa",3),rep("bb",3),rep("ab",3),rep("bba",3),
             rep("aab",3),rep("f",3),rep("S12",2),rep("S34",2),rep("SC_6h",22),
             rep("SC_4h",13),rep("KM_6h",16),rep("KM_32h",16))
df$group=group_list 
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')

#htclust
plot(hclust(dist(t(all_s_KM_fpkm))))
dist <- as.matrix(dist(t(all_s_KM_fpkm)))
write.csv(dist,file = "all_s_KM_dist.csv")





