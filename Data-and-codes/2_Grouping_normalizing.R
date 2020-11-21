# grouping and normalizing
rm(list = ls())

load("countdata_2.Rdata")
countdata2K <- countdata_2[,c(1,2,15)]
colnames(countdata2K) <- c("R19043246LR01-K-1","R19043247LR01-K-2","R19048100LR01-K-3")
countdata2aa <- countdata_2[,c(3,16,17)]
countdata2bb <- countdata_2[,c(4,5,18)]
countdata2ab <- countdata_2[,c(6,7,8)]
countdata2bba <- countdata_2[,c(9,10,11)]
countdata2f <- countdata_2[,c(12,13,14)]
countdata2aab <- countdata_2[,c(19,20,21)]
k_aa <- cbind(countdata2K,countdata2aa)
k_bb <- cbind(countdata2K,countdata2bb)
k_ab <- cbind(countdata2K,countdata2ab)
k_bba <- cbind(countdata2K,countdata2bba)
k_aab <- cbind(countdata2K,countdata2aab)
k_f <- cbind(countdata2K,countdata2f)
aa_bb <- cbind(countdata2aa,countdata2bb)
aa_ab <- cbind(countdata2aa,countdata2ab)
aa_bba <- cbind(countdata2aa,countdata2bba)
aa_aab <- cbind(countdata2aa,countdata2aab)
aa_f <- cbind(countdata2aa,countdata2f)
bb_ab <- cbind(countdata2bb,countdata2ab)
bb_bba <- cbind(countdata2bb,countdata2bba)
bb_aab <- cbind(countdata2bb,countdata2aab)
bb_f <- cbind(countdata2bb,countdata2f)
ab_bba <- cbind(countdata2ab,countdata2bba)
ab_aab <- cbind(countdata2ab,countdata2aab)
ab_f <- cbind(countdata2ab,countdata2f)
bba_f <- cbind(countdata2bba,countdata2f)
bba_aab <- cbind(countdata2bba,countdata2aab)
aab_f <- cbind(countdata2aab,countdata2f)
save(k_aa,k_bb,k_ab,k_bba,k_aab,k_f,aa_bb,aa_ab,aa_bba,aa_aab,aa_f,bb_ab,bb_bba,bb_aab,bb_f,ab_bba,
     ab_aab,ab_f,bba_f,bba_aab,aab_f,file = "exprSet.Rdata")

# k_aa
library(DESeq2)
colData <- data.frame(row.names=colnames(k_aa),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=k_aa, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
k_aa_vst <- as.data.frame(assay(vsd))
save(dds,k_aa_vst,file = "./DEseq_vst/k_aa_vst.Rdata")

# k_bb
colData <- data.frame(row.names=colnames(k_bb),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=k_bb, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
k_bb_vst <- as.data.frame(assay(vsd))
save(dds,k_bb_vst,vsd,file = "./DEseq_vst/k_bb_vst.Rdata")

# k_ab
colData <- data.frame(row.names=colnames(k_ab),
                      group_list=factor(c(rep(1,3),rep(2,3))))  
dds <-DESeqDataSetFromMatrix(countData=k_ab, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
k_ab_vst <- as.data.frame(assay(vsd))
save(dds,k_ab_vst,vsd,file = "./DEseq_vst/k_ab_vst.Rdata")

# k_bba
colData <- data.frame(row.names=colnames(k_bba),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=k_bba, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
k_bba_vst <- as.data.frame(assay(vsd))
save(dds,k_bba_vst,vsd,file = "./DEseq_vst/k_bba_vst.Rdata")

# k_aab
colData <- data.frame(row.names=colnames(k_aab),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=k_aab, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
k_aab_vst <- as.data.frame(assay(vsd))
save(dds,k_aab_vst,vsd,file = "./DEseq_vst/k_aab_vst.Rdata")

# k_f
colData <- data.frame(row.names=colnames(k_f),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=k_f, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
k_f_vst <- as.data.frame(assay(vsd))
save(dds,k_f_vst,vsd,file = "./DEseq_vst/k_f_vst.Rdata")

# aa_bb
colData <- data.frame(row.names=colnames(aa_bb),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=aa_bb, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
aa_bb_vst <- as.data.frame(assay(vsd))
save(dds,aa_bb_vst,vsd,file = "./DEseq_vst/aa_bb_vst.Rdata")

# aa_ab
colData <- data.frame(row.names=colnames(aa_ab),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=aa_ab, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
aa_ab_vst <- as.data.frame(assay(vsd))
save(dds,aa_ab_vst,vsd,file = "./DEseq_vst/aa_ab_vst.Rdata")

# aa_bba
colData <- data.frame(row.names=colnames(aa_bba),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=aa_bba, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
aa_bba_vst <- as.data.frame(assay(vsd))
save(dds,aa_bba_vst,vsd,file = "./DEseq_vst/aa_bba_vst.Rdata")

# aa_aab
colData <- data.frame(row.names=colnames(aa_aab),
                      group_list=factor(c(rep(1,3),rep(2,3))))  
dds <-DESeqDataSetFromMatrix(countData=aa_aab, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
aa_aab_vst <- as.data.frame(assay(vsd))
save(dds,aa_aab_vst,vsd,file = "./DEseq_vst/aa_aab_vst.Rdata")

# aa_f
colData <- data.frame(row.names=colnames(aa_f),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=aa_f, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
aa_f_vst <- as.data.frame(assay(vsd))
save(dds,aa_f_vst,vsd,file = "./DEseq_vst/aa_f_vst.Rdata")

# bb_ab
colData <- data.frame(row.names=colnames(bb_ab),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=bb_ab, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
bb_ab_vst <- as.data.frame(assay(vsd))
save(dds,bb_ab_vst,vsd,file = "./DEseq_vst/bb_ab_vst.Rdata")

# bb_bba
colData <- data.frame(row.names=colnames(bb_bba),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=bb_bba, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
bb_bba_vst <- as.data.frame(assay(vsd))
save(dds,bb_bba_vst,vsd,file = "./DEseq_vst/bb_bba_vst.Rdata")

# bb_aab
colData <- data.frame(row.names=colnames(bb_aab),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=bb_aab, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
bb_aab_vst <- as.data.frame(assay(vsd))
save(dds,bb_aab_vst,vsd,file = "./DEseq_vst/bb_aab_vst.Rdata")

# bb_f
colData <- data.frame(row.names=colnames(bb_f),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=bb_f, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
bb_f_vst <- as.data.frame(assay(vsd))
save(dds,bb_f_vst,vsd,file = "./DEseq_vst/bb_f_vst.Rdata")

# ab_bba
colData <- data.frame(row.names=colnames(ab_bba),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=ab_bba, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
ab_bba_vst <- as.data.frame(assay(vsd))
save(dds,ab_bba_vst,vsd,file = "./DEseq_vst/ab_bba_vst.Rdata")

# ab_aab
colData <- data.frame(row.names=colnames(ab_aab),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=ab_aab, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
ab_aab_vst <- as.data.frame(assay(vsd))
save(dds,ab_aab_vst,vsd,file = "./DEseq_vst/ab_aab_vst.Rdata")

# ab_f
colData <- data.frame(row.names=colnames(ab_f),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=ab_f, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
ab_f_vst <- as.data.frame(assay(vsd))
save(dds,ab_f_vst,vsd,file = "./DEseq_vst/ab_f_vst.Rdata")

# bba_f
colData <- data.frame(row.names=colnames(bba_f),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=bba_f, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
bba_f_vst <- as.data.frame(assay(vsd))
save(dds,bba_f_vst,vsd,file = "./DEseq_vst/bba_f_vst.Rdata")

# bba_aab
colData <- data.frame(row.names=colnames(bba_aab),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=bba_aab, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
bba_aab_vst <- as.data.frame(assay(vsd))
save(dds,bba_aab_vst,vsd,file = "./DEseq_vst/bba_aab_vst.Rdata")

# aab_f
colData <- data.frame(row.names=colnames(aab_f),
                      group_list=factor(c(rep(1,3),rep(2,3)))) 
dds <-DESeqDataSetFromMatrix(countData=aab_f, 
                             colData=colData, 
                             design=~group_list)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
aab_f_vst <- as.data.frame(assay(vsd))
save(dds,aab_f_vst,vsd,file = "./DEseq_vst/aab_f_vst.Rdata")
