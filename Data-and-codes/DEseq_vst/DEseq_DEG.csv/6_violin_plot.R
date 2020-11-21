library(ggplot2)
library(ggpubr)
library(ggsci)
# ??ĳGO???????л???????????????ͼ?????Ƶ?ͼ?ĺ???
# ????С????ͼ????
rm(list = ls())
p <<- list()
GO_GEP <- function(GO_ID){
  setwd("/Users/xufeng/Desktop/others/Kluyveromyces/countdata2/DEseq_vst/DEseq_DEG.csv/")
  load("KM_GO.Rdata")
  goterm_gene <- go_term2gene[go_term2gene$go_term == GO_ID,]
  k_aa <- read.csv(file = "k_aa_DEseq_DEG.csv")
  k_bb <- read.csv(file = "k_bb_DEseq_DEG.csv")
  ab_f <- read.csv(file = "ab_f_DEseq_DEG.csv")
  k_ab <- read.csv(file = "k_ab_DEseq_DEG.csv")
  k_aa_FC <- k_aa[k_aa$X %in% goterm_gene$gene,c("X","log2FoldChange")]
  k_bb_FC <- k_bb[k_bb$X %in% goterm_gene$gene,c("X","log2FoldChange")]
  ab_f_FC <- ab_f[ab_f$X %in% goterm_gene$gene,c("X","log2FoldChange")]
  k_ab_FC <- k_ab[k_ab$X %in% goterm_gene$gene,c("X","log2FoldChange")]
  
  k_aa_FC <- k_aa_FC[order(k_aa_FC$X),]
  k_bb_FC <- k_bb_FC[order(k_bb_FC$X),]
  ab_f_FC <- ab_f_FC[order(ab_f_FC$X),]
  k_ab_FC <- k_ab_FC[order(k_ab_FC$X),]
  
  goterm_gene_k_aa <- goterm_gene[goterm_gene$gene %in% k_aa_FC$X,]
  data.file_k_aa <- cbind(k_aa_FC,goterm_gene_k_aa)[,-1]
  data.file_k_aa$Group <- "k_aa"
  k_aa_median <- median(data.file_k_aa$log2FoldChange)
  
  goterm_gene_k_ab <- goterm_gene[goterm_gene$gene %in% k_ab_FC$X,]
  data.file_k_ab <- cbind(k_ab_FC,goterm_gene_k_ab)[,-1]
  data.file_k_ab$Group <- "k_ab"
  k_ab_median <- median(data.file_k_ab$log2FoldChange)
  
  goterm_gene_k_bb <- goterm_gene[goterm_gene$gene %in% k_bb_FC$X,]
  data.file_k_bb <- cbind(k_bb_FC,goterm_gene_k_bb)[,-1]
  data.file_k_bb$Group <- "k_bb"
  k_bb_median <- median(data.file_k_bb$log2FoldChange)
  
  goterm_gene_ab_f <- goterm_gene[goterm_gene$gene %in% ab_f_FC$X,]
  data.file_ab_f <- cbind(ab_f_FC,goterm_gene_ab_f)[,-1]
  data.file_ab_f$Group <- "ab_f"
  ab_f_median <- median(data.file_ab_f$log2FoldChange)
  
  data.file <- rbind(data.file_k_aa,data.file_k_bb,data.file_ab_f)
  data.file$Group <- factor(data.file$Group,levels = c("k_aa","k_bb","ab_f"))
  
  population <- data.file_k_aa[,1]
  n <- nrow(data.file_k_aa)
  reps <- 1000
  results <- numeric(reps)
  for (i in 1:1000){
    samp <- sample(population, n,replace = T)
    results[i] <- median(samp)
  }
  var <- var(results)
  x2 <- (median(population))^2/var
  p__kaa <- round(1-pchisq(x2,1),7)
  
  population <- data.file_k_bb[,1]
  n <- nrow(data.file_k_bb)
  reps <- 1000
  results <- numeric(reps)
  for (i in 1:1000){
    samp <- sample(population, n,replace = T)
    results[i] <- median(samp)
  }
  var <- var(results)
  x2 <- (median(population))^2/var
  p__kbb <- round(1-pchisq(x2,1),7)
  
  population <- data.file_ab_f[,1]
  n <- nrow(data.file_ab_f)
  reps <- 1000
  results <- numeric(reps)
  for (i in 1:1000){
    samp <- sample(population, n,replace = T)
    results[i] <- median(samp)
  }
  var <- var(results)
  x2 <- (median(population))^2/var
  p__abf <- round(1-pchisq(x2,1),7)
  
  data.file$p_value <- c(rep(p__kaa,nrow(data.file_k_aa)),rep(p__kbb,nrow(data.file_k_bb)),rep(p__abf,nrow(data.file_ab_f)))
  a <- data.file$p_value == 0
  data.file[a,"p_value"] <- 1e-32
#  for (i in nrow(data.file)){
#    if (data.file$p_value[i] == 0){
#      data.file$p_value[i] = 1}
#    }
#  write.csv(data.file,file = paste0(GO_ID,"_violin.csv"))
  p6 <<- ggviolin(data.file,x = "Group", y = "log2FoldChange", fill = "Group",
           palette = c("#00AFBB", "#E7B800", "#FC4E07"),
           title = GO[j], 
           font.title = c(20,"black"),
           font.y = c(20),
           xlab="",
           ylab="log2FoldChange",
           ylim=c(-3,3),
           add = "boxplot", 
           add.params = list(fill = "white")) +
    #  scale_fill_npg() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    font("axis.text",size = 25) + 
    geom_hline(aes(yintercept=0),color="black",linetype = "dashed")
}

GO <- c("GO:0006260","GO:0022857","GO:0007010","GO:0043436","GO:0000278")
for (j in 5) {
  GO_GEP(GO_ID = GO[j])
}

# 绘制all_gene的提琴图：
k_aa_all_FC <- k_aa[,c(1,4)]
k_aa_all_FC$group <- rep('k_aa',nrow(k_aa_all_FC))
k_bb_all_FC <- k_bb[,c(1,4)]
k_bb_all_FC$group <- rep('k_bb',nrow(k_bb_all_FC))
ab_f_all_FC <- k_aa[,c(1,4)]
ab_f_all_FC$group <- rep('ab_f',nrow(ab_f_all_FC))
all_data.file <- rbind(k_aa_all_FC,k_bb_all_FC,ab_f_all_FC)
p1 <<- ggviolin(all_data.file,x = "group", y = "log2FoldChange", fill = "group",
                palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                title = 'All genes', 
                font.title = c(20,"black"),
                font.y = c(20),
                xlab="",
                ylab="log2FoldChange",
                ylim=c(-3,3),
                add = "boxplot", 
                add.params = list(fill = "white")) +
  #  scale_fill_npg() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  font("axis.text",size = 25) + 
  geom_hline(aes(yintercept=0),color="black",linetype = "dashed")

# 将六张图合并为一张大图：
ggarrange(p1, p2, p3, p4, p5,p6,nrow = 2, ncol = 3, labels = c('a', 'b', 'c', 'd','e','f'), 
          font.label = list(size = 20,color = 'black',face = "bold"))


