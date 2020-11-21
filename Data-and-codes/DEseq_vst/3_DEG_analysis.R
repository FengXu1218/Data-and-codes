# DEG_analysis
rm(list = ls())
library(DESeq2)
library(ggpubr)
library(stringr)
library(gplots)
setwd("./DEseq_vst/")
load("../exprSet_vst.Rdata")
files <- (Sys.glob("*.Rdata"))
f<-length(files)
for (i in 1:f) {
  load(files[i])
  a <- str_split(files[i],"_vst",simplify = T)[,1]
  assay(vsd) <- get(a)
  exprSet <- as.data.frame(assay(vsd))
  res <- results(dds, tidy=TRUE)
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG=as.data.frame(resOrdered) 
  DEG=na.omit(DEG) 
  nrDEG=DEG
  rownames(nrDEG) <- nrDEG$row
  DEseq_DEG=nrDEG
  
  ## heatmap
  library(pheatmap)
  choose_gene=head(rownames(nrDEG),100) ## 100 maybe better
  choose_matrix= exprSet[choose_gene,]
  choose_matrix=t(scale(t(choose_matrix)))
  write.csv(choose_matrix,file = paste0(a,"_top100_DEG.csv"))

  # group_list=as.factor(c(rep("k",3),rep("aa",3)))
  # plot_color = c('orange','green')[group_list]
  pdf(file=paste0(a,"_DE_gene_heatmap.pdf"),width=60,height=90)
  par(oma=c(10,3,3,7))
  heatmapMat=choose_matrix
  heatmap.2(heatmapMat,col='greenred',trace="none",
            cexRow = 3, cexCol = 2 ,   
            srtCol = 45, offsetCol = -0.5
  )
  dev.off()
  
  # vocanol plot
  fdr = 0.05
  fold_change = 2
  DEseq_DEG$regulate=as.factor(ifelse(DEseq_DEG$pvalue<fdr&abs(DEseq_DEG$log2FoldChange) >=log2(fold_change), 
                                      ifelse(DEseq_DEG$log2FoldChange>log2(fold_change),'Up','Down'),'Normal'))
  write.csv(DEseq_DEG,file = paste0(a,"_DEseq_DEG.csv"))
  DEseq_DEG[6] <- -log10(DEseq_DEG[6])
  ggscatter(DEseq_DEG, 
            x = "log2FoldChange", 
            y = "pvalue", 
            ylab="-log10(P.value)", 
            size=0.6, 
            color = "regulate", 
            label = "row",
            #label.select = deg.data$Lable,
            label.select = DEseq_DEG$row[1:10],
            repel = TRUE,
            palette = c("#00AFBB", "#999999", "#FC4E07"))+
    geom_hline(yintercept = 1.30,linetype ="dashed")+
    geom_vline(xintercept = c(-1,1),linetype ="dashed")
  ggsave(paste0(a,"_vocanol.pdf"),width = 7.09, height =5.6,dpi = 300)
}
