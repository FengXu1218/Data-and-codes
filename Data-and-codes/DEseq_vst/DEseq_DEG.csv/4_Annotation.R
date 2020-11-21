# enrichment analysis
rm(list = ls())
setwd("./DEseq_vst/DEseq_DEG.csv/")
library(org.Sc.sgd.db)
library(clusterProfiler)
library(pathview)
library(topGO)
library(stringr)
library(ggplot2)
files <- (Sys.glob("*.csv"))
f<-length(files)
for (i in 1:f) {
  degs <- read.csv(file = files[i])
  a <- str_split(files[i],"_D",simplify = T)[,1]
  rownames(degs) <- degs[,1]
  degs <- degs[,-1]
  DEG_list <- as.character(degs$row[degs$regulate != "Normal"])
  ego <- enrichGO(gene          = DEG_list,
                  OrgDb         = org.Sc.sgd.db,
                  keyType = "GENENAME",
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  write.table(ego, file = paste0(a,"_GO_enrichment_stat.txt"),sep="\t", row.names =F, quote = F)
  ego_results<-summary(ego)
  ego_results

  barplot(ego, showCategory=10, x = "GeneRatio")
  ggsave(paste0(a,"_ego_barplot.pdf"),width = 13, height = 13)

  dotplot(ego,showCategory=10)
  ggsave(paste0(a,"_ego_dotplot.pdf"),width = 13, height = 13)
}

