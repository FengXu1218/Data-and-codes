# GO enrichment
rm(list = ls())
library(DOSE)
library(GOSemSim)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(pathview)
library(topGO)
library(stringr)
setwd("./DEseq_vst/DEseq_DEG.csv")
###############################################################################3
# Codes that do not need to be run
aaa <- read.csv(file = "dataTmp(1).csv",header = F,fill = T, row.names=1)

GOterms_genes <- read.csv(file = "GOterms_genes.csv",header = F,fill = T)
a <- as.matrix(GOterms_genes)
rownames(a) <- a[,1]
for (i in 1:nrow(a)) {
  b <- length(a[i,a[i,] %in% genelist])
  if(b == 0)
    a[i,1:ncol(a)] <- NA
  else
    a[i,1:b] <- a[i,a[i,] %in% genelist]
  a[i,(length(a[i,a[i,] %in% genelist])+2):ncol(a)] <- NA
}
#
get_GO_data <- function(OrgDb, ont, keytype) {
  GO_Env <- get_GO_Env()
  use_cached <- FALSE
  
  if (exists("organism", envir=GO_Env, inherits=FALSE) &&
      exists("keytype", envir=GO_Env, inherits=FALSE)) {
    
    org <- get("organism", envir=GO_Env)
    kt <- get("keytype", envir=GO_Env)
    
    if (org == DOSE:::get_organism(OrgDb) &&
        keytype == kt &&
        exists("goAnno", envir=GO_Env, inherits=FALSE)) {
      ## https://github.com/GuangchuangYu/clusterProfiler/issues/182
      ## && exists("GO2TERM", envir=GO_Env, inherits=FALSE)){
      
      use_cached <- TRUE
    }
  }
  
  if (use_cached) {
    goAnno <- get("goAnno", envir=GO_Env)
  } else {
    OrgDb <- GOSemSim:::load_OrgDb(OrgDb)
    kt <- keytypes(OrgDb)
    if (! keytype %in% kt) {
      stop("keytype is not supported...")
    }
    
    kk <- keys(OrgDb, keytype=keytype)
    goAnno <- suppressMessages(
      select(OrgDb, keys=kk, keytype=keytype,
             columns=c("GOALL", "ONTOLOGYALL")))
    
    goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])
    
    assign("goAnno", goAnno, envir=GO_Env)
    assign("keytype", keytype, envir=GO_Env)
    assign("organism", DOSE:::get_organism(OrgDb), envir=GO_Env)
  }
  
  if (ont == "ALL") {
    GO2GENE <- unique(goAnno[, c(2,1)])
  } else {
    GO2GENE <- unique(goAnno$ONTOLOGYALL == ont, c(2,1))
  }
  
  GO_DATA <- DOSE:::build_Anno(GO2GENE, get_GO2TERM_table())
  
  goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
  goOnt <- goOnt.df[,2]
  names(goOnt) <- goOnt.df[,1]
  assign("GO2ONT", goOnt, envir=GO_DATA)
  return(GO_DATA)
}

get_GO_Env <- function () {
  if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
  }
  get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}

get_GO2TERM_table <- function() {
  GOTERM.df <- get_GOTERM()
  GOTERM.df[, c("go_id", "Term")] %>% unique
}

get_GOTERM <- function() {
  pos <- 1
  envir <- as.environment(pos)
  if (!exists(".GOTERM_Env", envir=envir)) {
    assign(".GOTERM_Env", new.env(), envir)
  }
  GOTERM_Env <- get(".GOTERM_Env", envir = envir)
  if (exists("GOTERM.df", envir = GOTERM_Env)) {
    GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
  } else {
    GOTERM.df <- toTable(GOTERM)
    assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
  }
  return(GOTERM.df)
}
GO_DATA <- get_GO_data("org.Sc.sgd.db", "ALL", "GENENAME") 
######################################################################################
# ɾ??KM??û?еĻ?????
go <- goAnno[goAnno[,1] %in% genelist,][,c(2,3,1)]
# GOterms_list1 <- GO_DATA[["PATHID2EXTID"]]
# load("/Users/xufeng/Desktop/others/Kluyveromyces/countdata2/DEseq_vst/k_aa_vst.Rdata")
# genelist <- row.names(k_aa_vst)
# genelist <- matrix(nrow = 1,genelist)
# for (i in 1:5191) {
#  GOterms_list[[i]] <- GOterms_list[[i]][GOterms_list[[i]] %in% genelist]
# }

#######################################################################################
setwd("./DEseq_vst/DEseq_DEG.csv/")
library(org.Sc.sgd.db)
library(clusterProfiler)
library(pathview)
library(topGO)
library(stringr)
library(ggplot2)
library(tidyverse)
load("./DEseq_vst/k_aa_vst.Rdata")
genelist <- row.names(k_aa_vst)
# goAnno <- goAnno[,-3]
# goAnno <- goAnno[!duplicated(goAnno),]
# go <- goAnno[goAnno[,1] %in% genelist,][,c(2,3,1)]
# go <- go[order(go$GOALL),]
# rownames(go) <- 1:348343
# go_term <- go2term(go[,"GOALL"])
# rownames(go_term) <- go_term[,1]
# for (i in 1:348343) {
#   go[i,"Term"] <- go_term[as.character(go[,"GOALL"][i]),"Term"]
# }
# go_term2gene <- data.frame(go$GOALL,go$GENENAME)
# go_term2name <- data.frame(go$GOALL,go$Term)
# names(go_term2gene) <- c("go_term","gene")
# names(go_term2name) <- c("go_term","name")
# save(go_term2gene,go_term2name,go,file = "KM_GO.Rdata")
load("KM_GO.Rdata")
files <- (Sys.glob("k_aa_*.csv"))
f<-length(files)
for (i in 1:f) {
  degs <- read.csv(file = files[i],header = T)
  a <- str_split(files[i],"__",simplify = T)[,1]
  DEG_list <- as.character(degs$row[degs$regulate != "Normal"])
    ego = enricher(gene = DEG_list,
                   TERM2GENE = go_term2gene,
                   TERM2NAME = go_term2name,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
#  write.csv(ego, file = paste0("KM_GO2/",paste0(a,"_GO_enrichment_stat.csv")))
  #  ego_results<-summary(ego)
  #  ego_results
    
  #  barplot(ego, showCategory=10, x = "GeneRatio")
  #  ggsave(paste0(a,"_ego_barplot.pdf"))

  #  dotplot(ego,showCategory=10)
  #  ggsave(paste0(a,"_ego_dotplot.pdf"))
}
go_ID <- read.csv(file = "/Users/xufeng/Desktop/GO_k_aaORbb.csv",header = F)
a <- ego@result
a <- a[a$ID %in% go_ID[,1],]
Ontology <- go2ont(go_ID[,1])
a[,c(5,6,7)] <- 0.000001
a[c("GO:0044430","GO:0005856","GO:0008509","GO:0005342", "GO:0005777", 
    "GO:0005782", "GO:0015238", "GO:0022857", "GO:0031298", "GO:0031907", 
    "GO:0042579","GO:0140097"),c(5,6,7)] <- 1
ego@result <- a
ego@ontology <- "BP"
ego@organism <- "KM"

pdf(file = "k_aa_ORbb_BP.pdf")
plotGOgraph(ego,firstSigNodes = 48,
            useInfo = "all",
            sigForAll = TRUE,
            useFullNames = TRUE,)
dev.off()
pdf(file = "k_aa_ORbb_goplot_BP.pdf")
goplot(ego)
dev.off()
