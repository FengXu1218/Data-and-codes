library(ggpubr)
library(ggplot2)
setwd("/Users/xufeng/Desktop/k_f_Sce/")
coeff_km <- read.csv(file = "coeff_km.csv")
coeff_km <- coeff_km[,1:3]
group_list <- c("K","K","K","aa","aa","aa","bb","bb","bb","ab","ab","ab","f","f","f",
                rep("SC",39))
coeff_km$group <- group_list
means_k <- c(mean(coeff_km$PC1[1:3]),mean(coeff_km$PC2[1:3]))
means_aa_bb <- c(mean(coeff_km$PC1[4:9]),mean(coeff_km$PC2[4:9]))
means_ab <- c(mean(coeff_km$PC1[10:12]),mean(coeff_km$PC2[10:12]))
means_f <- c(mean(coeff_km$PC1[13:15]),mean(coeff_km$PC2[13:15]))
means_SC <- c(mean(coeff_km$PC1[16:54]),mean(coeff_km$PC2[16:54]))
coeff_km$X <- as.character(coeff_km$X)

means <- as.data.frame(t(data.frame(means_k,means_aa_bb,means_ab,means_f,means_SC)))
colnames(means) = c("PC1","PC2")
rownames(means) = c("means_k","means_aa_bb","means_ab","means_f","means_SC")
means$group <- c("K","aa","ab","f","SC")

ggplot() +
  geom_point(data = coeff_km,aes(x = PC1,y = PC2,color = group,shape = group),size = 3) +
  scale_shape_manual(values = c(15,16,17,18,13,8)) +
  geom_point(data = means,aes(x = PC1,y = PC2,color = group),shape = 3,size = 3)
ggsave()


