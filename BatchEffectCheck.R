setwd('/Users/mehrdad/Desktop/DataR')

library(Biobase)
library(pheatmap)
library(ggplot2)
library(gplots)
library(limma)
library(sva)


# Reading CAD  and Healthy data
CAD_Healthy <- readRDS('Datasets/Healthy_CAD_4dataset.Rds')
CAD_exp <- exprs(CAD_Healthy)
colnames(CAD_exp) <- CAD_Healthy@phenoData@data[["geo_accession"]]

Healthy_exp <- CAD_exp[,1:51]
write.csv(Healthy_exp, file="Results/BatchEffectCheck/Healthy.csv", quote = F, row.names = T, col.names = T)

Group_H = c(rep("75", 5), rep("09", 46))



CAD_exp <- CAD_exp[, 52:111]
write.csv(CAD_exp, file="Results/BatchEffectCheck/CAD.csv", quote = F, row.names = T, col.names = T)

Group_C = c(rep("46", 14), rep("67", 46))

# Reading MI data
MI_Healthy <- readRDS('Datasets/Healthy_MI_4dataset.Rds')
MI_exp <- exprs(MI_Healthy)
colnames(MI_exp) <- MI_Healthy@phenoData@data[["geo_accession"]]

MI_exp <- MI_exp[, 52:190]
write.csv(MI_exp, file="Results/BatchEffectCheck/MI.csv", quote = F, row.names = T, col.names = T)

Group_M = c(rep("46", 28), rep("67", 111))


##### Heatmap and PCA for native data

# Healthy Samples

pdf("Results/BatchEffectCheck/Healthy-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(Healthy_exp), labels_row = colnames(Healthy_exp), 
         labels_col = colnames(Healthy_exp), fontsize_row = 25, fontsize_col = 25,
         color=bluered(256), border_color = NA)
dev.off()

pdf("Results/BatchEffectCheck/Healthy-spearman.pdf", width = 25, height = 25)
pheatmap(cor(Healthy_exp, method = "spearman"), labels_row = colnames(Healthy_exp), 
         labels_col = colnames(Healthy_exp), fontsize_row = 25, fontsize_col = 25,
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(Healthy_exp)
pcr <- data.frame(pc$r[, 1:3] , Group = Group_H)

pdf("Results/BatchEffectCheck/Healthy-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group_H)) + geom_point(size = 7, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(Group_H)),
                     values = c("red", "blue"))
dev.off()


# CAD Samples

pdf("Results/BatchEffectCheck/CAD-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(CAD_exp), labels_row = colnames(CAD_exp), 
         labels_col = colnames(CAD_exp), fontsize_row = 25, fontsize_col = 25,
         color=bluered(256), border_color = NA)
dev.off()

pdf("Results/BatchEffectCheck/CAD-spearman.pdf", width = 25, height = 25)
pheatmap(cor(CAD_exp, method = "spearman"), labels_row = colnames(CAD_exp), 
         labels_col = colnames(CAD_exp), fontsize_row = 25, fontsize_col = 25,
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(CAD_exp)
pcr <- data.frame(pc$r[, 1:3] , Group = Group_C)

pdf("Results/BatchEffectCheck/CAD-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group_C)) + geom_point(size = 7, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(Group_C)),
                     values = c("red", "blue"))
dev.off()

# MI Samples

pdf("Results/BatchEffectCheck/MI-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(MI_exp), labels_row = colnames(MI_exp), 
         labels_col = colnames(MI_exp), fontsize_row = 15, fontsize_col = 15,
         color=bluered(256), border_color = NA)
dev.off()

pdf("Results/BatchEffectCheck/MI-spearman.pdf", width = 25, height = 25)
pheatmap(cor(MI_exp, method = "spearman"), labels_row = colnames(MI_exp), 
         labels_col = colnames(MI_exp), fontsize_row = 15, fontsize_col = 15,
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(MI_exp)
pcr <- data.frame(pc$r[, 1:3] , Group = Group_M)

pdf("Results/BatchEffectCheck/MI-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group_M)) + geom_point(size = 7, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(Group_M)),
                     values = c("red", "blue"))
dev.off()




### Batch effect removal by limma's removebatcheffect

Healthy_BFR <- removeBatchEffect(Healthy_exp, Group_H)
CAD_BFR <- removeBatchEffect(CAD_exp, Group_C)
MI_BFR <- removeBatchEffect(MI_exp, Group_M)

write.csv(Healthy_BFR, file="Results/BatchEffectCheck/Healthy_BFR.csv", quote = F, row.names = T, col.names = T)
write.csv(CAD_BFR, file="Results/BatchEffectCheck/CAD_BFR.csv", quote = F, row.names = T, col.names = T)
write.csv(MI_BFR, file="Results/BatchEffectCheck/MI_BFR.csv", quote = F, row.names = T, col.names = T)

# Healthy Samples

pdf("Results/BatchEffectCheck/Healthy-BFR-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(Healthy_BFR), labels_row = colnames(Healthy_BFR), 
         labels_col = colnames(Healthy_BFR), fontsize_row = 25, fontsize_col = 25,
         color=bluered(256), border_color = NA)
dev.off()

pdf("Results/BatchEffectCheck/Healthy-BFR-spearman.pdf", width = 25, height = 25)
pheatmap(cor(Healthy_BFR, method = "spearman"), labels_row = colnames(Healthy_BFR), 
         labels_col = colnames(Healthy_BFR), fontsize_row = 25, fontsize_col = 25,
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(Healthy_BFR)
pcr <- data.frame(pc$r[, 1:3] , Group = Group_H)

pdf("Results/BatchEffectCheck/Healthy-BFR-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group_H)) + geom_point(size = 7, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(Group_H)),
                     values = c("red", "blue"))
dev.off()


# CAD Samples

pdf("Results/BatchEffectCheck/CAD-BFR-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(CAD_BFR), labels_row = colnames(CAD_BFR), 
         labels_col = colnames(CAD_BFR), fontsize_row = 25, fontsize_col = 25,
         color=bluered(256), border_color = NA)
dev.off()

pdf("Results/BatchEffectCheck/CAD-BFR-spearman.pdf", width = 25, height = 25)
pheatmap(cor(CAD_BFR, method = "spearman"), labels_row = colnames(CAD_BFR), 
         labels_col = colnames(CAD_BFR), fontsize_row = 25, fontsize_col = 25,
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(CAD_BFR)
pcr <- data.frame(pc$r[, 1:3] , Group = Group_C)

pdf("Results/BatchEffectCheck/CAD_BFR-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group_C)) + geom_point(size = 7, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(Group_C)),
                     values = c("red", "blue"))
dev.off()

# MI Samples

pdf("Results/BatchEffectCheck/MI-BFR-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(MI_BFR), labels_row = colnames(MI_BFR), 
         labels_col = colnames(MI_BFR), fontsize_row = 15, fontsize_col = 15,
         color=bluered(256), border_color = NA)
dev.off()

pdf("Results/BatchEffectCheck/MI-BFR-spearman.pdf", width = 25, height = 25)
pheatmap(cor(MI_BFR, method = "spearman"), labels_row = colnames(MI_BFR), 
         labels_col = colnames(MI_BFR), fontsize_row = 15, fontsize_col = 15,
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(MI_BFR)
pcr <- data.frame(pc$r[, 1:3] , Group = Group_M)

pdf("Results/BatchEffectCheck/MI-BFR-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group_M)) + geom_point(size = 7, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(Group_M)), values = c("red", "blue"))
dev.off()


## Batch Effect Removal by sva's Combat

Healthy_Combat = ComBat(Healthy_exp, batch=Group_H)
CAD_Combat = ComBat(CAD_exp, batch=Group_C)
MI_Combat = ComBat(MI_exp, batch=Group_M)

write.csv(Healthy_Combat, file="Results/BatchEffectCheck/Healthy_Combat.csv", quote = F, row.names = T, col.names = T)
write.csv(CAD_Combat, file="Results/BatchEffectCheck/CAD_Combat.csv", quote = F, row.names = T, col.names = T)
write.csv(MI_Combat, file="Results/BatchEffectCheck/MI_Combat.csv", quote = F, row.names = T, col.names = T)


# Healthy Samples

pdf("Results/BatchEffectCheck/Healthy-Combat-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(Healthy_Combat), labels_row = colnames(Healthy_Combat), 
         labels_col = colnames(Healthy_Combat), fontsize_row = 25, fontsize_col = 25,
         color=bluered(256), border_color = NA)
dev.off()

pdf("Results/BatchEffectCheck/Healthy-Combat-spearman.pdf", width = 25, height = 25)
pheatmap(cor(Healthy_Combat, method = "spearman"), labels_row = colnames(Healthy_Combat), 
         labels_col = colnames(Healthy_Combat), fontsize_row = 25, fontsize_col = 25,
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(Healthy_Combat)
pcr <- data.frame(pc$r[, 1:3] , Group = Group_H)

pdf("Results/BatchEffectCheck/Healthy-Combat-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group_H)) + geom_point(size = 7, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(Group_H)),
                     values = c("red", "blue"))
dev.off()


# CAD Samples

pdf("Results/BatchEffectCheck/CAD-Combat-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(CAD_Combat), labels_row = colnames(CAD_Combat), 
         labels_col = colnames(CAD_Combat), fontsize_row = 25, fontsize_col = 25,
         color=bluered(256), border_color = NA)
dev.off()

pdf("Results/BatchEffectCheck/CAD-Combat-spearman.pdf", width = 25, height = 25)
pheatmap(cor(CAD_Combat, method = "spearman"), labels_row = colnames(CAD_Combat), 
         labels_col = colnames(CAD_Combat), fontsize_row = 25, fontsize_col = 25,
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(CAD_Combat)
pcr <- data.frame(pc$r[, 1:3] , Group = Group_C)

pdf("Results/BatchEffectCheck/CAD-Combat-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group_C)) + geom_point(size = 7, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(Group_C)),
                     values = c("red", "blue"))
dev.off()

# MI Samples

pdf("Results/BatchEffectCheck/MI-Combat-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(MI_Combat), labels_row = colnames(MI_Combat), 
         labels_col = colnames(MI_Combat), fontsize_row = 15, fontsize_col = 15,
         color=bluered(256), border_color = NA)
dev.off()

pdf("Results/BatchEffectCheck/MI-Combat-spearman.pdf", width = 25, height = 25)
pheatmap(cor(MI_Combat, method = "spearman"), labels_row = colnames(MI_Combat), 
         labels_col = colnames(MI_Combat), fontsize_row = 15, fontsize_col = 15,
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(MI_Combat)
pcr <- data.frame(pc$r[, 1:3] , Group = Group_M)

pdf("Results/BatchEffectCheck/MI-Combat-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group_M)) + geom_point(size = 7, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(Group_M)), values = c("red", "blue"))
dev.off()
