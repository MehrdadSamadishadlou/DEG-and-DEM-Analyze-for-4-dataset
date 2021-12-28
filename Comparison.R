setwd('/Users/mehrdad/Desktop/DataR')

library(ggVennDiagram)
library(pheatmap)
library(ggplot2)
library(gplots)

MI_Healthy <- read.csv("Results/MI-Healthy/DEGs.csv", header = TRUE, check.names = FALSE)
CAD_Healthy <- read.csv("Results/CAD-Healthy/DEGs_CAD_Healthy.csv", header = TRUE, check.names = FALSE)
MI_Healthy_GSE62646 <- read.csv("Results/Seperate Datasets/GSE62646/DEGs(GSE62646).csv", header = TRUE, check.names = FALSE)
MI_Healthy_GSE59867 <- read.csv("Results/Seperate Datasets/GSE59867/DEGs(GSE59867).csv", header = TRUE, check.names = FALSE)



probes <- list(MI_Healthy$ID, CAD_Healthy$ID)

pdf("Results/Comparison/DEGs in MI-Healthy and CAD-Healthy.pdf", width = 8, height = 6)
ggVennDiagram(probes, category.names = c("STEMI", "stable CAD")) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
dev.off()


probes <- list(MI_Healthy$ID, MI_Healthy_GSE62646$ID)

pdf("Results/Comparison/DEGs in MI-Healthy Tow Datasets and GSE62646.pdf", width = 8, height = 6)
ggVennDiagram(probes, category.names = c("STEMI (Two Datasets)", "STEMI GSE62646 Dataset")) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
dev.off()


probes <- list(MI_Healthy$ID, MI_Healthy_GSE59867$ID)

pdf("Results/Comparison/DEGs in MI-Healthy Tow Datasets and GSE59867.pdf", width = 8, height = 6)
ggVennDiagram(probes, category.names = c("STEMI (Two Datasets)", "STEMI GSE59867 Dataset")) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
dev.off()


########### Heatmap ############

superset <- readRDS('Datasets/Healthy_MI_4dataset.Rds')
exp <- exprs(superset)
colnames(exp) <- superset@phenoData@data[["geo_accession"]]

MI_Healthy <- MI_Healthy[order(abs(MI_Healthy$logFC)),]

top_50_logFC_ID <- as.character(MI_Healthy[4186:4235,1])
top_50_logFC_exp <- exp[row.names(exp) %in% top_50_logFC_ID,]
filler <- as.data.frame(matrix(rep(0, 190*140),nrow=140))
colnames(filler) <- colnames(top_50_logFC_exp)
top_50_logFC_exp <- rbind(top_50_logFC_exp, filler)

top_50_logFC_exp_trans <- t(top_50_logFC_exp)
filler2 <- as.data.frame(matrix(rep(0, 190*140),nrow=190))
top_50_logFC_exp_trans <- cbind(top_50_logFC_exp_trans, filler2)


pdf("Results/gene-sample_top-50-logFC.pdf", width = 25, height = 25)
pheatmap(b, labels_row = row.names(b), 
         labels_col = colnames(b),
         color=bluered(256), border_color = NA)
dev.off()


exp <- exprs(superset)
colnames(exp) <- superset@phenoData@data[["geo_accession"]]

MI_Healthy <- MI_Healthy[order(abs(MI_Healthy$adj.P.Val)),]

top_50_logFC_ID <- as.character(MI_Healthy[1:50,1])
top_50_logFC_exp <- exp[row.names(exp) %in% top_50_logFC_ID,]
filler <- as.data.frame(matrix(rep(0, 190*140),nrow=140))
colnames(filler) <- colnames(top_50_logFC_exp)
top_50_logFC_exp <- rbind(top_50_logFC_exp, filler)

top_50_logFC_exp_trans <- t(top_50_logFC_exp)
filler2 <- as.data.frame(matrix(rep(0, 190*140),nrow=190))
top_50_logFC_exp_trans <- cbind(top_50_logFC_exp_trans, filler2)

correlation <- cor(top_50_logFC_exp, top_50_logFC_exp_trans)
correlation <- correlation[,1:50]
pdf("Results/gene-sample_top-50-logFC.pdf", width = 25, height = 25)
pheatmap(correlation, labels_row = row.names(b), 
         labels_col = colnames(b),
         color=bluered(256), border_color = NA)
dev.off()
