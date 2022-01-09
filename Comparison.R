setwd('/Users/mehrdad/Desktop/DataR')

library(ggVennDiagram)
library(pheatmap)
library(ggplot2)
library(gplots)
library(ComplexHeatmap)
############################### Venn Diagram for MI_Healthy and CAD-Healthy DEGs ############################

MI_Healthy <- read.csv("Results/MI-Healthy/DEGs.csv", header = TRUE, check.names = FALSE)
CAD_Healthy <- read.csv("Results/CAD-Healthy/DEGs_CAD.csv", header = TRUE, check.names = FALSE)

probes <- list(MI_Healthy$ID, CAD_Healthy$ID)

pdf("Results/Comparison/DEGs in MI-Healthy and CAD-Healthy.pdf", width = 8, height = 6)
ggVennDiagram(probes, category.names = c("STEMI-Healthy", "stable CAD-Healthy")) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  labs(title="Differentially Expressed Genes") +
  theme(legend.position = "none", plot.title = element_text(color="black", size=18, face="bold.italic", hjust = 0.5))
dev.off()

############################### Venn Diagram for MI_Healthy Meta with single datasets ############################


MI_Healthy_GSE62646 <- read.csv("Results/Seperate Datasets/GSE62646/DEGs(GSE62646).csv", header = TRUE, check.names = FALSE)
MI_Healthy_GSE59867 <- read.csv("Results/Seperate Datasets/GSE59867/DEGs(GSE59867).csv", header = TRUE, check.names = FALSE)



probes <- list(MI_Healthy$ID, MI_Healthy_GSE59867$ID, MI_Healthy_GSE62646$ID)

pdf("Results/Comparison/DEGs in MI-Healthy meta-single.pdf", width = 8, height = 6)
ggVennDiagram(probes, category.names = c("Meta Analyze", "GSE59867", "GSE62646")) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  labs(title="DEGs in MI-Healthy meta analyze comparing with single datasets") +
  theme(legend.position = "none", plot.title = element_text(color="black", size=18, face="bold.italic", hjust = 0.5))
dev.off()


############################### Venn Diagram for CAD_Healthy Meta with single datasets ############################

CAD_Healthy_GSE62646 <- read.csv("Results/Seperate Datasets/GSE62646/CAD_Healthy/DEGs(GSE62646).csv", header = TRUE, check.names = FALSE)
CAD_Healthy_GSE59867 <- read.csv("Results/Seperate Datasets/GSE59867/CAD_Healthy/DEGs(GSE59867).csv", header = TRUE, check.names = FALSE)


probes <- list(CAD_Healthy$ID, CAD_Healthy_GSE59867$ID, CAD_Healthy_GSE62646$ID)

pdf("Results/Comparison/DEGs in CAD-Healthy meta-single.pdf", width = 9, height = 6)
ggVennDiagram(probes, category.names = c("Meta Analyze", "GSE59867", "GSE62646")) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  labs(title="DEGs in CAD-Healthy meta analyze comparing with single datasets") +
  theme(legend.position = "none", plot.title = element_text(color="black", size=18, face="bold.italic", hjust = 0.5))
dev.off()


################# Venn Diagram for MI_Healthy and CAD-Healthy differentially expressed microRNAs ###############

MI_Healthy_mir <- read.csv("Results/MI-Healthy/MIRS in MI DEGs.csv", header = TRUE, check.names = FALSE)
CAD_Healthy_mir <- read.csv("Results/CAD-Healthy/MIRS in CAD DEGs.csv", header = TRUE, check.names = FALSE)

probes <- list(MI_Healthy$ID, CAD_Healthy$ID)

pdf("Results/Comparison/DEMs in MI-Healthy and CAD-Healthy.pdf", width = 8, height = 6)
ggVennDiagram(probes, category.names = c("STEMI-Healthy", "stable CAD-Healthy")) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  labs(title="Differentially Expressed microRNAs") +
  theme(legend.position = "none", plot.title = element_text(color="black", size=18, face="bold.italic", hjust = 0.5))
dev.off()



########### Heatmap ############

superset <- readRDS('Datasets/Healthy_MI_4dataset.Rds')
exp <- exprs(superset)
colnames(exp) <- superset@phenoData@data[["geo_accession"]]

MI_Healthy <- MI_Healthy[order(abs(MI_Healthy$logFC)),]

top_50_logFC_ID <- as.character(MI_Healthy[3319:3368,1])
top_50_logFC_exp <- exp[row.names(exp) %in% top_50_logFC_ID,]

pdf("Results/Comparison/Heatmap.pdf", width = 12, height = 20)
Heatmap(t(top_50_logFC_exp))
dev.off()
