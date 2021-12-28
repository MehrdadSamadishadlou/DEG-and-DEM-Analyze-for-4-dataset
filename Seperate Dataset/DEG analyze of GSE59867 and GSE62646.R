setwd('/Users/mehrdad/Desktop/DataR')

library(GEOquery)
library(Biobase)
library(pheatmap)
library(ggplot2)
library(limma)
library(gplots)

superset <- readRDS('Datasets/Healthy_MI_4dataset.Rds')

MI46_Healthy <- superset[,superset@phenoData@data[["submission_date"]] != "Jul 29 2014"]
MI67_Healthy <- superset[,superset@phenoData@data[["submission_date"]] != "Oct 23 2014"]


##################################  GSE62646 ##############################

############################### Quality Control ############################  

exp <- exprs(MI46_Healthy)
colnames(exp) <- MI46_Healthy@phenoData@data[["geo_accession"]]

pdf("Results/Seperate Datasets/GSE62646/boxplot.pdf", width=64)
boxplot(exp)
dev.off()


pdf("Results/Seperate Datasets/GSE62646/CorHeatmap-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(exp), labels_row = superset@phenoData@data[["geo_accession"]], 
         labels_col = superset@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pdf("Results/Seperate Datasets/GSE62646/CorHeatmap-spearman.pdf", width = 25, height = 25)
pheatmap(cor(exp, method = "spearman"), labels_row = superset@phenoData@data[["geo_accession"]], 
         labels_col = superset@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(exp)

pdf("Results/Seperate Datasets/GSE62646/PC.pdf")
plot(pc)
dev.off()

pcr <- data.frame(pc$r[, 1:3] , Group = MI46_Healthy@phenoData@data[["title"]])

pdf("Results/Seperate Datasets/GSE62646/PCA-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group)) + geom_point(size = 4, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(MI46_Healthy@phenoData@data[["title"]])),
                     values = c("red", "green"))
dev.off()



############### Differential Expression Analysis with Limma #########################

sample_status <- factor(MI46_Healthy@phenoData@data[["title"]])
MI46_Healthy$description <- sample_status
Design <- model.matrix(~ description + 0, MI46_Healthy)
colnames(Design) <- levels(sample_status)

#colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(MI46_Healthy, Design)
cont.matrix <- makeContrasts(MI-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "adj.P.Val", "logFC"))

tT <- tT[order(tT$ID),]
sym <- read.csv("GPL6244_GeneSymbol.csv", header = TRUE, check.names = FALSE)
tT <- cbind(sym,tT[,-1])

write.csv(tT, "Results/Seperate Datasets/GSE62646/MI(GSE62646)-Healthy_statistical-data.csv", row.names=F, quote = F)

MI.up <- subset(tT, logFC>1 & adj.P.Val<0.05)
write.csv(MI.up, file="Results/Seperate Datasets/GSE62646/MI(GSE62646)_upper.csv", quote = F, row.names = F, col.names = F)

Normal.up <- subset(tT, logFC < -1 & adj.P.Val<0.05)
write.csv(Normal.up, file="Results/Seperate Datasets/GSE62646/Normal_upper(GSE62646).csv", quote = F, row.names = F, col.names = F)

DEGs <- rbind(MI.up,Normal.up)
write.csv(DEGs, file="Results/Seperate Datasets/GSE62646/DEGs(GSE62646).csv", quote = F, row.names = F, col.names = F)



################################## GSE59867 ##############################

############################### Quality Control ############################  

exp <- exprs(MI67_Healthy)
colnames(exp) <- MI67_Healthy@phenoData@data[["geo_accession"]]

pdf("Results/Seperate Datasets/GSE59867/boxplot.pdf", width=64)
boxplot(exp)
dev.off()


pdf("Results/Seperate Datasets/GSE59867/CorHeatmap-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(exp), labels_row = MI67_Healthy@phenoData@data[["geo_accession"]], 
         labels_col = MI67_Healthy@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pdf("Results/Seperate Datasets/GSE59867/CorHeatmap-spearman.pdf", width = 25, height = 25)
pheatmap(cor(exp, method = "spearman"), labels_row = MI67_Healthy@phenoData@data[["geo_accession"]], 
         labels_col = MI67_Healthy@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(exp)

pdf("Results/Seperate Datasets/GSE59867/PC.pdf")
plot(pc)
dev.off()

pcr <- data.frame(pc$r[, 1:3] , Group = MI67_Healthy@phenoData@data[["title"]])

pdf("Results/Seperate Datasets/GSE59867/PCA-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group)) + geom_point(size = 4, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(MI67_Healthy@phenoData@data[["title"]])),
                     values = c("red", "green"))
dev.off()



############### Differential Expression Analysis with Limma #########################

sample_status <- factor(MI67_Healthy@phenoData@data[["title"]])
MI67_Healthy$description <- sample_status
Design <- model.matrix(~ description + 0, MI67_Healthy)
colnames(Design) <- levels(sample_status)

#colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(MI67_Healthy, Design)
cont.matrix <- makeContrasts(MI-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "adj.P.Val", "logFC"))

tT <- tT[order(tT$ID),]
sym <- read.csv("GPL6244_GeneSymbol.csv", header = TRUE, check.names = FALSE)
tT <- cbind(sym,tT[,-1])

write.csv(tT, "Results/Seperate Datasets/GSE59867/MI(GSE59867)-Healthy_statistical-data.csv", row.names=F, quote = F)

MI.up <- subset(tT, logFC>1 & adj.P.Val<0.05)
write.csv(MI.up, file="Results/Seperate Datasets/GSE59867/MI(GSE59867)_upper.csv", quote = F, row.names = F, col.names = F)

Normal.up <- subset(tT, logFC < -1 & adj.P.Val<0.05)
write.csv(Normal.up, file="Results/Seperate Datasets/GSE59867/Normal_upper(GSE59867).csv", quote = F, row.names = F, col.names = F)

DEGs <- rbind(MI.up,Normal.up)
write.csv(DEGs, file="Results/Seperate Datasets/GSE59867/DEGs(GSE59867).csv", quote = F, row.names = F, col.names = F)
