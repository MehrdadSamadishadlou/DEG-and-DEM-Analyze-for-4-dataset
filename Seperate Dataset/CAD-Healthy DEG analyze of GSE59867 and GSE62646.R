setwd('/Users/mehrdad/Desktop/DataR')

library(GEOquery)
library(Biobase)
library(pheatmap)
library(ggplot2)
library(limma)
library(gplots)

superset <- readRDS('Datasets/Healthy_CAD_4dataset.Rds')

CAD46_Healthy <- superset[,superset@phenoData@data[["submission_date"]] != "Feb 27 2015"]
CAD67_Healthy <- superset[,superset@phenoData@data[["submission_date"]] != "Oct 23 2014"]


##################################  GSE62646 ##############################

############################### Quality Control ############################  

exp <- exprs(CAD46_Healthy)
colnames(exp) <- CAD46_Healthy@phenoData@data[["geo_accession"]]

pdf("Results/Seperate Datasets/GSE62646/CAD_Healthy/boxplot.pdf", width=64)
boxplot(exp)
dev.off()


pdf("Results/Seperate Datasets/GSE62646/CAD_Healthy/CorHeatmap-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(exp), labels_row = CAD46_Healthy@phenoData@data[["geo_accession"]], 
         labels_col = CAD46_Healthy@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pdf("Results/Seperate Datasets/GSE62646/CAD_Healthy/CorHeatmap-spearman.pdf", width = 25, height = 25)
pheatmap(cor(exp, method = "spearman"), labels_row = CAD46_Healthy@phenoData@data[["geo_accession"]], 
         labels_col = CAD46_Healthy@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(exp)

pdf("Results/Seperate Datasets/GSE62646/CAD_Healthy/PC.pdf")
plot(pc)
dev.off()

pcr <- data.frame(pc$r[, 1:3] , Group = CAD46_Healthy@phenoData@data[["title"]])

pdf("Results/Seperate Datasets/GSE62646/CAD_Healthy/PCA-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group)) + geom_point(size = 4, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(CAD46_Healthy@phenoData@data[["title"]])),
                     values = c("red", "green"))
dev.off()



############### Differential Expression Analysis with Limma #########################

sample_status <- factor(CAD46_Healthy@phenoData@data[["title"]])
CAD46_Healthy$description <- sample_status
Design <- model.matrix(~ description + 0, CAD46_Healthy)
colnames(Design) <- levels(sample_status)

#colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(CAD46_Healthy, Design)
cont.matrix <- makeContrasts(CAD-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "adj.P.Val", "logFC"))

tT <- tT[order(tT$ID),]
sym <- read.csv("GPL6244_GeneSymbol.csv", header = TRUE, check.names = FALSE)
tT <- cbind(sym,tT[,-1])

write.csv(tT, "Results/Seperate Datasets/GSE62646/CAD_Healthy/CAD(GSE62646)-Healthy_statistical-data.csv", row.names=F, quote = F)


CAD.up <- subset(tT, logFC>1 & adj.P.Val<0.05)

Normal.up <- subset(tT, logFC < -1 & adj.P.Val<0.05)

DEGs <- rbind(CAD.up,Normal.up)


#### Deleting Control probes

CAD.up <- CAD.up[CAD.up$address != "Control",]
Normal.up <- Normal.up[Normal.up$address != "Control",]
DEGs <- DEGs[DEGs$address != "Control",]


write.csv(CAD.up, file="Results/Seperate Datasets/GSE62646/CAD_Healthy/CAD(GSE62646)_upper.csv", quote = F, row.names = F, col.names = F)
write.csv(Normal.up, file="Results/Seperate Datasets/GSE62646/CAD_Healthy/Normal_upper(GSE62646).csv", quote = F, row.names = F, col.names = F)
write.csv(DEGs, file="Results/Seperate Datasets/GSE62646/CAD_Healthy/DEGs(GSE62646).csv", quote = F, row.names = F, col.names = F)

################################## GSE59867 ##############################

############################### Quality Control ############################  

exp <- exprs(CAD67_Healthy)
colnames(exp) <- CAD67_Healthy@phenoData@data[["geo_accession"]]

pdf("Results/Seperate Datasets/GSE59867/CAD_Healthy/boxplot.pdf", width=64)
boxplot(exp)
dev.off()


pdf("Results/Seperate Datasets/GSE59867/CAD_Healthy/CorHeatmap-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(exp), labels_row = CAD67_Healthy@phenoData@data[["geo_accession"]], 
         labels_col = CAD67_Healthy@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pdf("Results/Seperate Datasets/GSE59867/CAD_Healthy/CorHeatmap-spearman.pdf", width = 25, height = 25)
pheatmap(cor(exp, method = "spearman"), labels_row = CAD67_Healthy@phenoData@data[["geo_accession"]], 
         labels_col = CAD67_Healthy@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(exp)

pdf("Results/Seperate Datasets/GSE59867/CAD_Healthy/PC.pdf")
plot(pc)
dev.off()

pcr <- data.frame(pc$r[, 1:3] , Group = CAD67_Healthy@phenoData@data[["title"]])

pdf("Results/Seperate Datasets/GSE59867/CAD_Healthy/PCA-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group)) + geom_point(size = 4, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(CAD67_Healthy@phenoData@data[["title"]])),
                     values = c("red", "green"))
dev.off()



############### Differential Expression Analysis with Limma #########################

sample_status <- factor(CAD67_Healthy@phenoData@data[["title"]])
CAD67_Healthy$description <- sample_status
Design <- model.matrix(~ description + 0, CAD67_Healthy)
colnames(Design) <- levels(sample_status)

#colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(CAD67_Healthy, Design)
cont.matrix <- makeContrasts(CAD-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "adj.P.Val", "logFC"))

tT <- tT[order(tT$ID),]
sym <- read.csv("GPL6244_GeneSymbol.csv", header = TRUE, check.names = FALSE)
tT <- cbind(sym,tT[,-1])

write.csv(tT, "Results/Seperate Datasets/GSE59867/CAD_Healthy/CAD(GSE59867)-Healthy_statistical-data.csv", row.names=F, quote = F)


CAD.up <- subset(tT, logFC>1 & adj.P.Val<0.05)

Normal.up <- subset(tT, logFC < -1 & adj.P.Val<0.05)

DEGs <- rbind(CAD.up,Normal.up)


#### Deleting Control probes

CAD.up <- CAD.up[CAD.up$address != "Control",]
Normal.up <- Normal.up[Normal.up$address != "Control",]
DEGs <- DEGs[DEGs$address != "Control",]


write.csv(CAD.up, file="Results/Seperate Datasets/GSE59867/CAD_Healthy/CAD(GSE59867)_upper.csv", quote = F, row.names = F, col.names = F)
write.csv(Normal.up, file="Results/Seperate Datasets/GSE59867/CAD_Healthy/Normal_upper(GSE59867).csv", quote = F, row.names = F, col.names = F)
write.csv(DEGs, file="Results/Seperate Datasets/GSE59867/CAD_Healthy/DEGs(GSE59867).csv", quote = F, row.names = F, col.names = F)

