setwd('/Users/mehrdad/Desktop/DataR')

library(GEOquery)
library(Biobase)
library(pheatmap)
library(ggplot2)
library(limma)
library(gplots)


superset <- readRDS('Healthy_CAD_4dataset.Rds')
exp <- exprs(superset)
colnames(exp) <- superset@phenoData@data[["geo_accession"]]

## Quality Control

pdf("Results/CAD-Healthy/boxplot.pdf", width=64)
boxplot(exp)
dev.off()


pdf("Results/CAD-Healthy/CorHeatmap-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(exp), labels_row = superset@phenoData@data[["geo_accession"]], 
         labels_col = superset@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pdf("Results/CAD-Healthy/CorHeatmap-spearman.pdf", width = 25, height = 25)
pheatmap(cor(exp, method = "spearman"), labels_row = superset@phenoData@data[["geo_accession"]], 
         labels_col = superset@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(exp)

pdf("Results/CAD-Healthy/PC.pdf")
plot(pc)
dev.off()

pcr <- data.frame(pc$r[, 1:3] , Group = superset@phenoData@data[["title"]])

pdf("Results/CAD-Healthy/PCA-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group)) + geom_point(size = 4, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(superset@phenoData@data[["title"]])),
                     values = c("red", "green"))
dev.off()



############### Differential Expression Analysis with Limma #########################

sample_status <- factor(superset@phenoData@data[["title"]])
superset$description <- sample_status
Design <- model.matrix(~ description + 0, superset)
colnames(Design) <- levels(sample_status)

#colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(superset, Design)
cont.matrix <- makeContrasts(CAD-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "adj.P.Val", "logFC"))

tT <- tT[order(tT$ID),]
sym <- read.csv("GPL6244_GeneSymbol.csv", header = TRUE, check.names = FALSE)
tT <- cbind(sym,tT[,-1])

write.csv(tT, "Results/CAD-Healthy/CAD-Healthy_statistical-data.csv", row.names=F, quote = F)

CAD.up <- subset(tT, logFC>1 & adj.P.Val<0.05)
write.csv(MI.up, file="Results/CAD-Healthy/CAD_upper.csv", quote = F, row.names = F, col.names = F)

Normal.up <- subset(tT, logFC < -1 & adj.P.Val<0.05)
write.csv(Normal.up, file="Results/CAD-Healthy/Normal_upper.csv", quote = F, row.names = F, col.names = F)

DEGs <- rbind(MI.up,Normal.up)
write.csv(DEGs, file="Results/CAD-Healthy/DEGs.csv", quote = F, row.names = F, col.names = F)

#MI.up <- read.csv("Results/MI_upper.csv", header = TRUE, check.names = FALSE)
#Normal.up <- read.csv("Results/Normal_upper.csv", header = TRUE, check.names = FALSE)
#DEGs <- read.csv("Results/DEGs.csv", header = TRUE, check.names = FALSE)


##### Looking for control probes in DEGs

sum(CAD.up$address == "Control") #613

sum(Normal.up$address == "Control") #459

sum(DEGs$address == "Control") #1072

