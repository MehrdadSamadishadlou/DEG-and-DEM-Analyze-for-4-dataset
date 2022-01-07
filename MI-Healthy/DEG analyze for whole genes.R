setwd('/Users/mehrdad/Desktop/DataR')

library(GEOquery)
library(Biobase)
library(pheatmap)
library(ggplot2)
library(limma)
library(gplots)


superset <- readRDS('Super4set.Rds')
exp <- exprs(superset)
colnames(exp) <- superset@phenoData@data[["geo_accession"]]

## Quality Control

pdf("Results/boxplot.pdf", width=64)
boxplot(exp)
dev.off()


pdf("Results/CorHeatmap-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(exp), labels_row = superset@phenoData@data[["geo_accession"]], 
         labels_col = superset@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pdf("Results/CorHeatmap-spearman.pdf", width = 25, height = 25)
pheatmap(cor(exp, method = "spearman"), labels_row = superset@phenoData@data[["geo_accession"]], 
         labels_col = superset@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(exp)

pdf("Results/PC.pdf")
plot(pc)
dev.off()

pcr <- data.frame(pc$r[, 1:3] , Group = superset@phenoData@data[["title"]])

pdf("Results/PCA-Samples.pdf", width = 20, height = 13)        
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
cont.matrix <- makeContrasts(MI-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "adj.P.Val", "logFC"))

tT <- tT[order(tT$ID),]
sym <- read.csv("GPL6244_GeneSymbol.csv", header = TRUE, check.names = FALSE)
tT <- cbind(sym,tT[,-1])

write.csv(tT, "Results/MI-Healthy_statistical-data.csv", row.names=F, quote = F)

MI.up <- subset(tT, logFC>1 & adj.P.Val<0.05)

Normal.up <- subset(tT, logFC < -1 & adj.P.Val<0.05)

DEGs <- rbind(MI.up,Normal.up)


#MI.up <- read.csv("Results/MI-Healthy/MI_upper.csv", header = TRUE, check.names = FALSE)
#Normal.up <- read.csv("Results//MI-Healthy/Normal_upper.csv", header = TRUE, check.names = FALSE)
#DEGs <- read.csv("Results//MI-Healthy/DEGs.csv", header = TRUE, check.names = FALSE)


##### Looking for control probes in DEGs

sum(MI.up$address == "Control") #358

sum(Normal.up$address == "Control") #509

sum(DEGs$address == "Control") #867

#### Deleting Control probes

MI.up <- MI.up[MI.up$address != "Control",]
Normal.up <- Normal.up[Normal.up$address != "Control",]
DEGs <- DEGs[DEGs$address != "Control",]

write.csv(MI.up, file="Results/MI-Healthy/MI_upper.csv", quote = F, row.names = F, col.names = F)
write.csv(Normal.up, file="Results/MI-Healthy/Normal_upper.csv", quote = F, row.names = F, col.names = F)
write.csv(DEGs, file="Results/MI-Healthy/DEGs.csv", quote = F, row.names = F, col.names = F)
