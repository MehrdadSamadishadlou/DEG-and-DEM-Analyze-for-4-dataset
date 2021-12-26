setwd('/Users/mehrdad/Desktop/DataR')

library(GEOquery)
library(Biobase)
library(pheatmap)
library(ggplot2)
library(limma)
library(gplots)


mirset <- readRDS('mir4set.Rds')
exp <- exprs(mirset)
colnames(exp) <- mirset@phenoData@data[["geo_accession"]]

## Quality Control

pdf("Results/just microRNA/boxplot.pdf", width=64)
boxplot(exp)
dev.off()


pdf("Results/just microRNA/CorHeatmap-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(exp), labels_row = mirset@phenoData@data[["geo_accession"]], 
         labels_col = mirset@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pdf("Results/just microRNA/CorHeatmap-spearman.pdf", width = 25, height = 25)
pheatmap(cor(exp, method = "spearman"), labels_row = mirset@phenoData@data[["geo_accession"]], 
         labels_col = mirset@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(exp)

pdf("Results/just microRNA/PC.pdf")
plot(pc)
dev.off()

pcr <- data.frame(pc$r[, 1:3] , Group = mirset@phenoData@data[["title"]])

pdf("Results/just microRNA/PCA-Samples.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , color=Group)) + geom_point(size = 4, alpha = 0.5) +
  ggtitle("Samples PCA") +
  scale_color_manual(breaks = levels(as.factor(mirset@phenoData@data[["title"]])),
                     values = c("red", "green"))
dev.off()



############### Differential Expression Analysis with Limma #########################

sample_status <- factor(mirset@phenoData@data[["title"]])
mirset$description <- sample_status
Design <- model.matrix(~ description + 0, mirset)
colnames(Design) <- levels(sample_status)

#colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(mirset, Design)
cont.matrix <- makeContrasts(MI-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "adj.P.Val", "logFC"))

tT <- tT[order(tT$ID),]
sym <- read.csv("Results/just microRNA/mir_probes.csv", header = TRUE, check.names = FALSE)
sym <- sym[order(sym$ID),]
tT <- cbind(sym,tT[,-1])

write.csv(tT, "Results/just microRNA/MI-Healthy_statistical-data.csv", row.names=F, quote = F)

MI.up <- subset(tT, logFC>1 & adj.P.Val<0.05)
write.csv(MI.up, file="Results/just microRNA/MI_upper.csv", quote = F, row.names = F, col.names = F)

Normal.up <- subset(tT, logFC < -1 & adj.P.Val<0.05)
write.csv(Normal.up, file="Results/just microRNA/Normal_upper.csv", quote = F, row.names = F, col.names = F)

DEGs <- rbind(MI.up,Normal.up)
write.csv(DEGs, file="Results/just microRNA/DEGs.csv", quote = F, row.names = F, col.names = F)

#MI.up <- read.csv("Results/just microRNA/MI_upper.csv", header = TRUE, check.names = FALSE)
#Normal.up <- read.csv("Results/just microRNA/Normal_upper.csv", header = TRUE, check.names = FALSE)
#DEGs <- read.csv("Results/just microRNA/DEGs.csv", header = TRUE, check.names = FALSE)


##### Looking for control probes in DEGs

sum(MI.up$address == "Control") #0

sum(Normal.up$address == "Control") #20

sum(DEGs$address == "Control") #20

