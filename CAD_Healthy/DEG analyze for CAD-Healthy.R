setwd('/Users/mehrdad/Desktop/DataR')

library(GEOquery)
library(Biobase)
library(pheatmap)
library(ggplot2)
library(limma)
library(gplots)
library(stringr)



superset <- readRDS('Datasets/Healthy_CAD_4dataset.Rds')
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

Normal.up <- subset(tT, logFC < -1 & adj.P.Val<0.05)

DEGs <- rbind(MI.up,Normal.up)


#CAD.up <- read.csv("Results//CAD-Healthy/CAD_upper.csv", header = TRUE, check.names = FALSE)
#Normal.up <- read.csv("Results//CAD-Healthy/Normal_upperthan_CAD.csv", header = TRUE, check.names = FALSE)
#DEGs <- read.csv("Results//CAD-Healthy/DEGs_CAD.csv", header = TRUE, check.names = FALSE)


##### Looking for control probes in DEGs

sum(CAD.up$address == "Control") #613

sum(Normal.up$address == "Control") #459

sum(DEGs$address == "Control") #1072

#### Deleting Control probes

CAD.up <- CAD.up[CAD.up$address != "Control",]
Normal.up <- Normal.up[Normal.up$address != "Control",]
DEGs <- DEGs[DEGs$address != "Control",]


write.csv(MI.up, file="Results/CAD-Healthy/CAD_upper.csv", quote = F, row.names = F, col.names = F)
write.csv(Normal.up, file="Results/CAD-Healthy/Normal_upper_CAD.csv", quote = F, row.names = F, col.names = F)
write.csv(DEGs, file="Results/CAD-Healthy/DEGs_CAD.csv", quote = F, row.names = F, col.names = F)

#### Extracting miRs ######

data1 <- DEGs[str_detect(DEGs$`By Chr`, "MIR"), ]
data2 <- DEGs[str_detect(DEGs$`By ID`, "MIR"), ]
data3 <- DEGs[str_detect(DEGs$`By Platform`, "MIR"), ]
mirs <- rbind(data1, data2, data3)
mirs <- mirs[!duplicated(mirs$ID),]
mirs <- mirs[-1,]

write.csv(mirs, file="Results/CAD-Healthy/MIRS in CAD DEGs.csv", quote = F, row.names = F, col.names = F)


############## Extracting CAD samples expression profile for DEMs in MI-Healthy ##############

cadexp <- exp[,52:111]
DEGs <- read.csv("Results/MI-Healthy/DEGs.csv", header = TRUE, check.names = FALSE)
mirs <- DEGs[str_detect(DEGs$`By Platform`, "MIR"), ]
mirs <- mirs[!duplicated(mirs$ID),]
mirs <- mirs[-1,]

cadmirexp <- cadexp[rownames(cadexp) %in% mirs$ID,]
cadmirexp <- t(cadmirexp)
colnames(cadmirexp) <- names
cadmirexp <- cadmirexp[,c(-6, -8, -9)]
write.csv(cadmirexp, file = "Results/CAD-Healthy/MI-Healthy DEMs Expression of CAD samples.csv", quote = F, row.names = T, col.names = F)
