setwd('/Users/mehrdad/Desktop/DataR')

library(GEOquery)
library(Biobase)
library(pheatmap)
library(ggplot2)
library(limma)
library(gplots)


#Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 20)

GSE59867 <- getGEO("GSE59867", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Datasets/")
geneset67 <- GSE59867[[1]]

sub67MI <- geneset67[,geneset67@phenoData@data[["samples collection:ch1"]] == "on the 1st day of MI (admission)"]
sub67MI@phenoData@data[["geo_accession"]] <- paste0(rep("67.", 111), sub67MI@phenoData@data[["geo_accession"]], rep(".M", 111))

sub67C <- geneset67[,geneset67@phenoData@data[["samples collection:ch1"]] == "N/A"]
sub67C@phenoData@data[["geo_accession"]] <- paste0(rep("67.", 46), sub67C@phenoData@data[["geo_accession"]], rep(".C", 46))

GSE62646 <- getGEO("GSE62646", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Datasets/")
geneset46 <- GSE62646[[1]]

sub46MI <- geneset46[,geneset46@phenoData@data[["collection of blood samples:ch1"]] == "admission (on the 1st day of MI)"]
sub46MI@phenoData@data[["geo_accession"]] <- paste0(rep("46.", 28), sub46MI@phenoData@data[["geo_accession"]], rep(".M", 28))

sub46C <- geneset46[,geneset46@phenoData@data[["cardiovascular disease state:ch1"]] == "CAD"]
sub46C@phenoData@data[["geo_accession"]] <- paste0(rep("46.", 14), sub46C@phenoData@data[["geo_accession"]], rep(".C", 14))

superset <- combine(sub46C, sub67C, sub46MI, sub67MI)
sample_status <- c(rep('CAD', 60), rep('MI', 139))
superset@phenoData@data[["title"]] <- sample_status


saveRDS(superset, 'CAD_Healthy_2_Dataset.Rds')


exp <- exprs(superset)
colnames(exp) <- superset@phenoData@data[["geo_accession"]]


############################# Quality Control ############################

pdf("Results/MI-CAD/boxplot.pdf", width=64)
boxplot(exp)
dev.off()


pdf("Results/MI-CAD/CorHeatmap-euclidean.pdf", width = 25, height = 25)
pheatmap(cor(exp), labels_row = superset@phenoData@data[["geo_accession"]], 
         labels_col = superset@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pdf("Results/MI-CAD/CorHeatmap-spearman.pdf", width = 25, height = 25)
pheatmap(cor(exp, method = "spearman"), labels_row = superset@phenoData@data[["geo_accession"]], 
         labels_col = superset@phenoData@data[["geo_accession"]],
         color=bluered(256), border_color = NA)
dev.off()


pc <- prcomp(exp)

pdf("Results/MI-CAD/PC.pdf")
plot(pc)
dev.off()

pcr <- data.frame(pc$r[, 1:3] , Group = superset@phenoData@data[["title"]])

pdf("Results/MI-CAD/PCA-Samples.pdf", width = 20, height = 13)        
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
cont.matrix <- makeContrasts(MI-CAD, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "adj.P.Val", "logFC"))

tT <- tT[order(tT$ID),]
sym <- read.csv("GPL6244_GeneSymbol.csv", header = TRUE, check.names = FALSE)
tT <- cbind(sym,tT[,-1])

write.csv(tT, "Results/MI-CAD/MI_CAD_statistical-data.csv", row.names=F, quote = F)

MI.up <- subset(tT, logFC>1 & adj.P.Val<0.05)
write.csv(MI.up, file="Results/MI-CAD/MI_upper.csv", quote = F, row.names = F, col.names = F)

CAD.up <- subset(tT, logFC < -1 & adj.P.Val<0.05)
write.csv(CAD.up, file="Results/MI-CAD/CAD_upper.csv", quote = F, row.names = F, col.names = F)

DEGs <- rbind(MI.up,CAD.up)
write.csv(DEGs, file="Results/MI-CAD/DEGs.csv", quote = F, row.names = F, col.names = F)

#MI.up <- read.csv("Results/MI_upper.csv", header = TRUE, check.names = FALSE)
#Normal.up <- read.csv("Results/Normal_upper.csv", header = TRUE, check.names = FALSE)
#DEGs <- read.csv("Results/DEGs.csv", header = TRUE, check.names = FALSE)


##### Looking for control probes in DEGs

sum(CAD.up$address == "Control") #613

sum(Normal.up$address == "Control") #459

sum(DEGs$address == "Control") #1072