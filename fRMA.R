
library(frma)
library(hugene.1.0.st.v1frmavecs)
library(affy)
library(Biobase)
library(pheatmap)
library(ggplot2)
library(gplots)
library(limma)
library(oligo)
library(pd.hugene.1.0.st.v1)
library(tidyr)
library(dplyr)
library(stringr)
library(rgl)



data("hugene.1.0.st.v1frmavecs")
data("hugene.1.0.st.v1barcodevecs")
hugene10stv1frmavecs <- hugene.1.0.st.v1frmavecs

setwd("/Users/mehrdad/Desktop/affy")
setwd('75')
celFiles = list.celfiles()
affyRaw= read.celfiles(celFiles)
eset75 <- frma(affyRaw, input.vecs=hugene.1.0.st.v1frmavecs, target = "core", summarize="robust_weighted_average" )
b75 <- barcode(exprs(eset75), platform="GPL6244", mu=hugene.1.0.st.v1barcodevecs$mu,
          tau=hugene.1.0.st.v1barcodevecs$tau, output="binary")
exp75 <- exprs(eset75)


setwd("/Users/mehrdad/Desktop/affy")
setwd('09')
celFiles = list.celfiles()
affyRaw= read.celfiles(celFiles)
eset09 <- frma(affyRaw, input.vecs=hugene.1.0.st.v1frmavecs, target = "core", summarize="robust_weighted_average" )
b09 <- barcode(exprs(eset09), platform="GPL6244", mu=hugene.1.0.st.v1barcodevecs$mu,
               tau=hugene.1.0.st.v1barcodevecs$tau, output="binary")
exp09 <- exprs(eset09)


setwd("/Users/mehrdad/Desktop/affy")
setwd('46')
celFiles = list.celfiles()
affyRaw= read.celfiles(celFiles)
eset46 <- frma(affyRaw, input.vecs=hugene.1.0.st.v1frmavecs, target = "core", summarize="robust_weighted_average" )
b46 <- barcode(exprs(eset46), platform="GPL6244", mu=hugene.1.0.st.v1barcodevecs$mu,
               tau=hugene.1.0.st.v1barcodevecs$tau, output="binary")
exp46 <- exprs(eset46)


setwd("/Users/mehrdad/Desktop/affy")
setwd('67')
celFiles = list.celfiles()
affyRaw= read.celfiles(celFiles)
eset67 <- frma(affyRaw, input.vecs=hugene.1.0.st.v1frmavecs, target = "core", summarize="robust_weighted_average" )
b67 <- barcode(exprs(eset67), platform="GPL6244", mu=hugene.1.0.st.v1barcodevecs$mu,
               tau=hugene.1.0.st.v1barcodevecs$tau, output="binary")
exp67 <- exprs(eset67)


##############################################################################################################

setwd("/Users/mehrdad/Desktop/affy")

MH <- readRDS('/Users/mehrdad/Desktop/DataR/Datasets/Healthy_MI_4dataset.Rds')
CH <- readRDS('/Users/mehrdad/Desktop/DataR/Datasets/Healthy_CAD_4dataset.Rds')

neg_controls <- read.csv("/Users/mehrdad/Desktop/DataR/Platform/Negative Control genes IDs.csv", 
                         check.names = FALSE, header = FALSE)
neg_controls <- as.character(neg_controls[, 1])

pos_controls <- read.csv("/Users/mehrdad/Desktop/DataR/Platform/Positive Control genes IDs.csv", 
                         check.names = FALSE, header = FALSE)
pos_controls <- as.character(pos_controls[, 1])
pos_controls <- pos_controls[-1]


Group_H = c(rep("75", 5), rep("09", 46))
Group_M = c(rep("46", 28), rep("67", 111))
Group_C = c(rep("46", 14), rep("67", 46))

group_RLE <- c(rep("H", 51*33297), rep("C", 60*33297), rep("M", 139*33297))
group <- c(rep("H", 51), rep("C", 60), rep("M", 139))
dataset <- c(rep("75", 5), rep("09", 46), rep("46", 14), rep("67", 46), rep("46", 28), rep("67", 111))


all_fRMA_MI <- cbind(exp75, exp09, exp46[,15:42], exp67[,1:111])
colnames(all_fRMA_MI) <- MH@phenoData@data[["geo_accession"]]

all_fRMA_CAD <- cbind(exp75, exp09, exp46[,1:14], exp67[,112:157])
colnames(all_fRMA_CAD) <- CH@phenoData@data[["geo_accession"]]

all_fRMA <- cbind(all_fRMA_CAD, all_fRMA_MI[,52:190])

write.csv(all_fRMA, 'fRMA expression data.csv', quote = F)

all_MI_CAD <- all_fRMA[, 52:250]

pc <- prcomp(all_fRMA)
pcr <- data.frame(pc$r[, 1:3] , Group = group)
pdf("fRMA data.pdf", width = 20, height = 13)        
ggplot(pcr, aes(PC1 , PC2 , shape=dataset, color= group)) + geom_point(size = 7, alpha = 0.5) +
  ggtitle("PCA for fRMA data") +
  scale_color_manual(breaks = levels(as.factor(group)),
                     values = c("red", "darkgreen", "black")) +
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=25, face="bold"),
        legend.title = element_text(size = 18, face="bold"),
        legend.text = element_text(size = 14),
        plot.title = element_text(color="black", size=28, face="bold.italic", hjust = 0.5))
dev.off()

# RLE plot for primary data 
row_medians_assayData <- rowMedians(all_fRMA)
RLE_data <- sweep(all_fRMA, 1, row_medians_assayData)
RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- gather(RLE_data, patient_array, log2_expression_deviation)

pdf("RLE plot for fRMA data.pdf", width = 40, height = 15)
ggplot(RLE_data_gathered, aes(patient_array, log2_expression_deviation, color= group_RLE)) +
  geom_boxplot(outlier.shape = NA) + 
  ggtitle("RLE plot for fRMA data") +
  scale_color_manual(breaks = levels(as.factor(group_RLE)), values = c("red", "darkgreen", "black")) +
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "aquamarine4", angle = 90, size = 6.5, hjust = 1, face = "bold"),
        axis.text=element_text(size=24), 
        axis.title=element_text(size=25, face="bold"),
        legend.title = element_text(size = 22, face="bold"),
        legend.text = element_text(size = 18),
        plot.title = element_text(color="black", size=28, face="bold.italic", hjust = 0.5))
dev.off()


pdf("MI-Healthy/CorHeatmap-all-spearman.pdf", width = 25, height = 25)
pheatmap(cor(all_fRMA, method = "spearman"), labels_row = colnames(all_fRMA), 
         labels_col = colnames(all_fRMA), fontsize_row = 10, fontsize_col = 10,
         color=bluered(256), border_color = NA)
dev.off()

pdf("MI-Healthy/CorHeatmap-MI-Healthy-spearman.pdf", width = 25, height = 25)
pheatmap(cor(all_fRMA_MI, method = "spearman"), labels_row = colnames(all_fRMA_MI), 
         labels_col = colnames(all_fRMA_MI), fontsize_row = 10, fontsize_col = 10,
         color=bluered(256), border_color = NA)
dev.off()

pdf("MI-Healthy/CorHeatmap-CAD-Healthy-spearman.pdf", width = 25, height = 25)
pheatmap(cor(all_fRMA_CAD, method = "spearman"), labels_row = colnames(all_fRMA_CAD), 
         labels_col = colnames(all_fRMA_CAD), fontsize_row = 15, fontsize_col = 15,
         color=bluered(256), border_color = NA)
dev.off()


pdf("MI-Healthy/CorHeatmap-MI-CAD-spearman.pdf", width = 25, height = 25)
pheatmap(cor(all_MI_CAD, method = "spearman"), labels_row = colnames(all_MI_CAD), 
         labels_col = colnames(all_MI_CAD), fontsize_row = 11, fontsize_col = 11,
         color=bluered(256), border_color = NA)
dev.off()

##### Making an Expressionset ######

MH_fRMA <- all_fRMA[, -(52:111)]
CH_fRMA <- all_fRMA[, 1:111]

colnames(MH_fRMA) <- colnames(exprs(MH))
exprs(MH) <- MH_fRMA
expmh <- exprs(MH)
colnames(expmh) <- MH@phenoData@data[["geo_accession"]]

colnames(CH_fRMA) <- colnames(exprs(CH))
exprs(CH) <- CH_fRMA
expch <- exprs(CH)
colnames(expch) <- CH@phenoData@data[["geo_accession"]]

######## Differential expression Analysis with Limma using ComBat-SQN data ##########
### MI-Healthy ##############

sample_status <- factor(MH@phenoData@data[["title"]])
MH$description <- sample_status
Design <- model.matrix(~ description + 0, MH)
colnames(Design) <- levels(sample_status)


colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(MH, Design)
cont.matrix <- makeContrasts(MI-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "Gene.symbol","adj.P.Val", "logFC"))
tT <- tT[order(tT$ID),]
tt_MI <- tT

MI_up <- subset(tT, logFC>1 & adj.P.Val<0.05)

Healthy_up <- subset(tT, logFC < -1 & adj.P.Val<0.05)

DEGs <- rbind(MI_up, Healthy_up)

sum(DEGs$ID %in% pos_controls) 
# 83
sum(DEGs$ID %in% neg_controls) 
# 683

write.csv(tT, "MI-Healthy_statistical-data.csv", row.names=F, quote = F)
write.csv(MI_up, file="MI_upper.csv", quote = F, row.names = F, col.names = F)
write.csv(Healthy_up, file="Healthy_upper.csv", quote = F, row.names = F, col.names = F)
write.csv(DEGs, file="DEGs.csv", quote = F, row.names = F, col.names = F)

MI_up_no_cont <- MI_up[! MI_up$ID %in% c(neg_controls, pos_controls), ]
Healthy_up_no_cont <- Healthy_up[! Healthy_up$ID %in% c(neg_controls, pos_controls), ]
DEGs_no_cont <- rbind(MI_up_no_cont, Healthy_up_no_cont)

write.csv(MI_up_no_cont, file="MI_upper_no_cont.csv", quote = F, row.names = F, col.names = F)
write.csv(Healthy_up_no_cont, file="Healthy_upper_no_cont.csv", quote = F, row.names = F, col.names = F)
write.csv(DEGs_no_cont, file="DEGs_no_cont.csv", quote = F, row.names = F, col.names = F)


########## Finding mirs #########

mirs_id <- DEGs_no_cont[str_detect(string = DEGs_no_cont$Gene.symbol, pattern = 'MIR'), "ID"]
mirs_id <- mirs_id[! str_detect(mirs_id, "8063793")]
mirs_MI_exp <- expmh[rownames(expmh) %in% mirs_id, ]
mirs_MI_exp <- cbind(mirs_MI_exp, expch[rownames(expch) %in% mirs_id, 52:111])

write.csv(mirs_MI_exp, file="MI-Healthy/MI_Healthy_DEMs_Exp.csv", quote = F, row.names = T, col.names = F)

### CAD-Healthy ##############

sample_status <- factor(CH@phenoData@data[["title"]])
CH$description <- sample_status
Design <- model.matrix(~ description + 0, CH)
colnames(Design) <- levels(sample_status)


colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(CH, Design)
cont.matrix <- makeContrasts(CAD-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "Gene.symbol","adj.P.Val", "logFC"))
tT <- tT[order(tT$ID),]
tt_CAD <- tT

CAD_up <- subset(tT, logFC>1 & adj.P.Val<0.05)

Healthy_up <- subset(tT, logFC < -1 & adj.P.Val<0.05)

DEGs <- rbind(CAD_up, Healthy_up)

sum(DEGs$ID %in% pos_controls) 
# 91
sum(DEGs$ID %in% neg_controls) 
# 961

write.csv(tT, "CAD-Healthy/CAD-Healthy_statistical-data.csv", row.names=F, quote = F)
write.csv(CAD_up, file="CAD-Healthy/CAD_upper.csv", quote = F, row.names = F, col.names = F)
write.csv(Healthy_up, file="CAD-Healthy/Healthy_upper.csv", quote = F, row.names = F, col.names = F)
write.csv(DEGs, file="CAD-Healthy/DEGs.csv", quote = F, row.names = F, col.names = F)

CAD_up_no_cont <- CAD_up[! CAD_up$ID %in% c(neg_controls, pos_controls), ]
Healthy_up_no_cont <- Healthy_up[! Healthy_up$ID %in% c(neg_controls, pos_controls), ]
DEGs_no_cont <- rbind(CAD_up_no_cont, Healthy_up_no_cont)

write.csv(CAD_up_no_cont, file="CAD-Healthy/CAD_upper_no_cont.csv", quote = F, row.names = F, col.names = F)
write.csv(Healthy_up_no_cont, file="CAD-Healthy/Healthy_upper_no_cont.csv", quote = F, row.names = F, col.names = F)
write.csv(DEGs_no_cont, file="CAD-Healthy/DEGs_no_cont.csv", quote = F, row.names = F, col.names = F)

########## Finding mirs $$$$$$$$$$4

mirs_id <- DEGs_no_cont[str_detect(string = DEGs_no_cont$Gene.symbol, pattern = 'MIR'), "ID"]
mirs_id <- mirs_id[! str_detect(mirs_id, "8063793")]
mirs_CAD_exp <- expch[rownames(expch) %in% mirs_id, ]
mirs_CAD_exp <- cbind(mirs_CAD_exp, expmh[rownames(expmh) %in% mirs_id, 52:190])

write.csv(mirs_CAD_exp, file="CAD-Healthy/CAD_Healthy_DEMs_Exp.csv", quote = F, row.names = T, col.names = F)



## Finding the most consistent mir for normalization

consistent_mirs_MI <- tt_MI[str_detect(string = tt_MI$Gene.symbol, pattern = 'MIR'), ]
consistent_mirs_CAD <- tt_CAD[str_detect(string = tt_CAD$Gene.symbol, pattern = 'MIR'), ]

consistent_mirs_MI <- consistent_mirs_MI[abs(consistent_mirs_MI$logFC) < 0.05, ]
consistent_mirs_MI <- consistent_mirs_MI[abs(consistent_mirs_MI$adj.P.Val) > 0.5, ]
#consistent_mirs_MI <- consistent_mirs_MI[order(consistent_mirs_MI$adj.P.Val, decreasing = TRUE),]
consistent_mirs_CAD <- consistent_mirs_CAD[abs(consistent_mirs_CAD$logFC) < 0.05, ]
consistent_mirs_CAD <- consistent_mirs_CAD[abs(consistent_mirs_CAD$adj.P.Val) > 0.5, ]
#consistent_mirs_CAD <- consistent_mirs_CAD[order(consistent_mirs_CAD$adj.P.Val, decreasing = TRUE),]

common_mirs <- consistent_mirs_MI[consistent_mirs_MI$ID %in% consistent_mirs_CAD$ID, "ID"]

consistent_mirs_MI <- consistent_mirs_MI[consistent_mirs_MI$ID %in% common_mirs, ]
consistent_mirs_CAD <- consistent_mirs_CAD[consistent_mirs_CAD$ID %in% common_mirs, ]

# Candidate mir: 8087252  MIR191 

mir191_exp <- all_fRMA[rownames(all_fRMA) == '8087252', ]
summary(mir191_exp)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#6.049   6.273   6.356   6.372   6.472   6.877 
sd(mir191_exp)
#0.142

write.csv(all_fRMA[rownames(all_fRMA) == '8087252', ], file="consistant_191_Exp.csv", quote = F, row.names = T, col.names = F)

### MI_CAD ###

MC <- combine(MH[, 52:190], CH[, 52:111])

sample_status <- factor(MC@phenoData@data[["title"]])
MC$description <- sample_status
Design <- model.matrix(~ description + 0, MC)
colnames(Design) <- levels(sample_status)


colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(MC, Design)
cont.matrix <- makeContrasts(MI-CAD, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "Gene.symbol","adj.P.Val", "logFC"))
tT <- tT[order(tT$ID),]


MI_up <- subset(tT, logFC>1 & adj.P.Val<0.05)

CAD_up <- subset(tT, logFC < -1 & adj.P.Val<0.05)

DEGs <- rbind(MI_up, CAD_up)


write.csv(tT, "MI-CAD/MI-CAD-statistical-data.csv", row.names=F, quote = F)
write.csv(DEGs, file="MI-CAD/DEGs.csv", quote = F, row.names = F, col.names = F)


##### Working with barcode data

status_healthy <- rowMeans(cbind(b75, b09))
status_CAD <- rowMeans(cbind(b46[, 1:14], b67[, 112:157]))
status_MI <- rowMeans(cbind(b46[, 14:42], b67[, 1:111]))

barcode_status <- as.data.frame(cbind(status_healthy, status_CAD, status_MI))

attach(barcode_status)

CAD_up_barcode <- c()
MI_up_barcode <- c()
Healthy_upper_MI_barcode <- c()
Healthy_upper_CAD_barcode <- c()


for (i in 1:33297) {
  if ((status_healthy[i] >= 0.8 & status_MI[i] <= 0.5)) {
    Healthy_upper_MI_barcode <- c(Healthy_upper_MI_barcode, rownames(barcode_status)[i])
  } else if ((status_healthy[i] <= 0.5 & status_MI[i] >= 0.8)) {
    MI_up_barcode <- c(MI_up_barcode, rownames(barcode_status)[i])
  } else if ((status_healthy[i] >= 0.8 & status_CAD[i] <= 0.5)) {
    Healthy_upper_CAD_barcode <- c(Healthy_upper_CAD_barcode, rownames(barcode_status)[i])
  } else if ((status_healthy[i] <= 0.5 & status_CAD[i] >= 0.8)) {
    CAD_up_barcode <- c(CAD_up_barcode, rownames(barcode_status)[i])
  }
}


Healthy_upper_CAD_barcode <- Healthy_upper_CAD_barcode[! Healthy_upper_CAD_barcode %in% c(neg_controls, pos_controls)]
Healthy_upper_MI_barcode <- Healthy_upper_MI_barcode[! Healthy_upper_MI_barcode %in% c(neg_controls, pos_controls)]
MI_up_barcode <- MI_up_barcode[! MI_up_barcode %in% c(neg_controls, pos_controls)]
CAD_up_barcode <- CAD_up_barcode[! CAD_up_barcode %in% c(neg_controls, pos_controls)]


write.csv(MI_up_barcode, 'MI-Healthy/barcode/MI_up.csv', row.names = FALSE)
write.csv(CAD_up_barcode, 'MI-Healthy/barcode/CAD_up.csv', row.names = FALSE)
write.csv(Healthy_upper_CAD_barcode, 'MI-Healthy/barcode/Healthy_upper_CAD_barcode.csv', row.names = FALSE)
write.csv(Healthy_upper_MI_barcode, 'MI-Healthy/barcode/Healthy_upper_MI_barcode.csv', row.names = FALSE)






### 3D Plots

pc <- prcomp(all_fRMA_CAD[, 52:111])
pcr <- data.frame(pc$r[, 1:3] , Group = Group_C)

ggg <- c(Group_C, Group_M)
MI_color <- case_when(ggg == "46" ~ "red",
                      ggg == "67" ~ "gray",
                      TRUE ~ NA_character_)

mc_group <- c(rep("gray", 60), rep("red", 139))

plot3d(x = pcr$PC1, y = pcr$PC2, z = pcr$PC3, col = mc_group, type = "s", size = 1.5)
snapshot3d("MI3d.png")
movie3d(spin3d(axis = c(0,0,1), rpm = 4), duration = 15)

pc <- prcomp(all_MI_CAD)
pcr <- data.frame(pc$r[, 1:3] , Group = Group_M)

