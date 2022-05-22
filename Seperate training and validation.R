setwd('/Users/mehrdad/Desktop/affy')

library(Biobase)
library(pheatmap)
library(ggplot2)
library(gplots)
library(limma)
library(oligo)
library(tidyr)
library(dplyr)
library(stringr)
library(rgl)
library(affy)

neg_controls <- read.csv("/Users/mehrdad/Desktop/DataR/Platform/Negative Control genes IDs.csv", 
                         check.names = FALSE, header = FALSE)
neg_controls <- as.character(neg_controls[, 1])

pos_controls <- read.csv("/Users/mehrdad/Desktop/DataR/Platform/Positive Control genes IDs.csv", 
                         check.names = FALSE, header = FALSE)
pos_controls <- as.character(pos_controls[, 1])
pos_controls <- pos_controls[-1]

all_fRMA <- read.csv('fRMA expression data.csv', check.names = F)
rownames(all_fRMA) <- all_fRMA[, 1]
all_fRMA <- all_fRMA[, -1]

train_set <- all_fRMA[, 6:51]
train_set <- cbind(train_set, all_fRMA[66:111])
train_set <- cbind(train_set, all_fRMA[140:250])

write.csv(all_fRMA, 'fRMA train_set.csv', quote = F)

validation_set <- all_fRMA[, 1:5]
validation_set <- cbind(validation_set, all_fRMA[52:65])
validation_set <- cbind(validation_set, all_fRMA[112:139])

write.csv(all_fRMA, 'fRMA validation_set.csv', quote = F)

## Making an expression set

MH <- readRDS('/Users/mehrdad/Desktop/DataR/Datasets/Healthy_MI_4dataset.Rds')
CH <- readRDS('/Users/mehrdad/Desktop/DataR/Datasets/Healthy_CAD_4dataset.Rds')

MH_fRMA <- as.matrix(all_fRMA[, -(52:111)])
CH_fRMA <- as.matrix(all_fRMA[, 1:111])

colnames(MH_fRMA) <- colnames(exprs(MH))
exprs(MH) <- MH_fRMA
expmh <- exprs(MH)
colnames(expmh) <- MH@phenoData@data[["geo_accession"]]

colnames(CH_fRMA) <- colnames(exprs(CH))
exprs(CH) <- CH_fRMA
expch <- exprs(CH)
colnames(expch) <- CH@phenoData@data[["geo_accession"]]

validation_MI_expression_set <- Biobase::combine(MH[, 1:5], MH[, 52:79])
train_MI_expression_set <- Biobase::combine(MH[, 6:51], MH[, 80:190])
validation_CAD_expression_set <- Biobase::combine(CH[, 1:5], CH[, 52:65])
train_CAD_expression_set <- Biobase::combine(CH[, 6:51], CH[, 66:111])


## Differential expression Analysis with Limma 
### MI-Healthy 
#### Train Set

sample_status <- factor(train_MI_expression_set@phenoData@data[["title"]])
train_MI_expression_set$description <- sample_status
Design <- model.matrix(~ description + 0, train_MI_expression_set)
colnames(Design) <- levels(sample_status)


colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(train_MI_expression_set, Design)
cont.matrix <- makeContrasts(MI-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "Gene.symbol","adj.P.Val", "logFC"))
tT <- tT[order(tT$ID),]


MI_up <- subset(tT, logFC>1 & adj.P.Val<0.05)

Healthy_up <- subset(tT, logFC < -1 & adj.P.Val<0.05)

DEGs <- rbind(MI_up, Healthy_up)

sum(DEGs$ID %in% pos_controls) 
# 87
sum(DEGs$ID %in% neg_controls) 
# 733

MI_up_no_cont <- MI_up[! MI_up$ID %in% c(neg_controls, pos_controls), ]
Healthy_up_no_cont <- Healthy_up[! Healthy_up$ID %in% c(neg_controls, pos_controls), ]
DEGs_no_cont <- rbind(MI_up_no_cont, Healthy_up_no_cont)

write.csv(MI_up_no_cont, file="MI-Healthy/Train-Val/Train/MI_upper_train.csv", quote = F, row.names = F, col.names = F)
write.csv(Healthy_up_no_cont, file="MI-Healthy/Train-Val/Train/Healthy_upper_train.csv", quote = F, row.names = F, col.names = F)
write.csv(DEGs_no_cont, file="MI-Healthy/Train-Val/Train/DEGs_train.csv", quote = F, row.names = F, col.names = F)

#### Validation Set

sample_status <- factor(validation_MI_expression_set@phenoData@data[["title"]])
validation_MI_expression_set$description <- sample_status
Design <- model.matrix(~ description + 0, validation_MI_expression_set)
colnames(Design) <- levels(sample_status)


colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(validation_MI_expression_set, Design)
cont.matrix <- makeContrasts(MI-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "Gene.symbol","adj.P.Val", "logFC"))
tT <- tT[order(tT$ID),]


MI_up <- subset(tT, logFC>1 & adj.P.Val<0.05)

Healthy_up <- subset(tT, logFC < -1 & adj.P.Val<0.05)

DEGs <- rbind(MI_up, Healthy_up)

MI_up_no_cont <- MI_up[! MI_up$ID %in% c(neg_controls, pos_controls), ]
Healthy_up_no_cont <- Healthy_up[! Healthy_up$ID %in% c(neg_controls, pos_controls), ]
DEGs_no_cont <- rbind(MI_up_no_cont, Healthy_up_no_cont)

write.csv(MI_up_no_cont, file="MI-Healthy/Train-Val/Validation/MI_upper_val.csv", quote = F, row.names = F, col.names = F)
write.csv(Healthy_up_no_cont, file="MI-Healthy/Train-Val/Validation/Healthy_upper_val.csv", quote = F, row.names = F, col.names = F)
write.csv(DEGs_no_cont, file="MI-Healthy/Train-Val/Validation/DEGs_val.csv", quote = F, row.names = F, col.names = F)


### CAD-Healthy 
#### Train Set

sample_status <- factor(train_CAD_expression_set@phenoData@data[["title"]])
train_CAD_expression_set$description <- sample_status
Design <- model.matrix(~ description + 0, train_CAD_expression_set)
colnames(Design) <- levels(sample_status)


colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(train_CAD_expression_set, Design)
cont.matrix <- makeContrasts(CAD-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "Gene.symbol","adj.P.Val", "logFC"))
tT <- tT[order(tT$ID),]


CAD_up <- subset(tT, logFC>1 & adj.P.Val<0.05)

Healthy_up <- subset(tT, logFC < -1 & adj.P.Val<0.05)

DEGs <- rbind(CAD_up, Healthy_up)

sum(DEGs$ID %in% pos_controls) 
# 102
sum(DEGs$ID %in% neg_controls) 
# 1004

CAD_up_no_cont <- CAD_up[! CAD_up$ID %in% c(neg_controls, pos_controls), ]
Healthy_up_no_cont <- Healthy_up[! Healthy_up$ID %in% c(neg_controls, pos_controls), ]
DEGs_no_cont <- rbind(CAD_up_no_cont, Healthy_up_no_cont)

write.csv(CAD_up_no_cont, file="CAD-Healthy/Train-Val/Train/CAD_upper_train.csv", quote = F, row.names = F, col.names = F)
write.csv(Healthy_up_no_cont, file="CAD-Healthy/Train-Val/Train/Healthy_upper_train.csv", quote = F, row.names = F, col.names = F)
write.csv(DEGs_no_cont, file="CAD-Healthy/Train-Val/Train/DEGs_train.csv", quote = F, row.names = F, col.names = F)

#### Validation Set

sample_status <- factor(validation_CAD_expression_set@phenoData@data[["title"]])
validation_CAD_expression_set$description <- sample_status
Design <- model.matrix(~ description + 0, validation_CAD_expression_set)
colnames(Design) <- levels(sample_status)


colSums(Design)
#For checking correctness of the procedure.

fit <- lmFit(validation_CAD_expression_set, Design)
cont.matrix <- makeContrasts(CAD-Healthy, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID", "Gene.symbol","adj.P.Val", "logFC"))
tT <- tT[order(tT$ID),]


CAD_up <- subset(tT, logFC>1 & adj.P.Val<0.05)

Healthy_up <- subset(tT, logFC < -1 & adj.P.Val<0.05)

DEGs <- rbind(CAD_up, Healthy_up)

CAD_up_no_cont <- CAD_up[! CAD_up$ID %in% c(neg_controls, pos_controls), ]
Healthy_up_no_cont <- Healthy_up[! Healthy_up$ID %in% c(neg_controls, pos_controls), ]
DEGs_no_cont <- rbind(CAD_up_no_cont, Healthy_up_no_cont)

write.csv(CAD_up_no_cont, file="CAD-Healthy/Train-Val/Validation/CAD_upper_val.csv", quote = F, row.names = F, col.names = F)
write.csv(Healthy_up_no_cont, file="CAD-Healthy/Train-Val/Validation/Healthy_upper_val.csv", quote = F, row.names = F, col.names = F)
write.csv(DEGs_no_cont, file="CAD-Healthy/Train-Val/Validation/DEGs_val.csv", quote = F, row.names = F, col.names = F)
