setwd('/Users/mehrdad/Desktop/DataR')
library(Biobase)

superset <- readRDS('Super4set.Rds')
exp <- exprs(superset)
colnames(exp) <- superset@phenoData@data[["geo_accession"]]

MI.up <- read.csv("Results/MI_upper.csv", header = TRUE, check.names = FALSE)
Normal.up <- read.csv("Results/Normal_upper.csv", header = TRUE, check.names = FALSE)

ids <- as.character(MI.up$ID)
MI_up_expression <- exp[ids,]
write.csv(MI_up_expression, file = "Results/MI_up_expression.csv", quote = F)

ids <- as.character(Normal.up$ID)
Normal_up_expression <- exp[ids,]
write.csv(Normal_up_expression, file = "Results/Normal_up_expression.csv", quote = F)

DEGs_expression <- rbind(MI_up_expression, Normal_up_expression)
write.csv(DEGs_expression, file = "Results/DEGs_expression.csv", quote = F)

##################################### For mirset ##################################3

mirset <- readRDS('mir4set.Rds')
exp <- exprs(mirset)
colnames(exp) <- mirset@phenoData@data[["geo_accession"]]

MI.up <- read.csv("Results/just microRNA/MI_upper_mirs.csv", header = TRUE, check.names = FALSE)
Normal.up <- read.csv("Results/just microRNA/Normal_upper_mirs.csv", header = TRUE, check.names = FALSE)

ids <- as.character(MI.up$ID)
MI_up_expression <- exp[ids,]
write.csv(MI_up_expression, file = "Results/just microRNA/MI_up_mirs_expression.csv", quote = F)

ids <- as.character(Normal.up$ID)
Normal_up_expression <- exp[ids,]
write.csv(Normal_up_expression, file = "Results/just microRNA/Normal_up_mirs_expression.csv", quote = F)

DEGs_expression <- rbind(MI_up_expression, Normal_up_expression)
write.csv(DEGs_expression, file = "Results/just microRNA/DEGs_mirs_expression.csv", quote = F)
