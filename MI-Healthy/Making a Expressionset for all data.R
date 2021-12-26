setwd('/Users/mehrdad/Desktop/DataR')

library(GEOquery)
library(Biobase)


#Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 20)
GSE54475 <- getGEO("GSE54475", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Datasets/")
geneset75 <- GSE54475[[1]]

sub75 <- geneset75[,1:5]
sub75@phenoData@data[["geo_accession"]] <- paste0(rep("75.", 5), sub75@phenoData@data[["geo_accession"]], rep(".H", 5))


GSE56609 <- getGEO("GSE56609", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Datasets/")
geneset09 <- GSE56609[[1]]

sub09 <- geneset09[,geneset09@phenoData@data[["source_name_ch1"]] == "PBMC, fasting"]
sub09@phenoData@data[["geo_accession"]] <- paste0(rep("09.", 46), sub09@phenoData@data[["geo_accession"]], rep(".H", 46))


GSE59867 <- getGEO("GSE59867", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Datasets/")
geneset67 <- GSE59867[[1]]

sub67 <- geneset67[,geneset67@phenoData@data[["samples collection:ch1"]] == "on the 1st day of MI (admission)"]
sub67@phenoData@data[["geo_accession"]] <- paste0(rep("67.", 111), sub67@phenoData@data[["geo_accession"]], rep(".M", 111))

GSE62646 <- getGEO("GSE62646", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Datasets/")
geneset46 <- GSE62646[[1]]

sub46 <- geneset46[,geneset46@phenoData@data[["collection of blood samples:ch1"]] == "admission (on the 1st day of MI)"]
sub46@phenoData@data[["geo_accession"]] <- paste0(rep("46.", 28), sub46@phenoData@data[["geo_accession"]], rep(".M", 28))



superset <- combine(sub75, sub09, sub46, sub67)
sample_status <- c(rep('Healthy', 51), rep('MI', 139))
superset@phenoData@data[["title"]] <- sample_status


saveRDS(superset, 'Super4set.Rds')


#superset <- readRDS('Super4set.Rds')



#exp <- read.csv('ExpressionDataset.csv', header = TRUE, check.names = FALSE)
#rownames(exp) <- exp[,1]
#exp <- exp[,-1]

#sample_status <- c(rep('Healthy', 51), rep('MI', 139))
#p.data <- data.frame("Status" = sample_status, row.names = colnames(exp))
#phenoData <- new("AnnotatedDataFrame", data=p.data)


#total_geneset <- ExpressionSet(assayData = as.matrix(exp),
#                             phenoData = phenoData,
#                             experimentData = geneset75@experimentData,
#                             featureData = geneset75@featureData)
#

