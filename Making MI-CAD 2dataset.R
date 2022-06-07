setwd('/Users/mehrdad/Desktop/DataR')

library(GEOquery)
library(Biobase)


Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 20)

GSE59867 <- getGEO("GSE59867", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Datasets/")
geneset67 <- GSE59867[[1]]


CAD67 <- geneset67[,geneset67@phenoData@data[["samples collection:ch1"]] == "N/A"]
CAD67@phenoData@data[["geo_accession"]] <- paste0(rep("67.", 46), CAD67@phenoData@data[["geo_accession"]], rep(".C", 46))


MI67 <- geneset67[,geneset67@phenoData@data[["samples collection:ch1"]] == "on the 1st day of MI (admission)"]
MI67@phenoData@data[["geo_accession"]] <- paste0(rep("67.", 111), MI67@phenoData@data[["geo_accession"]], rep(".M", 111))

GSE62646 <- getGEO("GSE62646", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Datasets/")
geneset46 <- GSE62646[[1]]

CAD46 <- geneset46[,geneset46@phenoData@data[["collection of blood samples:ch1"]] == "admission" ]
CAD46@phenoData@data[["geo_accession"]] <- paste0(rep("46.", 14), CAD46@phenoData@data[["geo_accession"]], rep(".M", 14))

MI46 <- geneset46[,geneset46@phenoData@data[["collection of blood samples:ch1"]] == "admission (on the 1st day of MI)"]
MI46@phenoData@data[["geo_accession"]] <- paste0(rep("46.", 28), MI46@phenoData@data[["geo_accession"]], rep(".M", 28))



superset <- Biobase::combine(CAD46, MI46, CAD67, MI67)
sample_status <- c(rep('CAD', 14), rep('MI', 28), rep('CAD', 46), rep('MI', 111))
superset@phenoData@data[["title"]] <- sample_status


saveRDS(superset, 'MI_CAD_2dataset.Rds')

