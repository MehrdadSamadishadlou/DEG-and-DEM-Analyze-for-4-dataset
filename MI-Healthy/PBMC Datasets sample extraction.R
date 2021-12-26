setwd("/Users/mehrdad/Desktop/DataR")
library(GEOquery)

## Loading Data
gset <- getGEO("GSE54475", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
geneset <- gset[[1]]
Exp.set <- exprs(geneset)

## Extracting the samples
theset <- Exp.set[,1:5]

## Changing samples names and saving the dataset

colnames(theset) <- paste0(rep("75.", length(colnames(theset))), colnames(theset), 
                           rep(".H", length(colnames(theset))))

GSE75 <- theset

write.csv(GSE75, file="GSE54475_Healthy_Samples.csv", quote = F, row.names = T, col.names = T)

######## Repeating the above instruction for other datasets

## GSE56609
gset <- getGEO("GSE56609", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
geneset <- gset[[1]]
Exp.set <- exprs(geneset)

for (i in 1:92){
  if (as.integer(substr(colnames(Exp.set)[i],9,10)) %% 2 == 0){Exp.set <- Exp.set[,-i]}
}

theset <- Exp.set
colnames(theset) <- paste0(rep("09.", length(colnames(theset))), colnames(theset), 
                           rep(".H", length(colnames(theset))))

GSE09 <- theset

write.csv(GSE09, file="GSE56609_Healthy_Samples.csv", quote = F, row.names = T, col.names = T)


## GSE59867

gset <- getGEO("GSE59867", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
geneset <- gset[[1]]
Exp.set <- exprs(geneset)

colstodel=2
for (i in 3:436){
  title = geneset@phenoData@data[["title"]][i]
  if (substr(title, nchar(title) - 9, nchar(title)) != 'sampling 1'){
    colstodel <- cbind(colstodel, as.integer(i))
  }
}

Exp.set <- Exp.set[,-colstodel]

theset <- Exp.set

colnames(theset) <- paste0(rep("67.", length(colnames(theset))), colnames(theset), 
                           rep(".M", length(colnames(theset))))

GSE67 <- theset

write.csv(GSE6, file="GSE59867_MI_Samples.csv", quote = F, row.names = T, col.names = T)

##GSE62646

gset <- getGEO("GSE62646", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
geneset <- gset[[1]]
Exp.set <- exprs(geneset)

colstodel=1
for (i in 2:436){
  title = geneset@phenoData@data[["title"]][i]
  if (substr(title, 1, 43) != 'STEMI patient, blood collected on admission'){
    colstodel <- cbind(colstodel, as.integer(i))
  }
}

Exp.set <- Exp.set[,-colstodel]

theset <- Exp.set

colnames(theset) <- paste0(rep("46.", length(colnames(theset))), colnames(theset), 
                           rep(".M", length(colnames(theset))))

GSE46 <- theset

write.csv(GSE46, file="GSE62646_MI_Samples.csv", quote = F, row.names = T, col.names = T)


####Concatenating all datasets

TheFinal <- cbind(GSE75, GSE09, GSE46, GSE67)

write.csv(TheFinal, file="ExpressionDataset.csv", quote = F, row.names = T, col.names = T)



############ Combining  the Datasets and  the platform data

#TheFinal <- read.csv('ExpressionDataset.csv', header = T, check.names = F)

UltimateDataset <- cbind(Platform.final, TheFinal[,-1])

write.csv(UltimateDataset, file="UltimateDataset.csv", quote = F, row.names = F, col.names = T)




