################################ Extracting Platform Data (GPL6244) #######################################

setwd("/Users/mehrdad/Desktop/DataR")
library(GEOquery)
library(annotate)
library(quantmod)
library(hugene10sttranscriptcluster.db)
library(biomaRt)
library(reshape2)
library(AffyCompatible)


# Loading Data directly from the GEO

Platform <- getGEO("GPL6244", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")

## Finding number of control probes in the platform, which was 4468.
n = 0
for (i in 1:33297){
  if (Platform@dataTable@table[["Platform_SPOTID"]][i] == "control"){
    n = n + 1
  }
} # n = 4468

## Finding number of empty rows in 'Gene ID' column (which I think is representing Entrez ID) in the dataset.

n = 0
for (i in 1:33297){
  if (Platform@dataTable@table[["Gene ID"]][i] == ""){
    n = n + 1
  }
} ## n = 11102



## Ordering by ID and extracting the columns of interest.
colsofinterest <- c("ID", "Gene symbol", "Gene ID", "GI", "GenBank Accession")
Platform.data <- subset(Platform@dataTable@table, select=colsofinterest)
Platform.data <- Platform.data[order(Platform.data$ID),]

# Reading the downloaded .txt file from GEO, which has some difference with the directly downloaded data.
extraplatform <- read.csv('GPL6244-17930.txt', header = T, check.names = F, sep = "\t")

## Finding number of control probes in the this platform data, which was 4468.
n = 0
for (i in 1:33297){
    if (extraplatform$SPOT_ID[i] == "control"){
    n = n + 1
  }
} # n = 4468

## Ordering by ID and extracting the columns of interest.
extraplatform <- extraplatform[order(extraplatform$ID),]
extraplatform <- subset(extraplatform, select=c("ID", "gene_assignment"))


Platform.final <- cbind(Platform.data, extraplatform[,2])
colnames(Platform.final)[6] <- 'gene_assignment'

write.csv(Platform.final, file="Platform_data.csv", quote = F, row.names = F, col.names = T)


# Using "hugene10sttranscriptcluster.db" package for Data exteaction.

gset <- getGEO("GSE54475", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
geneset <- gset[[1]]
ID <- featureNames(geneset)

annodb <- "hugene10sttranscriptcluster.db"
Symbol <- as.character(lookUp(ID, annodb, "SYMBOL"))
Name   <- as.character(lookUp(ID, annodb, "GENENAME"))
Entrez <- as.character(lookUp(ID, annodb, "ENTREZID"))

#length(which(Name=="NA")) = 13491
#length(which(Symbol=="NA")) = 13491
#length(which(Entrez=="NA")) = 13491

plat <- cbind.data.frame(ID, Symbol, Name, Entrez)

### The empty Symbol and ID cells in the extracted data from the platform (Platform_data.csv) is 11102,
### which is lower than 13491 gained by hugene10sttranscriptcluster.db package. It should be kept in mind that
### 4468 of this empty cells are belong to "control" spots.







############################### Using BiomaRt for finding gene symbols #########################################


# listEnsembl() for seeing biomarts available for the first argument of the useEnsemble().

# With this we are using the speciy of interest:
## searchDatasets(mart = ensembl, pattern = 'Human')
## dataset              description    version
## 80 hsapiens_gene_ensembl Human genes (GRCh38.p13) GRCh38.p13

ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl')

##For selecting filters and attributes, following commands would be helpful.
#filters = listFilters(ensembl)
#attributes = listAttributes(ensembl)

## Example of using chromosomal_region
getBM(attributes = c('affy_hugene_1_0_st_v1', 'entrezgene_id', 'hgnc_symbol'), 
      filters = 'chromosomal_region', values = '1:54665840:54691137:+1', mart = ensembl)


## Example of using platform probe ID
getBM(attributes = c('affy_hugene_1_0_st_v1', 'entrezgene_id', 'hgnc_symbol'),
      filters = 'affy_hugene_1_0_st_v1', values = '7900488', mart = ensembl)


################ Using bioMart for finding unique gene symbols

#### Using Platform IDs

ids <- extraplatform$ID

Sym <- getBM(attributes = c('affy_hugene_1_0_st_v1', 'hgnc_symbol'),
                        filters = 'affy_hugene_1_0_st_v1', values = ids, mart = ensembl)


#write.csv(Sym, file="genesymbyidbiomartraw.csv", quote = F, row.names = F, col.names = T)
#Sym <- read.csv('GeneSymbolByID.csv', header = T, check.names = F)

attach(Sym)
GeneSymbolByID <- aggregate(Sym, by=list(affy_hugene_1_0_st_v1), FUN=toString, na.rm=TRUE)
GeneSymbolByID <- GeneSymbolByID[,-2]


# Finding missing probe ids in bioMart answer and adding them to the dataframe

absentprobes <- setdiff(ids, GeneSymbolByID$Group.1)

tmp <- data.frame(absentprobes, rep('NA', 4228))
colnames(tmp) <- colnames(GeneSymbolByID)
GeneSymbolByID <- rbind(GeneSymbolByID, tmp)
GeneSymbolByID <- GeneSymbolByID[order(GeneSymbolByID$Group.1),]


write.csv(GeneSymbolByID, file="GeneSymbolByID.csv", quote = F, row.names = F, col.names = T)



################### Using Chromosomal Region ############################


grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

getBM(attributes = c('affy_hugene_1_0_st_v1', 'entrezgene_id', 'hgnc_symbol'), 
      filters = 'chromosomal_region', values = '1:228327663:228336685:+1', mart = grch37)

# Preparing chromosomal region data in the format of needed for bioMart -> 1:69091:70008:1

attach(extraplatform)

chromosom <- data.frame('ID', 'address')

for (i in 1:33297) {
  if (SPOT_ID[i] == "control") {
    chromosom <- rbind(chromosom, c(ID[i], "Control"))
  } else {
    tmp <- paste(sub("chr", "", seqname[i]),":", RANGE_START[i], ":", RANGE_STOP[i], ":", RANGE_STRAND[i], "1", sep = "")
    chromosom <- rbind(chromosom, c(ID[i], tmp))
  }
}

chromosom <- chromosom[-1,]
colnames(chromosom) <- c("ID", "address")

write.csv(chromosom, file="chromosomal address.csv", quote = F, row.names = F, col.names = T)
#chromosom <- read.csv('chromosomal address.csv', header = T, check.names = F)

GeneSymbolByChr <- data.frame('ID', 'address', "hgnc_symbol")
colnames(GeneSymbolByChr) <- c('ID', 'address', "hgnc_symbol")
#GeneSymbolByChr <- read.csv('GeneSymbolByChr.csv', header = T, check.names = F)

attach(chromosom)
control <- address == "Control"

for (i in 1:33297) {
  if (!control[i]) {
    tmp <- getBM(attributes = c('hgnc_symbol'), filters = 'chromosomal_region', values = address[i], mart = grch37)
    GeneSymbolByChr <- rbind(GeneSymbolByChr, c(ID[i], address[i], as.character(tmp)))
  } else {GeneSymbolByChr <- rbind(GeneSymbolByChr, c(ID[i], address[i], 'NA'))}
  
}

GeneSymbolByChr <- GeneSymbolByChr[!duplicated(GeneSymbolByChr[,'ID']),]
GeneSymbolByChr <- GeneSymbolByChr[order(GeneSymbolByChr$ID),]

write.csv(GeneSymbolByChr, file="GeneSymbolByChr.csv", quote = F, row.names = F, col.names = T)


## Combining probe ID and chromosomal region data

GPL6244_GeneSymbol <- cbind(GeneSymbolByChr, GeneSymbolByID$hgnc_symbol)
write.csv(GPL6244_GeneSymbol, file="GPL6244_GeneSymbol.csv", quote = F, row.names = F, col.names = T)


####################################### NetAffx ########################################

rsrc <- NetAffxResource(user="samadishadlou@tbzmed.ac.ir", password='Mehrdad.69')

affxDescription(rsrc[["HuGene-1_0-st-v1"]])
# All Downloadable data for the platform using NetAffx package
#[1] "Probe Sequences, FASTA format"              "Probe Sequences, tabular format"           
#[3] "Probeset Annotations, CSV Format"           "Transcript Cluster Annotations, CSV"       
#[5] "Transcript Cluster Sequences, FASTA format" "BED File"                                  
#[7] "Background Probes File"                     "Intensity Layout File"                     
#[9] "cytoband csv"                               "EC Default Analysis Specifications"        
#[11] "Meta Probeset File"                         "Probe Group File"                          
#[13] "Probeset List File"                         "QC Control File"                           
#[15] "TAC Configuration File"                     "TAC 4.x Configuration file" 

annos <- rsrc[["HuGene-1_0-st-v1"]]

anno <- affxAnnotation(annos)[[11]]
df <- readAnnotation(rsrc, annotation=anno)

write.csv(df, file="df.csv", quote = F, row.names = F, col.names = T)

mpf <- read.table('Platform/NetAffx/Meta Probeset File.txt', header = T, check.names = F, sep = "\t")


