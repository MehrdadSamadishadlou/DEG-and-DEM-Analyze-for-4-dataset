## Combining All data from the platform

setwd("/Users/mehrdad/Desktop/DataR")
library(GEOquery)

newinfo <- read.csv('GPL6244_GeneSymbol.csv', header = TRUE, check.names = FALSE)
Platform <- getGEO("GPL6244", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Platform/")
plt <- Platform.data <- subset(Platform@dataTable@table, select=c('ID', 'Gene symbol'))
plt <- plt[order(plt$ID),]

finalinfo <- cbind(newinfo, plt$`Gene symbol`)
write.csv(finalinfo, file="GPL6244_GeneSymbol.csv", quote = F, row.names = F, col.names = T)


#finalinfo <- read.csv('GPL6244_GeneSymbol.csv', header = TRUE, check.names = FALSE)

#sum(is.na(finalinfo$`By Chr`)):  8107
#sum(is.na(finalinfo$`By ID`)): 5767
#sum(is.na(finalinfo$`By Platform`)): 11102
#Number of control probes: 4468

attach(finalinfo)
n = 0
for (i in 1:33297) {
  if (is.na(`By Chr`[i]) & is.na(`By ID`[i]) & is.na(`By Platform`[i])) {
    n <- n + 1
  }
}
