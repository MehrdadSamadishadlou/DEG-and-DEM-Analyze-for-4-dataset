setwd('/Users/mehrdad/Desktop/DataR')

library(Biobase)
library(stringr)


##################### Finding microRNA probes ##########################

sym <- read.csv('GPL6244_GeneSymbol.csv', header = TRUE, check.names = FALSE)

chor <- sym[str_detect(sym$`By Chr`, "MIR"), ]
chor <- chor[order(chor$ID),]
row.names(chor) <- 1:8909
chor <- chor[1:802,]

id <- sym[str_detect(sym$`By ID`, "MIR"), ]
id <- id[order(id$ID),]
row.names(id) <- 1:6077
id <- id[1:310,]

plat <- sym[str_detect(sym$`By Platform`, "MIR"), ]
plat <- plat[order(plat$ID),]
row.names(plat) <- 1:11320
plat <- plat[1:218,]

mir_probes <- rbind(chor, id, plat)
mir_probes <- mir_probes[!duplicated(mir_probes[,'ID']),]
write.csv(mir_probes, file="Results/just microRNA/mir_probes.csv", quote = F, row.names = F, col.names = F)
#mir_probes <- read.csv('Results/just microRNA/mir_probes.csv', header = TRUE, check.names = FALSE)



############# Extracting microRNA probes expression profile ###################

superset <- readRDS('Healthy_CAD_4dataset.Rds')

mirs_probe_id <- as.character(mir_probes$ID)
mirset <- superset[mirs_probe_id,]

saveRDS(mirset, 'Healthy_CAD_4dataset_microRNA.Rds')



#exp <- read.csv('ExpressionDataset.csv', header = TRUE, check.names = FALSE)
#rownames(exp) <- exp[,1]
#exp <- exp[,-1]

#mir_exp <- c()
#for (i in 1:894) {
#  mir_exp <- rbind(mir_exp, exp[match(mir_probes$ID[i],row.names(exp)),])
#}

#write.csv(mir_exp, file="Results/just microRNA/mir_exp.csv", quote = F, row.names = F, col.names = F)

