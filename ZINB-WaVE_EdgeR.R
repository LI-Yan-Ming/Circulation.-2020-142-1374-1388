rm(list=ls())
setwd("H:\\AHA_scRNAseq\\human\\TAA_Con_8vs3\\AllCells\\ZINB-WaVE_EdgeR\\")

library(zinbwave)
library(scRNAseq)
library(matrixStats)
library(ggplot2)
library(stringr)
library(edgeR)

meta_all <- read.csv("meta_addsubcluster.csv", header=T, row.names = "X", stringsAsFactors = F)
files <- list.files(path=".", pattern = "UMI_*")
names <- str_split_fixed(files,"_",3)[,3]
names <- str_split_fixed(names,"[.]",2)[,1]

for(i in 1:length(files)){
  counts <- as.matrix(read.csv(files[i], header=T, row.names = "X", stringsAsFactors = F))
  
  meta <- meta_all[meta_all$celltype2==names[i],]
  sce <- SingleCellExperiment(assays = list(counts = counts))
  
  ########## filtering #########
  filter <- rowSums(assay(sce)>0) > dim(counts)[2]/4
  sce <- sce[filter,]

  ########### fit ##############
  zinb <- zinbFit(sce, K=2, epsilon=1000)  ## time consuming ###
  sce_zinb <- zinbwave(sce, fitted_model = zinb, K = 2, epsilon=1000,
                       observationalWeights = TRUE)
  
  ######### differential expression #################
  weights <- assay(sce_zinb, "weights")
  
  dge <- DGEList(assay(sce_zinb))
  dge <- calcNormFactors(dge)
  
  design <- model.matrix(~meta$stim)
  dge$weights <- weights
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  
  lrt <- glmWeightedF(fit, coef = 2)
  topTags(lrt,20)
  
  write.csv(topTags(lrt, 2000), paste0("lrt2_", names[i], "_ConvsATAA.csv"), quote = F, row.names = T)
}


##################### Geneontology and KEGG analysis ##############################
files <- list.files(path=".", pattern = "lrt2_*")
dfs <- lapply(files, function(x) read.csv(x, header=T, row.names = "X", stringsAsFactors = F))

library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

data_background <- toTable(org.Hs.egSYMBOL)
for (i in 1: length(dfs)){
  each <- dfs[[i]]
  
  each_up <- each[which(each$logFC< (-0.3) & each$FDR<0.05),]
  test1 = bitr(row.names(each_up), fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
  each_up_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", 
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                         qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
  write.table(each_up_GO, paste0(substr(files[i], 6, (nchar(files[i])-14)),"_ATAAvsCON_up_GO.csv"), sep=",", quote=F)
  each_up_KEGG <- enrichKEGG(test1$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                             qvalueCutoff = 0.2)
  write.table(each_up_KEGG, paste0(substr(files[i], 6, (nchar(files[i])-14)),"_ATAAvsCON_up_KEGG.csv"), sep=",", quote=F)
  
  each_down <- each[which(each$logFC>0.3 & each$FDR<0.05),]
  test1 = bitr(row.names(each_down), fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
  each_down_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", 
                           pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                           qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
  write.table(each_down_GO, paste0(substr(files[i], 6, (nchar(files[i])-14)),"_ATAAvsCON_down_GO.csv"), sep=",", quote=F)
  each_down_KEGG <- enrichKEGG(test1$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                               qvalueCutoff = 0.2)
  write.table(each_down_KEGG, paste0(substr(files[i], 6, (nchar(files[i])-14)),"_ATAAvsCON_down_KEGG.csv"), sep=",", quote=F)
}


