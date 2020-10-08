rm(list=ls())
setwd("H:\\AHA_scRNAseq\\human\\TAA_Con_8vs3\\nonImmune\\")

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(DESeq2)
library(pheatmap)

combined <- readRDS("nonImmuneCombined.rds")

#################################### set up seperate seurat projects #########################################
files <- list.files(path=".", pattern = "*_fromIntegrative.csv")
set_up <- function(UMI,project,stim){
  data <- read.csv(UMI, header=T, row.names = "X", stringsAsFactors = F)
  each <- CreateSeuratObject(counts = data, project = project, min.cells = 5)
  each$stim <- stim
  each <- NormalizeData(each, verbose = FALSE)
  each <- FindVariableFeatures(each, selection.method = "vst", nfeatures = 2000)
}

Con4  <- set_up(files[9], "Con4", "Control")
Con6  <- set_up(files[10], "Con6", "Control")
Con9  <- set_up(files[11], "Con9", "Control")
TAA1  <- set_up(files[1], "ATAA1", "ATAA")
TAA2  <- set_up(files[2], "ATAA2", "ATAA")
TAA3  <- set_up(files[3], "ATAA3", "ATAA")
TAA4  <- set_up(files[4], "ATAA4", "ATAA")
TAA5  <- set_up(files[5], "ATAA5", "ATAA")
TAA6  <- set_up(files[6], "ATAA6", "ATAA")
TAA7  <- set_up(files[7], "ATAA7", "ATAA")
TAA8  <- set_up(files[8], "ATAA8", "ATAA")

###############################combined all#################################################
anchors <- FindIntegrationAnchors(object.list = list(Con4,Con6,Con9,TAA2,TAA3,TAA4,TAA5,TAA6,TAA7,TAA8), 
                                  dims = 1:20, k.filter = 150) ### take long time ####

combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

combined <- RunTSNE(combined, reduction = "pca", dims = 1:20) ##### option1 #########
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.7)

p1 <- DimPlot(combined, reduction = "tsne", group.by = "stim")
p2 <- DimPlot(combined, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)
DimPlot(combined, reduction = "tsne", split.by = "stim", label = TRUE)

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "AllMarkers_nonImmune_celltype.csv", quote = F)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(combined, features = top10$gene) #+ NoLegend()
ggsave("Fig2_top10Heatmap_nonImmune_res0.7.pdf", width=20, height=24)

saveRDS(combined, file = "nonImmuneCombined.rds")
write.csv(combined@meta.data, "meta_nonImmune.csv", row.names = T, quote = F)

############################## visulization #####################################
FeaturePlot(combined, features = c("MYH11", "ACTA2", "ATF3", "TPM4",
                                   "COL8A1", "DCN", "COL1A2", "PDGFRB", 
                                   "PECAM1", "SPARC", "MT1M", "MYC"),
            cols=c("darkblue", "darkorchid4", "goldenrod", "firebrick"))
ggsave("Fig3_FeaturePlot2.pdf", width=15, height=10)

DotPlot(combined, features = c("ERG"), split.by = "stim", cols = c("blue", "blue")) 

########################## cluster correlation by avg #########################
avg <- AverageExpression(combined)
avg <- avg$RNA
write.csv(avg, "avg_SeuratCluster.csv", row.names = T, quote = F)

res <- cor(avg, method = "spearman")
pheatmap(res, main="Correlation_SeuratCluster")

############################ rename cluster #############################################
combined <- RenameIdents(combined, `0` = "Contractile SMC", `1` = "Fibroblast2", `2` = "Contractile SMC", 
                         `3` = "Stressed SMC", `4` = "Proliferating SMC2", `5` = "Fibromyocyte", `6` = "Fibroblast1", 
                         `7` = "EC1", `8` = "Proliferating SMC1", `9` = "MSC1", 
                         `10` = "MSC2", `11` = "EC2", `12` = "Inflammatory2", 
                         `13` = "Inflammatory1", `14` = "Inflammatory3", `15` = "Schwann")
combined@meta.data$celltype <- Idents(combined)
DimPlot(combined, reduction = "tsne", label = TRUE) #

################################## average expression heatmap #######################################
avg <- AverageExpression(combined)
avg <- avg$RNA

list <- read.table("marker_top10.txt", header=F, sep = "\t", stringsAsFactors = F)
list <- list[,1]
list <- unique(list)

expr <- avg[row.names(avg) %in% list,]
library(pheatmap)
library(RColorBrewer)
breaksList = seq(-3, 3, by = 0.1)

pheatmap(t(expr[match(list, row.names(expr)),c(1,3,5,4,8,2,6,13,12,14,9,10,15,7,11)]), 
         scale = "column", main = "marker Genes", cluster_cols = F, cluster_rows = F,
         border_color=NA,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)

############################# calculate Module Score #############################

files <- list.files(path=".", pattern = "^marker_*")
for (i in 1:(length(files)-1)){
  markers <- read.table(files[i], header=F, sep="\t", stringsAsFactors = F)
  markers <- markers[,1]
  feature_name <- strsplit(strsplit(files[i],"_")[[1]][2],"[.]")[[1]][1]
  combined <- AddModuleScore(combined, features=list(markers),ctrl = 5,
                             name = feature_name)
}
meta <- combined@meta.data
module_score <- c()
for (i in unique(meta$celltype)){
  each <- meta[meta$celltype==i,]
  score <- apply(each[,9:dim(each)[2]],2,median)
  module_score <- rbind(module_score, score)
}
row.names(module_score) <- unique(meta$celltype)
colnames(module_score) <- substr(colnames(module_score), 1, (nchar(colnames(module_score))-1))
library(pheatmap)
library(RColorBrewer)
breaksList = seq(-2, 2, by = 0.01)
pheatmap(t(module_score[c(1,3,4,13,6,10,9,2,12,7,14,15,11,5,8),]), main="Scaled Module Score", scale = "row",cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList) #
write.csv(module_score,"ModuleScore_MedianOfCelltype.csv", row.names = T, quote = F)

############################## junction analysis between cell types ##############################
avg <- avg[,c(1,3,5,4,8,2,6,13,12,14,9,10,15,7,11)]
cell.cell <- read.table("input_cell_cell_junction.txt", header = T, sep="\t", stringsAsFactors = F)
cell.ECM <- read.table("input_cell_ECM_junction.txt", header = T, sep="\t", stringsAsFactors = F)

exp_cell.cell <- avg[row.names(avg) %in% cell.cell$Gene,]
score_all.all <- c()
for (i in 1:dim(exp_cell.cell)[2]){
  for (j in 1:dim(exp_cell.cell)[2]) {
    each <- sum(exp_cell.cell[,i]*exp_cell.cell[,j])
    score_all.all <- c(score_all.all, each)
  }
}
score_all.all <- matrix(score_all.all, nrow=15)
row.names(score_all.all) <- colnames(exp_cell.cell)
colnames(score_all.all) <- colnames(exp_cell.cell)
write.csv(score_all.all, "JunctionScore_cellcell.csv", row.names = T, quote = F)

library(ggplot2)
library(reshape2)
data <- melt(score_all.all)
p1 <- ggplot(data, aes(x=Var1, y=Var2, size=value)) + geom_point()+ theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title =" Cell-cell junction score")



exp_cell.ECM <- avg[row.names(avg) %in% cell.ECM$Gene,]
score_cell.ECM <- colSums(exp_cell.ECM)
score_cell.ECM <- matrix(score_cell.ECM, nrow=15)
row.names(score_cell.ECM) <- colnames(exp_cell.ECM)
write.csv(score_cell.ECM, "JunctionScore_cellECM.csv", row.names = T, quote = F)

data <- melt(score_cell.ECM)
p2 <- ggplot(data, aes(x=Var2, y=Var1, size=value)) + geom_point()+ theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title =" Cell-ECM junction score")

p1+p2

############################## avg heatmap of interested genes #########################################
target <- read.table("Input_Cyclins.txt", header = T, sep="\t", stringsAsFactors = F)
target_expr <- avg[row.names(avg) %in% target$Gene,]
target_expr <- target_expr[rowSums(target_expr) > 0.5,]

target <- target[target$Gene %in% row.names(target_expr),]
mycol <- data.frame("Group"=target$Group,
                    row.names = target$Gene)

breaksList = seq(-3, 3, by = 0.01)
pheatmap(target_expr, scale = "row", annotation_row = mycol,main = "Cyclin genes",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)

############################# ERG logFC ############################################################
avg <- AverageExpression(combined, add.ident = "stim")
avg <- avg$RNA
data <- avg[row.names(avg)=="ERG",]

Cons <- grep("*_Control", colnames(data))
FC <- c()
for (i in Cons){
  names <- substr(colnames(data)[i], 1, nchar(colnames(data)[i])-8)
  matched <- grep(paste0(names,"_ATAA"), colnames(data))
  each <- cbind(log10(data[,matched]/data[,i]), names)
  colnames(each) <- c("logFC", "CellType")
  FC <- rbind(FC, each)
}
FC <- as.data.frame(FC)
FC <- FC[FC$logFC != Inf, ]

barplot(as.numeric(as.character(FC$logFC)), ylab="ERG logFC", 
        names.arg = FC$CellType, las=1, las=2)

############################## ERG correlated genes in perticular cells ###################################
par(mfrow=c(2,3))
data <- c()
for(i in c("EC1", "Fibroblast1", "Proliferating SMC1", "Proliferating SMC2",  "Inflammatory1")){
  smc_like <- SubsetData(combined, ident.use = i)
  expr <- smc_like@assays$RNA@data
  ERG <- expr[row.names(expr)=="ERG",]
  
  res <- cor(ERG, y=as.data.frame(t(expr)), method = "spearman")
  names(res) <- row.names(expr)
  res <- na.omit(res[1,])
  
  plot(sort(res, decreasing = T), pch=20, ylim=c(-0.4,0.5),
       ylab="Correlation with ERG", xlab="Genes", main=i)
  
  each <- rbind(cbind(head(sort(res, decreasing = T)),  i),
                cbind(tail(sort(res, decreasing = T)),  i))
  colnames(each) <- c("Cor_to_ERG", "Cell_type")
  data <- rbind(data, each)
  
}
write.csv(data, "Cor_to_ERG.csv", row.names = T, quote = F)
