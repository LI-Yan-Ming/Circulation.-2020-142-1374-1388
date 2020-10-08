rm(list=ls())
setwd("H:\\AHA_scRNAseq\\human\\TAA_Con_8vs3\\AllCells\\")

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(DESeq2)
library(pheatmap)

combined <- readRDS("Combined_Con469TAA1to8.rds")

#################### set up seperate seurat projects #################################
files <- list.files(path=".", pattern = "*res0.6.rds")
set_up <- function(single_project,project,stim){
  data <- readRDS(single_project)
  UMI <- as.matrix(data@assays$RNA@counts)
  each <- CreateSeuratObject(counts = UMI, project = project, min.cells = 5)
  each$stim <- stim
  each <- NormalizeData(each, verbose = FALSE)
  each <- FindVariableFeatures(each, selection.method = "vst", nfeatures = 2000)
}

Con4  <- set_up(files[1], "Con4", "Control")
Con6  <- set_up(files[2], "Con6", "Control")
Con9  <- set_up(files[3], "Con9", "Control")
TAA1  <- set_up(files[4], "ATAA1", "ATAA")
TAA2  <- set_up(files[5], "ATAA2", "ATAA")
TAA3  <- set_up(files[6], "ATAA3", "ATAA")
TAA4  <- set_up(files[7], "ATAA4", "ATAA")
TAA5  <- set_up(files[8], "ATAA5", "ATAA")
TAA6  <- set_up(files[9], "ATAA6", "ATAA")
TAA7  <- set_up(files[10], "ATAA7", "ATAA")
TAA8  <- set_up(files[11], "ATAA8", "ATAA")

###############################combined all#################################################
anchors <- FindIntegrationAnchors(object.list = list(Con4,Con6,Con9,TAA1,TAA2,TAA3,TAA4,TAA5,TAA6,TAA7,TAA8), dims = 1:20) ### take long time ####
combined <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

combined <- RunTSNE(combined, reduction = "pca", dims = 1:20) ##### option1 #########
# combined <- RunUMAP(combined, reduction = "pca", dims = 1:20) ##### option2 #########
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.6)

p1 <- DimPlot(combined, reduction = "tsne", group.by = "stim")
p2 <- DimPlot(combined, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)
DimPlot(combined, reduction = "tsne", split.by = "stim")
ggsave("Fig1_tSNE_separate.pdf", width=30, height=4)

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "AllMarkers_combined_Con469TAA1to8_rename.csv", quote = F)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(combined, features = top10$gene) + NoLegend()
ggsave("Fig2_top10Heatmap_Con469TAA1to8_rename.pdf", width=40, height=40)

saveRDS(combined, file = "Combined_Con469TAA1to8.rds")

############################## visulization #####################################
FeaturePlot(combined, features = c("MYH11", "DCN", "PECAM1", "MT1M",
                                   "COL1A2", "PTPRC", "LYZ", "CD3D", 
                                   "KLRD1", "CD79A", "MZB1", "CPA3"),
            cols=c("darkblue", "goldenrod","firebrick"))
ggsave("Fig3_FeaturePlot1.pdf", width=16, height=10)


########################## cluster correlation by avg #########################
avg <- AverageExpression(combined)
avg <- avg$RNA
write.csv(avg, "avg_celltype.csv", row.names = T, quote = F)

res <- cor(avg, method = "spearman")
pheatmap(res, main="Correlation_CellType")

############################ rename cluster #############################################
combined <- RenameIdents(combined, `0` = "Tcell", `1` = "Tcell", `2` = "SMC1", 
                         `3` = "Tcell", `4` = "MonoMaphDC", `5` = "Tcell", `6` = "Fibroblast", 
                         `7` = "MonoMaphDC", `8` = "NK", `9` = "MonoMaphDC", 
                         `10` = "MonoMaphDC", `11` = "SMC2", `12` = "MonoMaphDC",
                         `13` = "MSC", `14` = "EC", `15` = "Plasma", `16` = "Tcell", 
                         `17` = "Mastcell", `18` = "Bcell", `19` = "MonoMaphDC", 
                         `20` = "MonoMaphDC", `21` = "MonoMaphDC", `22` = "Tcell",
                         `23` = "MonoMaphDC")
combined@meta.data$celltype <- Idents(combined)
DimPlot(combined, reduction = "tsne") #, label = TRUE

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

pheatmap(t(expr[match(list, row.names(expr)),c(1,3,2,5,4,6,8,7,9,11,10)]), 
         scale = "column", main = "marker Genes", cluster_cols = F, cluster_rows = F,
         border_color=NA,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)

###################### each type, their composition by stim ##############################
meta <- combined@meta.data 

actual <- c()
for (i in unique(meta$celltype)){
  each <- meta[meta$celltype==i,]
  for (j in unique(meta$stim)){
    sub <- each[each$stim==j,]
    frac <- dim(sub)[1]/dim(each)[1]
    frac <- cbind(i,j,frac)
    actual <- rbind(actual,frac)
  }
}
expected <- dim(meta[meta$stim=="Control",])[1]/dim(meta)[1]

colnames(actual) <- c("CellType", "Group",  "Fraction")
actual <- as.data.frame(actual)
actual$Fraction <- as.numeric(as.character(actual$Fraction))
ggplot(actual, aes(x=CellType,y=Fraction, fill=Group))  + 
  geom_bar(stat="identity")+ coord_flip() +
  geom_hline(yintercept=expected, linetype="dashed", color = "black") +
  theme_classic(base_size = 14.5)

###################### each individual, their composition by celltype ##############################
meta <- combined@meta.data 

actual <- c()
for (i in unique(meta$orig.ident)){
  each <- meta[meta$orig.ident==i,]
  for (j in unique(meta$celltype)){
    sub <- each[each$celltype==j,]
    frac <- dim(sub)[1]/dim(each)[1]
    frac <- cbind(i,j,frac)
    actual <- rbind(actual,frac)
  }
}

colnames(actual) <- c("Subject", "CellType",  "Fraction")
actual$CellType <- ifelse(actual$CellType %in% c("SMC1", "SMC2", "Fibroblast", "MSC", "EC"), "Non-Immune", "Immune")
actual <- as.data.frame(actual)
actual$Fraction <- as.numeric(as.character(actual$Fraction))
ggplot(actual, aes(x=Subject,y=Fraction, fill=CellType))  + 
  geom_bar(stat="identity") + coord_flip()+
  theme_classic()

#################### extract particular type of cells by subject ##########################
Iwant <- subset(combined, idents=c("MonoMaphDC"))

UMI <- Iwant@assays$RNA@counts
meta <- Iwant@meta.data

for(i in unique(meta$orig.ident)){
  each <- meta[meta$orig.ident==i,]
  UMI_each <- UMI[,colnames(UMI) %in% row.names(each)]
  write.csv(UMI_each, paste0("UMI_", i,"_MonoMaphDC_fromIntegrative.csv"), quote = F, row.names = T)
}

####################################### add subclusters ##############################################
meta <- combined@meta.data
meta_nonImm <- read.csv("meta_nonImmune.csv", header=T, row.names = "X", stringsAsFactors = F)
meta_maph <- read.csv("meta_MonoMaphDC.csv", header=T, row.names = "X", stringsAsFactors = F)
meta_Tcell <- read.csv("meta_Tcell.csv", header=T, row.names = "X", stringsAsFactors = F)
meta_sub <- rbind(meta_nonImm[,c(1:4,8)], meta_maph[,c(1:4,8)], meta_Tcell[,c(1:4,8)])

celltype2 <- c()
for(i in 1:dim(meta)[1]){
  each <- ifelse(row.names(meta)[i] %in% row.names(meta_sub),
                 meta_sub[row.names(meta_sub)==row.names(meta)[i],5],
                 as.character(meta$celltype[i]))
  celltype2 <- c(celltype2,each)
}
meta$celltype2 <- celltype2
Idents(combined) <- meta$celltype2
combined@meta.data$celltype2 <- celltype2
#DimPlot(combined, reduction = "tsne")
Iwant <- subset(combined, idents=unique(combined@meta.data$celltype2)[1:40]) # remove TAA1 nonImmune cells cause they did not submit to subcluster analysis

##overall correlation of subclusters
avg <- AverageExpression(Iwant)#, add.ident = "stim"
avg <- avg$RNA
write.csv(avg, "avg_addsubcluster.csv", row.names = T, quote = F)

data <- avg[row.names(avg) %in% VariableFeatures(Iwant),]
res <- cor(data, method = "spearman")
library(pheatmap)
pheatmap(res, main="Correlation_clusters")

## save avg by stim
avg <- AverageExpression(Iwant, add.ident = "stim")#
avg <- avg$RNA
write.csv(avg, "avg_addsubcluster_stim.csv", row.names = T, quote = F)

saveRDS(Iwant, file = "Combined_Con469TAA1to8_addSubCluster.rds")

################################## extract UMI count by celltype2 for DEG analysis ##################################
UMI <- Iwant@assays$RNA@counts
meta <- Iwant@meta.data
for (i in unique(meta$celltype2)){
  sub_meta <- meta[meta$celltype2==i,]
  sub_UMI <- UMI[,colnames(UMI) %in% row.names(sub_meta)]
  keep <- rowSums(sub_UMI)>5
  sub_UMI <- sub_UMI[keep,]
  write.csv(sub_UMI, paste0("UMI_subcluster_",i,".csv"), quote = F, row.names = T)
}
write.csv(meta, "meta_addsubcluster.csv", quote = F, row.names = T)


############### dimplot, selected or not selected ###################
meta <- combined@meta.data
selected <- ifelse(meta$celltype %in% c("MonoMaphDC"), "yes", "no")
Idents(combined) <- selected
DimPlot(combined, cols = c("grey", "blue")) + ggtitle("MonoMaphDC cells")

