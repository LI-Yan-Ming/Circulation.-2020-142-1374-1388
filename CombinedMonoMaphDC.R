rm(list=ls())
setwd("H:\\AHA_scRNAseq\\human\\TAA_Con_8vs3\\Maph\\")

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

combined <- readRDS( "MonoMaphDC_combined.rds")

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
anchors <- FindIntegrationAnchors(object.list = list(Con4,Con6,Con9,TAA1,TAA2,TAA3,TAA4,TAA5,TAA6,TAA7,TAA8), 
                                  dims = 1:20, k.filter = 200) ### take long time ####

combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

combined <- RunTSNE(combined, reduction = "pca", dims = 1:20) ##### option1 #########
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.4)

p1 <- DimPlot(combined, reduction = "tsne", group.by = "stim")
p2 <- DimPlot(combined, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)
DimPlot(combined, reduction = "tsne", split.by = "stim", label = TRUE)
ggsave("Fig1_tSNE_separate.pdf", width=26, height=4)

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "AllMarkers_MonoMaphDC_rename.csv", quote = F)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(combined, features = top10$gene) # + NoLegend()
ggsave("Fig2_top10Heatmap_MonoMaphDC_rename.pdf", width=14, height=20)

saveRDS(combined, file = "MonoMaphDC_combined.rds")
write.csv(combined@meta.data, "meta_MonoMaphDC.csv", row.names = T, quote = F)

############################## visulization #####################################
FeaturePlot(combined, features = c("CST3", "NFKB1", "IL1B", "TNF",
                                   "CD163", "STAB1", "MAF", "AREG", 
                                   "CD69", "IFI44L", "TAGLN", "TUBB"),
            cols=c("darkblue", "darkorchid4", "goldenrod", "firebrick"))
ggsave("Fig3_FeaturePlot3.pdf", width=15, height=10)

########################## cluster correlation by avg #########################
avg <- AverageExpression(combined)
avg <- avg$RNA
write.csv(avg, "avg_celltype.csv", row.names = T, quote = F)

res <- cor(avg[row.names(avg) %in% VariableFeatures(combined),], method = "spearman") #
pheatmap(res, main="Correlation_celltype")

############################ rename cluster #############################
combined <- RenameIdents(combined, `0` = "M1like1", `1` = "M2like1", `2` = "M1like1", 
                         `3` = "M1like1", `4` = "M2like2", `5` = "M2like1", `6` = "M1like2", 
                         `7` = "M1like2", `8` = "M1like3", `9` = "M_IFNresponse", 
                         `10` = "Monocyte", `11` = "M_Proliferating", `12` = "DC",`13` = "M_remodelling")
combined@meta.data$celltype <- Idents(combined)
DimPlot(combined, reduction = "tsne") #, label = TRUE

############################# avg expression heatmap ###################################
files <- list.files(path=".", pattern = "^Input_*")
list <- c()
role <- rep(c("Ligand", "Receptor"),4)
group <- c("Chemokine", "Chemokine","Interferon", "Interferon",
           "Interleukin", "Interleukin","TNFSF", "TNFSF")
for (i in 1:8){
  each <- read.table(files[i], header=F, sep = "\t", stringsAsFactors = F)
  each <- cbind(each,group[i], role[i])
  colnames(each) <- c("gene", "group", "role")
  list <- rbind(list, each)
}
expr <- avg[row.names(avg) %in% list$gene,]
expr <- expr[rowSums(expr) > 0.5,]

list <- list[list$gene %in% row.names(expr),]
mycol <- data.frame("Role"= list$role,
                    
                    row.names = list$gene)#"Group" = list$group,
library(pheatmap)
library(RColorBrewer)
breaksList = seq(-2, 2, by = 0.1)
p <- pheatmap(expr, 
              scale = "row", main = "Cytokine and receptors", cluster_cols = T, cluster_rows = T,
              show_rownames = F,annotation_row = mycol, border_color = F,
              color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
              breaks = breaksList)
ggsave(p, filename = "Fig5_heatmap_cytokines_showname.pdf", width = 6,height = 18)

##
list <- read.table("marker_proteases.txt", header=F, sep = "\t", stringsAsFactors = F)
list <- list[,1]

expr <- avg[row.names(avg) %in% list,]
expr <- expr[rowSums(expr)>0.5,]
library(pheatmap)
library(RColorBrewer)
breaksList = seq(-2, 2, by = 0.1)
pheatmap(expr, 
         scale = "row", main = "Protease genes", cluster_cols = T, cluster_rows = T,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)

##
list <- c("NFKB1", "IL1B", "TNF", "CD163", "STAB1","MRC1", "MERTK") #
expr <- avg[row.names(avg) %in% list,]
expr <- expr[match(list, row.names(expr)), sort(colnames(expr))]
breaksList = seq(-2, 2, by = 0.1)
pheatmap(t(expr), 
         scale = "column", main = "Marker genes", cluster_cols = F, cluster_rows = F,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)

################################## average top10 expression heatmap #######################################
avg <- read.csv("avg_celltype.csv", header=T, row.names = "X", stringsAsFactors = F)

list <- read.table("marker_top10.txt", header=F, sep = "\t", stringsAsFactors = F)
list <- list[,1]
list <- unique(list)

expr <- avg[row.names(avg) %in% list,]
library(pheatmap)
library(RColorBrewer)
breaksList = seq(-2.5, 2.5, by = 0.1)

pheatmap(t(expr[match(list, row.names(expr)),]), 
         scale = "column", main = "marker Genes", cluster_cols = F, cluster_rows = F,
         border_color=NA,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)