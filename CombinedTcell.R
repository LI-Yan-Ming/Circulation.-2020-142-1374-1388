rm(list=ls())
setwd("H:\\AHA_scRNAseq\\human\\TAA_Con_8vs3\\Tcell\\")

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(DESeq2)
library(pheatmap)

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
                                  dims = 1:20, k.filter = 90) ### take long time ####

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
write.csv(markers, "AllMarkers_Tcell_seuratCluster.csv", quote = F)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(combined, features = top10$gene) # + NoLegend()
ggsave("Fig2_top10Heatmap_Tcell_seuratCluster.pdf", width=14, height=20)

saveRDS(combined, file = "Tcell_combined.rds")
write.csv(combined@meta.data, "meta_Tcell.csv", row.names = T, quote = F)

############################## visulization #####################################
DefaultAssay(combined) <- "integrated"
FeaturePlot(combined, features = c("CD4", "CD8A", "CD8B", "HOPX",
                                   "HSPA1A", "TXNIP", "GIMAP4", "ISG15", 
                                   "STMN1", "C1QB", "KLF6", "DNAJA1"),
            cols=c("darkblue", "darkorchid4", "goldenrod", "firebrick"))
ggsave("Fig3_FeaturePlot2.pdf", width=15, height=10)

########################## cluster correlation by avg #########################
avg <- AverageExpression(combined)
avg <- avg$RNA
write.csv(avg, "avg_celltype.csv", row.names = T, quote = F)

res <- cor(avg[row.names(avg) %in% VariableFeatures(combined),], method = "spearman") #
pheatmap(res, main="Correlation_celltype")

############################### rename cluster #######################################
combined <- RenameIdents(combined, `0` = "CD4_active", `1` = "CD8_active", `2` = "CD4_EM_rest", 
                         `3` = "T_GIMAP", `4` = "T_HSP", `5` = "CD8_TEMRA", `6` = "T_stress", 
                         `7` = "CD4_Treg", `8` = "CD8_TEMRA", `9` = "T_proliferation",
                         `10` = "T_IFNresponse", `11` = "T_switched")
combined@meta.data$celltype <- Idents(combined)
DimPlot(combined, reduction = "tsne") #, label = TRUE


######################## avg expression heatmap ############################
list <- c("CD4","CD8A","CD8B",
          "CREM","CXCR6", "RGCC","NR3C1","GZMB",
          "CCR7", "IL7R", "CCL20", "KLRB1",
          "IL2RA", "TNFRSF18", "ID3", "LTB",  "CTLA4",
          "CCL4","GZMK", "CMC1", "CRTAM", 
          "GNLY", "HOPX", "KLRD1", "CCL5", "PRF1",
          "TXNIP","MALAT1", "GIMAP1","GIMAP4","GIMAP7",
          "HSPA1A","HSPA1B", "HSP90AA1","JUN",  
          "KLF6","DNAJA1", "FOS","TUBB4B",
          "ISG15", "LY6E", "ISG20", "MT2A", "IFI44","HERC5",
          "STMN1", "TUBA1B", "HMGN2", "CKS1B", "HMGB1",
           "TAGLN", "MYH11", "TPM1")
expr <- avg[row.names(avg) %in% list,]
expr <- expr[match(list, row.names(expr)),c(1,3,8,2,6,4,5,7,10,9,11)]

library(pheatmap)
library(RColorBrewer)
breaksList = seq(-3, 3, by = 0.1)
pheatmap(t(expr), scale = "column", cluster_cols = F, cluster_rows = F,
         border_color = NA, main="Marker gene",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)

##
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
mycol <- data.frame("Role"= list$role,"Group" = list$group,
                    row.names = list$gene)#
breaksList = seq(-3, 3, by = 0.1)
p <- pheatmap(t(expr), 
              scale = "column", main = "Cytokine and receptors", cluster_cols = T, cluster_rows = T,
              show_colnames = T,annotation_col= mycol, border_color = F,
              color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
              breaks = breaksList)
ggsave(p, filename = "Fig5_heatmap_cytokines_showname.pdf", width = 6,height = 18)

############################## CD4/CD8 each cell ######################################
data <- combined@assays$RNA@data
target <- data[row.names(data) %in% c("CD4", "CD8A"),]
target <- target+0.00001
target <- rbind(target,log2(target[2,]/target[1,]))

coordi <- as.data.frame(combined@reductions$tsne@cell.embeddings)
coordi$'log2(CD4/CD8A)' <- target[3,]

coordi$'SurfaceMarker'[coordi$'log2(CD4/CD8A)'<0] <- "CD8A"
coordi$'SurfaceMarker'[coordi$'log2(CD4/CD8A)'>0] <- "CD4"
coordi$'SurfaceMarker'[coordi$'log2(CD4/CD8A)'==0] <- "Not detected"
coordi$'SurfaceMarker' <- as.factor(coordi$'SurfaceMarker')

ggplot(coordi, aes(x=tSNE_1, y=tSNE_2, colour=SurfaceMarker)) + geom_point(size=0.2) +
   theme_classic()+scale_colour_manual(values = c("red", "blue", "grey")) 

