# =============================================================================
# Cachectic multiome (C26 and CTRL samples)
# Cell types annotation
# =============================================================================

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(scDblFinder)
library(scuttle)
library(scran)

set.seed(8569)

mul_comb <- readRDS("cachectic_mult_ctrl_c26.rds")

# Clustering
mul_comb <- FindNeighbors(mul_comb, reduction = "pca", dims = 1:30)
mul_comb <- FindClusters(mul_comb, resolution = 0.8)

# Markers useful for cluster annotation
markers <- read.table("cachectic_markers.txt", header = T)$Gene

p1 <- DotPlot(object = mul_comb, features = markers, scale = FALSE, assay = "SCT") + 
  RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")

p2 <- DotPlot(object = mul_comb, features = markers, assay = "SCT") + 
  RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")

patch <- p1/p2
ggsave(patch, file = "cachectic_dotplot_clust_markers.pdf", width = 20, height = 40, units = "cm")

mul_comb$clusters <- as.character(mul_comb$seurat_clusters)
mul_comb$clusters[c(mul_comb$clusters == "0"|mul_comb$clusters == "1"|mul_comb$clusters == "2"|
                      mul_comb$clusters == "4"|mul_comb$clusters == "5"|
                      mul_comb$clusters == "7"|mul_comb$clusters == "9"|
                      mul_comb$clusters == "10"|mul_comb$clusters == "11"|
                      mul_comb$clusters == "17"|mul_comb$clusters == "22")] <- "Myonuclei"

mul_comb$clusters[c(mul_comb$clusters == "3"|mul_comb$clusters == "12"|
                      mul_comb$clusters == "14"|mul_comb$clusters == "20")] <- "FAPs"

mul_comb$clusters[c(mul_comb$clusters == "6"|mul_comb$clusters == "23")] <- "Smooth_Muscular"

mul_comb$clusters[c(mul_comb$clusters == "8"|mul_comb$clusters == "19"|
                      mul_comb$clusters == "21")] <- "Endothelial"

mul_comb$clusters[mul_comb$clusters == "13"|mul_comb$clusters == "26"] <- "Immune_Cells"

mul_comb$clusters[mul_comb$clusters == "16"|mul_comb$clusters == "18"] <- "Tenocyte"
mul_comb$clusters[mul_comb$clusters == "15"] <- "MuSC"
mul_comb$clusters[mul_comb$clusters == "25"] <- "Nervous_System"
mul_comb$clusters[mul_comb$clusters == "24"] <- "Pericyte"

# Visualization
load("palette_color_ct.RData")
pltt_color <- palette_color[-c(6:11)]
names(pltt_color)[5] <- "Myonuclei"

p <- DimPlot(mul_comb, reduction = "umap", group.by = "clusters",label = TRUE, repel = TRUE) +
  NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") + scale_color_manual(values = palette_color)
ggsave(p, file = "cachetic_mult_umap_clu.pdf", width = 20, height = 20, units = "cm")

# Re-clustering on Myonuclei ####
myo_mult <- mul_comb[, mul_comb$clusters =="Myonuclei"]
DefaultAssay(myo_mult) <- "integratedRNA"

myo_mult <- RunPCA(myo_mult, npcs = 30, verbose = FALSE)
myo_mult <- FindNeighbors(myo_mult, reduction = "pca", dims = 1:30)
myo_mult <- FindClusters(myo_mult, resolution = 0.8)

# Find top genes in each cluster
myo_markers <- FindAllMarkers(myo_mult, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DefaultAssay(myo_mult) <- "SCT"
atrophy_genes <- c("Trim63", "Fbxo32", "Cdkal1", "Mybpc1", "Pde4d", "Lama2", "Nctc1", "Xirp2")

heat <- DoHeatmap(myo_mult, features=atrophy_genes, size = 3, angle = 90, draw.lines = TRUE,
                  lines.width=10) + NoLegend()+
  scale_fill_gradientn(colors = c("lightblue", "lightgrey", "red"))
ggsave(heat, file="cachectic_mult_myo_heatmap_atrophy.jpg", width = 30, height = 30, units = "cm")

myo_markers <- c("Ttn", "Col22a1", "Etv5", "Ano4", "Chrne","Myh1", "Myh2", "Myh4")

heat <- DoHeatmap(myo_mult, features=myo_markers, size = 3, angle = 90, draw.lines = TRUE,
                  lines.width=10) + NoLegend()+
  scale_fill_gradientn(colors = c("lightblue", "lightgrey", "red"))
ggsave(heat, file="cachectic_mult_myo_heatmap_markers.jpg", width = 30, height = 30, units = "cm")

p1 <- DotPlot(object = myo_mult, features = myo_markers, scale=FALSE, assay="SCT")+ RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")
p2 <- DotPlot(object = myo_mult, features = myo_markers, assay="SCT")+ RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")
patch <- p1/p2
ggsave(patch, file="multiome_myo_dotplot.jpg", width = 20, height = 40, units = "cm")

myo_mult$subcluster <- as.character(myo_mult$seurat_clusters)

myo_mult$subcluster[c(myo_mult$subcluster == "0"|myo_mult$subcluster == "1"|
                        myo_mult$subcluster == "4"|myo_mult$subcluster == "5"|myo_mult$subcluster == "12"|
                        myo_mult$subcluster == "16")] <- "Myonuclei_IIb"

myo_mult$subcluster[c(myo_mult$subcluster == "2"|myo_mult$subcluster == "7"|
                        myo_mult$subcluster == "11"|myo_mult$subcluster == "15")] <- "Myonuclei_IIx"

myo_mult$subcluster[c(myo_mult$subcluster == "3")] <- "Myonuclei_IIx_IIa"

myo_mult$subcluster[c(myo_mult$subcluster == "8"|myo_mult$subcluster == "9")] <- "Myonuclei_IIx_IIb"

myo_mult$subcluster[c(myo_mult$subcluster == "6")] <- "Myonuclei_Trim63" 

myo_mult$subcluster[c(myo_mult$subcluster == "10"|myo_mult$subcluster == "13")] <- "Myonuclei_MTJ"

myo_mult$subcluster[c(myo_mult$subcluster == "14")] <- "Myonuclei_NMJ" 

mul_comb$subclusters <- mul_comb$clusters
mul_comb$subclusters[mul_comb$subclusters =="Myonuclei"] <- myo_mult$subcluster

p <- DimPlot(mul_comb, reduction = "umap", label = TRUE, repel = TRUE, group.by = "subclusters") + 
  ggtitle("") + NoLegend() + scale_color_manual(values = palette_color)
ggsave(p, file = "cachectic_mult_umap_ct.pdf", width = 20, height = 20, units = "cm")

p <- DimPlot(mul_comb, reduction = "UMAP.wnn", group.by = "condition") + ggtitle("") +
  scale_color_manual(values = c("CTRL" = "blue", "C26" = "red"))
ggsave(p, file = "cachectic_mult_umap_cond.pdf", width = 20, height = 20, units = "cm")

saveRDS(mult_comb, file = "cachectic_mult_ctrl_c26.rds")