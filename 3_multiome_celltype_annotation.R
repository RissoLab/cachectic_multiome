.libPaths(new = "/mnt/ganimede/zuin/R/x86_64-pc-linux-gnu-library/4.3")

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(scDblFinder)
library(scuttle)
library(scran)

# For reproducibility
set.seed(8569)

# Load multiome data
mul_comb <- readRDS("multiome_comb_dblcalls.rds")
# Select singlet
sgl_mult <- mult_comb[, mult_comb$gex_dbl_calls == "singlet"]
DefaultAssay(sgl_mult) <- "integratedRNA"

# Re-compute dimensional reduction
sgl_mult <- RunPCA(sgl_mult, npcs = 30, verbose = FALSE)
sgl_mult <- RunUMAP(sgl_mult, reduction = "pca", dims = 1:30)

# Clustering
sgl_mult <- FindNeighbors(sgl_mult, reduction = "pca", dims = 1:30)
sgl_mult <- FindClusters(sgl_mult, resolution = 0.8)

# Cell-type specific markers
markers <- read.table("cachectic_markers.txt", header = T)
markers <- markers$Gene

p5 <- DotPlot(object = sgl_mult, features = markers, scale=FALSE, assay="SCT")+ RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")
p6 <- DotPlot(object = sgl_mult, features = markers, assay="SCT")+ RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")

patch <- p5/p6
ggsave(patch, file=paste("multiome_dotplot.jpg",sep = ""), width = 20, height = 40, units = "cm")

p <- DimPlot(sgl_mult, reduction = "umap", group.by = "seurat_clusters",label = TRUE, repel = TRUE)+
  NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") 
ggsave(p, file=paste("multiome_UMAP_clusters.jpg",sep = ""), width = 20, height = 20, units = "cm")

DefaultAssay(sgl_mult) <- "SCT"
heat <- DoHeatmap(sgl_mult, features=markers, size = 3, angle = 90, draw.lines = TRUE,
                  lines.width=10) + NoLegend()+
  scale_fill_gradientn(colors = c("lightblue", "lightgrey", "red"))
ggsave(heat, file="multiome_heatmap_markers.jpg", width = 30, height = 30, units = "cm")

sgl_mult$clusters <- as.character(sgl_mult$seurat_clusters)

sgl_mult$clusters[c(sgl_mult$clusters=="0"|sgl_mult$clusters=="1"|sgl_mult$clusters=="2"|
                      sgl_mult$clusters=="4"|sgl_mult$clusters=="5"|
                      sgl_mult$clusters=="7"|sgl_mult$clusters=="9"|
                      sgl_mult$clusters=="10"|sgl_mult$clusters=="11"|
                      sgl_mult$clusters=="17"|sgl_mult$clusters=="22")] <- "Myonuclei"

sgl_mult$clusters[c(sgl_mult$clusters=="3"|sgl_mult$clusters=="12"|
                      sgl_mult$clusters=="14"|sgl_mult$clusters=="20")] <- "FAPs"

sgl_mult$clusters[c(sgl_mult$clusters=="6"|sgl_mult$clusters=="23")] <- "Smooth_Muscular"

sgl_mult$clusters[c(sgl_mult$clusters=="8"|sgl_mult$clusters=="19"|
                      sgl_mult$clusters=="21")] <- "Endothelial"

sgl_mult$clusters[sgl_mult$clusters=="13"|sgl_mult$clusters=="26"] <- "Immune_Cells"

sgl_mult$clusters[sgl_mult$clusters=="16"|sgl_mult$clusters=="18"] <- "Tenocyte"
sgl_mult$clusters[sgl_mult$clusters=="15"] <- "MuSC"
sgl_mult$clusters[sgl_mult$clusters=="25"] <- "Nervous_System"
sgl_mult$clusters[sgl_mult$clusters=="24"] <- "Pericyte"

load("SubclusterPaletteColor.RData")

cluster_palette <- subcluster_palette[-c(6:10),]
cluster_palette$Label[cluster_palette$Label=="Myonuclei_IIb"] <- "Myonuclei"
palette_color <- cluster_palette$Color
names(palette_color) <- cluster_palette$Label

p <- DimPlot(sgl_mult, reduction = "umap", group.by = "clusters",label = TRUE, repel = TRUE)+
  NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") + scale_color_manual(values = palette_color)
ggsave(p, file=paste("multiome_UMAP_clusters.jpg",sep = ""), width = 20, height = 20, units = "cm")

# Re-clustering on Myonuclei ####
myo_mult <- sgl_mult[, sgl_mult$clusters=="Myonuclei"]

DefaultAssay(myo_mult) <- "integratedRNA"
# Run the standard workflow for visualization and clustering
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
ggsave(heat, file="multiome_myo_heatmap_atrophy.jpg", width = 30, height = 30, units = "cm")

myo_markers <- c("Ttn", "Col22a1", "Etv5", "Ano4", "Chrne","Myh1", "Myh2", "Myh4")

heat <- DoHeatmap(myo_mult, features=myo_markers, size = 3, angle = 90, draw.lines = TRUE,
                  lines.width=10) + NoLegend()+
  scale_fill_gradientn(colors = c("lightblue", "lightgrey", "red"))
ggsave(heat, file="multiome_myo_heatmap_markers.jpg", width = 30, height = 30, units = "cm")

p4 <- DotPlot(object = myo_mult, features = myo_markers, scale=FALSE, assay="SCT")+ RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")
p5 <- DotPlot(object = myo_mult, features = myo_markers, assay="SCT")+ RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")
patch <- p4/p5
ggsave(patch, file="multiome_myo_dotplot.jpg", width = 20, height = 40, units = "cm")

p <- DimPlot(myo_mult, reduction = "umap", group.by = "seurat_clusters",label = TRUE, repel = TRUE)+
  NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") 
ggsave(p, file=paste("multiome_myo_UMAP_clusters.jpg",sep = ""), width = 20, height = 20, units = "cm")

myo_mult$subcluster <- as.character(myo_mult$seurat_clusters)

myo_mult$subcluster[c(myo_mult$subcluster=="0"|myo_mult$subcluster=="1"|
                        myo_mult$subcluster=="4"|myo_mult$subcluster=="5"|myo_mult$subcluster=="12"|
                        myo_mult$subcluster=="16")] <- "Myonuclei_IIb"

myo_mult$subcluster[c(myo_mult$subcluster=="2"|myo_mult$subcluster=="7"|
                        myo_mult$subcluster=="11"|myo_mult$subcluster=="15")] <- "Myonuclei_IIx"

myo_mult$subcluster[c(myo_mult$subcluster=="3")] <- "Myonuclei_IIx_IIa"

myo_mult$subcluster[c(myo_mult$subcluster=="8"|myo_mult$subcluster=="9")] <- "Myonuclei_IIx_IIb"

myo_mult$subcluster[c(myo_mult$subcluster=="6")] <- "Myonuclei_Trim63" 

myo_mult$subcluster[c(myo_mult$subcluster=="10"|myo_mult$subcluster=="13")] <- "Myonuclei_MTJ"

myo_mult$subcluster[c(myo_mult$subcluster=="14")] <- "Myonuclei_NMJ" 

save(myo_mult, file = "multiome_myo_annot.RData")

sgl_mult$subclusters <- sgl_mult$clusters
sgl_mult$subclusters[sgl_mult$subclusters=="Myonuclei"] <- myo_mult$subcluster

saveRDS(sgl_mult, file = "multiome_sgl_comb.rds")

load("SubclusterPaletteColor.RData")
label_color <- subcluster_palette
label_color$Label[label_color$Label=="Myonuclei_MTJ_IIb"] <- "Myonuclei_MTJ"

add_color <- c( "#e19fc0")
names(add_color) <- c("Myonuclei_Trim63")

palette_color <- label_color$Color
names(palette_color) <- label_color$Label
palette_color <- c(palette_color, add_color)

p <- DimPlot(sgl_mult, reduction = "umap", label = TRUE, repel = TRUE, group.by="subclusters") + 
  ggtitle("") + NoLegend() + scale_color_manual(values=palette_color)
ggsave(p, file="multiome_UMAP_subclusters.jpg", width = 20, height = 20, units = "cm")

# UMAP plot per replicate ####
DefaultAssay(sgl_mult) <- "RNA"

ll <- lapply(unique(sgl_mult$dataset), function(x) sgl_mult[, sgl_mult$dataset==x])
names(ll) <- unique(sgl_mult$dataset)

ll <- lapply(ll, FUN = function(x) {
  x <- SCTransform(x)
})

ll <- lapply(ll, function(x) RunPCA(x, npcs = 30, verbose = FALSE))
ll <- lapply(ll, function(x) RunUMAP(x, reduction = "pca", dims = 1:30))

for (i in seq_along(ll)) {
  p <-  DimPlot(ll[[i]], reduction = "umap", label = TRUE, repel = TRUE, group.by="subclusters") + 
    ggtitle(names(ll)[i]) + NoLegend() + scale_color_manual(values=palette_color)
  ggsave(p, file=paste("multiome_UMAP_subclusters_",names(ll)[i],".jpg", sep = ""), width = 20, height = 20, units = "cm")
}

# UMAP plot per condition ####
DefaultAssay(sgl_mult) <- "RNA"

ll <- lapply(unique(sgl_mult$condition), function(x) sgl_mult[, sgl_mult$condition==x])
names(ll) <- unique(sgl_mult$condition)

ll <- lapply(ll, FUN = function(x) {
  x <- SCTransform(x)
})

ll <- lapply(ll, function(x) RunPCA(x, npcs = 30, verbose = FALSE))
ll <- lapply(ll, function(x) RunUMAP(x, reduction = "pca", dims = 1:30))

load("SubclusterPaletteColor.RData")
label_color <- subcluster_palette
label_color$Label[label_color$Label=="Myonuclei_MTJ_IIb"] <- "Myonuclei_MTJ"

add_color <- c( "#e19fc0")
names(add_color) <- c("Myonuclei_Trim63")

for (i in seq_along(ll)) {
  p <-  DimPlot(ll[[i]], reduction = "umap", label = TRUE, repel = TRUE, group.by="subclusters") + 
    ggtitle(names(ll)[i]) + NoLegend() + scale_color_manual(values=palette_color)
  ggsave(p, file=paste("multiome_UMAP_subclusters_",names(ll)[i],".jpg", sep = ""), width = 20, height = 20, units = "cm")
}

# UMAP colored by condition ####
sgl_mult$condition[sgl_mult$condition == "SHAM"] <- "CTRL"
sgl_mult$condition[sgl_mult$condition == "CAC"] <- "C26"
p <- DimPlot(sgl_mult, reduction = "UMAP.wnn", group.by = "condition") +
  NoLegend()+ labs(x = "WNN_UMAP1", y = "WNN_UMAP2") + ggtitle("") + scale_color_manual(values = c("CTRL" = "blue", "C26" = "red"))
ggsave(p, file=paste("multiome_umap_wnn_condition.pdf",sep = ""), width = 20, height = 20, units = "cm")
