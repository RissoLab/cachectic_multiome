.libPaths(new = "/mnt/ganimede/zuin/R/x86_64-pc-linux-gnu-library/4.3")

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(scDblFinder)
library(scuttle)
library(scran)

# for reproducibility
set.seed(9212)

# Load data
rna_comb <- readRDS("rna_comb_dblcalls.rds")
# Select singlet
sgl_rna <- rna_comb[, rna_comb$gex_dbl_calls == "singlet"]
DefaultAssay(sgl_rna) <- "integratedRNA"

sgl_rna <- RunPCA(sgl_rna, npcs = 30, verbose = FALSE)
sgl_rna <- RunUMAP(sgl_rna, reduction = "pca", dims = 1:30)

# Clustering
sgl_rna <- FindNeighbors(sgl_rna, reduction = "pca", dims = 1:30)
sgl_rna <- FindClusters(sgl_rna, resolution = 0.8)

# Cell-type specific markers
markers <- read.table("cachectic_markers.txt", header = T)
markers <- markers$Gene

p5 <- DotPlot(object = sgl_rna, features = markers, scale=FALSE, assay="SCT")+ RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")
p6 <- DotPlot(object = sgl_rna, features = markers, assay="SCT")+ RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")

patch <- p5/p6
ggsave(patch, file=paste("rna_dotplot.jpg",sep = ""), width = 20, height = 40, units = "cm")

p <- DimPlot(sgl_rna, reduction = "umap", group.by = "seurat_clusters",label = TRUE, repel = TRUE)+
  NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") 
ggsave(p, file=paste("rna_UMAP_clusters.jpg",sep = ""), width = 20, height = 20, units = "cm")

DefaultAssay(sgl_rna) <- "SCT"
heat <- DoHeatmap(sgl_rna, features=markers, size = 3, angle = 90, draw.lines = TRUE,
                  lines.width=10) + NoLegend()+
  scale_fill_gradientn(colors = c("lightblue", "lightgrey", "red"))
ggsave(heat, file="rna_heatmap_markers.jpg", width = 30, height = 30, units = "cm")

# Manual clustering
sgl_rna$clusters <- as.character(sgl_rna$seurat_clusters)
sgl_rna$clusters[c(sgl_rna$clusters=="0"|sgl_rna$clusters=="11"|
                     sgl_rna$clusters=="19")] <- "FAPs"

sgl_rna$clusters[c(sgl_rna$clusters=="1"|sgl_rna$clusters=="2"|
                     sgl_rna$clusters=="3"|sgl_rna$clusters=="5"|
                     sgl_rna$clusters=="6"|sgl_rna$clusters=="9"|
                     sgl_rna$clusters=="12"|sgl_rna$clusters=="13"|
                     sgl_rna$clusters=="16"|sgl_rna$clusters=="17"|
                     sgl_rna$clusters=="18"|sgl_rna$clusters=="24")] <- "Myonuclei"

sgl_rna$clusters[c(sgl_rna$clusters=="4"|sgl_rna$clusters=="15"|
                     sgl_rna$clusters=="21")] <- "Endothelial"

sgl_rna$clusters[c(sgl_rna$clusters=="7"|sgl_rna$clusters=="23")] <- "Smooth_Muscular"
sgl_rna$clusters[sgl_rna$clusters=="8"] <- "Immune_Cells"
sgl_rna$clusters[sgl_rna$clusters=="10"] <- "Tenocyte"
sgl_rna$clusters[sgl_rna$clusters=="14"] <- "MuSC"
sgl_rna$clusters[sgl_rna$clusters=="20"] <- "Nervous_System"
sgl_rna$clusters[sgl_rna$clusters=="22"] <- "Pericyte"

load("SubclusterPaletteColor.RData")

cluster_palette <- subcluster_palette[-c(6:10),]
cluster_palette$Label[cluster_palette$Label=="Myonuclei_IIb"] <- "Myonuclei"
palette_color <- cluster_palette$Color
names(palette_color) <- cluster_palette$Label

p <- DimPlot(sgl_rna, reduction = "umap", group.by = "clusters",label = TRUE, repel = TRUE)+
  NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") + scale_color_manual(values = palette_color)
ggsave(p, file=paste("rna_UMAP_clusters.jpg",sep = ""), width = 20, height = 20, units = "cm")

# Re-clustering on Myonuclei ####
myo_rna <- sgl_rna[, sgl_rna$clusters=="Myonuclei"]

DefaultAssay(myo_rna) <- "integratedRNA"
# Run the standard workflow for visualization and clustering
myo_rna <- RunPCA(myo_rna, npcs = 30, verbose = FALSE)

myo_rna <- FindNeighbors(myo_rna, reduction = "pca", dims = 1:30)
myo_rna <- FindClusters(myo_rna, resolution = 0.8)

myo_markers <- FindAllMarkers(myo_rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DefaultAssay(myo_rna) <- "SCT"
atrophy_genes <- c("Trim63", "Fbxo32", "Cdkal1", "Mybpc1", "Pde4d", "Lama2", "Nctc1", "Xirp2")

heat <- DoHeatmap(myo_rna, features=atrophy_genes, size = 3, angle = 90, draw.lines = TRUE,
                  lines.width=10) + NoLegend()+
  scale_fill_gradientn(colors = c("lightblue", "lightgrey", "red"))
ggsave(heat, file="rna_myo_heatmap_atrophy.jpg", width = 30, height = 30, units = "cm")

myo_markers <- c("Ttn", "Col22a1", "Etv5", "Ano4", "Chrne","Myh1", "Myh2", "Myh4")

heat <- DoHeatmap(myo_rna, features=myo_markers, size = 3, angle = 90, draw.lines = TRUE,
                  lines.width=10) + NoLegend()+
  scale_fill_gradientn(colors = c("lightblue", "lightgrey", "red"))
ggsave(heat, file="rna_myo_heatmap_markers.jpg", width = 30, height = 30, units = "cm")

p4 <- DotPlot(object = myo_rna, features = myo_markers, scale=FALSE, assay="SCT")+ RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")
p5 <- DotPlot(object = myo_rna, features = myo_markers, assay="SCT")+ RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")
patch <- p4/p5
ggsave(patch, file="rna_myo_dotplot.jpg", width = 20, height = 40, units = "cm")

p <- DimPlot(myo_rna, reduction = "umap", group.by = "seurat_clusters",label = TRUE, repel = TRUE)+
  NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") 
ggsave(p, file=paste("rna_myo_UMAP_clusters.jpg",sep = ""), width = 20, height = 20, units = "cm")

# Manual clustering
myo_rna$subcluster <- as.character(myo_rna$seurat_clusters)

myo_rna$subcluster[c(myo_rna$subcluster=="0")] <- "Myonuclei_Trim63" 
myo_rna$subcluster[c(myo_rna$subcluster=="1"|myo_rna$subcluster=="11")] <- "Myonuclei_IIx"

myo_rna$subcluster[c(myo_rna$subcluster=="2"|myo_rna$subcluster=="6"|
                       myo_rna$subcluster=="8"|myo_rna$subcluster=="13")] <- "Myonuclei_IIx_IIb"

myo_rna$subcluster[c(myo_rna$subcluster=="3"|
                       myo_rna$subcluster=="5"|myo_rna$subcluster=="9")] <- "Myonuclei_IIb"

myo_rna$subcluster[c(myo_rna$subcluster=="4")] <- "Myonuclei_Mybpc1" 

myo_rna$subcluster[c(myo_rna$subcluster=="7")] <- "Myonuclei_IIa"

myo_rna$subcluster[c(myo_rna$subcluster=="10")] <- "Myonuclei_IIx_IIa_IIb"

myo_rna$subcluster[c(myo_rna$subcluster=="12")] <- "Myonuclei_MTJ"

myo_rna$subcluster[c(myo_rna$subcluster=="14")] <- "Myonuclei_NMJ" 

save(myo_rna, file = "rna_myo_annot.RData")

sgl_rna$subclusters <- sgl_rna$clusters
sgl_rna$subclusters[sgl_rna$subclusters=="Myonuclei"] <- myo_rna$subcluster

label_color <- subcluster_palette
label_color$Label[label_color$Label=="Myonuclei_MTJ_IIb"] <- "Myonuclei_MTJ"
label_color <- label_color[!label_color$Label == "Myonuclei_IIx_IIa",]

add_color <- c("#9ef8f2","#bfffa5", "#e19fc0","#fff68f" )
names(add_color) <- c("Myonuclei_IIa", "Myonuclei_IIx_IIa_IIb", "Myonuclei_Trim63", "Myonuclei_Mybpc1")

palette_color <- label_color$Color
names(palette_color) <- label_color$Label
palette_color <- c(palette_color, add_color)

p <- DimPlot(sgl_rna, reduction = "umap", label = TRUE, repel = TRUE, group.by="subclusters") + 
  ggtitle("") + NoLegend() + scale_color_manual(values=palette_color)
ggsave(p, file="rna_UMAP_subclusters.jpg", width = 20, height = 20, units = "cm")

saveRDS(sgl_rna, file = "rna_sgl_annot.rds")

# UMAP plot per replicate ####
DefaultAssay(sgl_rna) <- "RNA"

ll <- lapply(unique(sgl_rna$dataset), function(x) sgl_rna[, sgl_rna$dataset==x])
names(ll) <- unique(sgl_rna$dataset)

ll <- lapply(ll, FUN = function(x) {
  x <- SCTransform(x)
})

ll <- lapply(ll, function(x) RunPCA(x, npcs = 30, verbose = FALSE))
ll <- lapply(ll, function(x) RunUMAP(x, reduction = "pca", dims = 1:30))

for (i in seq_along(ll)) {
  p <-  DimPlot(ll[[i]], reduction = "umap", label = TRUE, repel = TRUE, group.by="subclusters") + 
    ggtitle(names(ll)[i]) + NoLegend() + scale_color_manual(values=palette_color)
  ggsave(p, file=paste("rna_UMAP_subclusters_",names(ll)[i],".jpg", sep = ""), width = 20, height = 20, units = "cm")
}