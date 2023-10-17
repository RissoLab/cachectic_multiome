
.libPaths(new = "/mnt/ganimede/zuin/R/x86_64-pc-linux-gnu-library/4.3")

library(Seurat)
library(Signac)
library(SingleCellExperiment)
library(batchelor)
library(scuttle)
library(scran)
library(scater)
library(harmony)
library(pheatmap)

# For reproducibility
set.seed(8569)

# MNN ####
# Load multiome + RNA integrated with MNN
mnn <- readRDS("multiome_rna_mnn.rds")

# Clustering 
nn_cluster <- clusterCells(mnn, use.dimred="corrected")
mnn$clusters <- nn_cluster

p1 <- plotTSNE(mnn, colour_by="clusters", text_by="clusters",
               point_size = 0.5, add_legend = FALSE) 
ggsave(p1, file="multiome_rna_mnn_tsne_clusters.jpg")

# Heatmap 
markers <- read.table("cachectic_markers.txt", header = TRUE)
# Cell-type specific markers were selected inside the dataset
subset_mnn <- mnn[rownames(mnn) %in% markers$Gene, ]

# The log-normalized count average was computed in each group
label <- levels(factor(subset_mnn$clusters))
average <- matrix(nrow = length(rownames(subset_mnn)),
                  ncol = length(levels(factor(subset_mnn$clusters))))
for (j in seq_along(label)) {
  average[, j] <- rowMeans(logcounts(subset_mnn[, subset_mnn$clusters == label[j]]))
}
rownames(average) <- rownames(subset_mnn)
colnames(average) <- label

for (j in seq_along(label)) {
  sub_ct <- subset_mnn[, subset_mnn$clusters == label[j]]
  reconstr <- assays(sub_ct)$reconstructed
  average[, j] <- rowMeans(reconstr)
}
rownames(average) <- rownames(subset_mnn)
colnames(average) <- label

p <- pheatmap(average, fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE)
ggsave(p, file="multiome_rna_clustering_heatmap.jpg",  width = 10, height = 10, units = "cm")
# Heatmap scaled by row
p_scale <- pheatmap(average, scale = "row", fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE)
ggsave(p_scale, file="multiome_rna_clustering_heatmap_scaled.jpg", width = 10, height = 10, units = "cm")

mnn$subclusters_new <- as.character(mnn$clusters)
mnn$subclusters_new[c(mnn$subclusters_new=="12"|mnn$subclusters_new=="15"|
                        mnn$subclusters_new=="21")] <- "Endothelial"

mnn$subclusters_new[c(mnn$subclusters_new=="6"|mnn$subclusters_new=="23")] <- "FAPs"

mnn$subclusters_new[c(mnn$subclusters_new=="7"|mnn$subclusters_new=="18"|
                        mnn$subclusters_new=="29")] <- "Tenocyte"

mnn$subclusters_new[c(mnn$subclusters_new=="26"|mnn$subclusters_new=="27"|mnn$subclusters_new=="30")] <- "Nervous_System"

mnn$subclusters_new[c(mnn$subclusters_new=="1")] <- "MuSC"

mnn$subclusters_new[c(mnn$subclusters_new=="2"|mnn$subclusters_new=="13")] <- "Immune_Cells"

mnn$subclusters_new[c(mnn$subclusters_new=="5"|mnn$subclusters_new=="24"|
                        mnn$subclusters_new=="28")] <- "Smooth_Muscular"

mnn$subclusters_new[c(mnn$subclusters_new=="19")] <- "Pericyte"

mnn$subclusters_new[c(mnn$subclusters_new=="3"|mnn$subclusters_new=="4"|
                        mnn$subclusters_new=="8"|mnn$subclusters_new=="9"|
                        mnn$subclusters_new=="10"|mnn$subclusters_new=="11"|
                        mnn$subclusters_new=="14"|
                        mnn$subclusters_new=="16"|mnn$subclusters_new=="17"|
                        mnn$subclusters_new=="20"|mnn$subclusters_new=="22"|
                        mnn$subclusters_new=="25")] <- "Myonuclei"

load("SubclusterPaletteColor.RData")
label_color <- subcluster_palette
label_color$Label[label_color$Label=="Myonuclei_MTJ_IIb"] <- "Myonuclei_MTJ"
label_color$Label[label_color$Label=="Myonuclei_IIb"] <- "Myonuclei"

add_color <- c("#9ef8f2","#bfffa5", "#e19fc0","#fff68f" )
names(add_color) <- c("Myonuclei_IIa", "Myonuclei_IIx_IIa_IIb", "Myonuclei_Trim63", "Myonuclei_Mybpc1")

palette_color <- label_color$Color
names(palette_color) <- label_color$Label
palette_color <- c(palette_color, add_color)

p1 <- plotTSNE(mnn, colour_by="subclusters_new", text_by="subclusters_new",
               point_size = 0.8, add_legend = FALSE) +
  scale_color_manual(values=palette_color)
ggsave(p1, file="multiome_rna_mnn_tsne_subclusters_new.jpg")

## Re-clustering ####
subset_myo <- mnn[,mnn$subclusters_new =="Myonuclei"]

nn_cluster <- clusterCells(subset_myo, use.dimred="corrected")
subset_myo$clusters <- nn_cluster

markers_myo <- findMarkers(subset_myo, groups=subset_myo$clusters, assay.type="reconstructed")

p1 <- plotTSNE(subset_myo, colour_by="clusters", text_by="clusters",
               point_size = 0.5, add_legend = FALSE) 
ggsave(p1, file="multiome_rna_mnn_myo_tsne_clusters.jpg")

# Heatmap 
atrophy_genes <- c("Trim63", "Fbxo32", "Cdkal1", "Mybpc1", "Pde4d", "Lama2", "Nctc1", "Xirp2")
myo_atrophy <- subset_myo[rownames(subset_myo) %in% atrophy_genes, ]

myo_markers <- c("Ttn", "Col22a1", "Etv5", "Ano4", "Chrne","Myh1", "Myh2", "Myh4")

markers <- read.table("cachectic_markers.txt", header = TRUE)$Gene
markers <- c(markers, atrophy_genes, myo_markers)
# Cell-type specific markers were selected inside the dataset
myo_markers <- subset_myo[rownames(subset_myo) %in% markers, ]

# The log-normalized count average was computed in each group
label <- levels(factor(myo_markers$clusters))
average <- matrix(nrow = length(rownames(myo_markers)),
                  ncol = length(levels(factor(myo_markers$clusters))))
for (j in seq_along(label)) {
  sub_ct <- myo_markers[, myo_markers$clusters == label[j]]
  reconstr <- assays(sub_ct)$reconstructed
  average[, j] <- rowMeans(reconstr)
}
rownames(average) <- rownames(myo_markers)
colnames(average) <- label

myo_mnn_makers <- subset_myo[rownames(subset_myo) %in% myo_markers, ]

# The log-normalized count average was computed in each group
myo_atrophy <- subset_myo[rownames(subset_myo) %in% atrophy_genes, ]

label <- levels(factor(myo_atrophy$clusters))
average <- matrix(nrow = length(rownames(myo_atrophy)),
                  ncol = length(levels(factor(myo_atrophy$clusters))))
for (j in seq_along(label)) {
  average[, j] <- rowMeans(logcounts(myo_atrophy[, myo_atrophy$clusters == label[j]]))
}
rownames(average) <- rownames(myo_atrophy)
colnames(average) <- label

for (j in seq_along(label)) {
  sub_ct <- myo_atrophy[, myo_atrophy$clusters == label[j]]
  reconstr <- assays(sub_ct)$reconstructed
  average[,j] <- rowMeans(reconstr)
}
rownames(average) <- rownames(myo_atrophy)
colnames(average) <- label

p <- pheatmap(average, fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE)
ggsave(p, file="multiome_rna_myo_clustering_heatmap_atro.jpg",  width = 10, height = 10, units = "cm")
# Heatmap scaled by row
p_scale <- pheatmap(average, scale = "row", fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE)
ggsave(p_scale, file="multiome_rna_myo_clustering_heatmap_atro_scaled.jpg", width = 10, height = 10, units = "cm")

# The log-normalized count average was computed in each group
label <- levels(factor(myo_mnn_makers$clusters))
average <- matrix(nrow = length(rownames(myo_mnn_makers)),
                  ncol = length(levels(factor(myo_mnn_makers$clusters))))
for (j in seq_along(label)) {
  average[, j] <- rowMeans(logcounts(myo_mnn_makers[, myo_mnn_makers$clusters == label[j]]))
}
rownames(average) <- rownames(myo_mnn_makers)
colnames(average) <- label

for (j in seq_along(label)) {
  sub_ct <- myo_mnn_makers[, myo_mnn_makers$clusters == label[j]]
  reconstr <- assays(sub_ct)$reconstructed
  average[,j] <- rowMeans(reconstr)
}
rownames(average) <- rownames(myo_mnn_makers)
colnames(average) <- label

p <- pheatmap(average, fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE)
ggsave(p, file="multiome_rna_myo_clustering_heatmap.jpg",  width = 10, height = 10, units = "cm")
# Heatmap scaled by row
p_scale <- pheatmap(average, scale = "row", fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE)
ggsave(p_scale, file="multiome_rna_myo_clustering_heatmap_scaled.jpg", width = 10, height = 10, units = "cm")

subset_myo$subclusters_new <- subset_myo$clusters
subset_myo$subclusters_new <- as.character(subset_myo$subclusters_new)

subset_myo$subclusters_new[c(subset_myo$subclusters_new=="2")] <- "Myonuclei_NMJ"

subset_myo$subclusters_new[c(subset_myo$subclusters_new=="7"|subset_myo$subclusters_new=="9"|
                               subset_myo$subclusters_new=="13"|subset_myo$subclusters_new=="14"|
                               subset_myo$subclusters_new=="19")] <- "Myonuclei_MTJ"

subset_myo$subclusters_new[c(subset_myo$subclusters_new=="6"|subset_myo$subclusters_new=="8"|subset_myo$subclusters_new=="12")] <- "Myonuclei_IIb"

subset_myo$subclusters_new[c(subset_myo$subclusters_new=="10"|subset_myo$subclusters_new=="18")] <- "FAPs"
subset_myo$subclusters_new[c(subset_myo$subclusters_new=="11")] <- "Endothelial"

subset_myo$subclusters_new[c(subset_myo$subclusters_new=="1")] <- "Myonuclei_Trim63"
subset_myo$subclusters_new[c(subset_myo$subclusters_new=="16")] <- "Myonuclei_IIx_IIa_IIb"

subset_myo$subclusters_new[c(subset_myo$subclusters_new=="16")] <- "Myonuclei_IIx_IIa_IIb"
subset_myo$subclusters_new[c(subset_myo$subclusters_new=="4")] <- "Myonuclei_IIx"

subset_myo$subclusters_new[c(subset_myo$subclusters_new=="5")] <- "Myonuclei_IIx_IIb"

subset_myo$subclusters_new[c(subset_myo$subclusters_new=="3")] <- "Myonuclei_IIa"

subset_myo$subclusters_new[c(subset_myo$subclusters_new=="15"|subset_myo$subclusters_new=="17")] <- "Myonuclei_Xirp2"

mnn$subclusters_new[mnn$subclusters_new =="Myonuclei"] <- subset_myo$subclusters_new

table(subset_myo$subcluster_new)

save(mnn, file= "Cachectic_Multiome_RNA_MNN.RData")

load("SubclusterPaletteColor.RData")
label_color <- subcluster_palette
label_color$Label[label_color$Label=="Myonuclei_MTJ_IIb"] <- "Myonuclei_MTJ"

add_color <- c("#9ef8f2","#bfffa5", "#e19fc0","#fff68f" )
names(add_color) <- c("Myonuclei_IIa", "Myonuclei_IIx_IIa_IIb", "Myonuclei_Trim63", "Myonuclei_Xirp2")

palette_color <- label_color$Color
names(palette_color) <- label_color$Label
palette_color <- c(palette_color, add_color)

p1 <- plotTSNE(mnn, colour_by="subclusters_new", text_by="subclusters_new",
               point_size = 0.5, add_legend = FALSE) + scale_color_manual(values = palette_color)
ggsave(p1, file="multiome_rna_mnn_tsne_subclusters_new.jpg")
