.libPaths(new = "/mnt/ganimede/zuin/R/x86_64-pc-linux-gnu-library/4.3")

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(edgeR)
library(scuttle)
library(ggrepel)

# for reproducibility
set.seed(9212)

# Load RNA dataset
seu_rna <- readRDS("rna_list_filt.rds")
seu_rna[[1]]$condition <- "CAC"
seu_rna[[1]]$dataset <- "RNACAC"
seu_rna[[2]]$condition <- "SHAM"
seu_rna[[2]]$dataset <- "RNASHAM"

# RNA  Integration ####
# Normalization with SCT
seu_rna <- lapply(X = seu_rna, FUN = function(x) {
  x <- SCTransform(x)
})

# Set default assay to SCT for all objects in list 
for (i in seq_along(seu_rna)) {
  DefaultAssay(seu_rna[[i]]) <- "SCT"
}

# Choose the features to use when integrating multiple datasets.
features <- SelectIntegrationFeatures(seu_rna, nfeatures = 3000)

# Prepare an object list normalized with sctransform for integration.
seu_rna <- PrepSCTIntegration(seu_rna, anchor.features = features)

seu_rna <- lapply(X = seu_rna, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = seu_rna, normalization.method = "SCT",anchor.features = features, 
                                        dims = 1:30, reduction = "rpca", k.anchor = 20)

rna_comb <- IntegrateData(anchorset = anchors, normalization.method = "SCT", new.assay.name = "integratedRNA", dims = 1:30)

rna_comb <- RunPCA(rna_comb, npcs = 30, verbose = FALSE)
rna_comb <- RunUMAP(rna_comb, reduction = "pca", dims = 1:30)

# Doublets identification #####
DefaultAssay(rna_comb) <- "RNA"

# Create a list object
ll <- lapply(unique(rna_comb$dataset), function(x) rna_comb[, rna_comb$dataset==x])
names(ll) <- unique(rna_comb$dataset)

# From Seurat to SingleCellExperiment object
scelist <- lapply(ll, function(x) as.SingleCellExperiment(x))

# HVGs were calculated for each sample
topgs <- lapply(scelist, function(x) getTopHVGs(x, prop=0.1))

scores <- mapply(function(x,y) computeDoubletDensity(x, subset.row=y) , scelist, topgs)
dbl_calls <- lapply(scores, function(x) doubletThresholding(data.frame(score=x), 
                                                            method="griffiths", returnType="call")) 

save(scores, dbl_calls, file = "RNA_Doublets.RData")

rna_comb$gex_scores <- rna_comb$dataset
rna_comb$gex_dbl_calls <- rna_comb$dataset
for (i in seq_along(unique(rna_comb$dataset))) {
  rna_comb$gex_scores[rna_comb$gex_scores==unique(rna_comb$dataset)[i]] <- scores[[i]]
  rna_comb$gex_dbl_calls[rna_comb$gex_dbl_calls==unique(rna_comb$dataset)[i]] <- dbl_calls[[i]]
}
rna_comb$gex_dbl_calls[rna_comb$gex_dbl_calls=="1"] <- "singlet"
rna_comb$gex_dbl_calls[rna_comb$gex_dbl_calls=="2"] <- "doublet"

color <- c("grey","red")
names(color) <- c("singlet","doublet")

p <- DimPlot(rna_comb, reduction = "umap", group.by = "gex_dbl_calls") + 
  NoLegend() + scale_color_manual(values=color) +ggtitle("Doublets")
ggsave(p, file=paste("rna_UMAP_doublets.jpg",sep = ""), width = 20, height = 20, units = "cm")

saveRDS(rna_comb, file = "rna_comb_dblcalls.rds")