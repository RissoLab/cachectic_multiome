library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(patchwork)

# For reproducibility
set.seed(8569)

# This file contains RNA and ATAC integration.
# Load dataset after QC
seu1_sham <- readRDS("multiome_filt_rep1_sham.rds")
seu2_sham <- readRDS("multiome_filt_rep2_sham.rds")
seu1_cac <- readRDS("multiome_filt_rep1_cac.rds")
seu2_cac <- readRDS("multiome_filt_rep2_cac.rds")

seu1_sham$dataset <- "REP1SHAM"
seu1_sham$condition <- "SHAM"
seu2_sham$dataset <- "REP2SHAM"
seu2_sham$condition <- "SHAM"
seu1_cac$dataset <- "REP1CAC"
seu1_cac$condition <- "CAC"
seu2_cac$dataset <- "REP2CAC"
seu2_cac$condition <- "CAC"

# Create a list object
ll <- list(seu1SHAM = seu1_sham, seu2SHAM = seu2_sham, seu1CAC = seu1_cac, seu2CAC = seu2_cac)

# Data was normalized with SCT
ll <- lapply(X = ll, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x)
})

# Set default assay to SCT for all objects in list
for (i in seq_along(ll)) {
  DefaultAssay(ll[[i]]) <- "SCT"
}
# Choose the features to use when integrating multiple datasets.
features <- SelectIntegrationFeatures(ll, nfeatures = 3000)

# Prepare an object list normalized with sctransform for integration.
ll <- PrepSCTIntegration(ll, anchor.features = features)

ll <- lapply(X = ll, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Integration ####
## RNA integration ####
# Instead of utilizing canonical correlation analysis (CCA) to identify anchors,
# we instead utilize reciprocal PCA (‘RPCA’).
# The anchors of RNA integration was computed, using k.anchor = 20
anchors <- FindIntegrationAnchors(object.list = ll, normalization.method = "SCT",
                                  anchor.features = features, dims = 1:30,
                                  reduction = "rpca", k.anchor = 20)

mult_comb <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                          new.assay.name = "integratedRNA", dims = 1:30)

mult_comb <- RunPCA(mult_comb, npcs = 30, verbose = FALSE)
mult_comb <- RunUMAP(mult_comb, reduction = "pca", dims = 1:30)

## ATAC integration ####
# run LSI on new seurat object with integrated RNA assay
DefaultAssay(mult_comb) <- "ATAC"
mult_comb <- RunTFIDF(mult_comb)
mult_comb <- FindTopFeatures(mult_comb, min.cutoff = "q25")
mult_comb <- RunSVD(mult_comb)
# The RNA anchors were used for ATAC integration
so_atac <- IntegrateEmbeddings(
  anchorset = anchors,
  new.reduction.name = "integratedLSI",
  reductions = mult_comb@reductions$lsi
)

mult_comb@reductions$integratedLSI <- so_atac@reductions$integratedLSI

mult_comb <- RunUMAP(mult_comb, reduction = "integratedLSI", dims = 2:30,
                    reduction.name = "UMAP.atac",  reduction.key = "atacUMAP_")

## WNN ####
mult_comb <- FindMultiModalNeighbors(mult_comb, reduction.list = list("pca", "integratedLSI"),
                                    dims.list = list(1:30, 2:30))
mult_comb <- RunUMAP(mult_comb, nn.name = "weighted.nn", reduction.name = "UMAP.wnn",
                    reduction.key = "wnnUMAP_")

# Doublets removal ####
# We identified doublets in each replicate
DefaultAssay(mult_comb) <- "RNA"

# Create a list object
ll <- lapply(unique(mult_comb$dataset), function(x) mult_comb[, mult_comb$dataset == x])
names(ll) <- unique(mult_comb$dataset)

ll <- lapply(ll, function(x){DefaultAssay(x) <- "RNA"; x})
# From Seurat to SingleCellExperiment object
scelist <- lapply(ll, function(x) as.SingleCellExperiment(x))

# HVGs were calculated for each sample
topgs <- lapply(scelist, function(x) getTopHVGs(x, prop = 0.1))

scores <- mapply(function(x,y) computeDoubletDensity(x, subset.row = y) , scelist, topgs)
dbl_calls <- lapply(scores, function(x) doubletThresholding(data.frame(score = x), 
                                                            method = "griffiths", returnType = "call")) 

save(scores, dbl_calls, file = "Multiome_Doublets.RData")

mult_comb$gex_scores <- mult_comb$dataset
mult_comb$gex_dbl_calls <- mult_comb$dataset
for (i in seq_along(unique(mult_comb$dataset))) {
  mult_comb$gex_scores[mult_comb$gex_scores == unique(mult_comb$dataset)[i]] <- scores[[i]]
  mult_comb$gex_dbl_calls[mult_comb$gex_dbl_calls == unique(mult_comb$dataset)[i]] <- dbl_calls[[i]]
}
mult_comb$gex_dbl_calls[mult_comb$gex_dbl_calls == "1"] <- "singlet"
mult_comb$gex_dbl_calls[mult_comb$gex_dbl_calls == "2"] <- "doublet"

color <- c("grey","red")
names(color) <- c("singlet","doublet")

p <- DimPlot(mult_comb, reduction = "umap", group.by = "gex_dbl_calls") + 
  NoLegend() + scale_color_manual(values = color) + ggtitle("Doublets")
ggsave(p, file=paste("multiome_UMAP_doublets.jpg",sep = ""), width = 20, height = 20, units = "cm")

saveRDS(mult_comb, file = "multiome_comb_dblcalls.rds")