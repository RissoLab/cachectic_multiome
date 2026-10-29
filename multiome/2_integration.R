# =============================================================================
# Cachectic multiome (C26 and CTRL samples)
# Multi-omics integration between RNA-seq and ATAC-seq
# =============================================================================

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(8569)

ll <- readRDS("cachectic_mult_filt_ll.rds")

ll <- lapply(ll, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x)
  DefaultAssay(x) <- "RNA"
  x
})

features <- SelectIntegrationFeatures(ll, nfeatures = 3000)

ll <- PrepSCTIntegration(ll, anchor.features = features)

ll <- lapply(ll, function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Doublets removal ####
# From Seurat to SingleCellExperiment object
scelist <- lapply(ll, function(x) as.SingleCellExperiment(x))
topgs <- lapply(scelist, function(x) getTopHVGs(x, prop = 0.1))

scores <- mapply(function(x,y) computeDoubletDensity(x, subset.row = y) , scelist, topgs)
dbl_calls <- lapply(scores, function(x) doubletThresholding(data.frame(score = x), 
                                                            method = "griffiths", returnType = "call")) 

ll <- lapply(seq_along(ll), function(x){
  ll[[x]]$gex_scores <- scores[[i]]
  ll[[x]]$gex_dbl_calls <- dbl_calls[[i]]
  
  ll[[x]]$gex_dbl_calls[ll[[x]]$gex_dbl_calls == "1"] <- "singlet"
  ll[[x]]$gex_dbl_calls[ll[[x]]$gex_dbl_calls == "2"] <- "doublet"
  ll[[x]] <- ll[[x]][, ll[[x]]$gex_dbl_calls == "singlet"]
})

# Integration ####
# Multi-omics data was integrated using the weighted nearest neighbor (WNN) method.
# Reciprocal PCA (RPCA) approach was used to identify anchors
anchors <- FindIntegrationAnchors(object.list = ll, normalization.method = "SCT",
                                  anchor.features = features, dims = 1:30,
                                  reduction = "rpca", k.anchor = 20)

mult_comb <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                           new.assay.name = "integratedRNA", dims = 1:30)
# We computed PCA and UMAP on integratedRNA asssay
mult_comb <- RunPCA(mult_comb, npcs = 30, verbose = FALSE)
mult_comb <- RunUMAP(mult_comb, reduction = "pca", dims = 1:30)

DefaultAssay(mult_comb) <- "ATAC"
mult_comb <- RunTFIDF(mult_comb)
mult_comb <- FindTopFeatures(mult_comb, min.cutoff = "q25")
mult_comb <- RunSVD(mult_comb)

# The RNA anchors were used for ATAC integration
so_atac <- IntegrateEmbeddings(anchorset = anchors, new.reduction.name = "integratedLSI",
                               reductions = mult_comb@reductions$lsi)
mult_comb@reductions$integratedLSI <- so_atac@reductions$integratedLSI

mult_comb <- RunUMAP(mult_comb, reduction = "integratedLSI", dims = 2:30,
                     reduction.name = "UMAP.atac",  reduction.key = "atacUMAP_")

# WNN computation
mult_comb <- FindMultiModalNeighbors(mult_comb, reduction.list = list("pca", "integratedLSI"),
                                     dims.list = list(1:30, 2:30))
mult_comb <- RunUMAP(mult_comb, nn.name = "weighted.nn", reduction.name = "UMAP.wnn",
                     reduction.key = "wnnUMAP_")

saveRDS(mult_comb, file = "cachectic_mult_ctrl_c26.rds")
