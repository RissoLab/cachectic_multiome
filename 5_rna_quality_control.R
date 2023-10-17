library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(ggplot2)

# for reproducibility
set.seed(9212)
# 10X Multiome data processing: REP CAC and REP SHAM. Two replicate without ATAC ####
# REP CAC Import ####
path <- "/mnt/tebe/data/2023_Sandri/multiome_cachectic/"
cond <- "libraries/RNA_CAC"
# RNA matrix
genemat <- "outs/filtered_feature_bc_matrix.h5"
counts <- Read10X_h5(file.path(path, cond, genemat))

seu_cac <- CreateSeuratObject(counts = counts, assay = "RNA", project = "693")
seu_cac[["percent.mt"]] <- PercentageFeatureSet(seu_cac, pattern = "^mt-")

seu_cac <- subset(
  x = seu_cac,
  subset = nFeature_RNA > 200 & nFeature_RNA< 3500 & percent.mt < 5 &
    nCount_RNA < 20000 & nCount_RNA > 500
)

p1 <- VlnPlot(
  object = seu_cac,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file=paste("rna_CAC_QC.jpg",sep = ""))

# REP SHAM Import ####
path <- "/mnt/tebe/data/2023_Sandri/multiome_cachectic/"
cond <- "libraries/RNA_SHA"
# RNA matrix
genemat <- "outs/filtered_feature_bc_matrix.h5"
counts <- Read10X_h5(file.path(path, cond, genemat))

seu_sham <- CreateSeuratObject(counts = counts, assay = "RNA", project= "692")
seu_sham[["percent.mt"]] <- PercentageFeatureSet(seu_sham, pattern = "^mt-")

seu_sham <- subset(
  x = seu_sham,
  subset = nFeature_RNA > 200 & nFeature_RNA< 3500 & percent.mt < 5 &
    nCount_RNA < 20000 & nCount_RNA >  500
)

p1 <- VlnPlot(
  object = seu_sham,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file=paste("rna_SHAM_QC.jpg",sep = ""))

# Save replicates into a list ####
seu_rna <- list(CAC=seu_cac, SHAM=seu_sham)
saveRDS(seu_rna, file = "rna_list_filt.rds")