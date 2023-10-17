.libPaths(new = "/mnt/ganimede/zuin/R/x86_64-pc-linux-gnu-library/4.3")

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)

# for reproducibility
set.seed(9212)

# This file contains the quality control computed for each multiome replicate.
# REP1 CAC: Replicate 1 CACHECTIC condition ####
# Import data
path <- "/mnt/tebe/data/2023_Sandri/multiome_cachectic/"
cond <- "libraries/REP1_CAC"
# RNA matrix
genemat <- "outs/filtered_feature_bc_matrix.h5"
# ATAC matrix
fragmat <- "outs/atac_fragments.tsv.gz"

fragpath <- file.path(path, cond, fragmat)

counts <- Read10X_h5(file.path(path, cond, genemat))

# Get gene annotation and create Seurat object.
ensdb <- EnsDb.Mmusculus.v79
seqlevelsStyle(ensdb) <- "UCSC"
annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotation) <- "mm10"

# Metadata can be found in the cellranger-arc output: outs/per_barcode_metrics.csv
metadata <- read.csv(
  file = "/mnt/tebe/data/2023_Sandri/multiome_cachectic/libraries/REP1_CAC/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

# Create a Seurat object containing the RNA matrix
seu1CAC <- CreateSeuratObject(counts = counts$'Gene Expression', assay = "RNA", project = "749", meta.data = metadata)
# create ATAC assay and add it to the object
seu1CAC[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)

## QC ####
# For RNA
seu1CAC[["percent.mt"]] <- PercentageFeatureSet(seu1CAC, pattern = "^mt-")

# Violin plot for RNA QC metrics
p1 <- VlnPlot(
  object = seu1CAC,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP1CAC_QC_RNA_unfiltered.jpg")

# For ATAC
DefaultAssay(seu1CAC) <- "ATAC"

seu1CAC <- NucleosomeSignal(seu1CAC)
seu1CAC <- TSSEnrichment(seu1CAC)

seu1CAC$pct_reads_in_peaks <- seu1CAC$atac_peak_region_fragments / seu1CAC$atac_fragments * 100
seu1CAC$blacklist_ratio <- FractionCountsInRegion(object = seu1CAC, assay = 'ATAC', regions = blacklist_mm10)

# Violin plot for ATAC QC metrics
p2 <- VlnPlot(
  object = seu1CAC,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP1CAC_QC_ATAC_unfiltered.jpg", width = 25, height = 25, units = "cm")

# Save the Seurat object
saveRDS(seu1CAC, file = "multiome_rep1_cac.rds")

# Remove the low quality droplets
seu1CAC <- subset(
  x = seu1CAC,
  subset = atac_peak_region_fragments > 3000 &
    atac_peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5  & nCount_RNA < 20000 & nCount_RNA > 500 & nCount_ATAC < 100000 &  nCount_ATAC > 1000 
)

# For RNA
p1 <- VlnPlot(
  object = seu1CAC,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP1CAC_QC_RNA.jpg")

# For ATAC
p2 <- VlnPlot(
  object = seu1CAC,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP1CAC_QC_ATAC.jpg", width = 25, height = 25, units = "cm")

# Save filtered Seurat object into a .rds
saveRDS(seu1CAC, file = "multiome_filt_rep1_cac.rds")

# REP1 SHAM: Replicate 1 SHAM/control condition ####
# Import data
path <- "/mnt/tebe/data/2023_Sandri/multiome_cachectic/"
cond <- "libraries/REP1_SHAM"
# RNA matrix
genemat <- "outs/filtered_feature_bc_matrix.h5"
# ATAC matrix
fragmat <- "outs/atac_fragments.tsv.gz"

fragpath <- file.path(path, cond, fragmat)

counts <- Read10X_h5(file.path(path, cond, genemat))

# Get gene annotation and create Seurat object.
ensdb <- EnsDb.Mmusculus.v79
seqlevelsStyle(ensdb) <- "UCSC"
annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotation) <- "mm10"

metadata <- read.csv(
  file = "/mnt/tebe/data/2023_Sandri/multiome_cachectic/libraries/REP1_SHAM/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

# Create a Seurat object containing the RNA data
seu1SHAM <- CreateSeuratObject(counts = counts$'Gene Expression', assay = "RNA", project = "768", meta.data = metadata)
# create ATAC assay and add it to the object
seu1SHAM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)

## QC ####
# For RNA
seu1SHAM[["percent.mt"]] <- PercentageFeatureSet(seu1SHAM, pattern = "^mt-")

p1 <- VlnPlot(
  object = seu1SHAM,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP1SHAM_QC_RNA_unfiltered.jpg")

# For ATAC
DefaultAssay(seu1SHAM) <- "ATAC"

seu1SHAM <- NucleosomeSignal(seu1SHAM)
seu1SHAM <- TSSEnrichment(seu1SHAM)

seu1SHAM$pct_reads_in_peaks <- seu1SHAM$atac_peak_region_fragments/seu1SHAM$atac_fragments * 100
seu1SHAM$blacklist_ratio <- FractionCountsInRegion(object = seu1SHAM, assay = 'ATAC', regions = blacklist_mm10)

p2 <- VlnPlot(
  object = seu1SHAM,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP1SHAM_QC_ATAC_unfiltered.jpg", width = 25, height = 25, units = "cm")
# Save Seurat object
saveRDS(seu1SHAM, file = "multiome_rep1_sham.rds")

# Remove low-quality droplets
seu1SHAM <- subset(
  x = seu1SHAM,
  subset = atac_peak_region_fragments > 3000 &
    atac_peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5  & nCount_RNA < 20000 & nCount_RNA > 500 & nCount_ATAC < 100000 &  nCount_ATAC > 1000 
)

p1 <- VlnPlot(
  object = seu1SHAM,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP1SHAM_QC_RNA.jpg")

p2 <- VlnPlot(
  object = seu1SHAM,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP1SHAM_QC_ATAC.jpg", width = 25, height = 25, units = "cm")
# Save filtered Seurat object into a .rds
saveRDS(seu1SHAM, file = "multiome_filt_rep1_sham.rds")

# REP2 CAC: Replicate 2 CACHECTIC condition ####
# Import data
path <- "/mnt/tebe/data/2023_Sandri/multiome_cachectic/"
cond <- "libraries/REP2_CAC"
# RNA matrix
genemat <- "outs/filtered_feature_bc_matrix.h5"
# ATAC matrix
fragmat <- "outs/atac_fragments.tsv.gz"

fragpath <- file.path(path, cond, fragmat)

counts <- Read10X_h5(file.path(path, cond, genemat))

# Get gene annotation and create Seurat object.
ensdb <- EnsDb.Mmusculus.v79
seqlevelsStyle(ensdb) <- "UCSC"
annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotation) <- "mm10"

metadata <- read.csv(
  file = "/mnt/tebe/data/2023_Sandri/multiome_cachectic/libraries/REP2_CAC/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

# Create a Seurat object containing the RNA data
seu2CAC <- CreateSeuratObject(counts = counts$'Gene Expression', assay = "RNA", project = "750", meta.data = metadata)
# create ATAC assay and add it to the object
seu2CAC[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)

## QC ####
# For RNA
seu2CAC[["percent.mt"]] <- PercentageFeatureSet(seu2CAC, pattern = "^mt-")

p1 <- VlnPlot(
  object = seu2CAC,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP2CAC_QC_RNA.jpg")

# For ATAC
DefaultAssay(seu2CAC) <- "ATAC"

seu2CAC <- NucleosomeSignal(seu2CAC)
seu2CAC <- TSSEnrichment(seu2CAC)

seu2CAC$pct_reads_in_peaks <- seu2CAC$atac_peak_region_fragments/seu2CAC$atac_fragments * 100
seu2CAC$blacklist_ratio <- FractionCountsInRegion(object = seu2CAC, assay = 'ATAC', regions = blacklist_mm10)

p2 <- VlnPlot(
  object = seu2CAC,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP2CAC_QC_ATAC.jpg", width = 25, height = 25, units = "cm")
# Save Seurat object
saveRDS(seu2CAC, file = "multiome_rep2_cac.rds")

# Remove low-quality droplets
seu2CAC <- subset(
  x = seu2CAC,
  subset = atac_peak_region_fragments > 3000 &
    atac_peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5  & nCount_RNA < 20000 & nCount_RNA > 500 & nCount_ATAC < 100000 &  nCount_ATAC > 1000 
)

p1 <- VlnPlot(
  object = seu2CAC,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP2CAC_QC_RNA.jpg")

p2 <- VlnPlot(
  object = seu2CAC,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP2CAC_QC_ATAC.jpg", width = 25, height = 25, units = "cm")

saveRDS(seu2CAC, file = "multiome_filt_rep2_cac.rds")

# REP2 SHAM : Replicate 2 SHAM condition ####
# Import data
path <- "/mnt/tebe/data/2023_Sandri/multiome_cachectic/"
cond <- "libraries/REP2_SHAM"
# RNA matrix
genemat <- "outs/filtered_feature_bc_matrix.h5"
# ATAC matrix
fragmat <- "outs/atac_fragments.tsv.gz"

fragpath <- file.path(path, cond, fragmat)

counts <- Read10X_h5(file.path(path, cond, genemat))

# Get gene annotation and create Seurat object.
ensdb <- EnsDb.Mmusculus.v79
seqlevelsStyle(ensdb) <- "UCSC"
annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotation) <- "mm10"

metadata <- read.csv(
  file = "/mnt/tebe/data/2023_Sandri/multiome_cachectic/libraries/REP2_SHAM/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

# Create a Seurat object containing the RNA data
seu2SHAM <- CreateSeuratObject(counts = counts$'Gene Expression', assay = "RNA", project = "769", meta.data = metadata)
# create ATAC assay and add it to the object
seu2SHAM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)

## QC ####
# For RNA
seu2SHAM[["percent.mt"]] <- PercentageFeatureSet(seu2SHAM, pattern = "^mt-")

p1 <- VlnPlot(
  object = seu2SHAM,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP2SHAM_QC_RNA.jpg")

# For ATAC
DefaultAssay(seu2SHAM) <- "ATAC"

seu2SHAM <- NucleosomeSignal(seu2SHAM)
seu2SHAM <- TSSEnrichment(seu2SHAM)

seu2SHAM$pct_reads_in_peaks <- seu2SHAM$atac_peak_region_fragments/seu2SHAM$atac_fragments * 100
seu2SHAM$blacklist_ratio <- FractionCountsInRegion(object = seu2SHAM, assay = 'ATAC', regions = blacklist_mm10)

p2 <- VlnPlot(
  object = seu2SHAM,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP2SHAM_QC_ATAC.jpg", width = 25, height = 25, units = "cm")
# Save Seurat object
saveRDS(seu2SHAM, file = "multiome_rep2_sham.rds")

# Remove low-quality droplets
seu2SHAM <- subset(
  x = seu2SHAM,
  subset = atac_peak_region_fragments > 3000 &
    atac_peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5  & nCount_RNA < 20000 & nCount_RNA > 500 & nCount_ATAC < 100000 &  nCount_ATAC > 1000 
)

p1 <- VlnPlot(
  object = seu2SHAM,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP2SHAM_QC_RNA.jpg")

p2 <- VlnPlot(
  object = seu2SHAM,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP2SHAM_QC_ATAC.jpg", width = 25, height = 25, units = "cm")
# Save filtered Seurat object
saveRDS(seu2SHAM, file="multiome_filt_rep2_sham.rds")