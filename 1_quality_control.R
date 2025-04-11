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
seu1_cac <- CreateSeuratObject(counts = counts$'Gene Expression', assay = "RNA", project = "749", meta.data = metadata)
# create ATAC assay and add it to the object
seu1_cac[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)

## QC ####
# For RNA
seu1_cac[["percent.mt"]] <- PercentageFeatureSet(seu1_cac, pattern = "^mt-")

# Violin plot for RNA QC metrics
p1 <- VlnPlot(
  object = seu1_cac,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP1CAC_QC_RNA_unfiltered.jpg")

# For ATAC
DefaultAssay(seu1_cac) <- "ATAC"

seu1_cac <- NucleosomeSignal(seu1_cac)
seu1_cac <- TSSEnrichment(seu1_cac)

seu1_cac$pct_reads_in_peaks <- seu1_cac$atac_peak_region_fragments / seu1_cac$atac_fragments * 100
seu1_cac$blacklist_ratio <- FractionCountsInRegion(object = seu1_cac, assay = 'ATAC', regions = blacklist_mm10)

# Violin plot for ATAC QC metrics
p2 <- VlnPlot(
  object = seu1_cac,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP1CAC_QC_ATAC_unfiltered.jpg", width = 25, height = 25, units = "cm")

# Save the Seurat object
saveRDS(seu1_cac, file = "multiome_rep1_cac.rds")

# Remove the low quality droplets
seu1_cac <- subset(
  x = seu1_cac,
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
  object = seu1_cac,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP1CAC_QC_RNA.jpg")

# For ATAC
p2 <- VlnPlot(
  object = seu1_cac,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP1CAC_QC_ATAC.jpg", width = 25, height = 25, units = "cm")

# Save filtered Seurat object into a .rds
saveRDS(seu1_cac, file = "multiome_filt_rep1_cac.rds")

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
seu1_sham <- CreateSeuratObject(counts = counts$'Gene Expression', assay = "RNA", project = "768", meta.data = metadata)
# create ATAC assay and add it to the object
seu1_sham[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)

## QC ####
# For RNA
seu1_sham[["percent.mt"]] <- PercentageFeatureSet(seu1_sham, pattern = "^mt-")

p1 <- VlnPlot(
  object = seu1_sham,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP1SHAM_QC_RNA_unfiltered.jpg")

# For ATAC
DefaultAssay(seu1_sham) <- "ATAC"

seu1_sham <- NucleosomeSignal(seu1_sham)
seu1_sham <- TSSEnrichment(seu1_sham)

seu1_sham$pct_reads_in_peaks <- seu1_sham$atac_peak_region_fragments/seu1_sham$atac_fragments * 100
seu1_sham$blacklist_ratio <- FractionCountsInRegion(object = seu1_sham, assay = 'ATAC', regions = blacklist_mm10)

p2 <- VlnPlot(
  object = seu1_sham,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP1SHAM_QC_ATAC_unfiltered.jpg", width = 25, height = 25, units = "cm")
# Save Seurat object
saveRDS(seu1_sham, file = "multiome_rep1_sham.rds")

# Remove low-quality droplets
seu1_sham <- subset(
  x = seu1_sham,
  subset = atac_peak_region_fragments > 3000 &
    atac_peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5  & nCount_RNA < 20000 & nCount_RNA > 500 & nCount_ATAC < 100000 &  nCount_ATAC > 1000 
)

p1 <- VlnPlot(
  object = seu1_sham,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP1SHAM_QC_RNA.jpg")

p2 <- VlnPlot(
  object = seu1_sham,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP1SHAM_QC_ATAC.jpg", width = 25, height = 25, units = "cm")
# Save filtered Seurat object into a .rds
saveRDS(seu1_sham, file = "multiome_filt_rep1_sham.rds")

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
seu2_cac <- CreateSeuratObject(counts = counts$'Gene Expression', assay = "RNA", project = "750", meta.data = metadata)
# create ATAC assay and add it to the object
seu2_cac[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)

## QC ####
# For RNA
seu2_cac[["percent.mt"]] <- PercentageFeatureSet(seu2_cac, pattern = "^mt-")

p1 <- VlnPlot(
  object = seu2_cac,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP2CAC_QC_RNA.jpg")

# For ATAC
DefaultAssay(seu2_cac) <- "ATAC"

seu2_cac <- NucleosomeSignal(seu2_cac)
seu2_cac <- TSSEnrichment(seu2_cac)

seu2_cac$pct_reads_in_peaks <- seu2_cac$atac_peak_region_fragments/seu2_cac$atac_fragments * 100
seu2_cac$blacklist_ratio <- FractionCountsInRegion(object = seu2_cac, assay = 'ATAC', regions = blacklist_mm10)

p2 <- VlnPlot(
  object = seu2_cac,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP2CAC_QC_ATAC.jpg", width = 25, height = 25, units = "cm")
# Save Seurat object
saveRDS(seu2_cac, file = "multiome_rep2_cac.rds")

# Remove low-quality droplets
seu2_cac <- subset(
  x = seu2_cac,
  subset = atac_peak_region_fragments > 3000 &
    atac_peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5  & nCount_RNA < 20000 & nCount_RNA > 500 & nCount_ATAC < 100000 &  nCount_ATAC > 1000 
)

p1 <- VlnPlot(
  object = seu2_cac,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP2CAC_QC_RNA.jpg")

p2 <- VlnPlot(
  object = seu2_cac,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP2CAC_QC_ATAC.jpg", width = 25, height = 25, units = "cm")

saveRDS(seu2_cac, file = "multiome_filt_rep2_cac.rds")

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
seu2_sham <- CreateSeuratObject(counts = counts$'Gene Expression', assay = "RNA", project = "769", meta.data = metadata)
# create ATAC assay and add it to the object
seu2_sham[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)

## QC ####
# For RNA
seu2_sham[["percent.mt"]] <- PercentageFeatureSet(seu2_sham, pattern = "^mt-")

p1 <- VlnPlot(
  object = seu2_sham,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP2SHAM_QC_RNA.jpg")

# For ATAC
DefaultAssay(seu2_sham) <- "ATAC"

seu2_sham <- NucleosomeSignal(seu2_sham)
seu2_sham <- TSSEnrichment(seu2_sham)

seu2_sham$pct_reads_in_peaks <- seu2_sham$atac_peak_region_fragments/seu2_sham$atac_fragments * 100
seu2_sham$blacklist_ratio <- FractionCountsInRegion(object = seu2_sham, assay = 'ATAC', regions = blacklist_mm10)

p2 <- VlnPlot(
  object = seu2_sham,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP2SHAM_QC_ATAC.jpg", width = 25, height = 25, units = "cm")
# Save Seurat object
saveRDS(seu2_sham, file = "multiome_rep2_sham.rds")

# Remove low-quality droplets
seu2_sham <- subset(
  x = seu2_sham,
  subset = atac_peak_region_fragments > 3000 &
    atac_peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5  & nCount_RNA < 20000 & nCount_RNA > 500 & nCount_ATAC < 100000 &  nCount_ATAC > 1000 
)

p1 <- VlnPlot(
  object = seu2_sham,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(p1, file = "multiome_REP2SHAM_QC_RNA.jpg")

p2 <- VlnPlot(
  object = seu2_sham,
  features = c("pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
  ncol = 3
)
ggsave(p2, file = "multiome_REP2SHAM_QC_ATAC.jpg", width = 25, height = 25, units = "cm")
# Save filtered Seurat object
saveRDS(seu2_sham, file = "multiome_filt_rep2_sham.rds")
