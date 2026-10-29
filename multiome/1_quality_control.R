# =============================================================================
# Cachectic multiome (C26 and CTRL samples)
# Quality control
# =============================================================================

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)

set.seed(9212)

# Import data #####
# Get gene annotation
ensdb <- EnsDb.Mmusculus.v79
seqlevelsStyle(ensdb) <- "UCSC"
annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotation) <- "mm10"

# C26 Replicate 1 ###
# RNA expression matrix file path
genepath <- "libraries/REP1_CAC/outs/filtered_feature_bc_matrix.h5"

# ATAC expression matrix file path
fragpath <- "libraries/REP1_CAC/outs/atac_fragments.tsv.gz"

counts <- Read10X_h5(genepath)
# Metadata can be found in the cellranger-arc output: outs/per_barcode_metrics.csv
metadata <- read.csv(file = "libraries/REP1_CAC/outs/per_barcode_metrics.csv", header = TRUE, row.names = 1)

seu1_cac <- CreateSeuratObject(counts = counts$'Gene Expression', 
                               assay = "RNA", project = "749", meta.data = metadata)
# Add ATAC assay
seu1_cac[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, 
                                           sep = c(":", "-"), fragments = fragpath, 
                                           annotation = annotation)

# C26 Replicate 2 ####
# RNA expression matrix file path
genepath <- "libraries/REP2_CAC/outs/filtered_feature_bc_matrix.h5"

# ATAC expression matrix file path
fragpath <- "libraries/REP2_CAC/outs/atac_fragments.tsv.gz"

counts <- Read10X_h5(genepath)
# Metadata can be found in the cellranger-arc output: outs/per_barcode_metrics.csv
metadata <- read.csv(file = "libraries/REP2_CAC/outs/per_barcode_metrics.csv", header = TRUE, row.names = 1)

seu2_cac <- CreateSeuratObject(counts = counts$'Gene Expression', 
                               assay = "RNA", project = "750", meta.data = metadata)
# Add ATAC assay
seu2_cac[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, 
                                           sep = c(":", "-"), fragments = fragpath, 
                                           annotation = annotation)

# CTRL Replicate 1 ####
# RNA expression matrix file path
genepath <- "libraries/REP1_SHAM/outs/filtered_feature_bc_matrix.h5"

# ATAC expression matrix file path
fragpath <- "libraries/REP1_SHAM/outs/atac_fragments.tsv.gz"

counts <- Read10X_h5(genepath)
# Metadata can be found in the cellranger-arc output: outs/per_barcode_metrics.csv
metadata <- read.csv(file = "libraries/REP1_SHAM/outs/per_barcode_metrics.csv", header = TRUE, row.names = 1)

seu1_sham <- CreateSeuratObject(counts = counts$'Gene Expression', 
                               assay = "RNA", project = "768", meta.data = metadata)
# Add ATAC assay
seu1_sham[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, 
                                           sep = c(":", "-"), fragments = fragpath, 
                                           annotation = annotation)

# CTRL Replicate 2 ####
# RNA expression matrix file path
genepath <- "libraries/REP2_SHAM/outs/filtered_feature_bc_matrix.h5"

# ATAC expression matrix file path
fragpath <- "libraries/REP2_SHAM/outs/atac_fragments.tsv.gz"

counts <- Read10X_h5(genepath)
# Metadata can be found in the cellranger-arc output: outs/per_barcode_metrics.csv
metadata <- read.csv(file = "libraries/REP2_SHAM/outs/per_barcode_metrics.csv", header = TRUE, row.names = 1)

seu2_sham <- CreateSeuratObject(counts = counts$'Gene Expression', 
                                assay = "RNA", project = "769", meta.data = metadata)
# Add ATAC assay
seu2_sham[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, 
                                            sep = c(":", "-"), fragments = fragpath, 
                                            annotation = annotation)

ll <- list(C26_1 = seu1_cac, C26_2 = seu2_cac, CTRL_1 = seu1_sham, CTRL_2 = seu2_sham)
lables <- names(ll) 

# Quality control ####
ll <- lapply(ll, function(x) {
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  
  DefaultAssay(x) <- "ATAC"
  x <- NucleosomeSignal(x)
  x <- TSSEnrichment(x)
  x$pct_reads_in_peaks <- x$atac_peak_region_fragments / x$atac_fragments * 100
  x$blacklist_ratio <- FractionCountsInRegion(object = x, assay = 'ATAC', regions = blacklist_mm10)
  x
})

# Remove low quality droplets
ll <- lapply(ll, function(x){
  x <- subset(x = x,
              subset = atac_peak_region_fragments > 3000 &
                atac_peak_region_fragments < 20000 &
                pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 &
                nucleosome_signal < 4 & TSS.enrichment > 2 & nFeature_RNA > 200 & 
                nFeature_RNA < 3500 & percent.mt < 5  & nCount_RNA < 20000 & 
                nCount_RNA > 500 & nCount_ATAC < 100000 &  nCount_ATAC > 1000 
  )
})

ll <- lapply(seq_along(ll), function(x) {
  ll[[x]]$sample_id <- labels[x]
  ll[[x]]$condition <- sub("\\_.*", "", labels[x])
  ll[[x]]
})

saveRDS(ll, file = "cachectic_mult_filt_ll.rds")

for (i in seq_along(ll)) {
  p1 <- VlnPlot(object = ll[[i]],
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(p1, file = paste0("cachectic_mult_qc_rna_", lables[i], ".pdf"), width = 25, height = 25, units = "cm")
  
  p2 <- VlnPlot(object = ll[[i]],
                features = c("pct_reads_in_peaks", "atac_peak_region_fragments", 
                             "TSS.enrichment", "blacklist_ratio","nucleosome_signal"),
                ncol = 3)
  ggsave(p1, file = paste0("cachectic_mult_qc_atac_", lables[i], ".pdf"), width = 25, height = 25, units = "cm")
}
