.libPaths(new = "/mnt/ganimede/zuin/R/x86_64-pc-linux-gnu-library/4.3")

library(Seurat)
library(Signac)
library(SingleCellExperiment)
library(batchelor)
library(scuttle)
library(scran)
library(scater)
library(harmony)

# For reproducibility
set.seed(8569)

# Load multiome and RNA samples
sgl_mult <- readRDS("multiome_sgl_comb.rds")
DefaultAssay(sgl_mult) <- "RNA"

sgl_rna <- readRDS("rna_sgl_annot.rds")
DefaultAssay(sgl_rna) <- "RNA"

# Convert into a SingleCellExperiment
sce_mult <- as.SingleCellExperiment(sgl_mult)
colData(sce_mult) <- colData(sce_mult)[,c(43,44,55)]
altExp(sce_mult) <- NULL  

sce_rna <- as.SingleCellExperiment(sgl_rna)
colData(sce_rna) <- colData(sce_rna)[,c(5,6,17)]

altExp(sce_rna) <- NULL

scelist_mult <- lapply(unique(sce_mult$dataset), function(x)
  sce_mult[, sce_mult$dataset==x])
names(scelist_mult) <- unique(sce_mult$dataset)

scelist_gex <- lapply(unique(sce_rna$dataset), function(x)
  sce_rna[, sce_rna$dataset==x])
names(scelist_gex) <- unique(sce_rna$dataset)

ll <- c(scelist_mult, scelist_gex)

ll <- lapply(ll, logNormCounts)
all.dec <- lapply(ll, modelGeneVar)

combined.dec <- combineVar(all.dec[[1]], all.dec[[2]], 
                           all.dec[[3]], all.dec[[4]], 
                           all.dec[[5]], all.dec[[6]])
chosen.hvgs <- combined.dec$bio > 0
sum(chosen.hvgs)

## MNN ####
mnn <- fastMNN(ll[[1]], ll[[2]], ll[[3]], ll[[4]],
               ll[[5]], ll[[6]], d=50, k=20, subset.row=chosen.hvgs,
               BSPARAM=BiocSingular::RandomParam(deferred=TRUE))

uncorrected <- cbind(ll[[1]], ll[[2]], ll[[3]], ll[[4]], ll[[5]], ll[[6]])
uncorrected <- uncorrected[rownames(uncorrected) %in% rownames(mnn), ]

assay(mnn, "counts") <- counts(uncorrected)
assay(mnn, "logcounts") <- logcounts(uncorrected)

mnn <- runTSNE(mnn, dimred="corrected")
mnn <- runUMAP(mnn, dimred="corrected")
mnn$subclusters <- uncorrected$subclusters
mnn$sample_id <- uncorrected$dataset
mnn$condition <- uncorrected$condition

load("SubclusterPaletteColor.RData")
label_color <- subcluster_palette
label_color$Label[label_color$Label=="Myonuclei_MTJ_IIb"] <- "Myonuclei_MTJ"

add_color <- c("#9ef8f2","#bfffa5", "#e19fc0","#fff68f" )
names(add_color) <- c("Myonuclei_IIa", "Myonuclei_IIx_IIa_IIb", "Myonuclei_Trim63", "Myonuclei_Mybpc1")

palette_color <- label_color$Color
names(palette_color) <- label_color$Label
palette_color <- c(palette_color, add_color)

p1<- plotTSNE(mnn, colour_by="sample_id") 
ggsave(p1, file="multiome_rna_mnn_tsne_replicate.jpg")

p1 <- plotTSNE(mnn, colour_by="subclusters", text_by="subclusters",
               point_size = 0.8, add_legend = FALSE) +
  scale_color_manual(values=palette_color)
ggsave(p1, file="multiome_rna_mnn_tsne_subclusters.jpg")

p1<- plotTSNE(mnn, colour_by="condition") 
ggsave(p1, file="multiome_rna_mnn_tsne_condition.jpg")

p1<- plotUMAP(mnn, colour_by="sample_id") 
ggsave(p1, file="multiome_rna_mnn_umap_replicate.jpg")

p1 <- plotUMAP(mnn, colour_by="subclusters", text_by="subclusters",
               point_size = 0.8, add_legend = FALSE) +
  scale_color_manual(values=palette_color)
ggsave(p1, file="multiome_rna_mnn_umap_subclusters.jpg")

p1<- plotUMAP(mnn, colour_by="condition") 
ggsave(p1, file="multiome_rna_mnn_umap_condition.jpg")

# Save MNN object
saveRDS(mnn, file = "multiome_rna_mnn.rds")

## Harmony ####
mult_counts <- sgl_mult@assays$RNA@counts
meta_mult <- cbind(sgl_mult$dataset, sgl_mult$condition, sgl_mult$subclusters)

rna_counts <- sgl_rna@assays$RNA@counts
meta_rna <- cbind(sgl_rna$dataset, sgl_rna$condition, sgl_rna$subclusters)

mt_rna_counts <- cbind(mult_counts, rna_counts)
metadata <- as.data.frame(rbind(meta_mult, meta_rna))
colnames(metadata) <- c("dataset", "condition", "subclusters")

mult_rna <- CreateSeuratObject(counts = mt_rna_counts, assay = "RNA", meta.data=metadata)

mult_rna <- mult_rna %>%
  NormalizeData(mult_rna = FALSE)

VariableFeatures(mult_rna) <- split(row.names(mult_rna@meta.data), 
                                    mult_rna@meta.data$dataset) %>% lapply(function(cells_use) {
                                      mult_rna[,cells_use] %>%
                                        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
                                        VariableFeatures()
                                    }) %>% unlist %>% unique

mult_rna <- mult_rna %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(features = VariableFeatures(mult_rna), npcs = 20, verbose = FALSE)

mult_rna <- RunHarmony(mult_rna, "dataset")

mult_rna <- mult_rna %>%
  RunTSNE(reduction = "harmony")

p <- DimPlot(mult_rna, reduction = "tsne", group.by = "dataset", pt.size = .5) +
  ggtitle("") 
ggsave("multiome_rna_harmony_tsne_dataset.jpg")

p1 <- DimPlot(object = mult_rna, reduction = "tsne", pt.size = .5, 
              label = TRUE, repel = TRUE, group.by = "subclusters") +
  ggtitle("") +scale_color_manual(values=palette_color) + NoLegend() 
ggsave("multiome_rna_harmony_tsne_subclusters.jpg")

p1 <- DimPlot(object = mult_rna, reduction = "tsne", pt.size = .5, 
              group.by = "condition") +
  ggtitle("") 
ggsave("multiome_rna_harmony_tsne_condition.jpg")

mult_rna <- mult_rna %>%
  RunUMAP(reduction = "harmony",  dims = 1:20)

p <- DimPlot(mult_rna, reduction = "umap", group.by = "dataset", pt.size = .5) +
  ggtitle("") 
ggsave("multiome_rna_harmony_umap_dataset.jpg")

p1 <- DimPlot(object = mult_rna, reduction = "umap", pt.size = .5, 
              label = TRUE, repel = TRUE, group.by = "subclusters") +
  ggtitle("") +scale_color_manual(values=palette_color) + NoLegend()
ggsave("multiome_rna_harmony_umap_subclusters.jpg")

p1 <- DimPlot(object = mult_rna, reduction = "umap", pt.size = .5, 
              group.by = "condition") +
  ggtitle("") 
ggsave("multiome_rna_harmony_umap_condition.jpg")

saveRDS(mult_rna, file = "multiome_rna_harmony.rds")