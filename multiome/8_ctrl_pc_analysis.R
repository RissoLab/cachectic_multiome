# =============================================================================
# Cachectic multiome (PC and CTRL samples)
# Complete analysis
# =============================================================================

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scDblFinder)

set.seed(9212)
# Quality control ######

# Get gene annotation and create Seurat object.
ensdb <- EnsDb.Mmusculus.v79
seqlevelsStyle(ensdb) <- "UCSC"
annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotation) <- "mm10"

genemat <- list.files(pattern = ".h5")
genemat <- genemat[!grepl("CAC_1", genemat)]
fragmat <- list.files(pattern = ".tsv.gz")
fragmat <- fragmat[!grepl("CAC_1", fragmat)]

counts <- lapply(genemat, function(i) Read10X_h5(i))

files <- list.files(pattern = "metrics.csv")
metadata <- lapply(files, function(i) read.csv(file = i, header = TRUE, row.names = 1))

for(i in seq_along(fragmat)) {
  fragments <- fragmat[[i]] 
  system(paste("tabix -p bed", fragments))
}

sample_id <- c("PC_1", "PC_2", "PC_3", "CTRL_3")

ll <- lapply(seq_along(counts), function(x) {
  so <- CreateSeuratObject(counts = counts[[x]]$'Gene Expression', assay = "RNA",  meta.data = metadata[[x]])
  so[["ATAC"]] <- CreateChromatinAssay(counts = counts[[x]]$Peaks, sep = c(":", "-"), 
                                       fragments = fragmat[[x]], annotation = annotation)
  so$sample_id <- sample_id[x]
  so
})

ll <- lapply(ll, function(x) {
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  DefaultAssay(x) <- "ATAC"
  x <- NucleosomeSignal(x)
  x <- TSSEnrichment(x)
  x$pct_reads_in_peaks <- x$atac_peak_region_fragments / x$atac_fragments * 100
  x$blacklist_ratio <- FractionCountsInRegion(object = x, assay = 'ATAC', regions = blacklist_mm10)
  x
})

ll <- lapply(ll, function(x) {
  x <- subset(
    x = x,
    subset = atac_peak_region_fragments > 3000 &
      atac_peak_region_fragments < 20000 & pct_reads_in_peaks > 15 &
      blacklist_ratio < 0.05 & nucleosome_signal < 4 &
      TSS.enrichment > 2 & nFeature_RNA > 200 &
      nFeature_RNA < 3500 & percent.mt < 5  &
      nCount_RNA < 20000 & nCount_RNA > 500 & 
      nCount_ATAC < 100000 & nCount_ATAC > 1000 
  )
})

solist <- readRDS("cachectic_mult_filt_ll.rds")
solist <- c(solist, CTRL_3 = ll[[4]], PC_1 = ll[[1]], PC_2 = ll[[2]], PC_3 = ll[[3]])

saveRDS(solist, file = "cachectic_mult_filt_ll.rds")

# Doublets removal #####
solist <- solist[5:8]

scelist <- lapply(solist, function(x){ 
  DefaultAssay(x) <- "RNA"
  x <- as.SingleCellExperiment(x)
  x <- logNormCounts(x)
  })

topgs <- lapply(scelist, function(x) getTopHVGs(x, prop = 0.1))

scores <- mapply(function(x,y) computeDoubletDensity(x, subset.row = y) , scelist, topgs)
dbl_calls <- lapply(scores, function(x) doubletThresholding(data.frame(score = x), 
                                                            method = "griffiths", returnType = "call")) 

solist <- lapply(seq_along(solist), function(x) {
  solist[[x]]$gex_dbl <- dbl_calls[[x]]
  solist[[x]] <- solist[[x]][, solist[[x]]$gex_dbl == "singlet"]
})

# Cell type annotation ####
solist <- lapply(solist, function(x){
  DefaultAssay(x) <- "RNA"
  x[["ATAC"]] <- NULL
  x <- SCTransform(x)
  DefaultAssay(x) <- "SCT"
  x
})

## For PC samples  #####
ll <- solist[2:4]

features <- SelectIntegrationFeatures(ll, nfeatures = 3000)
ll <- PrepSCTIntegration(ll, anchor.features = features)

ll <- lapply(X = ll, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = ll, normalization.method = "SCT",
                                  anchor.features = features, dims = 1:30,
                                  reduction = "rpca", k.anchor = 20)

mult_comb <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                           new.assay.name = "integratedRNA", dims = 1:30)

mult_comb <- RunPCA(mult_comb, npcs = 30, verbose = FALSE)
mult_comb <- RunUMAP(mult_comb, reduction = "pca", dims = 1:30)

mult_comb <- FindNeighbors(mult_comb, reduction = "pca", dims = 1:30)
mult_comb <- FindClusters(mult_comb, resolution = 0.5)

markers <- read.table("cachectic_markers.txt", header = T)
markers <- markers$Gene

DotPlot(object = mult_comb, features = markers, assay = "SCT")+ RotatedAxis() +
  scale_color_continuous(low = "grey", high = "red")


mult_comb$clusters <- as.character(mult_comb$seurat_clusters)
mult_comb$clusters[mult_comb$clusters %in% as.character(c(0:3, 5:8, 10))] <- "Myonuclei"
mult_comb$clusters[mult_comb$clusters %in% as.character(c(4))] <- "FAPs"
mult_comb$clusters[mult_comb$clusters %in% as.character(c(11,14))] <- "Smooth_Muscular"
mult_comb$clusters[mult_comb$clusters %in% as.character(c(9))] <- "Endothelial"
mult_comb$clusters[mult_comb$clusters %in% as.character(c(13))] <- "Immune_Cells"
mult_comb$clusters[mult_comb$clusters %in% as.character(c(12))] <- "MuSC"

# Re-clustering on Myonuclei
myo_mult <- mult_comb[, mult_comb$clusters == "Myonuclei"]
DefaultAssay(myo_mult) <- "integratedRNA"
myo_mult <- RunPCA(myo_mult, npcs = 30, verbose = FALSE)

myo_mult <- FindNeighbors(myo_mult, reduction = "pca", dims = 1:30)
myo_mult <- FindClusters(myo_mult, resolution = 0.5)

myo_markers <- c("Ttn", "Col22a1", "Etv5", "Ano4", "Chrne","Myh1", "Myh2", "Myh4")

p2 <- DotPlot(object = myo_mult, features = myo_markers, assay="SCT") + 
  RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")
ggsave(p2, file="plots/cachectic_mult_dotplot_myo_markers_clust.jpg", width = 20, height = 20, units = "cm")

myo_mult$subcluster <- as.character(myo_mult$seurat_clusters)
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(0, 2, 3, 6))] <- "Myonuclei_IIb"
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(1, 8, 11))] <- "Myonuclei_IIx_IIb"
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(4))] <- "Myonuclei_IIx"
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(5))] <- "Myonuclei_IIx_IIa"
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(7))] <- "Myonuclei_MTJ"
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(10))] <- "Myonuclei_NMJ"
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(9, 12))] <- "Myonuclei_Trim63"

DimPlot(myo_mult, reduction = "umap", group.by = "subcluster", 
        label = TRUE, repel = TRUE)

mult_comb$subclusters <- mult_comb$clusters
mult_comb$subclusters[mult_comb$subclusters =="Myonuclei"] <- myo_mult$subcluster

p <- DimPlot(mult_comb, reduction = "umap", group.by = "subclusters",label = TRUE, repel = TRUE)+
  NoLegend() + ggtitle("WNN") + scale_color_manual(values = palette_color)
ggsave(p, file=paste("plots/cachectic_mult_wnn_umap_subclust_pc.jpg", sep = ""),
       width = 20, height = 20, units = "cm")

# For CTRL Replicate 3
ctrl <- solist[[1]]
ctrl <- FindNeighbors(ctrl, reduction = "pca", dims = 1:30)
ctrl <- FindClusters(ctrl, resolution = 0.8)

p2 <- DotPlot(object = ctrl, features = markers, assay = "SCT") + RotatedAxis() +
  scale_color_continuous(low = "grey", high = "red")
ggsave(p2, file=paste("plots/cachectic_mult_dotplot_markers_clust.jpg", sep = ""), 
       width = 20, height = 20, units = "cm")

ctrl$clusters <- as.character(ctrl$seurat_clusters)
ctrl$clusters[ctrl$clusters %in% as.character(c(0:5, 7, 13, 14))] <- "Myonuclei"
ctrl$clusters[ctrl$clusters %in% as.character(c(6, 8))] <- "FAPs"
ctrl$clusters[ctrl$clusters %in% as.character(c(9, 16))] <- "Smooth_Muscular"
ctrl$clusters[ctrl$clusters %in% as.character(c(10))] <- "Endothelial"
ctrl$clusters[ctrl$clusters %in% as.character(c(15))] <- "Immune_Cells"
ctrl$clusters[ctrl$clusters %in% as.character(c(11))] <- "MuSC"
ctrl$clusters[ctrl$clusters %in% as.character(c(12))] <- "Tenocyte"

myo_mult <- ctrl[, ctrl$clusters == "Myonuclei"]
myo_mult <- RunPCA(myo_mult, npcs = 30, verbose = FALSE)

myo_mult <- FindNeighbors(myo_mult, reduction = "pca", dims = 1:30)
myo_mult <- FindClusters(myo_mult, resolution = 0.8)

myo_markers <- c("Ttn", "Col22a1", "Etv5", "Ano4", "Chrne","Myh1", "Myh2", "Myh4")

p2 <- DotPlot(object = myo_mult, features = myo_markers, assay="SCT") + 
  RotatedAxis() +  scale_color_continuous(low = "grey", high = "red")
ggsave(p2, file="plots/cachectic_mult_dotplot_markers_clust.jpg", width = 20, height = 20, units = "cm")

p <- DimPlot(myo_mult, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5,label = TRUE, repel = TRUE)+
  NoLegend() 
ggsave(p, file=paste("plots/cachectic_mult_umap_ctrl.jpg",sep = ""), width = 20, height = 20, units = "cm")

myo_mult$subcluster <- as.character(myo_mult$seurat_clusters)
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(2, 3, 5, 6))] <- "Myonuclei_IIb"
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(1))] <- "Myonuclei_IIx_IIb"
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(0, 4, 10))] <- "Myonuclei_IIx"
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(8))] <- "Myonuclei_IIx_IIa"
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(7))] <- "Myonuclei_MTJ"
myo_mult$subcluster[myo_mult$subcluster %in% as.character(c(9))] <- "Myonuclei_NMJ"

ctrl$subclusters <- ctrl$clusters
ctrl$subclusters[ctrl$subclusters =="Myonuclei"] <- myo_mult$subcluster

mult_ctrl <- readRDS("cachecti_mult_ctrl_c26.rds")
mult_ctrl <- mult_ctrl[, mult_ctrl$condition == "CTRL"]
mult_ctrl <- lapply(unique(mult_ctrl$sample_id), function(x) mult_ctrl[, mult_ctrl$sample_id == x])

pcs <- lapply(unique(mult_pc$sample_id), function(x) mult_pc[, mult_pc$sample_id == x])

solist <- list(CTRL_1 = mult_ctrl[[1]], CTRL_2 = mult_ctrl[[2]], CTRL_3 = ctrl,
               PC_1 = pcs[[1]], PC_2 = pcs[[2]], PC_3 = pcs[[3]])


int <- Reduce(intersect, lapply(solist, rownames))
solist <- lapply(solist, function(x) x[int, ])
solist <- lapply(solist, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x)
  DefaultAssay(x) <- "SCT"
  x
})

features <- SelectIntegrationFeatures(solist, nfeatures = 3000)
solist <- PrepSCTIntegration(solist, anchor.features = features)

solist <- lapply(solist, function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = solist, normalization.method = "SCT",
                                  anchor.features = features, dims = 1:30,
                                  reduction = "rpca", k.anchor = 20)
mult_comb <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                           new.assay.name = "integratedRNA", dims = 1:30)

mult_comb <- RunPCA(mult_comb, npcs = 30, verbose = FALSE)
mult_comb <- RunUMAP(mult_comb, reduction = "pca", dims = 1:30)
mult_comb$subclusters <- gsub("Trim63", "IIb", mult_comb$subclusters)

load("palette_color_ct.RData")
p <- DimPlot(mult_comb, reduction = "umap", group.by = "subclusters", label = T, repel = T) + 
  scale_color_manual(values = palette_color) + ggtitle("RNA") + NoLegend()
ggsave(p, file=paste("plots/cachectic_mult_umap_ctrl_pc_sample.jpg",sep = ""), width = 20, height = 20, units = "cm")

# ATAC integration ###
DefaultAssay(mult_comb) <- "ATAC"
mult_comb <- RunTFIDF(mult_comb)
mult_comb <- FindTopFeatures(mult_comb, min.cutoff = "q25")
mult_comb <- RunSVD(mult_comb)
mult_comb <- RunUMAP(mult_comb, reduction = "integratedLSI", dims = 2:30,
                     reduction.name = "UMAP.atac",  reduction.key = "atacUMAP_")

## WNN ####
mult_comb <- FindMultiModalNeighbors(mult_comb, reduction.list = list("pca", "integratedLSI"),
                                     dims.list = list(1:30, 2:30))
mult_comb <- RunUMAP(mult_comb, nn.name = "weighted.nn", reduction.name = "UMAP.wnn",
                     reduction.key = "wnnUMAP_")

p <- DimPlot(mult_comb, reduction = "UMAP.wnn", group.by = "subclusters",label = TRUE, repel = TRUE)+
  NoLegend() + ggtitle("WNN") + scale_color_manual(values = palette_color)

saveRDS(mult_comb, file = "cachectic_mult_ctrl_pc.rds")

# DEA #####
DefaultAssay(mult_comb) <- "RNA"

so <- JoinLayers(mult_comb)

sce <- SingleCellExperiment(assays=list(counts = so@assays$RNA$count))
colData(sce) <- DataFrame(so@meta.data)

# Genes filtering: keep the genes that have an expression average greater than zero in at least of 5% cells
# for each replicate.
keep <- apply(counts(sce), 1, function(y) all(tapply(y, sce$sample_id, function(x) mean(x > 0)) >= 0.05))
sce <- sce[keep, ]

# Pseudo-bulk ####
pb <- aggregateAcrossCells(sce, use.assay.type = "counts", 
                           id=DataFrame(label = sce$subclusters, sample = sce$sample_id))
pb <- logNormCounts(pb)
colnames(pb) <- paste(pb$subclusters, pb$sample_id, sep = "_")

ct_pb <- pb[, pb$subclusters == unique(pb$subclusters)[2]]
groups <- makeGroups(ct_pb$condition)
set <- newSeqExpressionSet(as.matrix(counts(ct_pb)),
                           phenoData = data.frame(colData(ct_pb)))
uq <- betweenLaneNormalization(set, which="upper")

ruv1 <- RUVs(uq, cIdx = empirical, scIdx = groups, k=1)
col <- c("red", "blue")
x <- as.factor(rep(c("CTRL", "PC"), each=3))

jpeg(filename = "plots/pca_ruv1.jpg")
plotPCA(ruv1, col = col[x], cex=1.2)
dev.off()

## Negative control genes #####
pb_cond <- aggregateAcrossCells(sce, use.assay.type = "counts", 
                                id=DataFrame(sample = sce$sample_id))
pb_cond <- logNormCounts(pb_cond)

design <- model.matrix(~0+condition, data=colData(pb_cond))
colnames(design) <- c("CTRL", "PC")
y <- DGEList(counts=counts(pb_cond), group=pb_cond$condition)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

contrast <- makeContrasts(PC - CTRL,levels = design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast = contrast)

empirical <- topTags(lrt, n=nrow(pb_cond))$table
empirical <- empirical[empirical$FDR > 0.80,]
empirical <- rownames(empirical)

# DEA #####
diff_exp_analysis <- function(pb, i) {
  print(unique(pb$subclusters)[i])
  ct_pb <- pb[, pb$subclusters == unique(pb$subclusters)[i]]
  groups <- makeGroups(ct_pb$condition)
  set <- newSeqExpressionSet(as.matrix(counts(ct_pb)),
                             phenoData = data.frame(colData(ct_pb)))
  uq <- betweenLaneNormalization(set, which="upper")
  
  ruv1 <- RUVs(uq, cIdx = empirical, scIdx = groups, k=1)
  
  y <- DGEList(counts(ct_pb), samples = colData(ct_pb))
  y <- calcNormFactors(y, method = "upperquartile")
  
  design <- model.matrix(~0 + y$samples$condition + ruv1$W_1, y$samples)
  colnames(design) <- c("CTRL", "PC", "W_1")
  
  contrast <- makeContrasts(PC - CTRL,levels = design)
  
  y <- estimateDisp(y, design)
  
  fit <- glmFit(y, design)
  
  res <- glmLRT(fit, contrast = contrast)
  
  x <- as.data.frame(topTags(res, n = length(rownames(pb))))
  
  x <- x[order(x$logFC, decreasing = TRUE), ]
}

df_mult <- lapply(c(1:4, 6:11, 14), function(i) diff_exp_analysis(pb, i))
names(df_mult) <- unique(pb$subclusters)[c(1:4, 6:11, 14)]
saveRDS(df_mult, file = "cachectic_mult_ctrl_pc_degs.rds")
