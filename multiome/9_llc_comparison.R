# =============================================================================
# Comparison Sandri (PC vs CTRL) and Kir datasets 
# =============================================================================

library(Seurat)
library(Signac)
library(SingleCellExperiment)
library(scuttle)
library(edgeR)
library(ggvenn)
library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)
library(tidyverse)

kir <- readRDS("llc_dataset.rds")
kir$subclusters <- kir$celltype
kir$sample_id <- "kir"
kir$condition <- kir$orig.ident
kir <- lapply(unique(kir$condition), function(x) kir[, kir$condition == x])
sandri <- readRDS("cachectic_mult_ctrl_pc.rds")
sandri <- lapply(unique(sandri$sample_id), function(x) sandri[, sandri$sample_id == x])

# Batch effects correction ##########
ll <- c(kir, sandri)

ll <- lapply(ll, function(x) {
  DefaultAssay(x) <- "RNA"
  x[["integratedRNA"]] <- NULL
  x[["integrated"]] <- NULL
  x[["ATAC"]] <- NULL
  x[["SCT"]] <- NULL
  x@meta.data <- x@meta.data[, c("sample_id", "condition", "subclusters")]
  x
})

inter <- Reduce(intersect, lapply(ll, rownames))
ll <- lapply(ll, function(x) x[inter, ])

ll <- lapply(ll, function(x) {
  x <- SCTransform(x)
  DefaultAssay(x) <- "SCT"
  x
})

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
saveRDS(mult_comb, file = "cachectic_mult_kir_pc.rds")

mult_comb <-readRDS("cachectic_mult_kir.rds")
mult_comb$celltype <- mult_comb$subclusters
mult_comb$celltype[mult_comb$sample_id == "kir"] <- paste0(mult_comb$celltype[mult_comb$sample_id == "kir"], "_Kir")

mult_comb$celltype[grepl("Smooth M.", mult_comb$celltype)] <- "Smooth_Muscular_Kir"
mult_comb$celltype[grepl("Macrophage", mult_comb$celltype)] <- "Immune_Cells_Kir"
mult_comb$celltype[grepl("MTJ_K", mult_comb$celltype)]  <- "Myonuclei_MTJ_Kir"
mult_comb$celltype[grepl("NMJ_K", mult_comb$celltype)]  <- "Myonuclei_NMJ_Kir"
mult_comb$celltype[grepl("Macrophage", mult_comb$celltype)]  <- "Smooth_Muscular_Kir"
mult_comb$celltype[grepl("Type IIa_", mult_comb$celltype)] <- "Myonuclei_IIa_Kir"
mult_comb$celltype[grepl("Type IIb_", mult_comb$celltype)] <- "Myonuclei_IIb_Kir"
mult_comb$celltype[grepl("Type IIx_", mult_comb$celltype)] <- "Myonuclei_IIx_Kir"
mult_comb$celltype[grepl("Type IIx/b", mult_comb$celltype)] <- "Myonuclei_IIx_IIb_Kir"
mult_comb$celltype[grepl("Myonuclei_Gsk3b", mult_comb$celltype)] <- "Myonuclei_Trim63"

load("palette_color_ct.RData")
pltt_kir <- palette_color
names(pltt_kir) <- paste0(names(pltt_kir), "_Kir")   
pltt_kir <- c(pltt_kir, "Myonuclei_IIa_Kir" = "#3F8A4B")
palette_color <- c(palette_color, pltt_kir)

DimPlot(mult_comb, reduction = "umap", group.by = "celltype", pt.size = 0.5, label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("RNA") + scale_color_manual(values = palette_color)


# DEGs ####
DefaultAssay(kir) <- "RNA"
kir$condition <- kir$orig.ident
Idents(object = kir) <- "condition"

kir <- JoinLayers(kir)
sce <- as.SingleCellExperiment(kir)


keep <- apply(counts(sce), 1, function(y) all(tapply(y, sce$condition, function(x) mean(x > 0)) >= 0.1))
sce <- sce[keep, ]                                                

sce <- sce[, c(grepl("Type", sce$celltype) | sce$celltype %in% c("NMJ", "MTJ"))]
sce$celltype <- as.character(sce$celltype)

pb_mult <- aggregateAcrossCells(sce, use.assay.type = "counts",
                                id = DataFrame(sample  = sce$condition))
pb_mult$celltype <- "Myonuclei"
pb_mult <- logNormCounts(pb_mult)
colnames(pb_mult) <- paste(pb_mult$celltype, pb_mult$condition, sep = "_")

y <- DGEList(counts(pb_mult))
y <- calcNormFactors(y, method = "upperquartile")

design <- model.matrix(~0 + condition, data = colData(pb_mult))
colnames(design) <- c("C26", "CTRL")

contrast <- makeContrasts(C26 - CTRL, levels = design)

bcv <- 0.4
y$common.dispersion <- bcv^2

fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast = contrast)

df <- as.data.frame(topTags(lrt, n = Inf))
df <- cbind(Gene_Name_ID = rownames(df), df)
saveRDS(df, file = "cachectic_kir_degs_myonuclei.rds")
# writexl::write_xlsx(x = df, path = "cachectic_kir_degs_myonuclei.xlsx")

sandri <- readRDS("degs/cachectic_mult_ctrl_pc_degs_overmyo.rds")
sandri <- sandri[sandri$FDR < 0.05, ]
sandri_up <- dplyr::filter(sandri, logFC >= 0)
sandri_down <- dplyr::filter(sandri, logFC <= 0)

kir_up <- dplyr::filter(df, logFC > 0)
kir_down <- dplyr::filter(df, logFC <= 0)

venn_list <- list(Sandri_Up = rownames(sandri_up), Kir_Up = kir_up$Gene_Name_ID)
inter <- intersect(rownames(sandri_up), kir_up$Gene_Name_ID)
inter <- data.frame(Gene_Name_ID = inter)
writexl::write_xlsx(inter, path = "degs/common_degs_up_kir_sandri.xlsx")

ggvenn(venn_list, fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
       stroke_size = 0.5, set_name_size = 4, show_percentage = T)

venn_list <- list(Sandri_Down = rownames(sandri_down),Kir_Down = kir_down$Gene_Name_ID)
inter <- intersect(rownames(sandri_down), kir_down$Gene_Name_ID)
inter <- data.frame(Gene_Name_ID = inter)
writexl::write_xlsx(inter, path = "common_degs_down_kir_sandri.xlsx")

ggvenn(venn_list, fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
       stroke_size = 0.5, set_name_size = 4, show_percentage = T)

# ORA ####
genes <- readxl::read_xlsx("common_degs_up_kir_sandri.xlsx")$Gene_Name_ID

sandri <- readRDS("cachectic_mult_ctrl_pc_degs_overmyo.rds")
sig_genes_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
universe <- bitr(rownames(sandri), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

universe_entrez <- universe$ENTREZID

enrich_KEGG <- enrichKEGG(
  gene          = sig_genes_entrez$ENTREZID,
  organism      = 'mmu',                     
  keyType       = "kegg",                   
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)     

enrich_KEGG <- enrich_KEGG@result
enrich_KEGG$geneSymbol <- sapply(enrich_KEGG$geneID, function(x) {
  entrez_ids <- unlist(strsplit(x, "/"))
  symbols <- sig_genes_entrez %>%
    filter(ENTREZID %in% entrez_ids) %>%
    pull(SYMBOL)
  paste(symbols, collapse = "/")
})

enrich_KEGG$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", enrich_KEGG$Description)
enrich_KEGG <- enrich_KEGG[order(enrich_KEGG$p.adjust), ]
enrich_KEGG <- enrich_KEGG[enrich_KEGG$p.adjust < 0.05, ]

enrich_KEGG <- enrich_KEGG %>%
  separate(GeneRatio, into = c("num", "den"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_num = num / den)
writexl::write_xlsx(enrich_KEGG, path = "cachectic_kir_sandri_gsea_common_up.xlsx")

ggplot(head(enrich_KEGG,15), aes(-log10(p.adjust),
                                 y=reorder(Description,-log10(p.adjust)),
                                 color= GeneRatio_num, 
                                 size= Count)) + geom_point() +
  scale_color_gradient(low = "red", high = "blue")+
  theme_minimal() + ylab("Description") + xlab("-log10(p.adjust)") + ggtitle("Top 15KEGG")

