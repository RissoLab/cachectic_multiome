library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(scDblFinder)
library(scuttle)
library(scran)
library(edgeR)
library(ggrepel)
library(EDASeq)
library(UpSetR)
library(ComplexUpset)
library(patchwork)

set.seed(8569)

# This file contains the differential expression analysis (DEA) and differential abundance analysis (DAA)
sgl_mult <- readRDS("multiome_sgl_comb.rds")
DefaultAssay(sgl_mult) <- "RNA"
# From Seurat to SingleCellExperiment object
sce_mult <- as.SingleCellExperiment(sgl_mult)

# Genes filtering: keep the genes that have an expression average greater than zero in at least of 10% for each replicate.
keep <- apply(counts(sce_mult), 1, function(y) all(tapply(y, sce_mult$dataset, function(x) mean(x > 0)) >= 0.1))

filtered_sce <- sce_mult[keep, ]                                                

# Create pseudo-bulk 
pb_mult <- aggregateAcrossCells(filtered_sce, use.assay.type = "counts", 
                                id=DataFrame(label = filtered_sce$subclusters, sample = filtered_sce$dataset))
pb_mult <- logNormCounts(pb_mult)

colnames(pb_mult) <- paste(pb_mult$subclusters, pb_mult$dataset, sep = "_")
# Save pseudo-bulk
saveRDS(pb_mult, file = "multiome_pb.rds")

# DEA ####
## Negative binomial generalized model (GLM) with Upper Quantile normalization ####
res_mult <- df_mult <- list()
for (i in seq_along(unique(pb_mult$subclusters))) {
  ct_pb <- pb_mult[, pb_mult$subclusters == unique(pb_mult$subclusters)[i]]
  
  y <- DGEList(counts(ct_pb), samples = colData(ct_pb))
  y <- calcNormFactors(y, method = "upperquartile")
  
  design <- model.matrix(~0 + y$samples$condition, y$samples)
  colnames(design) <- c("C26", "CTRL")
  
  contrast <- makeContrasts(C26 - CTRL,levels = design)
  
  y <- estimateDisp(y, design)
  
  fit <- glmFit(y, design)
  
  res_mult[[i]] <- glmLRT(fit, contrast = contrast)
  
  df_mult[[i]] <- as.data.frame(topTags(res_mult[[i]], n = length(rownames(pb_mult))))

  df_mult[[i]] <- df_mult[[i]][order(df_mult[[i]]$logFC, decreasing = TRUE), ]
}
names(df_mult) <- unique(pb_mult$subclusters)

saveRDS(res_mult, file = "multiome_UQ_edgeRResults.rds")
saveRDS(df_mult, file = "multiome_UQ_DFResults.rds")

df_degs <- lapply(df_mult, function(u) u[u$FDR < 0.05, ]) 

df_up_degs <- lapply(df_degs, function(u) u[u$logFC >= 0, ])

df_down_degs <- lapply(df_degs, function(u) u[u$logFC < 0, ])

## Myonuclei overall ####
myo <- filtered_sce[, filtered_sce$clusters == "Myonuclei"]

pb_myo <- aggregateAcrossCells(myo, use.assay.type = "counts", 
                                id = DataFrame(label = myo$clusters, sample = myo$dataset))
pb_myo <- logNormCounts(pb_myo)

y <- DGEList(counts(pb_myo), samples = colData(pb_myo))
y <- calcNormFactors(y, method = "upperquartile")
  
design <- model.matrix(~0 + y$samples$condition, y$samples)
colnames(design) <- c("C26", "CTRL")
  
contrast <- makeContrasts(C26 - CTRL, levels = design)
  
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
res_myo <- glmLRT(fit, contrast = contrast)
  
df_myo <- as.data.frame(topTags(res_myo, n = length(rownames(pb_myo))))
df_myo <- df_myo[order(df_myo$logFC, decreasing = TRUE), ]

saveRDS(df_myo, file = "multiome_myo_overall_UQ_DFResults.rds")
  
## C26 Trim63 vs CTRL Myonuclei overall ####
subset_myo <- myo[, c(myo)]
subset_myo1 <- subset_myo[, c(subset_myo$subclusters == "Myonuclei_Trim63" & subset_myo$condition == "CAC")]
subset_myo2 <- subset_myo[, c(subset_myo$subclusters != "Myonuclei_Trim63" & subset_myo$condition == "SHAM")]
subset_myo <- cbind(subset_myo1, subset_myo2)

subset_myo$clusters[subset_myo$subclusters == "Myonuclei_Trim63"] <- "Myonuclei_Trim63"

# Create pseudo-bulk 
pb_mult <- aggregateAcrossCells(subset_myo, use.assay.type = "counts", 
                                id = DataFrame(label = subset_myo$clusters, sample = subset_myo$dataset))
pb_mult <- logNormCounts(pb_mult)
colnames(pb_mult) <- paste(pb_mult$clusters, pb_mult$dataset, sep="_")

y <- DGEList(counts(pb_mult), samples = colData(pb_mult))
y <- calcNormFactors(y, method = "upperquartile")
design <- model.matrix(~0 + y$samples$condition, y$samples)
colnames(design) <- c("C26", "CTRL")
contrast <- makeContrasts(C26 - CTRL, levels = design)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
res_mult <- glmLRT(fit, contrast = contrast)
df_trim63 <- as.data.frame(topTags(res_mult, n = length(rownames(pb_mult))))
df_trim63 <- df_trim63[order(df_trim63$logFC, decreasing = TRUE), ]

df_trim63 <- cbind(rownames(df_trim63), df_trim63)
colnames(df_trim63)[1] <- "Gene_Name_ID"

writexl::write_xlsx(df_trim63, path = "multiome_DEGlist_c26trim63_ctrlMyo.xlsx")

## Upset plot ####
dea_out <- df_mult
dea_out <- dea_out[grepl("Myonuclei", names(dea_out))]

dea_out <- lapply(dea_out, function(x) {
  x$Gene_Name_ID <- rownames(x)
  x
})

dea_out$Myonuclei_Trim63 <- NULL
dea_out$Myonuclei_C26Trim63_CTRLOverall <- df_trim63

degs <- lapply(dea_out, function(x) x[x$FDR < 0.05, ]$Gene_Name_ID)

df.edgeRList <- fromList(degs)

p1 <- upset(df.edgeRList, colnames(df.edgeRList)[1:7],
            sort_intersections_by = c("degree", "cardinality"),
            sort_intersections = "ascending",
            name = "", width_ratio = 0.1, keep_empty_groups = TRUE,
            queries = list(
              upset_query(set = "Myonuclei_IIb", fill = "#cbe7f5"),
              upset_query(set = "Myonuclei_IIx", fill = "#60adff"),
              upset_query(set = "Myonuclei_IIx_IIa", fill = "#ff7068"),
              upset_query(set = "Myonuclei_IIx_IIb", fill = "#ffce99"),
              upset_query(set = "Myonuclei_MTJ", fill = "#0000eb"),
              upset_query(set = "Myonuclei_NMJ", fill = "#fa983d"),
              upset_query(set = "Myonuclei_C26Trim63_CTRLOverall", fill = "#e19fc0")
            ),
            intersections = list(
              # Unique DEGs were visualized for each cell-type
              "Myonuclei_IIb", "Myonuclei_IIx", "Myonuclei_IIx_IIa", "Myonuclei_IIx_IIb",
              "Myonuclei_MTJ", "Myonuclei_NMJ", "Myonuclei_C26Trim63_CTRLOverall", 
              
              # Intersection of all cell-types was visualized
              c("Myonuclei_C26Trim63_CTRLOverall", "Myonuclei_IIb"),
              c("Myonuclei_IIx", "Myonuclei_C26Trim63_CTRLOverall"),
              c("Myonuclei_IIx_IIa", "Myonuclei_C26Trim63_CTRLOverall"),
              c("Myonuclei_IIx_IIb","Myonuclei_C26Trim63_CTRLOverall"),
              c("Myonuclei_MTJ", "Myonuclei_C26Trim63_CTRLOverall"),
              c("Myonuclei_NMJ", "Myonuclei_C26Trim63_CTRLOverall")
            ),
            base_annotations = list(
              "Intersection size" = (
                intersection_size(
                  bar_number_threshold = 1, # show all numbers on top of bars
                  width = 0.4,  # reduce width of the bars
                  text = list(size = 4)
                )
                # add some space on the top of the bars
                + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
                + theme(
                  # hide grid lines
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  # show axis lines
                  axis.line = element_line(colour = "black")
                )
              )
            ),
            stripes = upset_stripes(
              geom = geom_segment(size = 3),
              colors = c("white", "white")
            ),
            matrix = intersection_matrix(
              geom = geom_point(
                shape = "circle filled",
                size = 3,
                stroke = 0
              )
            ),
            set_sizes = (
              upset_set_size(geom = geom_bar(width = 0.5), filter_intersections = TRUE)
              + theme(
                axis.line.x = element_line(colour = "black"),
                axis.ticks.x = element_line()
              )
            ),
            themes = upset_default_themes(text = element_text(size = 10, color = "black"))) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

ggsave(p1,  file = "multiome_upsetplot_myo.pdf", width = 14, height = 7)

# Differential Abundance analysis####
abundances <- table(sce_mult$subclusters, sce_mult$dataset)
abundances <- unclass(abundances) 
abundances

df <- as.data.frame(abundances)
df <- cbind(rowSums(df[,c(1,3)]), rowSums(df[,c(2,4)]))
colnames(df) <- c("C26", "CTRL")

# Create a df for barplot
ct <- rep(rownames(df), 2)
cond <- c(rep("C26", 15), rep("CTRL",15))
value <- c(df[,1], df[,2])

df_barplot <- data.frame(ct, cond, value)

# Barplot
p <- ggplot(df_barplot, aes(y=value, x=cond, fill=ct, color=ct)) + 
  geom_bar(position="fill", stat="identity") + scale_y_continuous(labels = scales::percent) + 
  scale_color_manual(values=palette_color)+
  scale_fill_manual(values=palette_color) + labs(y="Nuclei proportion(%)", x="", fill=NULL, color=NULL) +theme_classic()  #+theme(legend.position="none")
ggsave(p, file="multiome_barplot.jpg")

### Performing the DA analysis ####
extra_info <- as.data.frame(cbind(condition=c("C26","CTRL", "C26", "CTRL"), 
                                  sample.id=c("REP1C26", "REP1CTRL", "REP2C26", "REP2CTRL")))
rownames(extra_info) <- extra_info$sample.id

y_ab <- DGEList(abundances, samples = extra_info)
# We filter out low-abundance labels 
keep <- filterByExpr(y_ab, group = extra_info$condition)
y_ab <- y_ab[keep, ]

design <- model.matrix(~0 + factor(condition), y_ab$samples)
colnames(design) <- c("C26", "CTRL")

y_ab <- estimateDisp(y_ab, design, trend = "none")

fit_ab <- glmQLFit(y_ab, design,robust = TRUE, abundance.trend = FALSE)

contrast <- makeContrasts(C26 - CTRL, levels = design)
res1_contr <- glmQLFTest(fit_ab, contrast = contrast)
summary(decideTests(res1_contr))
topTags(res1_contr)

#### Assuming most labels do not change #####
y_ab2 <- calcNormFactors(y_ab)

y_ab2 <- estimateDisp(y_ab2, design, trend = "none")
fit_ab2 <- glmQLFit(y_ab2, design, robust = TRUE, abundance.trend = FALSE)
res2_contr <- glmQLFTest(fit_ab2, contrast = contrast)
summary(decideTests(res2_contr))
topTags(res2_contr, n = 15)
