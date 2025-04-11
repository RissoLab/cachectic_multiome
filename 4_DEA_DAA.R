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
library(gprofiler2)
library(UpSetR)
library(ComplexUpset)

# For reproducibility
set.seed(8569)

# This file contains the differential expression analysis (DEA) and differential abundance analysis (DAA)
sgl_mult <- readRDS("multiome_sgl_comb.rds")
DefaultAssay(sgl_mult) <- "RNA"
# From Seurat to SingleCellExperiment object
sce_mult <- as.SingleCellExperiment(sgl_mult)

# genes filtering: 10% for each replicate
keep <- apply(counts(sce_mult), 1, function(y) all(tapply(y, sce_mult$dataset, function(x) mean(x>0))>=0.1))

filtered_sce <- sce_mult[keep, ]                                                

# Create pseudo-bulk 
pb_mult <- aggregateAcrossCells(filtered_sce, use.assay.type = "counts", 
                                id=DataFrame(label=filtered_sce$subclusters, sample=filtered_sce$dataset))
pb_mult <- logNormCounts(pb_mult)

colnames(pb_mult) <- paste(pb_mult$subclusters, pb_mult$dataset, sep="_")
# Save pseudo-bulk
saveRDS(pb_mult, file = "multiome_pb.rds")

# DEA ####
# Negative binomial generalized model (GLM) with Upper Quantile normalization ####
res_mult <- df_mult <- list()
for (i in seq_along(unique(pb_mult$subclusters))) {
  ct_pb <- pb_mult[, pb_mult$subclusters == unique(pb_mult$subclusters)[i]]
  
  y <- DGEList(counts(ct_pb), samples = colData(ct_pb))
  # keep <- filterByExpr(y, group = ct_pb$condition)
  # y <- y[keep, ]
  y <- calcNormFactors(y, method = "upperquartile")
  
  design <- model.matrix(~0 + y$samples$condition, y$samples)
  colnames(design) <- c("C26", "CTRL")
  
  contrast <- makeContrasts(C26 - CTRL,levels = design)
  
  y <- estimateDisp(y, design)
  
  # fit <- glmQLFit(y, design)
  fit <- glmFit(y, design)
  
  res_mult[[i]] <- glmLRT(fit, contrast = contrast)
  
  # df_mult[[i]] <- as.data.frame(res_mult[[i]]$table)
  
  # df_mult[[i]] <- as.data.frame(df_mult[[i]])[order(as.data.frame(df_mult[[i]])$PValue), ]
  
  df_mult[[i]] <- as.data.frame(topTags(res_mult[[i]], n = length(rownames(pb_mult))))

  df_mult[[i]] <- df_mult[[i]][order(df_mult[[i]]$logFC, decreasing = TRUE), ]
}
names(df_mult) <- unique(pb_mult$subclusters)

saveRDS(res_mult, file = "multiome_UQ_edgeRResults.rds")
saveRDS(df_mult, file = "multiome_UQ_DFResults.rds")

# Results #####
n_genes <- lapply(df_mult, function(u) length(rownames(u))) 
df_degs <- lapply(df_mult, function(u) u[u$FDR<0.05,]) 
n_degs <- lapply(df_degs, function(u) length(rownames(u))) 

df_up_degs <- lapply(df_degs, function(u) u[u$logFC>=0,])
n_up_degs <-  lapply(df_up_degs, function(u) length(rownames(u))) 

df_down_degs <- lapply(df_degs, function(u) u[u$logFC<0,])
n_down_degs <-  lapply(df_down_degs, function(u) length(rownames(u)))

n_nuclei <- matrix(pb_mult$ncells,15,4,byrow=TRUE) 
n_nuclei <- rowSums(n_nuclei) # N. of nuclei

tab <- cbind(n_nuclei, as.matrix(n_genes), as.matrix(n_degs), as.matrix(n_up_degs), as.matrix(n_down_degs))
colnames(tab) <- c("nuclei","genes","n_DEGs", "Up","Down")
# colnames(tab) <- c("N. of nuclei","N. of detected genes","N. of DEGs", "N. of Up DEGs","N. of Down DEGs")
rownames(tab) <- levels(factor(pb_mult$subclusters))

# Save DE list
df_save <- lapply(df_degs, function(x){
  x$Gene_Name_ID <- rownames(x); x
}  )

df_save <- lapply(df_save, function(x) cbind(x[,6], x[,1:5]))
df_save <- lapply(df_save, function(x) {
  colnames(x) <- c("Gene_Name_ID", "logFC","logCPM","LR", "PValue", "FDR");x
})
writexl::write_xlsx(df_save, "multiome_DEGlist_after_filtering.xlsx")

df_save <- lapply(df_up_degs, function(x){
  x$Gene_Name_ID <- rownames(x); x
}  )

df_save <- lapply(df_save, function(x) cbind(x[,6], x[,1:5]))
df_save <- lapply(df_save, function(x) {
  colnames(x) <- c("Gene_Name_ID", "logFC","logCPM","LR", "PValue", "FDR");x
})
writexl::write_xlsx(df_save, "multiome_UP_DEList.xlsx")

df_save <- lapply(df_down_degs, function(x){
  x$Gene_Name_ID <- rownames(x); x
}  )

df_save <- lapply(df_save, function(x) cbind(x[,6], x[,1:5]))
df_save <- lapply(df_save, function(x) {
  colnames(x) <- c("Gene_Name_ID", "logFC","logCPM","LR", "PValue", "FDR");x
})
writexl::write_xlsx(df_save, "multiome_Down_DEList.xlsx")

# Myonuclei overall ###
myo <- filtered_sce[, filtered_sce$clusters=="Myonuclei"]

pb_myo <- aggregateAcrossCells(myo, use.assay.type = "counts", 
                                id=DataFrame(label=myo$clusters, sample=myo$dataset))
pb_myo <- logNormCounts(pb_myo)

  y <- DGEList(counts(pb_myo), samples = colData(pb_myo))
  # keep <- filterByExpr(y, group = ct_pb$condition)
  # y <- y[keep, ]
  y <- calcNormFactors(y, method = "upperquartile")
  
  design <- model.matrix(~0 + y$samples$condition, y$samples)
  colnames(design) <- c("C26", "CTRL")
  
  contrast <- makeContrasts(C26 - CTRL,levels = design)
  
  y <- estimateDisp(y, design)
  
  # fit <- glmQLFit(y, design)
  fit <- glmFit(y, design)
  
  res_myo <- glmLRT(fit, contrast = contrast)
  

  df_myo <- as.data.frame(topTags(res_myo, n = length(rownames(pb_myo))))
  
  df_myo <- df_myo[order(df_myo$logFC, decreasing = TRUE), ]

  saveRDS(df_myo, file = "multiome_myo_overall_UQ_DFResults.rds")
  
## Scatter plot####
tab <- as.data.frame(tab)
tab$Label <- factor(rownames(tab))

tab$n_DEGs <- as.numeric(tab$n_DEGs)
tab$nuclei <- as.numeric(tab$nuclei)

load("SubclusterPaletteColor.RData")
label_color <- subcluster_palette
label_color$Label[label_color$Label=="Myonuclei_MTJ_IIb"] <- "Myonuclei_MTJ"

add_color <- c( "#e19fc0")
names(add_color) <- c("Myonuclei_Trim63")

palette_color <- label_color$Color
names(palette_color) <- label_color$Label
palette_color <- c(palette_color, add_color)
palette_color <- palette_color[order(names(palette_color))]

p1 <- ggplot(tab, aes(x = log(nuclei), y = log(n_DEGs))) + 
  geom_point(aes(color = Label), size = 7.5) +
  scale_color_manual(values = palette_color) + 
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, col = "grey") +
  theme_classic(base_size = 18) + theme(legend.position = "none") +
  geom_label_repel(aes(label = Label), size=5,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', colour = palette_color)
ggsave(p1, file="multiome_scatterplot.jpg",  width = 20, height = 20, units = "cm")

# PCA ####
label <- levels(factor(pb_mult$subclusters))

ct_pb <- set_mult<-  list()
for (i in seq_along(unique(pb_mult$subclusters))) {
  ct_pb[[i]] <- pb_mult[,pb_mult$subclusters == unique(pb_mult$subclusters)[i]]
  set_mult[[i]] <- newSeqExpressionSet(as.matrix(round(counts(ct_pb[[i]]))), 
                                       phenoData = data.frame(colData(ct_pb[[i]]),
                                                              row.names = colnames(ct_pb[[i]])))
  set_mult[[i]] <- betweenLaneNormalization(set_mult[[i]], which = "upper")
} 

file_name <- paste("multiome_PCA_part1.jpg", sep = "")
jpeg(file = file_name,  width = 600, height = 900)
par(mfrow = c(3,2))
plotPCA(set_mult[[1]], xlim = c(-1, 1),label=FALSE, pch=20, 
        col=as.numeric(as.factor(set_mult[[1]]$condition)), cex=2, 
        main= paste(label[1],sep = ""))
plotPCA(set_mult[[2]], xlim = c(-1, 1), label=FALSE, pch=20, 
        col=as.numeric(as.factor(set_mult[[2]]$condition)), cex=2, 
        main= paste(label[2],sep = ""))
plotPCA(set_mult[[3]], xlim = c(-1, 1),label=FALSE,  pch=20, 
        col=as.numeric(as.factor(set_mult[[3]]$condition)), cex=2, 
        main= paste(label[3],sep = ""))

plotPCA(set_mult[[4]], xlim = c(-1, 1),label=FALSE,  pch=20, 
        col=as.numeric(as.factor(set_mult[[4]]$condition)), cex=2, 
        main= paste(label[4],sep = ""))
plotPCA(set_mult[[5]], xlim = c(-1, 1),label=FALSE,  pch=20, 
        col=as.numeric(as.factor(set_mult[[5]]$condition)), cex=2, 
        main= paste(label[5],sep = ""))
plotPCA(set_mult[[6]], xlim = c(-1, 1),label=FALSE,  pch=20, 
        col=as.numeric(as.factor(set_mult[[6]]$condition)), cex=2, 
        main= paste(label[6],sep = ""))
legend("right", legend=levels(factor(set_mult[[1]]$condition)),  pch=20, cex=1,
       col= levels(as.factor(as.numeric(as.factor(set_mult[[1]]$condition)))), xpd = TRUE, 
       horiz = FALSE,  bty = "n",title="Condition")
dev.off()

file_name <- paste("multiome_PCA_part2.jpg", sep = "")
jpeg(file = file_name,  width = 600, height = 900)
par(mfrow = c(3,2))
plotPCA(set_mult[[7]], xlim = c(-1, 1),label=FALSE, pch=20, 
        col=as.numeric(as.factor(set_mult[[7]]$condition)), cex=2, 
        main= paste(label[7],sep = ""))
plotPCA(set_mult[[8]], xlim = c(-1, 1), label=FALSE, pch=20, 
        col=as.numeric(as.factor(set_mult[[8]]$condition)), cex=2, 
        main= paste(label[8],sep = ""))
plotPCA(set_mult[[9]], xlim = c(-1, 1),label=FALSE,  pch=20, 
        col=as.numeric(as.factor(set_mult[[9]]$condition)), cex=2, 
        main= paste(label[9],sep = ""))

plotPCA(set_mult[[10]], xlim = c(-1, 1),label=FALSE,  pch=20, 
        col=as.numeric(as.factor(set_mult[[10]]$condition)), cex=2, 
        main= paste(label[10],sep = ""))
plotPCA(set_mult[[11]], xlim = c(-1, 1),label=FALSE,  pch=20, 
        col=as.numeric(as.factor(set_mult[[11]]$condition)), cex=2, 
        main= paste(label[11],sep = ""))
plotPCA(set_mult[[12]], xlim = c(-1, 1),label=FALSE,  pch=20, 
        col=as.numeric(as.factor(set_mult[[12]]$condition)), cex=2, 
        main= paste(label[12],sep = ""))
legend("right", legend=levels(factor(set_mult[[7]]$condition)),  pch=20, cex=1,
       col= levels(as.factor(as.numeric(as.factor(set_mult[[7]]$condition)))), xpd = TRUE, 
       horiz = FALSE,  bty = "n",title="Condition")
dev.off()

file_name <- paste("multiome_PCA_part3.jpg", sep = "")
jpeg(file = file_name,  width = 600, height = 600)
par(mfrow = c(2,2))
plotPCA(set_mult[[13]], xlim = c(-1, 1),label=FALSE, pch=20, 
        col=as.numeric(as.factor(set_mult[[13]]$condition)), cex=2, 
        main= paste(label[13],sep = ""))
plotPCA(set_mult[[14]], xlim = c(-1, 1), label=FALSE, pch=20, 
        col=as.numeric(as.factor(set_mult[[14]]$condition)), cex=2, 
        main= paste(label[14],sep = ""))
plotPCA(set_mult[[15]], xlim = c(-1, 1),label=FALSE,  pch=20, 
        col=as.numeric(as.factor(set_mult[[15]]$condition)), cex=2, 
        main= paste(label[15],sep = ""))

legend("right", legend=levels(factor(set_mult[[13]]$condition)),  pch=20, cex=1,
       col= levels(as.factor(as.numeric(as.factor(set_mult[[13]]$condition)))), xpd = TRUE, 
       horiz = FALSE,  bty = "n",title="Condition")
dev.off()

# Error in svd(Y) : infinite or missing values in 'x'

# Differential Abundance analysis####
abundances <- table(sce_mult$subclusters, sce_mult$dataset)
abundances <- unclass(abundances) 
abundances

df <- as.data.frame(abundances)
df <- cbind(rowSums(df[,c(1,3)]), rowSums(df[,c(2,4)]))
colnames(df) <- c("CAC", "SHAM")

# Create a df for barplot
ct <- rep(rownames(df), 2)
cond <- c(rep("CAC", 15), rep("SHAM",15))
value <- c(df[,1], df[,2])

df_barplot <- data.frame(ct, cond, value)

# Barplot
p <- ggplot(df_barplot, aes(y=value, x=cond, fill=ct, color=ct)) + 
  geom_bar(position="fill", stat="identity") + scale_y_continuous(labels = scales::percent) + 
  scale_color_manual(values=palette_color)+
  scale_fill_manual(values=palette_color) + labs(y="Nuclei proportion(%)", x="", fill=NULL, color=NULL) +theme_classic()  #+theme(legend.position="none")
ggsave(p, file="multiome_barplot.jpg")

### Performing the DA analysis ####
extra_info <- as.data.frame(cbind(condition=c("CAC","SHAM", "CAC", "SHAM"), 
                                  sample.id=c("REP1CAC", "REP1SHAM", "REP2CAC", "REP2SHAM")))
rownames(extra_info) <- extra_info$sample.id

y_ab <- DGEList(abundances, samples=extra_info)
# We filter out low-abundance labels 
keep <- filterByExpr(y_ab, group=extra_info$condition)
y_ab <- y_ab[keep,]

# design <- model.matrix(~factor(condition), y_ab$samples) # Without create contrast
# colnames(design) <- c("Intercept", "SHAM")

design <- model.matrix(~0+factor(condition), y_ab$samples)
colnames(design) <- c("CAC", "SHAM")

y_ab <- estimateDisp(y_ab, design, trend="none")

fit_ab <- glmQLFit(y_ab, design,robust=TRUE, abundance.trend=FALSE)

# res <- glmQLFTest(fit_ab, coef=2)
contrast <- makeContrasts(CAC-SHAM,levels=design)
res1_contr <- glmQLFTest(fit_ab, contrast =contrast)
summary(decideTests(res1_contr))
topTags(res1_contr)

#### Assuming most labels do not change #####
y_ab2 <- calcNormFactors(y_ab)

y_ab2 <- estimateDisp(y_ab2, design, trend="none")
fit_ab2 <- glmQLFit(y_ab2, design, robust=TRUE, abundance.trend=FALSE)
res2_contr <- glmQLFTest(fit_ab2, contrast=contrast)
summary(decideTests(res2_contr))
topTags(res2_contr, n = 15)

# Gene Set Enrichmente analysis####
background <- lapply(df_mult, function(x) rownames(x))

up_degs <- lapply(df_up_degs, function(x) rownames(x))
down_degs <- lapply(df_down_degs, function(x) rownames(x))

# Load the GMT file coming from Muscle Gene Sets
upload_GMT_file(gmtfile = "MuscleGeneSets_v3_FINAL1517.gmt")
# Note: every time the file was uploaded, the custom annotations ID changed

gostres_up <- gostres_down <- list()
for (i in seq_along(up_degs)) {
  gostres_up[[i]] <- gost(query = up_degs[[i]], organism = "gp__67Aw_iaWe_s94", custom_bg = background[[i]])
  gostres_down[[i]] <- gost(query = down_degs[[i]], organism = "gp__67Aw_iaWe_s94", custom_bg = background[[i]])
}
save(gostres_down, gostres_up, file = "multiome_GSEA.RData")

gsea_up <- lapply(gostres_up, function(x) x$result)
names(gsea_up) <- names(up_degs)

writexl::write_xlsx(gsea_up, "multiome_GSEA_Up.xlsx")

sum_up <- lapply(gsea_up, function(x) x[,c(3,9:11)])
sum_up <- lapply(sum_up, function(x) as.data.frame(cbind(x$source, x$term_id, x$term_name, x$p_value)))
sum_up <- lapply(sum_up, function(x){colnames(x)=c("source", "term_id", "term_name","p_value");x})
writexl::write_xlsx(sum_up, "multiome_GSEA_Up_sum.xlsx")

gsea_down <- lapply(gostres_down, function(x) x$result)
names(gsea_down) <- names(up_degs)

writexl::write_xlsx(gsea_down, "multiome_GSEA_Down.xlsx")

sum_down <- lapply(gsea_down, function(x) x[,c(3,9:11)])
sum_down <- lapply(sum_down, function(x) as.data.frame(cbind(x$source, x$term_id, x$term_name, x$p_value)))
sum_down <- lapply(sum_down, function(x){colnames(x)=c("source", "term_id", "term_name","p_value");x})
writexl::write_xlsx(sum_down, "multiome_GSEA_Down_sum.xlsx")
