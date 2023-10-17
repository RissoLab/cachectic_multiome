.libPaths(new = "/mnt/ganimede/zuin/R/x86_64-pc-linux-gnu-library/4.3")

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
library(muscat)

# For reproducibility
set.seed(8569)

mnn <- readRDS("multiome_rna_mnn.rds")

# Create pseudo-bulk 
pb_mult_rna <- aggregateAcrossCells(mnn, use.assay.type = "counts", 
                                    id=DataFrame(label=mnn$subclusters, sample=mnn$sample_id))
pb_mult_rna <- logNormCounts(pb_mult_rna)

colnames(pb_mult_rna) <- paste(pb_mult_rna$subclusters, pb_mult_rna$sample_id, sep="_")
# Save pseudo-bulk
saveRDS(pb_mult_rna, file = "multiome_rna_pb.rds")

## MDS plot ####
source("pbMDSfunction.R")

load("SubclusterPaletteColor.RData")
label_color <- subcluster_palette[-7,]
label_color$Label[label_color$Label=="Myonuclei_MTJ_IIb"] <- "Myonuclei_MTJ"

add_color <- c("#9ef8f2","#bfffa5", "#e19fc0","#fff68f" )
names(add_color) <- c("Myonuclei_IIa", "Myonuclei_IIx_IIa_IIb", "Myonuclei_Trim63", "Myonuclei_Xirp2")

palette_color <- label_color$Color
names(palette_color) <- label_color$Label
palette_color <- c(palette_color, add_color)

palette_color <- palette_color[order(names(palette_color))]

# MDS plot according to the Glutamatergic labels
prep_sce <- prepSCE(pb_mult_rna,
                    kid = "subclusters_new", # subpopulation assignments
                    gid = "condition",  # group IDs
                    sid = "sample_id",   # sample_id IDs
                    drop = TRUE)
pb <-  aggregateData(prep_sce, assay = "logcounts",
                     by = c("cluster_id", "sample_id"))

mds_plot <- pbMDS(pb) + scale_color_manual(values = palette_color) +
  labs(col = "Cell-types", shape = "Condition") +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 7))
ggsave(mds_plot, file="multiome_rna_mdsplot.jpg")

prep_sce <- prepSCE(pb_mult_rna[,grep("Myo", pb_mult_rna$subclusters_new)],
                    kid = "subclusters_new", # subpopulation assignments
                    gid = "condition",  # group IDs
                    sid = "sample_id",   # sample_id IDs
                    drop = TRUE)
pb <-  aggregateData(prep_sce, assay = "logcounts",
                     by = c("cluster_id", "sample_id"))
myo_palette <- palette_color[5:13]

mds_plot <- pbMDS(pb) + scale_color_manual(values = myo_palette) +
  labs(col = "Cell-types", shape = "Condition") +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 7))
ggsave(mds_plot, file="multiome_rna_myo_mdsplot.jpg")

# DEA ####
# Negative binomial generalized model (GLM) with Upper Quantile normalization ####
res_mult_rna <- df_mult_rna <- list()
for (i in seq_along(unique(pb_mult_rna$subclusters_new))) {
  ct_pb <- pb_mult_rna[, pb_mult_rna$subclusters_new == unique(pb_mult_rna$subclusters_new)[i]]
  
  y <- DGEList(counts(ct_pb), samples = colData(ct_pb))
  keep <- filterByExpr(y, group = ct_pb$condition)
  y <- y[keep, ]
  y <- calcNormFactors(y, method = "upperquartile")
  
  design <- model.matrix(~0 + y$samples$condition, y$samples)
  colnames(design) <- c("CAC", "SHAM")
  
  contrast <- makeContrasts(CAC - SHAM,levels = design)
  
  y <- estimateDisp(y, design)
  
  fit <- glmFit(y, design)
  
  res_mult_rna[[i]] <- glmLRT(fit, contrast = contrast)
  
  df_mult_rna[[i]] <- as.data.frame(res_mult_rna[[i]]$table)
  
  df_mult_rna[[i]] <- as.data.frame(df_mult_rna[[i]])[order(as.data.frame(df_mult_rna[[i]])$PValue), ]
  
  FDR <- as.data.frame(topTags(res_mult_rna[[i]], n = length(rownames(pb_mult_rna)))[, 5])
  
  df_mult_rna[[i]] <- cbind(df_mult_rna[[i]], FDR)
  df_mult_rna[[i]] <- df_mult_rna[[i]][order(df_mult_rna[[i]]$FDR), ]
}
names(df_mult_rna) <- unique(pb_mult_rna$subclusters_new)

saveRDS(res_mult_rna, file = "multiome_rna_UQ_edgeRResults.rds")
saveRDS(df_mult_rna, file = "multiome_rna_UQ_DFResults.rds")

# Results #####
n_genes <- lapply(df_mult_rna, function(u) length(rownames(u))) 
df_degs <- lapply(df_mult_rna, function(u) u[u$FDR<0.05,]) 
n_degs <- lapply(df_degs, function(u) length(rownames(u))) 

df_up_degs <- lapply(df_degs, function(u) u[u$logFC>=0,])
n_up_degs <-  lapply(df_up_degs, function(u) length(rownames(u))) 

df_down_degs <- lapply(df_degs, function(u) u[u$logFC<0,])
n_down_degs <-  lapply(df_down_degs, function(u) length(rownames(u)))

n_nuclei <- pb_mult_rna$ncells
n_nuclei <- c(n_nuclei[1:75],0, n_nuclei[76:101])
n_nuclei <- matrix(n_nuclei, 17,6,byrow=TRUE) 
n_nuclei <- rowSums(n_nuclei) # N. of nuclei

tab <- cbind(n_nuclei, as.matrix(n_genes), as.matrix(n_degs), as.matrix(n_up_degs), as.matrix(n_down_degs))
colnames(tab) <- c("nuclei","genes","n_DEGs", "Up","Down")
# colnames(tab) <- c("N. of nuclei","N. of detected genes","N. of DEGs", "N. of Up DEGs","N. of Down DEGs")
rownames(tab) <- levels(factor(pb_mult_rna$subclusters_new))

# Save DE list
df_save <- lapply(df_degs, function(x){
  x$Gene_Name_ID <- rownames(x); x
}  )

df_save <- lapply(df_save, function(x)
  as.data.frame(cbind(x$Gene_Name_ID,
                      x$logFC, x$logCPM, x$LR, x$PValue, x$FDR))
)

df_save <- lapply(df_save, function(x) {
  colnames(x) <- c("Gene_Name_ID", "logFC","logCPM","LR", "PValue", "FDR");x
})
writexl::write_xlsx(df_save, "multiome_rna_DEList.xlsx")


df_save <- lapply(df_up_degs, function(x){
  x$Gene_Name_ID <- rownames(x); x
}  )

df_save <- lapply(df_save, function(x)
  as.data.frame(cbind(x$Gene_Name_ID,
                      x$logFC, x$logCPM, x$LR, x$PValue, x$FDR))
)
df_save <- lapply(df_save, function(x) {
  colnames(x) <- c("Gene_Name_ID", "logFC","logCPM","LR", "PValue", "FDR");x
})
writexl::write_xlsx(df_save, "multiome_rna_Up_DEList.xlsx")

df_save <- lapply(df_down_degs, function(x){
  x$Gene_Name_ID <- rownames(x); x
}  )

df_save <- lapply(df_save, function(x)
  as.data.frame(cbind(x$Gene_Name_ID,
                      x$logFC, x$logCPM, x$LR, x$PValue, x$FDR))
)
df_save <- lapply(df_save, function(x) {
  colnames(x) <- c("Gene_Name_ID", "logFC","logCPM","LR", "PValue", "FDR");x
})
writexl::write_xlsx(df_save, "multiome_rna_Down_DEList.xlsx")

## Scatter plot####
tab <- as.data.frame(tab)
tab$Label <- factor(rownames(tab))

tab$n_DEGs <- as.numeric(tab$n_DEGs)
tab$nuclei <- as.numeric(tab$nuclei)

load("SubclusterPaletteColor.RData")
label_color <- subcluster_palette[-7,]
label_color$Label[label_color$Label=="Myonuclei_MTJ_IIb"] <- "Myonuclei_MTJ"

add_color <- c("#9ef8f2","#bfffa5", "#e19fc0","#fff68f" )
names(add_color) <- c("Myonuclei_IIa", "Myonuclei_IIx_IIa_IIb", "Myonuclei_Trim63", "Myonuclei_Xirp2")

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
ggsave(p1, file="multiome_rna_scatterplot.jpg",  width = 20, height = 20, units = "cm")

# PCA ####
label <- levels(factor(pb_mult_rna$subclusters_new))

ct_pb <- set_mult<-  list()
for (i in seq_along(unique(pb_mult_rna$subclusters_new))) {
  ct_pb[[i]] <- pb_mult_rna[,pb_mult_rna$subclusters_new == unique(pb_mult_rna$subclusters_new)[i]]
  set_mult[[i]] <- newSeqExpressionSet(as.matrix(round(counts(ct_pb[[i]]))), 
                                       phenoData = data.frame(colData(ct_pb[[i]]),
                                                              row.names = colnames(ct_pb[[i]])))
  set_mult[[i]] <- betweenLaneNormalization(set_mult[[i]], which = "upper")
} 

file_name <- paste("multiome_rna_PCA_part1.jpg", sep = "")
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

file_name <- paste("multiome_rna_PCA_part2.jpg", sep = "")
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

file_name <- paste("multiome_rna_PCA_part3.jpg", sep = "")
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
plotPCA(set_mult[[16]], xlim = c(-1, 1), label=FALSE, pch=20, 
        col=as.numeric(as.factor(set_mult[[16]]$condition)), cex=2, 
        main= paste(label[16],sep = ""))
plotPCA(set_mult[[17]], xlim = c(-1, 1),label=FALSE,  pch=20, 
        col=as.numeric(as.factor(set_mult[[17]]$condition)), cex=2, 
        main= paste(label[17],sep = ""))

legend("right", legend=levels(factor(set_mult[[13]]$condition)),  pch=20, cex=1,
       col= levels(as.factor(as.numeric(as.factor(set_mult[[13]]$condition)))), xpd = TRUE, 
       horiz = FALSE,  bty = "n",title="Condition")
dev.off()

# Error in svd(Y) : infinite or missing values in 'x'

# Differential Abundance ####
abundances <- table(mnn$subclusters_new, mnn$sample_id)
abundances <- unclass(abundances) 
abundances

df <- as.data.frame(abundances)
df <- cbind(rowSums(df[,c(1,3,5)]), rowSums(df[,c(2,4,6)]))
colnames(df) <- c("CAC", "SHAM")

# Create a df for barplot
ct <- rep(rownames(df), 2)
cond <- c(rep("CAC", 17), rep("SHAM",17))
value <- c(df[,1], df[,2])

df_barplot <- data.frame(ct, cond, value)

# Barplot
p <- ggplot(df_barplot, aes(y=value, x=cond, fill=ct, color=ct)) + 
  geom_bar(position="fill", stat="identity") + scale_y_continuous(labels = scales::percent) + 
  scale_color_manual(values=palette_color)+
  scale_fill_manual(values=palette_color) + labs(y="Nuclei proportion(%)", x="", fill=NULL, color=NULL) +theme_classic()  #+theme(legend.position="none")
ggsave(p, file="multiome_rna_barplot.jpg")

### Performing the DA analysis ####
extra_info <- as.data.frame(cbind(condition=c("CAC","SHAM", "CAC", "SHAM", "CAC", "SHAM"), 
                                  sample.id=c("REP1CAC", "REP1SHAM", "REP2CAC", "REP2SHAM",
                                              "RNACAC", "RNASHAM")))
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
topTags(res2_contr)

# Gene Set Enrichmente analysis####
df_mult_rna <- readRDS("multiome_rna_UQ_DFResults.rds")

background <- lapply(df_mult_rna, function(x) rownames(x))

up_degs <- lapply(df_up_degs, function(x) rownames(x))
down_degs <- lapply(df_down_degs, function(x) rownames(x))

# Load the GMT file coming from Muscle Gene Sets
upload_GMT_file(gmtfile = "MuscleGeneSets_v3_FINAL1517.gmt")
# Note: every time the file was uploaded, the custom annotations ID changed

gostres_up <- gostres_down <- list()
for (i in seq_along(up_degs)) {
  gostres_up[[i]] <- gost(query = up_degs[[i]], organism = "gp__tWEE_oRjg_PqE", custom_bg = background[[i]])
  gostres_down[[i]] <- gost(query = down_degs[[i]], organism = "gp__tWEE_oRjg_PqE", custom_bg = background[[i]])
}
save(gostres_down, gostres_up, file = "multiome_rna_GSEA.RData")

gsea_up <- lapply(gostres_up, function(x) x$result)
names(gsea_up) <- names(up_degs)

writexl::write_xlsx(gsea_up, "multiome_rna_GSEA_Up_output.xlsx")

gsea_up_sum <- lapply(gsea_up, function(x) x[,c(3,9:11)])
gsea_up_sum <- lapply(gsea_up_sum, function(x) as.data.frame(cbind(x$source, x$term_id, x$term_name, x$p_value)))
gsea_up_sum <- lapply(gsea_up_sum, function(x){
  colnames(x)=c("source", "term_id", "term_name","p_value");x})

writexl::write_xlsx(gsea_up_sum, "multiome_rna_GSEA_Up.xlsx")

gsea_down <- lapply(gostres_down, function(x) as.data.frame(x$result))
names(gsea_down) <- names(up_degs)

writexl::write_xlsx(gsea_down, "multiome_rna_GSEA_Down_output.xlsx")

sum_down <- lapply(gsea_down, function(x) x[,c(3,9:11)])
sum_down <- lapply(sum_down, function(x) as.data.frame(cbind(x$source, x$term_id, x$term_name, x$p_value)))
sum_down <- lapply(sum_down, function(x){colnames(x)=c("source", "term_id", "term_name","p_value");x})
writexl::write_xlsx(sum_down, "multiome_rna_GSEA_Down.xlsx")