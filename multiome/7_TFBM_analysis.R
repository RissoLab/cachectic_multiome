# =============================================================================
# Cachectic multiome (C26 and CTRL samples)
# Transcription Factor Binding Motifs
# =============================================================================

library(TENxMultiomeTools)
library(SingleCellExperiment)
library(DEScan2)
library(GenomicRanges)
library(rtracklayer)
library(ChIPpeakAnno)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(RColorBrewer)
library(GenomicFeatures)
library(ChIPseeker)
library(edgeR)
library(biomaRt)
library(DEScan2)
library(RUVSeq)
library(motifTestR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(MotifDb)
library(scuttle)

bams <- list.files("./libraries", pattern="*bam$", recursive = TRUE, full.names = TRUE)
bams <- unique(dirname(bams)[-c(9,10)])

scelist <- readRDS("cachectic_mult_scelist_annot.rds")

label <- levels(factor(scelist[[1]]$subclusters))
ct_scelist <- scelist

# For each cell types
for (j in seq_along(label)) {
  for (i in seq_along(ct_scelist)) {
    ct_scelist[[i]]$dirbam <- bams[i]
    ct_scelist[[i]] <- ct_scelist[[i]][, ct_scelist[[i]]$subclusters == label[j]]
  }
  
  lapply(ct_scelist, function(sce){
    sce$Barcode <- sce$gex_barcode
    sce$Sample <- sce$sample_id
    lapply(unique(sce$subclusters), function(ct){
      createBamCt(sce, cellType = ct,
                  cellTypesCol = "subclusters",
                  bamdir = unique(sce$dirbam), bamType = "ATAC",
                  outdir = paste0("/mnt/ganimede/zuin/SandriProj/bams/", tolower(label[j])), ncores = 10)
    })
  })
}

# Peak calling (MACS2) ####
label_lower <- tolower(label)
bamfiles <- lapply(label_lower, function(x) list.files(label_lower, pattern = "*_sorted.bam$", full.names = TRUE, recursive=TRUE))

outdir <- "./MACS2_peaks"
additional.params <- "--nomodel "
extsize <- 200
shift  <- -(extsize/2)

for (i in seq_along(bamfiles)) {
  bamfile <- bamfiles[[i]]
  for (file in bamfile) {
    cmd <- paste0("macs2 callpeak --treatment ", file,
                  " --name ", basename(file),
                  " --outdir ", outdir,
                  " --format BAM --gsize mm ",
                  additional.params,
                  " --extsize ", extsize,
                  " --shift ", shift,
                  " --broad")
    
    print(cmd)
    system(cmd)
  }
}

# Peaks alignment ####
# Align peaks to form common regions then filter regions for presence in multiple
# replicates taking in input a GRangesList where each element is a sample of called
# peaks.
bpeaks <- list.files(paste0("./MACS2_peaks/", label[1]), pattern=".broadPeak$", full.names=TRUE)

grl <- lapply(bpeaks, import, format = "broadPeak")
names(grl) <- gsub("_sorted.bam_peaks.broadPeak", "", basename(bpeaks))

ccs <- c("C26","CTRL")
grl_ll <- lapply(ccs, function(c) GRangesList(grl[grep(c, names(grl))]))
names(grl_ll) <- ccs

rrlistcond <- lapply(grl_ll, function(c)
  finalRegions(c,
               zThreshold=0, minCarriers=1, saveFlag=FALSE,
               scorecolname="score", verbose=TRUE,
               BPPARAM=BiocParallel::MulticoreParam(workers=1))
)

# finalRegions lists
c26 <- rrlistcond[[1]]
ctrl <- rrlistcond[[2]]

# Peak filtering
c26 <- c26[c26$`k-carriers` >=2, ]
ctrl <- ctrl[ctrl$`k-carriers` >=2, ]

region_ll <- findOverlapsOfPeaks(c26, ctrl)

# unique peaks per condition
unique_c26 <- region_ll$uniquePeaks[grep("c26", names(region_ll$uniquePeaks)), ]
unique_c26$details <- "C26"

unique_ctrl <- region_ll$uniquePeaks[grep("ctrl", names(region_ll$uniquePeaks)), ]
unique_ctrl$details <- "CTRL"

# overlapped peaks
merged <- region_ll$mergedPeaks
rr <- as.data.frame(merged)
names(merged) <- paste0(seqnames(merged), ":", rr$start, "-", rr$end)
mcols(merged) <- NULL
merged$details <- "merged"

# consensus peak list
allpeaks <- sort(c(unique_c26, unique_ctrl, merged))

gtftxdb <- makeTxDbFromGFF(file="./refdata-gex-mm10-2020-A/genes/genes.gtf", format="gtf")

annotPeak <- annotatePeak(allpeaks, tssRegion=c(-5000, 5000), TxDb=gtftxdb, annoDb="org.Mm.eg.db")
annotPeak <- annotPeak@anno

save(region_ll, annotPeak, file = paste0("peaks_FAPs.RData"))

# TFBM analysis #####
sce_mult <- readRDS("cachectic_mult_sce.rds")

dea_out <- readRDS("cachectic_mult_degs_ctrl_c26.rds")
up_degs <- lapply(dea_out, function(x) rownames(x[c(x$FDR < 0.05 & x$logFC >=0), ]))

trim63 <- readxl::read_xlsx("cachectic_mult_degs_c26trim63_ctrlMyo.xlsx")
trim63 <- trim63[c(trim63$FDR < 0.05 & trim63$logFC >=0), ]
up_degs$Myonuclei_Trim63 <- trim63$Gene_Name_ID

myo_overall <- readRDS("cachectic_mult_degs_myo_ctrl_c26.rds")
myo_overall <- myo_overall[c(myo_overall$FDR < 0.05 & myo_overall$logFC >=0), ]
up_degs$Myonuclei_Overall <- rownames(myo_overall)

# Reference 
ex_pfm <- MotifDb |>
  subset(organism == "Mmusculus") |>
  query("jaspar2022") |>
  as.list() 
names(ex_pfm) <- gsub(".+jaspar2022-(.+)-.+", "\\1", names(ex_pfm))

# Peaks 
ct <- "Tenocyte"

load(paste0("peaks_DARs_", ct, ".RData"))

select_degs <- up_degs$Tenocyte 

annotPro_5kb <- annotPeak[annotPeak$details != "CTRL", ]
annotPro_5kb <- annotPro_5kb[grepl("Promoter", annotPro_5kb$annotation), ]
annotPro_1kb <- annotPro_5kb[annotPro_5kb$annotation == "Promoter (<=1kb)", ]

annotPro_1kb <- annotPro_1kb[annotPro_1kb$SYMBOL %in% select_degs, ]
annotPro_5kb <- annotPro_5kb[annotPro_5kb$SYMBOL %in% select_degs, ]

length(unique(annotPro_1kb$SYMBOL))
length(unique(annotPro_5kb$SYMBOL))

length(annotPro_1kb)
length(annotPro_5kb)

ll <- list("1kb" = annotPro_1kb, "5kb" = annotPro_5kb)
test <- lapply(ll, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

match_tf <- lapply(test, function(x) {
  out <- getPwmMatches(ex_pfm, x, best_only = TRUE, break_ties = "all")
  names(out) <- paste(names(out), seq(1, length(names(out))), sep = "_")
  out
})

remove_zero <- function(obj) {
  n <- lapply(obj, function(x) length(rownames(as.data.frame(x))))
  n <- do.call(rbind, n)
  n <- n[n>0, ] 
  obj <- obj[names(n)]
  names(obj) <- sub("\\_.*", "", names(obj))
  return(obj)
}

match_tf <- lapply(match_tf, function(x) remove_zero(x))
match_tf <- lapply(match_tf, function(x) {
  names(x) <- sub("Mmusculus-jaspar2022-", "", names(x))
  names(x) <- do.call(rbind, stringr::str_split(names(x), "-MA"))[,1]
  x
})

label <- sce_mult[, sce_mult$subclusters == ct]
pb_mult <- aggregateAcrossCells(label, use.assay.type = "counts",
                                id=DataFrame(label=label$subclusters, sample=label$dataset))
pb_mult <- logNormCounts(pb_mult)
colnames(pb_mult) <- paste(pb_mult$subclusters, pb_mult$dataset, sep = "_")

expressed_tf <- lapply(match_tf, function(x) {
  ex <- pb_mult[, pb_mult$subclusters == ct]
  ex <- ex[rownames(ex) %in% names(x), ]
  ex <- ex[rowSums(logcounts(ex))>10, ]
  ex <- rownames(ex)
})

pb_ll <- lapply(expressed_tf, function(x) pb_mult[rownames(pb_mult) %in% x, ])

sapply(pb_ll, function(x) length(rownames(x)))

# ## For C26 Myonuclei_Trim63 vs CTRL Myonuclei overall #####
# label <- sce_mult[, sce_mult$clusters == "Myonuclei"]
# label$subclusters[label$subclusters != "Myonuclei_Trim63"] <- "Myonuclei"
# label <- label[, !(label$condition != "CAC" & label$subclusters == "Myonuclei_Trim63")]
# pb_mult <- aggregateAcrossCells(label, use.assay.type = "counts", 
#                                 id=DataFrame(label=label$clusters, sample=label$dataset))

# colnames(pb_mult) <- paste(pb_mult$clusters, pb_mult$dataset, sep = "_")

expressed_tf <- lapply(match_tf, function(x) {
  ex <- pb_mult[, pb_mult$subclusters == ct]
  ex <- ex[rownames(ex) %in% names(x), ]
  ex <- ex[rowSums(logcounts(ex))>10, ]
  ex <- rownames(ex)
})

df_freq <- lapply(1:2, function(i) {
  tf_freq <- sort(sapply(match_tf[[i]], function(x) nrow(x)))
  tf_freq <- data.frame(Var1 = names(tf_freq), Freq = tf_freq)
  tf_freq <- tf_freq[order(tf_freq$Freq, decreasing = T), ]
  tf_freq <- tf_freq[tf_freq$Var1 %in% expressed_tf[[i]], ]
})

p <- ggplot(data = df_freq[[1]], aes(x = reorder(Var1, +Freq), y = Freq)) +
  geom_bar(stat = "identity")+ coord_flip() +theme_classic() +
  labs(y = "N. of targets", x = "TF", fill = NULL, color = NULL) +
  ggtitle(ct) + theme(text=element_text(size=13))
ggsave(p, file = paste0("cachectic_barplot_tf_ntarget_", tolower(ct),"_1kb.pdf"))

p <- ggplot(data = df_freq[[2]], aes(x = reorder(Var1, +Freq), y = Freq)) +
  geom_bar(stat = "identity")+ coord_flip() +theme_classic() +
  labs(y = "N. of targets", x = "TF", fill = NULL, color = NULL) +
  ggtitle(ct) + theme(text=element_text(size=13))
ggsave(p, file = paste0("cachectic_barplot_tf_ntarget_", tolower(ct),"_5kb.pdf"))

up_reg <- lapply(expressed_tf, function(x) {
  pb <- pb_mult[rownames(pb_mult) %in% x, ]
  y <- DGEList(counts(pb), samples = colData(pb))
  y <- calcNormFactors(y, method = "upperquartile")
  
  # design <- model.matrix(~0 + y$samples$subclusters, y$samples)
  # colnames(design) <- c("Myonuclei", "Myonuclei_Trim63")
  # contrast <- makeContrasts(Myonuclei_Trim63 - Myonuclei,levels = design)
  
  design <- model.matrix(~0 + y$samples$condition, y$samples)
  colnames(design) <- c("C26", "CTRL")
  contrast <- makeContrasts(C26 - CTRL,levels = design)
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  res <- glmLRT(fit, contrast = contrast)
  
  up_reg <- as.data.frame(topTags(res, n = length(rownames(pb))))
  up_reg <- up_reg[order(up_reg$logFC, decreasing = TRUE), ]
  
  up_reg <- up_reg[c(up_reg$logFC >= 0), ]
  up_reg <- rownames(up_reg)
})

sapply(up_reg, function(x) length(x))

df_freq <- lapply(1:2, function(i) df_freq[[i]][df_freq[[i]]$Var1 %in% up_reg[[i]], ])

p <- ggplot(data = df_freq[[1]], aes(x = reorder(Var1, +Freq), y = Freq)) +
  geom_bar(stat = "identity")+ coord_flip() +theme_classic() +
  labs(y = "N. of targets", x = "TF", fill = NULL, color = NULL) +
  ggtitle(ct) + theme(text=element_text(size=13))
ggsave(p, file = paste0("cachectic_barplot_uptf_ntarget_", tolower(ct),"_1kb.pdf"))

p <- ggplot(data = df_freq[[2]], aes(x = reorder(Var1, +Freq), y = Freq)) +
  geom_bar(stat = "identity")+ coord_flip() +theme_classic() +
  labs(y = "N. of targets", x = "TF", fill = NULL, color = NULL) +
  ggtitle(ct) + theme(text=element_text(size=13))
ggsave(p, file = paste0("cachectic_barplot_uptf_ntarget_", tolower(ct),"_5kb.pdf"))

match_expressed <- lapply(1:2, function(i) match_tf[[i]][expressed_tf[[i]]])

df_peaks <- lapply(ll, function(x) {
  df <- as.data.frame(x)
  rownames(df) <- paste(df$seqnames, paste(df$start, df$end, sep = "-"), sep = ":")
  df
})

genes_TF <- lapply(match_expressed[[1]], function(x) {
  inter <- intersect(x$seq, rownames(df_peaks[[1]]))
  genes <- df_peaks[[1]][inter, ]$SYMBOL
})

genes_TF <- lapply(genes_TF, function(x) data.frame(DEGs = x))
writexl::write_xlsx(genes_TF, path = paste0("cachectic_DEGs_TF_", tolower(ct), "_ctrloverall_1kb.xlsx"))

genes_TF <- lapply(match_expressed[[2]], function(x) {
  inter <- intersect(x$seq, rownames(df_peaks[[2]]))
  genes <- df_peaks[[2]][inter, ]$SYMBOL
})

genes_TF <- lapply(genes_TF, function(x) data.frame(DEGs = x))
writexl::write_xlsx(genes_TF, path = paste0("cachectic_DEGs_TF_", tolower(ct), "_ctrloverall_5kb.xlsx"))

