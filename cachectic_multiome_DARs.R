.libPaths(new = "/mnt/ganimede/zuin/R/x86_64-pc-linux-gnu-library/4.3")

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

bams <- list.files("/mnt/tebe/data/2023_Sandri/multiome_cachectic/libraries", pattern="*bam$",
                   recursive=TRUE, full.names=TRUE)
bams <- unique(dirname(bams)[-c(9,10)])

scelist <- readRDS("SCElist_multiome_annot.rds")

label <- levels(factor(scelist[[1]]$subclusters))
ct_scelist <- scelist

# FAPs
for (i in seq_along(ct_scelist)) {
  ct_scelist[[i]]$dirbam <- bams[i]
  ct_scelist[[i]] <- ct_scelist[[i]][, ct_scelist[[i]]$subclusters == label[15]]
}

lapply(ct_scelist, function(sce){
  sce$Barcode <- sce$gex_barcode
  sce$Sample <- sce$dataset
  lapply(unique(sce$subclusters), function(ct){
    createBamCt(sce, cellType = ct,
                cellTypesCol = "subclusters",
                bamdir = unique(sce$dirbam), bamType = "ATAC",
                outdir = paste0("/mnt/ganimede/zuin/SandriProj/bams/", tolower(label[15])), ncores = 10)
  })
})

# Peak calling (MACS2) ####
bamfiles <- list.files("/mnt/ganimede/zuin/SandriProj/bams/tenocyte", pattern="*_sorted.bam$", full.names=TRUE, recursive=TRUE)

outdir <- "/mnt/ganimede/zuin/SandriProj/MACS2_peaks"
additional.params <- "--nomodel "
extsize <- 200
shift  <- -(extsize/2)
for(file in bamfiles)
{
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

# Peaks alignment ####
# Align peaks to form common regions then filter regions for presence in multiple
# replicates taking in input a GRangesList where each element is a sample of called
# peaks.
bpeaks <- list.files("/mnt/ganimede/zuin/SandriProj/MACS2_peaks/Tenocyte", pattern=".broadPeak$", full.names=TRUE)

grl <- lapply(bpeaks, import, format="broadPeak")
names(grl) <- gsub("_sorted.bam_peaks.broadPeak", "", basename(bpeaks))

ccs <- c("CAC","SHAM")
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

gtftxdb <- makeTxDbFromGFF(file="/mnt/ganimede/zuin/SandriProj/refdata-gex-mm10-2020-A/genes/genes.gtf", format="gtf")

annotPeak <- annotatePeak(allpeaks, tssRegion=c(-5000, 5000), TxDb=gtftxdb, annoDb="org.Mm.eg.db")
annotPeak <- annotPeak@anno

save(region_ll, annotPeak, file = paste0("peaks_DARs_Tenocyte.RData"))


# Negative control genes
# df <- readRDS("multiome_UQ_DFResults.rds")
# df <- lapply(df, function(x) x[x$FDR>=0.05,])
# df <- lapply(df, function(x) rownames(x))
#
# neg_ctrl <- Reduce(intersect, df)
#
# # We converted the gene symbol to Ensembl ID
# ensembl_biomaRt <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
# neg_ctrl <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
#                   filters = "external_gene_name", values = neg_ctrl, mart = ensembl_biomaRt)
#
# # assign peaks to negative control genes
# allpeaks$negCtrl <- FALSE
# annotPeak$negCtrl <- FALSE
#
# int <- intersect(annotPeak$geneId, neg_ctrl$ensembl_gene_id)
# annotPeak$negCtrl[annotPeak$geneId %in% int] <- TRUE
#
# neg_peak <- names(annotPeak[annotPeak$negCtrl == TRUE, ])
# allpeaks$negCtrl[names(allpeaks) %in% neg_peak] <- TRUE

# DARs #####
bamdir <- paste0("/mnt/ganimede/zuin/SandriProj/bams/bams_sorted/FAPs")
count_final_regions <- countFinalRegions(allpeaks, readsFilePath = bamdir,
                                         fileType="bam", minCarriers=1)
colnames(count_final_regions) <- c("REP1CAC", "REP1SHAM","REP2CAC", "REP2SHAM")

counts <- SummarizedExperiment::assay(count_final_regions)
counts <- counts[, order(colnames(counts))]
head(counts)

# regions <- SummarizedExperiment::rowRanges(count_final_regions)
# negpeak <- rownames(count_final_regions[rowData(count_final_regions)$negCtrl==TRUE,])

colors <- brewer.pal(3, "Set2")
trt <- factor(c("CAC", "SHAM", "CAC", "SHAM"))
groups <- makeGroups(trt)

# s <- RUVSeq::RUVs(counts, cIdx = negpeak, scIdx = groups, k = 1)

EDASeq::plotPCA(counts, col = colors[trt], main = "PCA",
                labels = T, pch = 19)

design <- model.matrix(~0 + trt)
colnames(design) <- c(levels(trt))

y <- edgeR::DGEList(counts = counts, group = trt)
y <- edgeR::estimateDisp(y, design)

fit <- edgeR::glmQLFit(y, design, robust=TRUE)

con <- limma::makeContrasts(CAC - SHAM, levels = design)

qlf <- edgeR::glmQLFTest(fit, contrast=con)
res <- edgeR::topTags(qlf, n=Inf, p.value=0.05)
res <- as.data.frame(res)
res <- res[order(res$logFC, decreasing = T), ]

head(res)
dars <- res[res$FDR < 0.05, ]
dim(dars)

up_dars <- dars[dars$logFC >= 0.5, ]
down_dars <- dars[dars$logFC <= -0.5, ]

save(region_ll, annotPeak, annotPeak5kb, file = paste0("peaks_DARs_Myonuclei_Trim63.RData"))
