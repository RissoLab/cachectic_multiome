# ==============================================================
# Cachectic visium
# Geneset enrichment analysis
# ==============================================================

library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(org.Mm.eg.db)

degs <- read.table(file = 'c26_treated_GFPvsBlack.tsv', sep = '\t', header = TRUE)
degs <- degs[degs$FDR < 0.05, ]

gene <- degs$logFC
names(gene) <- rownames(degs)
ids <- bitr(names(gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]), ]

df2 <- degs[rownames(degs) %in% dedup_ids$SYMBOL, ]
df2$Y <- dedup_ids$ENTREZID

kegg_gene <- df2$logFC
names(kegg_gene) <- df2$Y
kegg_gene < -na.omit(kegg_gene)
kegg_gene <- sort(kegg_gene, decreasing = TRUE)

kk2 <- gseKEGG(geneList = kegg_gene,
               organism = "mmu",
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType = "kegg")
kk2@result$Description <- sub(' - Mu.*', '', kk2@result$Description)

kk2 <- as.data.frame(kk2)

writexl::write_xlsx(kk2, path = "cachectic_visium_gsea.xlsx")
