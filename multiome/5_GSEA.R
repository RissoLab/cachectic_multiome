# =============================================================================
# Cachectic multiome (C26 and CTRL samples)
# Geneset enrichment analysis (GSEA)
# =============================================================================

library(openxlsx)
library(clusterProfiler)
library(ggplot2)
library(org.Mm.eg.db)

df_mult <- readRDS("cachectic_mult_degs_ctrl_c26.rds")
df_mult <- df_mult[c(2, 5:10)]

df_mult <- lapply(df_mult, function(x) x[x$FDR < 0.05, ])

# we want the log2 fold change 
kegg_gene_list <- lapply(df_mult, function(x) {
  gene <- x$logFC
  names(gene) <- rownames(x)
  ids <- bitr(names(gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]), ]
  
  df2 <- x[rownames(x) %in% dedup_ids$SYMBOL, ]
  
  df2$Y <- dedup_ids$ENTREZID

  kegg_gene <- df2$logFC
  names(kegg_gene) <- df2$Y
  kegg_gene < -na.omit(kegg_gene)

  kegg_gene <- sort(kegg_gene, decreasing = TRUE)
})


kk1 <- lapply(kegg_gene_list, function(x) {
  out <- gseKEGG(geneList = x,
               organism = "mmu",
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType = "kegg")
  out@result$Description <- sub(' - Mu.*', '', out@result$Description)
  out
  })

options(enrichplot.colours = c("red","blue"))
lapply(1:9, function(i) {
  p <- dotplot(kk1[[i]], showCategory = 30, title = "Enriched Pathways", split=".sign") + facet_grid(.~.sign) 
  ggsave(p, filename = paste0("multiome_bubble_plot_", names(kk2)[i], "KEGG.pdf"), height = 30, width = 30, units = "cm")
         })

# Myonuclei overall
myo_overall <- readRDS("cachectic_mult_degs_myo_ctrl_c26.rds")
myo_overall <- myo_overall[myo_overall$FDR < 0.05, ]

gene <- myo_overall$logFC
names(gene) <- rownames(myo_overall)
ids <- bitr(names(gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]), ]
  
df2 <- myo_overall[rownames(myo_overall) %in% dedup_ids$SYMBOL, ]
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

options(enrichplot.colours = c("red","blue"))

p <- dotplot(kk2, showCategory = 30, title = "Enriched Pathways", split=".sign") + facet_grid(.~.sign) 
ggsave(p, filename = paste0("multiome_bubble_plot_myo_KEGG.pdf"), height = 30, width = 30, units = "cm")

labels <- c(names(kk1), "Myonuclei_Overall")
kk <- c(kk1, kk2)
names(kk) <- labels
writexl::write_xlsx(x = kk, path = "cachectic_mult_gsea_ctrl_c26.xlsx")