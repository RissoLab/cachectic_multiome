library(devtools)
library(CellChat)
library(scran)
library(Seurat)
library(dplyr)

# Load data
sgl_mult <- readRDS("multiome_sgl_comb.rds")
DefaultAssay(sgl_mult) <- "RNA"
sce <- as.SingleCellExperiment(sgl_mult)

# Genes filtering: keep the genes that have an expression average greater than zero in at least of 10% for each replicate.
keep <- apply(counts(sce), 1, function(y) all(tapply(y, sce$dataset, function(x) mean(x>0))>=0.1))
sce <- sce[keep, ]  

sce_ll <- lapply(unique(sce$condition), function(i) sce[, sce$condition == i])
sce_ll <- lapply(sce_ll, function(x) logNormCounts(x))
names(sce_ll) <- c("CTRL", "C26")

# Ligand-receptor interactions in CellChat database for mouse
CellChatDB <- CellChatDB.mouse 
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 

# Create a new CellChat object for each condition
cellchat <- lapply(sce_ll, function(x) 
  createCellChat(object = x, group.by = "subclusters"))

cellchat <- lapply(cellchat, function(x) {
  x <- setIdent(x, ident.use = "subclusters")
  x@DB <- CellChatDB.use
  x <- subsetData(x) # This step is necessary even if using the whole database
  x <- identifyOverExpressedGenes(x)
  x <- identifyOverExpressedInteractions(x)
  x <- computeCommunProb(x, type = "triMean")
  x <- filterCommunication(x, min.cells = 10)
})  

cellchat <- lapply(cellchat, function(x) {
  # Communication probability on signaling pathway level  
  x <- computeCommunProbPathway(x)
  #  calculates the aggregated cell-cell communication network by counting the number of links
  x <- aggregateNet(x, thresh = 1); x})

# Extract the inferred communications at the level of signaling pathways
ctrl_all <- cellchat$SHAM
df_ctrl <- subsetCommunication(ctrl_all, slot.name = "netP")
# Communication probability
prob <- df_ctrl %>% group_by(pathway_name) %>% summarise(sum = sum(prob))
prob <- prob[order(prob$sum, decreasing = T), ]

c26_all <- cellchat$CAC
df_c26 <- subsetCommunication(c26_all, slot.name = "netP")
prob <- df_c26 %>% group_by(pathway_name) %>% summarise(sum = sum(prob))
prob <- prob[order(prob$sum, decreasing = T), ]

save(ctrl_all, file = "cellchat_ctrl_all.RData")
save(c26_all, file = "cellchat_c26_all.RData")

load("SubclusterPaletteColor.RData")
palette_color <- color_palette$Color
names(palette_color) <- color_palette$Label

pathways.show <- cellchat[[1]]@netP$pathways
for (i in seq_along(pathways.show)) {
  pdf(file = paste0("multiome_chord_ctrl_", pathways.show[i], ".pdf"), width = 9, height = 9)
  p <- netVisual_aggregate(cellchat[[1]], signaling = pathways.show[i], 
                           layout = "chord",
                           color.use = palette_color,
                           signaling.name = paste0("CTRL - ", pathways.show[i]),
                           vertex.label.cex = 0.8, title.space = 10)
  print(p)
  dev.off()
}

pathways.show <- cellchat[[2]]@netP$pathways
for (i in seq_along(pathways.show)) {
  pdf(file = paste0("multiome_chord_c26_", pathways.show[i], ".pdf"), width = 9, height = 9)
  p <- netVisual_aggregate(cellchat[[2]], signaling = pathways.show[i], 
                           layout = "chord",
                           color.use = palette_color,
                           signaling.name = paste0("C26 - ", pathways.show[i]),
                           vertex.label.cex = 0.8, title.space = 10)
  print(p)
  dev.off()
}
save(cellchat, file = "multiome_cellchat_cond.RData")

ctrl <- cellchat[[1]]
c26 <- cellchat[[2]]

pdf("chord_ctrl_fgf.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(ctrl, signaling = c("FGF"), color.use = palette_color, lab.cex = 0.6, thresh = 0.05, slot.name = "netP",
                          legend.pos.y = 10, title.name = "FGF - CTRL (84)")
print(p)
dev.off()


pdf("chord_c26_fgf.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(c26, signaling = c("FGF"), color.use = palette_color, lab.cex = 0.6, thresh = 0.05, 
                          slot.name = "netP", legend.pos.y = 10, title.name = "FGF - C26 (232)")
print(p)
dev.off()

pdf("chord_ctrl_tgfb.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(ctrl, signaling = c("TGFb"), color.use = palette_color, lab.cex = 0.6, thresh = 0.05, slot.name = "netP",
                          legend.pos.y = 10, title.name = "TGFb - CTRL (64)")
print(p)
dev.off()


pdf("chord_ctrl_fgf.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(ctrl, signaling = c("FGF"), color.use = palette_color, lab.cex = 0.6, thresh = 0.05, slot.name = "netP", 
                          legend.pos.y = 10, title.name = "FGF - CTRL (84)")
print(p)
dev.off()

pdf("chord_ctrl_vegf.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(ctrl, signaling = c("VEGF"), color.use = palette_color, lab.cex = 0.6, thresh = 0.05, slot.name = "netP",
                          legend.pos.y = 10, title.name = "VEGF - CTRL (60)")
print(p)
dev.off()

pdf("chord_ctrl_sema.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(ctrl, signaling = c("SEMA3"), color.use = palette_color, lab.cex = 0.6, thresh = 0.05, slot.name = "netP",
                          legend.pos.y = 10, title.name = "SEMA3 - CTRL (102)")
print(p)
dev.off()

pdf("chord_ctrl_pdgf.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(ctrl, signaling = c("PDGF"), color.use = palette_color, lab.cex = 0.6, thresh = 0.05, slot.name = "netP",
                          legend.pos.y = 10, title.name = "PDGF - CTRL (8)")
print(p)
dev.off()

pdf("chord_c26_bmp.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(c26,
                          signaling = c("BMP"),
                          color.use = palette_color,
                          lab.cex = 0.6, thresh = 0.05, slot.name = "netP", legend.pos.y = 10,
                          title.name = "BMP - C26 (546)")
print(p)
dev.off()

pdf("chord_c26_tgfb.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(c26,
                          signaling = c("TGFb"),
                          color.use = palette_color,
                          lab.cex = 0.6, thresh = 0.05, slot.name = "netP",
                          legend.pos.y = 10,
                          title.name = "TGFb - C26 (64)")
print(p)
dev.off()

pdf("chord_c26_fgf.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(ctrl,
                          signaling = c("FGF"),
                          color.use = palette_color,
                          lab.cex = 0.6, thresh = 0.05, slot.name = "netP",
                          legend.pos.y = 10,
                          title.name = "FGF - C26 (232)")
print(p)
dev.off()

pdf("chord_c26_vegf.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(c26,
                          signaling = c("VEGF"), color.use = palette_color,
                          lab.cex = 0.6, thresh = 0.05, slot.name = "netP",
                          legend.pos.y = 10,
                          title.name = "VEGF - C26 (76)")
print(p)
dev.off()

pdf("chord_c26_sema.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(c26,
                          signaling = c("SEMA3"), color.use = palette_color,
                          lab.cex = 0.6, thresh = 0.05, slot.name = "netP",
                          legend.pos.y = 10,
                          title.name = "SEMA3 - C26 (124)")
print(p)
dev.off()

pdf("chord_c26_pdgf.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(c26,
                          signaling = c("PDGF"), color.use = palette_color,
                          lab.cex = 0.6, thresh = 0.05, slot.name = "netP",
                          legend.pos.y = 10,
                          title.name = "PDGF - C26 (25)")
print(p)
dev.off()

pdf("chord_c26_bmp_tgfb_fgf.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(c26,
                          signaling = c("BMP", "TGFb", "FGF"), color.use = palette_color,
                          lab.cex = 0.6, thresh = 0.05, slot.name = "netP",
                          title.name = "BMP, TGFb, FGF - C26 (842)")
print(p)
dev.off()

pdf("chord_c26_vegf_sema_pdgf.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(c26,
                          signaling = c("VEGF", "SEMA3", "PDGF"), color.use = palette_color,
                          lab.cex = 0.6, thresh = 0.05, slot.name = "netP",
                          title.name = "VEGF, SEMA3, PDGF - C26 (225)")
print(p)
dev.off()

# select ligands from TGFb family
pos <- which(grepl("Igf", df_net_c26$ligand) == "TRUE")  
length(pos)
pdf("chord_c26_igf.pdf", width = 20, height = 20, bg = "white", paper = "A4") 
p <- netVisual_chord_gene(c26, signaling = "IGF", color.use = palette_color, 
                          lab.cex = 0.6, title.name = "IGF - C26 (17)")
print(p)
dev.off()
