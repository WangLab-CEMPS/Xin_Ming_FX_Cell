library(Seurat)
library(eulerr)
library(dplyr)

# read leaf scRNA-seq data
obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public_anno.Rds")
obj <- JoinLayers(obj)
table(obj$orig.ident)

# DEG between JC and LC, as the gene set affected by protoplasted effect
Idents(obj) <- "orig.ident"
deg <- FindMarkers(obj, ident.1 = "JC", ident.2 = "LC", logfc.threshold = 0.25)
head(deg)

sel <- deg$p_val_adj < 0.05 & deg$avg_log2FC > 1.58
sum(sel)
pseudo_protolasted_genes <- deg[sel, ]

write.csv(deg, "./dataLib/athLeaf/sc_protolasted_leaf.csv")

# read in real protoplasted genes
protolast <- read.csv("./dataLib/athLeaf/Wolf_B_Frommer_protolasted.csv")

intersect(rownames(pseudo_protolasted_genes), protolast$Gene_ID) %>% length()

# venn diagram
gl <- list(
  SC_protolasted = rownames(pseudo_protolasted_genes),
  Bulk_protoplasted = protolast$Gene_ID
)

fit <- euler(gl)

p <- plot(
  fit,
  edges = list(col = "lightgray", lex = 3),
  labels = list(font = 1.2, cex = 1.2),
  quantities = list(col = "black", font = 1.5),
  fill = c("#c9d9f4", "#e7ffd2")
)

pdf("./figures/venn_sc_protolasted.pdf", width = 7, height = 7)
p
dev.off()


gs1 <- setdiff(protolast$Gene_ID, rownames(pseudo_protolasted_genes))
gs2 <- intersect(protolast$Gene_ID, rownames(pseudo_protolasted_genes))
gs3 <- setdiff(rownames(pseudo_protolasted_genes), protolast$Gene_ID)

library(clusterProfiler)
library(org.At.tair.db)
library(ggplot2)

goall <- compareCluster(
  geneCluster = list(
    SC_unique = gs3,
    Bulk_unque = gs1,
    intersect = gs2
  ),
  fun = "enrichGO",
  OrgDb = org.At.tair.db,
  ont = "BP",
  keyType = "TAIR"
)

goall <- clusterProfiler::simplify(goall,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

pp <- dotplot(goall, showCategory = 10, label_format = 80) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("./figures/go_protolasted_sc_bulk_wound_leaf.pdf", width = 7, height = 7)
pp
dev.off()
