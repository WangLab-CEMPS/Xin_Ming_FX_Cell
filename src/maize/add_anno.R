library(Seurat)
library(dplyr)
library(ggplot2)

obj <- readRDS("./dataLib/maize/maize_pub.Rds")

obj1 <- readRDS("./dataLib/maize/maize.rds")
obj2 <- readRDS("./dataLib/maize/science_maize.rds")

meta1 <- obj1@meta.data
rownames(meta1) <- paste0("maize_", rownames(meta1))

meta2 <- obj2@meta.data
rownames(meta2) <- paste0("science_maize_", rownames(meta2))

meta1 <- meta1[, c("orig.ident", "celltype")]
meta2 <- meta2[, c("orig.ident", "celltype")]

table(rownames(meta1) %in% colnames(obj))
table(rownames(meta2) %in% colnames(obj))

meta <- rbind(meta1, meta2)
meta <- mutate(meta, cell = rownames(meta))
meta <- meta[, c("celltype", "cell")]
head(meta)

obj@meta.data <- obj@meta.data %>%
  mutate(cell = rownames(.)) %>%
  left_join(meta)
rownames(obj@meta.data) <- obj$cell

col <- pigRutils::select_colors("col_33")[9:22]


p_umap <- DimPlot(obj, group.by = "celltype", label = TRUE, pt.size = 0.7, cols = col) &
  tidydr::theme_dr() &
  theme(panel.grid = element_blank(), title = element_blank()) &
  guides(
    color = guide_legend(
      keywidth = 0.2,
      override.aes = list(size = 2, stroke = 2)
    )
  )

ggsave("./Plots/maize/maize_pub_anno_umap.pdf", p_umap, height = 6, width = 7)


p_umap_split <- DimPlot(obj, group.by = "celltype", label = TRUE, pt.size = 0.7, cols = col, split.by = "orig.ident") &
  tidydr::theme_dr() &
  theme(panel.grid = element_blank(), title = element_blank()) &
  guides(
    color = guide_legend(
      keywidth = 0.2,
      override.aes = list(size = 2, stroke = 2)
    )
  )

ggsave("./Plots/maize/maize_pub_anno_umap_split.pdf", p_umap_split, height = 6, width = 13)

saveRDS(obj, "./dataLib/maize/maize_pub.Rds")
