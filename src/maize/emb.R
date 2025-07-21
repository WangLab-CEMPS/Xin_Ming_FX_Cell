library(Seurat)
library(ggplot2)

obj <- readRDS("./dataLib/maize/science_maize.rds")

col <- pigRutils::select_colors("col_33")[9:22]

p_umap <- DimPlot(
  obj,
  group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.9,
  cols = col
) &
  tidydr::theme_dr() &
  theme(panel.grid = element_blank(), title = element_blank()) &
  guides(
    color = guide_legend(
      keywidth = 0.2,
      override.aes = list(size = 2, stroke = 2)
    )
  )

p_umap

ggsave("./Plots/maize/science_maize_anno_umap.pdf", height = 6, width = 7)
