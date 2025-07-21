library(Seurat)
library(dplyr)


obj <- readRDS("./dataLib/athLeaf/LC_LPC_normal.Rds")

drop_gene <- read.table("./dataLib/athLeaf/drop_gene.txt", header = TRUE)[[1]]
table(drop_gene %in% VariableFeatures(obj))

obj <- RunPCA(obj, features = VariableFeatures(obj)[!VariableFeatures(obj) %in% drop_gene], npcs = 100)

ElbowPlot(object = obj, ndims = 100)

dimension <- 1:30

obj %<>% FindNeighbors(dims = dimension) %>%
  FindClusters(resolution = 0.4) %>%
  RunTSNE(dims = dimension) %>% 
  RunUMAP(dims = dimension, metric = "correlation")

library(ggplot2)

p_umap <- DimPlot(
  obj,
  group.by = "orig.ident",
  reduction = "umap", pt.size = .6
) +
  DimPlot(
    obj,
    reduction = "umap",
    label = TRUE,
    repel = TRUE, pt.size = .6,
    cols = pigRutils::select_colors("col_33")
  ) &
  tidydr::theme_dr() &
  theme(
    panel.grid = element_blank(),
    legend.position = "top", title = element_blank()
  ) &
  guides(
    color = guide_legend(
      keywidth = 0.4,
      nrow = 2,
      override.aes = list(size = 2, stroke = 2)
    )
  )

p_umap


ggsave(filename = "./Plots/LC_LPZ_drop_gene_merge.pdf", p_umap, width = 16, height = 8)

saveRDS(obj, "./dataLib/athLeaf/LC_LPC_drop_gene.Rds")
