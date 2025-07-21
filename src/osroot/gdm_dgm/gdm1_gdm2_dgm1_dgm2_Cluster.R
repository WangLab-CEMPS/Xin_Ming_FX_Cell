library(Seurat)
library(magrittr)
library(ggplot2)

nn <- c("gdm1", "gdm2", "dgm1", "dgm2")

for (i in nn) {
  on <- paste0("./dataLib/osRoot/gdm_dgm/", i, ".rds")
  message(on)
  obj <- readRDS(on)

  obj %<>%
    NormalizeData() %>%
    FindVariableFeatures(
      selection.method = "vst",
      nfeatures = round(0.1 * length(rownames(.)))
    ) %>%
    ScaleData(features = rownames(.), var.to.regress = "nCount_RNA") %>%
    RunPCA(features = VariableFeatures(.), npcs = 100, verbose = FALSE)

  ElbowPlot(obj, reduction = "pca", ndims = 100)


  dims <- 1:50
  obj <- FindNeighbors(obj, dims = dims)

  obj %<>%
    FindClusters(resolution = 0.3) %>%
    RunTSNE(dims = dims) %>%
    RunUMAP(dims = dims, metric = "correlation")


  pp <- DimPlot(
    obj,
    label = TRUE,
    reduction = "tsne",
    pt.size = .3,
    repel = TRUE,
  ) +
    DimPlot(
      obj,
      reduction = "umap",
      label = TRUE,
      repel = TRUE,
      pt.size = .3,
    ) &
    tidydr::theme_dr() &
    scale_color_manual(values = pigRutils::select_colors("col_33")) &
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
  pp

  ggsave(filename = paste0("./Plots/", i, "_tsne_umap.pdf"), pp, width = 11.5, height = 6)

  saveRDS(obj, on)
}
