library(Seurat)
library(ggplot2)
library(magrittr)
library(reticulate)


fl <- list.files("./dataLib/osRoot/gdm_dgm", full.names = TRUE, pattern = "Rds$")
# fl <- list.files("./dataLib/osRoot/gdm_dgm", full.names = TRUE, pattern = "^gdm")
# fl <- list.files("./dataLib/osRoot/gdm_dgm", full.names = TRUE, pattern = "^dgm")

n <- gsub("12.Rds", "", basename(fl))

obj_list <- lapply(fl, function(x) {
  readRDS(x) %>%
    GetAssayData(layer = "counts") %>%
    CreateSeuratObject()
})

(obj <- merge(
  x = obj_list[[1]],
  y = obj_list[2:length(fl)],
  add.cell.ids = n
))

obj <- JoinLayers(obj)
obj

rm(obj_list)
gc()

# ----
obj %<>%
  NormalizeData() %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = round(0.1 * length(rownames(.)))
  ) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(features = VariableFeatures(.), npcs = 100, verbose = FALSE)

ElbowPlot(obj, reduction = "pca", ndims = 100)


dims <- 1:50
obj <- FindNeighbors(obj, dims = dims)

obj %<>%
  FindClusters(resolution = 0.4) %>%
  RunTSNE(dims = dims) %>%
  RunUMAP(dims = dims, metric = "correlation")


p_umap <- DimPlot(
  obj,
  group.by = "orig.ident",
  reduction = "umap", pt.size = .5
) +
  DimPlot(
    obj,
    reduction = "umap",
    label = TRUE,
    repel = TRUE, pt.size = .5,
    cols = pigRutils::select_colors("col29")
  ) &
  tidydr::theme_dr() &
  theme(
    panel.grid = element_blank(),
    title = element_blank()
  ) &
  guides(
    color = guide_legend(
      keywidth = 0.4,
      ncol = 2,
      override.aes = list(size = 2, stroke = 2)
    )
  )


p_umap2 <- DimPlot(
  obj,
  split.by = "orig.ident",
  reduction = "umap",
  label = TRUE,
  repel = TRUE, pt.size = .5,
  cols = pigRutils::select_colors("col29")
) &
  tidydr::theme_dr() &
  theme(
    panel.grid = element_blank(),
    title = element_blank()
  ) &
  guides(
    color = guide_legend(
      keywidth = 0.4,
      ncol = 2,
      override.aes = list(size = 2, stroke = 2)
    )
  )

# ggsave(filename = "./Plots/osroot/gdm_dgm/dgm12_umap.pdf", p_umap, width = 13.4, height = 6)
# ggsave(filename = "./Plots/osroot/gdm_dgm/dgm12_umap_split.pdf", p_umap2, width = 13, height = 6)
# 
# saveRDS(obj, "./dataLib/osRoot/gdm_dgm/dgm12.Rds")

# ggsave(filename = "./Plots/osroot/gdm_dgm/gdm12_umap.pdf", p_umap, width = 13.4, height = 6)
# ggsave(filename = "./Plots/osroot/gdm_dgm/gdm12_umap_split.pdf", p_umap2, width = 13, height = 6)
# 
# saveRDS(obj, "./dataLib/osRoot/gdm_dgm/gdm12.Rds")


ggsave(filename = "./Plots/osroot/gdm_dgm/dgm12_gdm12_umap.pdf", p_umap, width = 16, height = 8)
ggsave(filename = "./Plots/osroot/gdm_dgm/dgm12_gdm12_umap_split.pdf", p_umap2, width = 18, height = 6)

saveRDS(obj, "./dataLib/osRoot/gdm_dgm/dgm12_gdm12.Rds")
