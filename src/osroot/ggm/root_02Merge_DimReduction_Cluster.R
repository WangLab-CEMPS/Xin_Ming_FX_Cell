library(Seurat)
library(ggplot2)
library(magrittr)
library(reticulate)


fl <- list.files("./dataLib/osRoot", full.names = TRUE, pattern = "rds$")
n <- gsub(".rds", "", basename(fl))

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
  FindClusters(resolution = 0.3) %>%
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


p_umap2 <- DimPlot(
  obj,
  split.by = "orig.ident",
  reduction = "umap",
  label = TRUE,
  repel = TRUE, pt.size = .5,
  cols = pigRutils::select_colors("col_33")
)  &
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


ggsave(filename = "./Plots/osroot/ggm12_normal_umap.pdf", p_umap, width = 16, height = 8)
ggsave(filename = "./Plots/osroot/ggm12_normal_umap_split.pdf", p_umap2, width = 16, height = 8)

saveRDS(obj, "./dataLib/osRoot/ggm12_normal.Rds")

