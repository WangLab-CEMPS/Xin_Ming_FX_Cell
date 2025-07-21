library(Seurat)
library(ggplot2)
library(magrittr)
library(reticulate)


fl <- c(
  "./dataLib/osRoot/ggm/ggm12_normal.Rds",
  "./dataLib/osRoot/gdm_dgm/gdm12.Rds",
  "./dataLib/osRoot/gdm_dgm/dgm12.Rds"
)

n <- c("ggm", "gdm", "dgm")

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

table(obj$orig.ident)
obj[["tech"]] <- gsub("[1-2]$", "", obj$orig.ident)
table(obj$tech)

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

obj$tech <- factor(obj$tech, levels = c("ggm", "gdm","dgm"))

p_umap2 <- DimPlot(
  obj,
  split.by = "tech",
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

ggsave(filename = "./Plots/osroot/ggm_dgm_gdm_umap.pdf", p_umap, width = 13.5, height = 6)
ggsave(filename = "./Plots/osroot/ggm_dgm_gdm_umap_split.pdf", p_umap2, width = 18.8, height = 6)

saveRDS(obj, "./dataLib/osRoot/ggm_dgm_gdm.Rds")
