library(Seurat)
library(ggplot2)
library(magrittr)
library(reticulate)


fl <- c("./dataLib/osRoot/ggm2.rds", "./dataLib/osRoot/root_tq_strata.rds")

obj_list <- lapply(fl, function(x) {
  readRDS(x) %>%
    GetAssayData(layer = "counts") %>%
    CreateSeuratObject()
})


(obj <- merge(
  x = obj_list[[1]],
  y = obj_list[[2]]
))

# obj <- JoinLayers(obj)

rm(obj_list)
gc()

# obj <- subset(obj, subset = orig.ident != "osRoot1")
# obj$orig.ident[obj$orig.ident == "ggm1"] <- "ggm"
# obj$orig.ident[obj$orig.ident == "ggm2"] <- "ggm"

table(obj$orig.ident)

# ----
obj %<>%
  NormalizeData() %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = round(0.1 * length(rownames(.)))
  ) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(features = VariableFeatures(.), npcs = 100, verbose = FALSE)


obj <- IntegrateLayers(
  object = obj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = TRUE,
  max.iter.cluster = 40,
  theta = 4
)

dims <- 1:40
mm <- "harmony"

obj %<>%
  FindNeighbors(obj, dims = dims, reduction = mm) %>%
  FindClusters(resolution = 0.3, cluster.name = "harmony_clusters") %>%
  RunUMAP(dims = dims, metric = "correlation", reduction.name = "harmony_umap", reduction = mm)


# -----------
p_umap <- DimPlot(
  obj,
  group.by = "orig.ident",
  reduction = "harmony_umap", pt.size = .3
) +
  DimPlot(
    obj,
    reduction = "harmony_umap",
    label = TRUE,
    repel = TRUE, pt.size = .3,
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

p_split <- DimPlot(obj, split.by = "orig.ident", cols = pigRutils::select_colors("col_33"), pt.size = 0.3)


ggsave(filename = "./Plots/osroot/ggm2_TQ2_harmony_umap_strata.pdf", p_umap, width = 12.4, height = 7)
ggsave(filename = "./Plots/osroot/ggm2_TQ2_harmony_split_umap_strata.pdf", p_split, width = 12.7, height = 6)


saveRDS(obj, "./dataLib/osRoot/ggm2_TQ2_strata.Rds")
