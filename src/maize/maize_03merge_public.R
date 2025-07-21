library(Seurat)
library(ggplot2)
library(magrittr)
library(reticulate)


fl <- c("./dataLib/maize/maize.rds", "./dataLib/maize/science_maize.rds")
n <- gsub("maize", "science")


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

obj

rm(obj_list)
gc()

obj$orig.ident[obj$orig.ident == "SeuratProject"] <- "science"
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
  verbose = FALSE,
  theta = 5
)

dims <- 1:40
mm <- "harmony"

obj %<>%
  FindNeighbors(obj, dims = dims, reduction = mm) %>%
  FindClusters(resolution = 0.4, cluster.name = "harmony_clusters") %>%
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
    title = element_blank()
  ) &
  guides(
    color = guide_legend(
      keywidth = 0.4,
      ncol = 2,
      override.aes = list(size = 2, stroke = 2)
    )
  )

p_split <- DimPlot(obj, split.by = "orig.ident", group.by = "orig.ident")


ggsave(filename = "./Plots/maize/maize_public_harmony_umap.pdf", p_umap, width = 13, height = 6)
ggsave(filename = "./Plots/maize/maize_public_harmony_split_umap.pdf", p_split, width = 12, height = 6)


saveRDS(obj, "./dataLib/maize/maize_pub.Rds")
