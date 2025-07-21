library(patchwork)
library(ggrastr)
library(ggplot2)
library(Seurat)
library(dplyr)
library(pigRutils)

# rhizome ---------
obj <- readRDS("./dataLib/tillering_node_rhizome/rhizome_S0.rds")
table(obj$seurat_clusters)

col <- c(pigRutils::select_colors("col28"), pigRutils::select_colors("col_6"))

embSCdim(
  obj, "umap",
  group_by = "seurat_clusters",
  colors = col[c(29:34, 10:22, 2:9)],
  add_label = TRUE,
  save = TRUE,
  width = 7, height = 6,
  save_path = "./figures/tillering_node_rhizome/rhizome_S0_umap_with_Cluster.pdf"
)

embSCdim(
  obj, "umap",
  group_by = "celltype",
  colors = col[c(29:33, 10:17)],
  save = TRUE,
  width = 11.5, height = 6,
  save_path = "./figures/tillering_node_rhizome/rhizome_S0_umap_with_anno.pdf"
)

# tillering -------
obj <- readRDS("./dataLib/tillering_node_rhizome/tillering.rds")
table(obj$seurat_clusters)

embSCdim(
  obj, "umap",
  group_by = "seurat_clusters",
  colors = col[c(29:34, 10:17)],
  add_label = TRUE,
  point_size = 0.2,
  save = TRUE,
  width = 7, height = 6,
  save_path = "./figures/tillering_node_rhizome/tillering_umap_with_Cluster.pdf"
)

embSCdim(
  obj, "umap",
  group_by = "celltype",
  colors = col[c(9:16, 30:33)],
  save = TRUE,
  width = 8, height = 6,
  save_path = "./figures/tillering_node_rhizome/tillering_umap_with_anno.pdf"
)
