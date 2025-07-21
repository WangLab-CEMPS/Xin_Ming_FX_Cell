library(patchwork)
library(ggrastr)
library(ggplot2)
library(Seurat)
library(dplyr)
library(pigRutils)

# only maize -----------------
if (TRUE) {
  # read in data
  obj <- readRDS("./dataLib/maize/maize.rds")
  DimPlot(obj)
  colnames(obj@meta.data)

  # Anno UMAP -----------------
  embSCdim(
    obj, "umap",
    group_by = "celltype",
    colors = pigRutils::select_colors("col_33")[9:22],
    save = TRUE,
    width = 7, height = 5.4,
    save_path = "./figures/maize/maize_umap_with_anno.pdf"
  )

  # Cluster UMAP --------------
  embSCdim(
    obj, "umap",
    group_by = "seurat_clusters",
    colors = pigRutils::select_colors("col29"),
    save = TRUE,
    save_path = "./figures/maize/maize_umap_with_Cluster.pdf",
    width = 7,
    height = 6
  )

  embSCdim(
    obj, "umap",
    group_by = "seurat_clusters",
    colors = pigRutils::select_colors("col29"),
    save = TRUE,
    add_label = TRUE,
    save_path = "./figures/maize/maize_umap_with_Cluster_label.pdf",
    width = 7,
    height = 6
  )
}

# only science ---------------
if (TRUE) {
  # read in data
  obj <- readRDS("./dataLib/maize/science_maize.rds")
  DimPlot(obj)
  colnames(obj@meta.data)

  # Anno UMAP -----------------
  embSCdim(
    obj, "umap",
    group_by = "celltype",
    colors = pigRutils::select_colors("col_33")[9:22],
    save = TRUE,
    save_path = "./figures/maize/science_maize_umap_with_anno.pdf", width = 7, height = 5.4
  )

  # Cluster UMAP --------------
  embSCdim(
    obj, "umap",
    group_by = "seurat_clusters",
    colors = pigRutils::select_colors("col29"),
    save = TRUE,
    save_path = "./figures/maize/science_maize_umap_with_Cluster.pdf", width = 7, height = 6
  )

  embSCdim(
    obj, "umap",
    group_by = "seurat_clusters",
    colors = pigRutils::select_colors("col29"),
    add_label = TRUE,
    save = TRUE,
    save_path = "./figures/maize/science_maize_umap_with_Cluster_label.pdf", width = 7, height = 6
  )
}

# maize and science public -----------------
if (TRUE) {
  # read in data
  obj <- readRDS("./dataLib/maize/maize_pub.Rds")
  DimPlot(obj)
  colnames(obj@meta.data)
  obj@meta.data <- obj@meta.data[, -6]
  Reductions(obj)
  # samples UMAP --------
  embSCdim(
    obj, "harmony_umap",
    group_by = "orig.ident",
    colors = c("#f97148", "#00b6e4"),
    save = TRUE,
    save_path = "./figures/maize/maize_pub_umap_with_samples.pdf",
    width = 7, height = 5.4
  )
  # Anno UMAP merge/split -----------------
  # merge
  embSCdim(
    obj, "harmony_umap",
    group_by = "celltype",
    colors = pigRutils::select_colors("col_33")[9:22],
    save = TRUE,
    save_path = "./figures/maize/maize_pub_umap_with_anno.pdf", width = 7, height = 5.4
  )

  # split
  embSCdim(
    obj, "harmony_umap",
    group_by = "celltype",
    split_by = "orig.ident",
    colors = pigRutils::select_colors("col_33")[9:22],
    save = TRUE,
    save_path = "./figures/maize/maize_pub_umap_with_anno_split.pdf", width = 12, height = 5.5
  )

  # Cluster UMAP merge/split ------------------
  # merge
  embSCdim(
    obj, "harmony_umap",
    group_by = "seurat_clusters",
    colors = pigRutils::select_colors("col29"),
    save = TRUE,
    save_path = "./figures/maize/maize_pub_umap_with_Cluster.pdf", width = 7, height = 6
  )
  embSCdim(
    obj, "harmony_umap",
    group_by = "seurat_clusters",
    colors = pigRutils::select_colors("col29"),
    add_label = TRUE,
    save = TRUE,
    save_path = "./figures/maize/maize_pub_umap_with_Cluster_label.pdf", width = 7, height = 6
  )
  # split
  embSCdim(
    obj, "harmony_umap",
    group_by = "seurat_clusters",
    split_by = "orig.ident",
    colors = pigRutils::select_colors("col29"),
    save = TRUE,
    save_path = "./figures/maize/maize_pub_umap_with_Cluster_split.pdf", width = 12, height = 6
  )
}
