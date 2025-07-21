library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(patchwork)

# gdm12
if (TRUE) {
  # 读取、准备画图数据 --------------
  obj <- readRDS("./dataLib/osRoot/gdm_dgm/gdm12.Rds")
  DimPlot(obj)
  colnames(obj@meta.data)
  Reductions(obj)

  pd <- SeuratObject::Embeddings(obj, reduction = "umap") %>%
    as_tibble(rownames = "cell") %>%
    mutate(c = obj$seurat_clusters, sample = obj$orig.ident)

  colnames(pd)[2:3] <- c("UMAP_1", "UMAP_2")

  # Cluster UMAP merge/split -----------------
  # merge
  p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = c)) +
    ggrastr::geom_point_rast(size = 0.1) +
    scale_color_manual(values = pigRutils::select_colors("col29")) +
    tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(
        keywidth = 0.4,
        ncol = 2,
        override.aes = list(size = 2, stroke = 2)
      )
    )

  ggsave("./figures/fig2s2/gdm_umap_with_cluster.pdf", width = 6.5, height = 5)

  p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
    ggrastr::geom_point_rast(size = 0.1) +
    scale_color_manual(values = c("#f9ca55", "#3d63ed")) +
    tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(
        keywidth = 0.4,
        ncol = 1,
        override.aes = list(size = 2, stroke = 2)
      )
    )

  ggsave("./figures/fig2s2/gdm_umap_with_sample.pdf", width = 6.5, height = 5)

  # split
  p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = c)) +
    ggrastr::geom_point_rast(size = 0.1) +
    scale_color_manual(values = pigRutils::select_colors("col29")) +
    facet_wrap(~sample) +
    tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(
        keywidth = 0.4,
        ncol = 2,
        override.aes = list(size = 2, stroke = 2)
      )
    )

  ggsave("./figures/fig2s2/gdm_umap_with_cluster_split.pdf", width = 10.2, height = 5)
}
# dgm12
if (TRUE) {
  # 读取、准备画图数据 --------------
  obj <- readRDS("./dataLib/osRoot/gdm_dgm/dgm12.Rds")
  DimPlot(obj)
  colnames(obj@meta.data)
  Reductions(obj)

  pd <- SeuratObject::Embeddings(obj, reduction = "umap") %>%
    as_tibble(rownames = "cell") %>%
    mutate(c = obj$seurat_clusters, sample = obj$orig.ident)

  colnames(pd)[2:3] <- c("UMAP_1", "UMAP_2")

  # Cluster UMAP merge/split -----------------
  # merge
  p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = c)) +
    ggrastr::geom_point_rast(size = 0.1) +
    scale_color_manual(values = pigRutils::select_colors("col29")) +
    tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(
        keywidth = 0.4,
        ncol = 2,
        override.aes = list(size = 2, stroke = 2)
      )
    )

  ggsave("./figures/fig2s2/dgm_umap_with_cluster.pdf", width = 6.5, height = 5)

  p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
    ggrastr::geom_point_rast(size = 0.1) +
    scale_color_manual(values = c("#9dd462", "#d34eb8")) +
    tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(
        keywidth = 0.4,
        ncol = 1,
        override.aes = list(size = 2, stroke = 2)
      )
    )

  ggsave("./figures/fig2s2/dgm_umap_with_sample.pdf", width = 6.5, height = 5)

  # split
  p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = c)) +
    ggrastr::geom_point_rast(size = 0.1) +
    scale_color_manual(values = pigRutils::select_colors("col29")) +
    facet_wrap(~sample) +
    tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(
        keywidth = 0.4,
        ncol = 2,
        override.aes = list(size = 2, stroke = 2)
      )
    )

  ggsave("./figures/fig2s2/dgm_umap_with_cluster_split.pdf", width = 10.2, height = 5)
}

# gdm12 dgm12 ggm12
if (TRUE) {
  # 读取数据、准备画图数据 -----------
  obj <- readRDS("./dataLib/osRoot/ggm_dgm_gdm.Rds")
  DimPlot(obj)
  colnames(obj@meta.data)
  Reductions(obj)

  pd <- SeuratObject::Embeddings(obj, reduction = "umap") %>%
    as_tibble(rownames = "cell") %>%
    mutate(c = obj$seurat_clusters, cc = obj$celltype, sample = obj$orig.ident, con = obj$tech)

  colnames(pd)[2:3] <- c("UMAP_1", "UMAP_2")

  # Cluster UMAP merge/split -----------------
  # merge
  p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = c)) +
    ggrastr::geom_point_rast(size = 0.1) +
    scale_color_manual(values = pigRutils::select_colors("col29")) +
    tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(
        keywidth = 0.4,
        ncol = 2,
        override.aes = list(size = 2, stroke = 2)
      )
    )

  ggsave("./figures/fig2s2/ggm_dgm_gdm_umap_with_cluster.pdf", width = 6.5, height = 5)

  p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = con)) +
    ggrastr::geom_point_rast(size = 0.1) +
    scale_color_manual(values = c("#eb9aeb", "#0d9eff", "#e8b964")) +
    tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(
        keywidth = 0.4,
        ncol = 1,
        override.aes = list(size = 2, stroke = 2)
      )
    )

  ggsave("./figures/fig2s2/ggm_dgm_gdm_umap_with_condition.pdf", width = 6.5, height = 5)

  # split
  p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = c)) +
    ggrastr::geom_point_rast(size = 0.1) +
    scale_color_manual(values = pigRutils::select_colors("col29")) +
    facet_wrap(~con) +
    tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(
        keywidth = 0.4,
        ncol = 2,
        override.aes = list(size = 2, stroke = 2)
      )
    )

  ggsave("./figures/fig2s2/ggm_dgm_gdm_umap_with_cluster_split.pdf", width = 15.2, height = 5)

  # Anno UMAP merge/split---------------
  # merge
  p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = cc)) +
    ggrastr::geom_point_rast(size = 0.08) +
    scale_color_manual(values = pigRutils::select_colors("col29")) +
    tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(
        keywidth = 0.4,
        ncol = 1,
        override.aes = list(size = 2, stroke = 2)
      )
    )

  ggsave("./figures/fig2s2/ggm_dgm_gdm_umap_with_anno.pdf", width = 7.2, height = 5)

  # split
  p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = cc)) +
    ggrastr::geom_point_rast(size = 0.08) +
    scale_color_manual(values = pigRutils::select_colors("col29")) +
    facet_wrap(~con) +
    tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(
        keywidth = 0.4,
        ncol = 1,
        override.aes = list(size = 2, stroke = 2)
      )
    )

  ggsave("./figures/fig2s2/ggm_dgm_gdm_umap_with_anno_split.pdf", width = 16, height = 5)

  # marker genes UMAP -------------------
  mm <- c(
    "Os02g0662000", "Os03g0279200", "Os12g0555600", "Os04g0423800",
    "Os03g0850900", "Os03g0216700", "Os10g0467800", "Os03g0720800",
    "Os07g0104100", "Os04g0554500", "Os03g0103200", "Os02g0581200"
  )

  p <- FeaturePlot(obj, mm)

  p_list <- lapply(p, function(x) {
    pd <- x$data %>%
      as_tibble(rownames = "cell") %>%
      left_join(pd)
    t <- x$labels$title

    ggplot(pd, aes(x = umap_1, y = umap_2, color = !!sym(t))) +
      geom_point_rast(size = 0.1) +
      facet_wrap(~con) +
      ggtitle(t)
  })


  pp <- wrap_plots(p_list, ncol = 2) &
    scale_color_gradientn(colours = c("#c8c6c3", "#f6ee6c", "#f9a432", "#eb212f", "#88181d")) &
    theme_linedraw() &
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank()
    )

  pdf("./figures/fig2s2/ggm_gdm_dgm_marker_umap.pdf", width = 16.4, height = 20)
  pp
  dev.off()
}
