library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(patchwork)

# ggm2 + TQ2 样本UMPA 注释UMAP Cluster UMAP -------------
# 准备画图数据
obj <- readRDS("./dataLib/osRoot/ggm/ggm2_TQ2_strata.Rds")
colnames(obj@meta.data)
Reductions(obj)

pd <- SeuratObject::Embeddings(obj, reduction = "harmony_umap") %>%
  as_tibble(rownames = "cell") %>%
  mutate(c = obj$seurat_clusters, ct = obj$celltype, sample = obj$orig.ident)

colnames(pd)[2:3] <- c("UMAP_1", "UMAP_2")

# 注释 UMAP -----------------
p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = cc)) +
  ggrastr::geom_point_rast(size = 0.1) +
  scale_color_manual(values = pigRutils::select_colors("col_33")) +
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

ggsave("./figures/fig1s1/ggm_tqStrata_umap_with_annotation.pdf", width = 7.2, height = 5)
# Cluster UMAP merge/split -----------
# merge
p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = c)) +
  ggrastr::geom_point_rast(size = 0.1) +
  scale_color_manual(values = pigRutils::select_colors("col_33")) +
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

ggsave("./figures/fig1s1/ggm_tqStrata_umap_with_Cluster.pdf", width = 7, height = 5.6)

# split
p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = c)) +
  ggrastr::geom_point_rast(size = 0.1) +
  scale_color_manual(values = pigRutils::select_colors("col_33")) +
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

ggsave("./figures/fig1s1/ggm_tqStrata_umap_with_Cluster_split.pdf", width = 10.1, height = 5)

# Sample UMAP merge --------------
p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
  ggrastr::geom_point_rast(size = 0.1) +
  scale_color_manual(values = c("#f56b62", "#02a5aa")) +
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

ggsave("./figures/fig1s1/ggm_tqStrata_umap_with_sample.pdf", width = 7, height = 5.3)

# marker umap plot -----------
mm <- c(
  "Os10g0122600", "Os03g0279200", "Os08g0512600", "Os04g0423800",
  "Os03g0850900", "Os03g0216700", "Os10g0467800", "Os04g0684300",
  "Os07g0104100", "Os04g0554500", "Os03g0103200", "Os02g0581200"
)

p <- FeaturePlot(obj, mm)

p_list <- lapply(p, function(x) {
  dd <- x$data %>%
    as_tibble(rownames = "cell") %>%
    left_join(pd)
  t <- x$labels$title

  ggplot(dd, aes(x = harmonyumap_1, y = harmonyumap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
    facet_wrap(~sample) +
    ggtitle(t)
})


pp <- wrap_plots(p_list, ncol = 3) &
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

pdf("./figures/fig1s1/ggm2_tqStrata_marker_umap_split.pdf", width = 18.6, height = 14.2)
pp
dev.off()

# clear--------
rm(list = ls())
gc()

# ggm1 ggm2----------
# 准备画图数据

obj <- readRDS("./dataLib/osRoot/ggm/ggm12_normal.Rds")
colnames(obj@meta.data)
Reductions(obj)

pd <- SeuratObject::Embeddings(obj, reduction = "umap") %>%
  as_tibble(rownames = "cell") %>%
  mutate(c = obj$seurat_clusters, sample = obj$orig.ident)

colnames(pd)[2:3] <- c("UMAP_1", "UMAP_2")

# Cluster UMAP merge/split -----------
# merge
p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = c)) +
  ggrastr::geom_point_rast(size = 0.1) +
  scale_color_manual(values = pigRutils::select_colors("col_33")) +
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

ggsave("./figures/fig1s1/ggm12_umap_with_Cluster.pdf", width = 7, height = 5.6)

# split
p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = c)) +
  ggrastr::geom_point_rast(size = 0.1) +
  scale_color_manual(values = pigRutils::select_colors("col_33")) +
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

ggsave("./figures/fig1s1/ggm12_umap_with_Cluster_split.pdf", width = 10.3, height = 5)

# marker umap plot -----------
mm <- c(
  "Os02g0662000", "Os03g0279200", "Os12g0555600", "Os04g0423800",
  "Os03g0850900", "Os03g0216700", "Os10g0467800", "Os03g0720800",
  "Os07g0104100", "Os04g0554500", "Os03g0103200", "Os02g0581200"
)

p <- FeaturePlot(obj, mm)

p_list <- lapply(p, function(x) {
  pd <- x$data
  t <- x$labels$title

  ggplot(pd, aes(x = umap_1, y = umap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
    ggtitle(t)
})


pp <- wrap_plots(p_list, ncol = 4) &
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

pdf("./figures/fig1s1/ggm12_marker_umap.pdf", width = 12.4, height = 9.9)
pp
dev.off()
