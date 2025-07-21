library(patchwork)
library(ggplot2)
library(ggrastr)
library(Seurat)
library(dplyr)

obj <- readRDS("./dataLib/tillering_node_rhizome/tillering.rds")
DimPlot(obj, label = TRUE)

rc <- c("Os04g0552000", "Os01g0248000")

p <- FeaturePlot(obj, rc, order = TRUE)

p_list <- lapply(p, function(x) {
  pd <- x$data
  t <- x$labels$title

  ggplot(pd, aes(x = umap_1, y = umap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
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

pdf("./figures/tillering_node_rhizome/tillering_rc_marker_umap.pdf", width = 6.6, height = 3.6)
pp
dev.off()


p_list_sub <- lapply(p, function(x) {
  pd <- x$data
  pd <- pd[pd$umap_1 < -2.5 & pd$umap_2 < -0.5 & pd$umap_2 > -5, ]
  t <- x$labels$title

  ggplot(pd, aes(x = umap_1, y = umap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.5) +
    ggtitle(t)
})

pp_sub <- wrap_plots(p_list_sub, ncol = 2) &
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

pdf("./figures/tillering_node_rhizome/tillering_rc_marker_umap_sub.pdf", width = 6.6, height = 3.6)
pp_sub
dev.off()
