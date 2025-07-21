library(patchwork)
library(ggplot2)
library(ggrastr)
library(Seurat)
library(dplyr)

obj <- readRDS("./dataLib/tillering_node_rhizome/rhizome_S0.rds")
DimPlot(obj, label = TRUE)

QC <- c("Bochr01G042780", "Bochr02G072640") # "Bochr06G206520"
PLT2 <- c("Bochr06G213020")
PLT1 <- c("Bochr04G157850")

QC <- c(QC, PLT2, PLT1)

p <- FeaturePlot(obj, QC, order = TRUE)

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

pdf("./figures/tillering_node_rhizome/rhizome_qc_marker_umap.pdf", width = 14.7, height = 4)
pp
dev.off()
