library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(patchwork)

obj <- readRDS("./dataLib/osRoot/ggm_dgm_gdm.Rds")

mm <- c(
  "Os02g0662000", "Os03g0279200", "Os12g0555600", "Os04g0423800",
  "Os03g0850900", "Os03g0216700", "Os10g0467800", "Os03g0720800",
  "Os07g0104100", "Os04g0554500", "Os03g0103200", "Os02g0581200"
)

p <- FeaturePlot(obj, mm)

meta <- obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  select(cell, tech)

p_list <- lapply(p, function(x) {
  pd <- x$data %>%
    as_tibble(rownames = "cell") %>%
    left_join(meta)
  t <- x$labels$title

  ggplot(pd, aes(x = umap_1, y = umap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
    facet_wrap(~tech) +
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

pdf("./Plots/osroot/ggm_gdm_dmg_marker_umap.pdf", width = 16.3, height = 20)
pp
dev.off()
