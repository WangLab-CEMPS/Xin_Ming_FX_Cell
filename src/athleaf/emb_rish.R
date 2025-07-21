library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(patchwork)


mm <- c(
  "AT3G17609",
  "AT1G56060",
  "AT2G02990",
  "AT3G04300",
  "AT2G35290",
  "AT5G54300",
  "AT1G22810",
  "AT4G34410",
  "AT1G65390"
)

obj <- readRDS("./dataLib/athLeaf/LC_LPC_public.Rds")

p <- FeaturePlot(obj, features = mm, reduction = "harmony_umap")

meta <- obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  select(cell, orig.ident)

p_list <- lapply(seq.int(p), function(x) {
  xx <- p[[x]]
  pd <- xx$data %>%
    as_tibble(rownames = "cell") %>%
    left_join(meta)
  t <- xx$labels$title
  print(t)
  ggplot(pd, aes(x = harmonyumap_1, y = harmonyumap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
    facet_wrap(~orig.ident) +
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

pdf("./Plots/leaf_rish.pdf", width = 17.2, height = 16)
pp
dev.off()
