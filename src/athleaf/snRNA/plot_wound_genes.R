library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(pigRutils)

obj <- readRDS("./dataLib/athLeaf/scsnRNA_harmony.Rds")


genes <- c(
  "AT5G64750", "AT4G34410", "AT3G59220", "AT3G28210", "AT3G25780",
  "AT3G23250", "AT3G02840", "AT2G38530", "AT2G29500", "AT2G02990",
  "AT1G53540", "AT1G07400", "AT1G02920", "AT1G02400"
)

p <- FeaturePlot(obj, features = genes, reduction = "harmony_umap")

meta <- obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  select(cell, orig.ident)

sapply(seq.int(p), function(x) {
  xx <- p[[x]]
  pd <- xx$data %>%
    as_tibble(rownames = "cell") %>%
    left_join(meta)
  t <- xx$labels$title
  print(t)
  p <- ggplot(pd, aes(x = harmonyumap_1, y = harmonyumap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
    facet_wrap(~orig.ident, nrow = 1) +
    scale_color_gradientn(colours = c("#c8c6c3", "#f6ee6c", "#f9a432", "#eb212f", "#88181d")) +
    ggtitle(t) +
    theme_linedraw() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank()
    )

  ggsave(paste0("./Plots/athleaf/scsn_wound/", t, ".pdf"), p, height = 5, width = 22.2)
})
