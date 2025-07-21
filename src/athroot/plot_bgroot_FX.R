library(pigRutils)
library(ggplot2)
library(ggrastr)
library(dplyr)

# read data
pd <- read.csv("./dataLib/athRoot/bigroot_FX_umap.csv", row.names = 1)
head(pd)

pd$cell_cluster <- ifelse(is.na(pd$cell_cluster), "FX-Cell", pd$cell_cluster)

cc <- c(
  "#f670cb", "#b1cfda", "#ead054", "#da9f91", "#d3ef9b",
  "#e7555b", "#e97617", "#3ea6bb", "#b489ac", "#adba26",
  "#299e92", "#f3dc5b", "#f09679", "#4dc684", "grey60"
)

names(cc) <- unique(pd$cell_cluster)

p <- ggplot() +
  # 先绘制非FX-Cell的细胞（底层）
  geom_point_rast(
    data = pd[pd$cell_cluster != "FX-Cell", ],
    aes(x = UMAP_1, y = UMAP_2, color = cell_cluster),
    size = 0.1
  ) +
  # 再绘制FX-Cell（上层）
  geom_point_rast(
    data = pd[pd$cell_cluster == "FX-Cell", ],
    aes(x = UMAP_1, y = UMAP_2, color = "FX-Cell"),
    size = 0.1, alpha = 0.7
  ) +
  # 设置颜色，非FX-Cell使用cc中的颜色，FX-Cell使用灰色
  scale_color_manual(values = cc) +
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

ggsave("./figures/athroot/bigroot_FX_umap.pdf", p, width = 15, height = 12)


p <- ggplot() +
  # 先绘制非FX-Cell的细胞（底层）
  geom_point_rast(
    data = pd[pd$cell_cluster != "FX-Cell", ],
    aes(x = UMAP_1, y = UMAP_2, color = "grey75"),
    size = 0.1
  ) +
  # 再绘制FX-Cell（上层）
  geom_point_rast(
    data = pd[pd$cell_cluster == "FX-Cell", ],
    aes(x = UMAP_1, y = UMAP_2, color = "#f67c7a"),
    size = 0.2, alpha = 0.8
  ) +
  scale_color_manual(
    values = c("grey75" = "grey75", "#f67c7a" = "#f67c7a"),
    labels = c("grey75" = "Other cells", "#f67c7a" = "FX-Cell")
  ) +
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

ggsave("./figures/athroot/bigroot_FX_umap_1.pdf", p, width = 15, height = 12)



head(pd)
cc <- c(
  "#f670cb", "#b1cfda", "#ead054", "#da9f91", "#d3ef9b",
  "#e7555b", "#e97617", "#3ea6bb", "#b489ac", "#adba26",
  "#299e92", "#f3dc5b", "#f09679", "#4dc684", "grey70"
)

names(cc) <- unique(pd$cell_cluster)


p2 <- ggplot(
  data = pd,
  aes(x = UMAP_1, y = UMAP_2, color = cell_cluster)
) +
  geom_point_rast(
    size = 0.1
  ) +
  scale_color_manual(values = cc) +
  facet_wrap(~t) +
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

ggsave("./figures/athroot/bigroot_FX_umap_split.pdf", p2, width = 20, height = 8)
