library(ggplot2)
library(gtools)
library(tidyr)
library(dplyr)

obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public_anno.Rds")
table(obj$celltype)
obj$celltype <- gsub("Mesophyll\\*", "Mesophyll", obj$celltype)
meta <- obj@meta.data %>% as_tibble(rownames = "cell")
colnames(meta)

# bar -------------
props <- meta %>%
  group_by(celltype, orig.ident) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype) %>%
  mutate(proportion = n / sum(n) * 100)

props$orig.ident <- factor(props$orig.ident, levels = c("LC", "LPZ", "JC"))
props$cc <- as.numeric(props$orig.ident)

p <- ggplot(props, aes(x = celltype, y = proportion, fill = orig.ident)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("#1dc1ab", "#56acfc", "#fc877f")) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12)
  ) +
  scale_y_continuous(expand = c(0, 0))

pdf("./figures/leaf/leaf_LC_LPC_pub_bar.pdf", width = 5.6, height = 3.8)
p
dev.off()


props <- meta %>%
  group_by(seurat_clusters, orig.ident) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = n / sum(n) * 100)

props$orig.ident <- factor(props$orig.ident, levels = c("LC", "LPZ", "JC"))
props$cc <- as.numeric(props$orig.ident)

p <- ggplot(props, aes(x = seurat_clusters, y = proportion, fill = orig.ident)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("#1dc1ab", "#56acfc", "#fc877f")) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12)
  ) +
  scale_y_continuous(expand = c(0, 0))

pdf("./figures/leaf/leaf_LC_LPC_pub_bar_clusters.pdf", width = 3, height = 4)
p
dev.off()

# Roe-----------------
library(RColorBrewer)
library(pigRutils)
library(ComplexHeatmap)

# cluster
dd <- calcConDist(
  dat.tb = meta,
  colname.cluster = "seurat_clusters",
  colname.condition = "orig.ident",
  method = "chisq"
)

pd <- dd %>%
  as.data.frame() %>%
  pivot_wider(names_from = "Var2", values_from = "Freq") %>%
  as.data.frame()

pd$Var1 <- paste0("C", pd$Var1)
rownames(pd) <- pd$Var1
pd <- arrange(pd, desc(LC))
pd <- pd[, c(3, 4, 2)]

p <- Heatmap(
  as.matrix(pd),
  col = brewer.pal(9, "YlGnBu"),
  name = "Ro/e",
  cluster_columns = FALSE,
  rect_gp = gpar(col = "white", lwd = 1),
  layer_fun = function(j, i, x, y, width, height, fill) {
    v <- pindex(as.matrix(pd), i, j)
    l <- v > 1
    grid.text(sprintf("%.2f", v[l]), x[l], y[l], gp = gpar(fontsize = 12, col = "white"))
  }
)

pdf("./figures/leaf/leaf_LC_LPC_pub_Roe_cluster.pdf", height = 6.4, width = 3)
p
dev.off()


# celltype
ddc <- calcConDist(
  dat.tb = meta,
  colname.cluster = "celltype",
  colname.condition = "orig.ident",
  method = "chisq"
)

pdc <- ddc %>%
  as.data.frame() %>%
  pivot_wider(names_from = "Var2", values_from = "Freq") %>%
  as.data.frame()

rownames(pdc) <- pdc$Var1
pdc <- arrange(pdc, desc(LC))
pdc <- pdc[, c(3, 4, 2)]
# pdc <- pdc - min(pdc)

quantile(unlist(pdc), probs = c(0, 0.2, 0.4, 0.6, 0.8, 1))

p <- Heatmap(
  as.matrix(pdc),
  col = brewer.pal(9, "YlGnBu"),
  name = "Ro/e",
  cluster_columns = FALSE,
  rect_gp = gpar(col = "white", lwd = 1),
  layer_fun = function(j, i, x, y, width, height, fill) {
    v <- pindex(as.matrix(pd), i, j)
    l <- v > 1
    grid.text(sprintf("%.2f", v[l]), x[l], y[l], gp = gpar(fontsize = 12, col = "white"))
  }
  # layer_fun = function(j, i, x, y, width, height, fill) {
  #   v <- pindex(as.matrix(pdc), i, j)
  #   # Get symbols based on value ranges
  #   symbols <- ifelse(v > 1, "+++",
  #     ifelse(v > 0.8, "++",
  #       ifelse(v >= 0.2, "+",
  #         ifelse(v > 0, "+/-", "-")
  #       )
  #     )
  #   )
  #   # Add symbols to all cells
  #   grid.text(symbols, x, y, gp = gpar(fontsize = 12, col = "#7c7c7c"))
  # }
)

pdf("./figures/leaf/leaf_LC_LPC_pub_Roe_celltype.pdf", height = 5.4, width = 5)
p
dev.off()
