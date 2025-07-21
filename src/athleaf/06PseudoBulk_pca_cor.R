library(dplyr)
library(Seurat)
library(ggplot2)
library(pigRutils)

obj <- readRDS("./dataLib/athLeaf/LC_LPC_public.Rds")

obj[["tmp"]] <- paste0(obj$orig.ident, "_", obj$harmony_clusters)

cts <- AggregateExpression(obj, group.by = "tmp")
cts <- cts$RNA %>% as.matrix()
head(cts)

coldat <- data.frame(
  groups = colnames(cts) %>% stringr::str_split("-", simplify = TRUE) %>% .[,1],
  samples = colnames(cts)
)

pd <- doPCA(cts, coldata = coldat, ntop = 4000)
head(pd)

percentVar <- round(100 * attr(pd, "percentVar"))


pd2 <- pd %>%
  group_by(group) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))


p <- ggplot(pd, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = group), stat = "unique", size = 4, alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = pd2,
    aes(label = group),
    size = 3.5,
    bg.r = 0.25,
    color = "#000103",
    bg.color = "#f5f0f0",
    seed = 1,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(1, "lines")
  ) +
  labs(title = "PCA analysis") +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  coord_fixed()
p


library(ComplexHeatmap)

cc <- cor(cts, method = "spearman")
cc <- cc[1:18, 19:53]

Heatmap(cc, cluster_columns = FALSE)
