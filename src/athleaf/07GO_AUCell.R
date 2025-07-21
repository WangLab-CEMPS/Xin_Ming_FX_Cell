library(AUCell)
library(Seurat)
library(scales)
library(dplyr)
library(ggplot2)


go <- read.csv("./dataLib/athLeaf/GO_response_to_wounding.tsv", sep = "\t", header = FALSE)
go <- go[[1]]

obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public_anno.Rds")
obj <- JoinLayers(obj)
Layers(obj)
cts <- LayerData(obj, layer = "counts")

GeneSets <- list(
  wound = unique(go)
)

CellsAUC <- AUCell_run(cts, GeneSets)

pd <- Embeddings(obj, reduction = "harmony_umap") %>% as.data.frame()
pd$auc <- pigRutils::mf_minmax(getAUC(CellsAUC)["wound", ])
pd$celltype <- obj$celltype
pd$group <- obj$orig.ident

head(pd)

p <- ggplot(pd, aes(x = harmonyumap_1, y = harmonyumap_2, color = auc)) +
  ggrastr::geom_point_rast(size = 0.3) +
  facet_wrap(~group) +
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) +
  theme_linedraw() +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  ggtitle("AUCell Score: Response to wounding")

ggsave("./figures/leaf/goID_wound_AUCell_split.pdf", p, height = 6.2, width = 18)

table(pd$celltype)
pd$celltype <- gsub("Mesophyll\\*", "Mesophyll", pd$celltype)
table(pd$celltype)

p3 <- ggplot(pd, aes(x = celltype, y = auc)) +
  geom_violin(aes(fill = group)) +
  geom_boxplot(width = 0.1) +
  geom_hline(aes(yintercept = round(mean(auc), 1)), linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~group) +
  ggtitle("AUCell Score: Response to wounding") +
  ylab("AUCell Score") +
  xlab("CellType") +
  ggprism::theme_prism(base_size = 11) +
  scale_fill_manual(values = c("#f3827e", "#22bba5", "#59a1d7")) +
  coord_flip() +
  theme(
    legend.title = element_blank(), axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.text = element_text(face = "plain")
  )

pdf("./figures/leaf/goID_wound_AUCell_violin.pdf", width = 10, height = 5)
p3
dev.off()
