library(clusterProfiler)
library(org.At.tair.db)
library(Seurat)
library(dplyr)

alias <- read.csv("~/ref/Athaliana/Araport11/gene_aliases_20220331_tidy.csv")

obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public.Rds")
obj <- JoinLayers(obj)

obj <- subset(obj, subset = orig.ident %in% c("LPZ"))

deg <- FindMarkers(
  object = obj,
  group.by = "celltype",
  ident.1 = "Mesophyll",
  ident.2 = "Mesophyll*"
)

deg %>%
  as_tibble(rownames = "geneId") %>%
  write.csv("./dataLib/athLeaf/LPZ_Mesophyll_de.csv", row.names = FALSE)

dd <- subset(deg, p_val < 0.05 & p_val_adj < 0.05 & abs(avg_log2FC) > 2) %>%
  as_tibble(rownames = "geneId") %>%
  left_join(alias, by = c("geneId" = "name"))


up <- dd$geneId[dd$avg_log2FC > 0]
down <- dd$geneId[dd$avg_log2FC < 0]

orgdb <- org.At.tair.db

go_up <- enrichGO(
  gene = up,
  OrgDb = orgdb,
  keyType = "TAIR",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
) %>% clusterProfiler::simplify(
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

write.csv(
  as.data.frame(go_up),
  "./dataLib/athLeaf/LPZ_Mesophyll_go_up.csv",
  row.names = FALSE
)

barplot(go_up, showCategory = 20)

go_down <- enrichGO(
  gene = down,
  OrgDb = orgdb,
  keyType = "TAIR",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
) %>% clusterProfiler::simplify(
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

write.csv(
  as.data.frame(go_down),
  "./dataLib/athLeaf/LPZ_Mesophyll_go_down.csv",
  row.names = FALSE
)

barplot(go_down, showCategory = 20)


library(patchwork)
library(ggplot2)

p1 <- barplot(go_up, showCategory = 7, label_format = 100) + ggtitle("Up")
p2 <- barplot(go_down, showCategory = 15, label_format = 100) + ggtitle("Down")

p <- p1 / p2 + plot_layout(height = c(1, 2))

pdf("./figures/leaf/LPZ_Mesophyll_de_go.pdf", width = 8, height = 8)
p
dev.off()
