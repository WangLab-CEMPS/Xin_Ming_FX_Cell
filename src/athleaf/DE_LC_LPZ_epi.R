library(dplyr)
library(Seurat)
library(ggplot2)
library(magrittr)
library(patchwork)
library(org.At.tair.db)
library(clusterProfiler)

# load data
obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public.Rds")
obj <- subset(obj, subset = orig.ident %in% c("LC", "LPZ"))

DimPlot(obj, label = TRUE, label.box = TRUE)

epi <- subset(obj, subset = seurat_clusters %in% c("0", "8"))
epi <- JoinLayers(epi)

DimPlot(epi, group.by = "orig.ident", label = TRUE, label.box = TRUE)

epi <- FindClusters(epi, resolution = 0.2)
DimPlot(epi, label = TRUE, label.box = TRUE)

epi[["wound"]] <- ifelse(epi$seurat_clusters == "2", "Y", "N")
DimPlot(epi, group.by = "wound")

deg <- FindMarkers(
  object = epi,
  group.by = "wound",
  ident.1 = "N",
  ident.2 = "Y"
)

alias <- read.csv("~/ref/Athaliana/Araport11/gene_aliases_20220331_tidy.csv")

deg %>%
  as_tibble(rownames = "geneId") %>%
  left_join(alias, by = c("geneId" = "name")) %T>%
  write.csv("./dataLib/athLeaf/LC_LPZ_epi_de.csv", row.names = FALSE) %>%
  subset(p_val < 0.05 & p_val_adj < 0.05 & abs(avg_log2FC) > 2) -> dd

dd %>%
  mutate(sig = ifelse(avg_log2FC > 0, "Up", "Down")) %>%
  write.csv("./dataLib/athLeaf/LC_LPZ_epi_de_sig.csv", row.names = FALSE)

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
  "./dataLib/athLeaf/LC_LPZ_epi_go_up.csv",
  row.names = FALSE
)

dotplot(go_up, showCategory = 20)

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
  "./dataLib/athLeaf/LC_LPZ_epi_go_down.csv",
  row.names = FALSE
)

dotplot(go_down, showCategory = 20)


p1 <- dotplot(go_up, showCategory = 20, label_format = 100) + ggtitle("Up")
p2 <- dotplot(go_down, showCategory = 20, label_format = 100) + ggtitle("Down")

p <- p1 / p2 + plot_layout(height = c(1, 1.5))

pdf("./figures/leaf/LC_LPZ_epi_de_go.pdf", width = 8, height = 10)
p
dev.off()
