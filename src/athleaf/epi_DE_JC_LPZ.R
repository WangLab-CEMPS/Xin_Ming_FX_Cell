library(clusterProfiler)
library(org.At.tair.db)
library(Seurat)
library(dplyr)

orgdb <- org.At.tair.db

alias <- read.csv("~/ref/Arabidopsis_thaliana/Araport11/gene_aliases_20220331_tidy.csv")

obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public.Rds")
obj <- subset(obj, subset = orig.ident %in% c("JC", "LPZ"))
obj <- JoinLayers(obj)

DimPlot(obj, label = TRUE, label.box = TRUE)


epi <- subset(obj, subset = seurat_clusters %in% c("0", "8"))
epi_jc <- subset(epi, subset = orig.ident == "JC")
epi_lpz <- subset(epi, subset = orig.ident == "LPZ")


epi_jc <- FindClusters(epi_jc, resolution = 0.2)
DimPlot(epi_jc, label = TRUE, label.box = TRUE)
epi_jc[["wound"]] <- ifelse(epi_jc$seurat_clusters == "1", "Y", "N")


epi_lpz <- FindClusters(epi_lpz, resolution = 0.3)
DimPlot(epi_lpz, label = TRUE, label.box = TRUE)
epi_lpz[["wound"]] <- ifelse(epi_lpz$seurat_clusters %in% c("5", "3", "0"), "Y", "N")

epi_tmp <- epi_lpz

deg <- FindMarkers(
  object = epi_tmp,
  group.by = "wound",
  ident.1 = "N",
  ident.2 = "Y"
)

deg %>%
  as_tibble(rownames = "geneId") %>%
  left_join(alias, by = c("geneId" = "name")) %T>%
  write.csv("./dataLib/athLeaf/LPZ_Epidermis_de.csv", row.names = FALSE) %>%
  subset(p_val < 0.05 & p_val_adj < 0.05 & abs(avg_log2FC) > 2) %>%
  mutate(sig = ifelse(avg_log2FC > 0, "Up", "Down")) -> dd


up <- dd$geneId[dd$avg_log2FC > 0]
down <- dd$geneId[dd$avg_log2FC < 0]


gores <- compareCluster(
  gene = list(up = up, down = down),
  fun = "enrichGO",
  OrgDb = orgdb,
  keyType = "TAIR",
  ont = "BP"
) %>% clusterProfiler::simplify(
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

write.csv(
  as.data.frame(gores),
  "./dataLib/athLeaf/LPZ_Epidermis_go.csv",
  row.names = FALSE
)
