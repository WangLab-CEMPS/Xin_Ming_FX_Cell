library(clusterProfiler)
library(org.At.tair.db)
library(Seurat)
library(dplyr)


orgdb <- org.At.tair.db


alias <- read.csv("~/ref/Arabidopsis_thaliana/Araport11/gene_aliases_20220331_tidy.csv")

obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public.Rds")
obj <- JoinLayers(obj)
table(obj$orig.ident)

obj_JCLPZ <- subset(obj, subset = orig.ident %in% c("JC", "LPZ"))
obj_JC <- subset(obj, subset = orig.ident %in% c("JC"))
obj_LPZ <- subset(obj, subset = orig.ident %in% c("LPZ"))


obj_tmp <- obj_LPZ

deg <- FindMarkers(
  object = obj_tmp,
  group.by = "celltype",
  ident.1 = "Mesophyll",
  ident.2 = "Mesophyll*"
)

deg %>%
  as_tibble(rownames = "geneId") %>%
  left_join(alias, by = c("geneId" = "name")) %T>%
  write.csv("./dataLib/athLeaf/LPZ_Mesophyll_de.csv", row.names = FALSE) %>%
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
  "./dataLib/athLeaf/LPZ_Mesophyll_go.csv",
  row.names = FALSE
)
