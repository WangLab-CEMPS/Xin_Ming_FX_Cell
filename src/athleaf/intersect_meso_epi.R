library(clusterProfiler)
library(org.At.tair.db)
library(ggplot2)
library(eulerr)
library(Seurat)
library(dplyr)


# read data
alias <- read.csv("~/ref/Athaliana/Araport11/gene_aliases_20220331_tidy.csv")

epi_sig <- read.csv("./dataLib/athLeaf/LC_LPZ_epi_de_sig.csv")
dim(epi_sig)

meso_sig <- read.csv("./dataLib/athLeaf/LC_LPZ_Mesophyll_de.csv") %>%
  subset(p_val < 0.05 & p_val_adj < 0.05 & abs(avg_log2FC) > 2)
dim(meso_sig)

epi_up <- epi_sig$geneId[epi_sig$avg_log2FC > 0]
epi_down <- epi_sig$geneId[epi_sig$avg_log2FC < 0]
length(epi_up)
length(epi_down)

meso_up <- meso_sig$geneId[meso_sig$avg_log2FC > 0]
meso_down <- meso_sig$geneId[meso_sig$avg_log2FC < 0]
length(meso_up)
length(meso_down)

# Venn diagram

plot(
  euler(list(epi = epi_up, meso = meso_up), shape = "ellipse"),
  quantities = list(type = "percent", font = 11),
  col = c("lightgray", "lightgray"),
  fill = c("#a1edc7", "#9dd9e8")
)

plot(
  euler(list(epi = epi_down, meso = meso_down), shape = "ellipse"),
  quantities = list(type = "percent", font = 11),
  col = c("lightgray", "lightgray"),
  fill = c("#FF9999", "#FFCC66")
)

# intersect
up_intersect <- intersect(epi_up, meso_up)
length(up_intersect)
down_intersect <- intersect(epi_down, meso_down)
length(down_intersect)


epi_down_specific <- setdiff(epi_down, meso_down)
length(epi_down_specific)

alias[alias$name %in% epi_down_specific, ]

meso_down_specific <- setdiff(meso_down, epi_down)
length(meso_down_specific)

orgdb <- org.At.tair.db

epi_down_specific_go <- enrichGO(
  epi_down_specific,
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

dotplot(epi_down_specific_go, showCategory = 20)

meso_down_specific_go <- enrichGO(
  meso_down_specific,
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

dotplot(meso_down_specific_go, showCategory = 20)


#
obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public.Rds")
obj <- JoinLayers(obj)
obj <- subset(obj, subset = orig.ident %in% c("LC", "LPZ"))

for (g in epi_down_specific) {
  print(g)
  p <- FeaturePlot(object = obj, features = g)
  ggsave(paste0("./Plots/athleaf/epi_down/", g, ".png"), p)
}
