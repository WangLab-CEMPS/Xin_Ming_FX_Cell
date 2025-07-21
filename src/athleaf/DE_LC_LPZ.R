library(Seurat)
library(magrittr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(org.At.tair.db)
library(clusterProfiler)

# load data
obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public_anno.Rds")
obj <- JoinLayers(obj)
obj <- subset(obj, subset = orig.ident %in% c("LC", "LPZ"))

alias <- read.csv("~/ref/Athaliana/Araport11/gene_aliases_20220331_tidy.csv")

# DE
table(obj$celltype)
obj$celltype <- gsub("Mesophyll\\*", "Mesophyll", obj$celltype)
selected_ct <- c("Bundle-Sheet", "Phloem", "Epidermis", "Hydathode", "Mesophyll")
orgdb <- org.At.tair.db

for (ct in selected_ct) {
  print(ct)
  sub_obj <- subset(obj, subset = celltype == ct)
  deg <- FindMarkers(
    object = sub_obj,
    group.by = "orig.ident",
    ident.1 = "LC",
    ident.2 = "LPZ"
  )
  deg %>%
    as_tibble(rownames = "geneId") %>%
    left_join(alias, by = c("geneId" = "name")) %T>%
    write.csv(paste0("./dataLib/athLeaf/de/LC_LPZ_", ct, "_de.csv"), row.names = FALSE) %>%
    subset(p_val < 0.05 & p_val_adj < 0.05 & abs(avg_log2FC) > 1.58) -> dd

  dd %>%
    mutate(sig = ifelse(avg_log2FC > 0, "Up", "Down")) %>%
    write.csv(paste0("./dataLib/athLeaf/de/LC_LPZ_", ct, "_de_sig.csv"), row.names = FALSE)

  up <- dd$geneId[dd$avg_log2FC > 0]
  down <- dd$geneId[dd$avg_log2FC < 0]

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
    paste0("./dataLib/athLeaf/go/LC_LPZ_", ct, "_go_up.csv"),
    row.names = FALSE
  )

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
    paste0("./dataLib/athLeaf/go/LC_LPZ_", ct, "_go_down.csv"),
    row.names = FALSE
  )
}
