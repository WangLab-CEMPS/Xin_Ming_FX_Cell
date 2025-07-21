library(Seurat)
library(magrittr)
library(dplyr)
library(org.At.tair.db)
library(clusterProfiler)

# load data
obj <- readRDS("./dataLib/athLeaf/scsnRNA_harmony.Rds")
obj <- JoinLayers(obj)
table(obj$orig.ident)
obj <- subset(obj, subset = orig.ident %in% c("snLC", "snLPZ"))

table(obj$celltype)
obj$celltype <- gsub("mesophyll\\([123]\\)", "mesophyll", obj$celltype)
obj$celltype <- gsub(" ", "-", obj$celltype)
table(obj$celltype)

selected_ct <- c("bundle-sheath-cell", "phloem", "epidermis", "hydathode", "mesophyll")


# DE
alias <- read.csv("~/ref/Arabidopsis_thaliana/Araport11/gene_aliases_20220331_tidy.csv")

orgdb <- org.At.tair.db

for (ct in selected_ct) {
  print(ct)
  sub_obj <- subset(obj, subset = celltype == ct)
  deg <- FindMarkers(
    object = sub_obj,
    group.by = "orig.ident",
    ident.1 = "snLC",
    ident.2 = "snLPZ"
  )
  deg %>%
    as_tibble(rownames = "geneId") %>%
    left_join(alias, by = c("geneId" = "name")) %T>%
    write.csv(paste0("./dataLib/athLeaf/snRNA/de/snLC_snLPZ_", ct, "_de.csv"), row.names = FALSE) %>%
    subset(p_val < 0.05 & p_val_adj < 0.05 & abs(avg_log2FC) > 1.58) -> dd

  dd %>%
    mutate(sig = ifelse(avg_log2FC > 0, "Up", "Down")) %>%
    write.csv(paste0("./dataLib/athLeaf/snRNA/de/snLC_snLPZ_", ct, "_de_sig.csv"), row.names = FALSE)

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
    paste0("./dataLib/athLeaf/snRNA/go/snLC_snLPZ_", ct, "_go_up.csv"),
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
    paste0("./dataLib/athLeaf/snRNA/go/snLC_snLPZ_", ct, "_go_down.csv"),
    row.names = FALSE
  )
}
