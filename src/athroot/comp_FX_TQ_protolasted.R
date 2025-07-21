library(Seurat)
library(dplyr)

obj <- readRDS("./dataLib/athRoot/FX_TQ_v5.Rds")
colnames(obj@meta.data)

Idents(obj) <- "celltype"
markers <- FindAllMarkers(obj, logfc.threshold = 1)
markers <- markers[markers$p_val_adj <= 0.05, ]

ct <- unique(obj$celltype)
table(obj$tech)

# a named list
vs <- lapply(ct, function(x) {
  subobj <- subset(obj, celltype == x)
  message(x)
  if (length(unique(subobj$tech)) < 2) {
    return(NULL)
  }
  mm <- FindMarkers(
    subobj,
    group.by = "tech",
    ident.1 = "tq",
    ident.2 = "atroot"
  )
  submarkers <- markers[markers$cluster == x, ]
  mm[mm$p_val_adj <= 0.05 & rownames(mm) %in% submarkers$gene & abs(mm$avg_log2FC) > 1, ] %>%
    mutate(gene = rownames(.))
})

names(vs) <- ct

sapply(vs, nrow)

openxlsx::write.xlsx(vs, file = "./dataLib/athRoot/FX_TQ_protolasted_only_pos_TRUE.xlsx")


table(obj$tech)
Idents(obj) <- "tech"

mm <- FindMarkers(
  obj,
  group.by = "tech",
  ident.1 = "tq",
  ident.2 = "atroot"
)

mm <- mm[mm$p_val_adj < 0.05 & abs(mm$avg_log2FC) >= 2, ]
dim(mm)

write.csv(mm, file = "./dataLib/athRoot/TQ_vs_FX_protolasted.csv")
