library(Seurat)
library(dplyr)

obj <- readRDS("./dataLib/osRoot/ggm/ggm2_TQ2_strata.Rds")
obj <- JoinLayers(obj)
colnames(obj@meta.data)

Idents(obj) <- "celltype"
markers <- FindAllMarkers(obj, logfc.threshold = 1)
markers <- markers[markers$p_val_adj <= 0.05, ]

ct <- unique(obj$celltype)
table(obj$orig.ident)
obj[["tech"]] <- obj[["orig.ident"]]
obj$tech <- ifelse(obj$tech == "ggm2", "FX", "TQ")

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
    ident.1 = "TQ",
    ident.2 = "FX"
  )
  submarkers <- markers[markers$cluster == x, ]
  mm[mm$p_val_adj <= 0.05 & rownames(mm) %in% submarkers$gene & abs(mm$avg_log2FC) > 1, ] %>%
    mutate(gene = rownames(.))
})

names(vs) <- ct

sapply(vs, nrow)

openxlsx::write.xlsx(vs, file = "./dataLib/osRoot/TQ_FX_protolasted_celltype.xlsx")


table(obj$tech)
Idents(obj) <- "tech"

mm <- FindMarkers(
  obj,
  group.by = "tech",
  ident.1 = "TQ",
  ident.2 = "FX"
)

mm <- mm[mm$p_val_adj < 0.05 & abs(mm$avg_log2FC) >= 2, ]
dim(mm)

write.csv(mm, file = "./dataLib/osRoot/TQ_vs_FX_protolasted.csv")
