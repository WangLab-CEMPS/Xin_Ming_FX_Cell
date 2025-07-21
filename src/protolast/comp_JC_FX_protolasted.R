library(Seurat)
library(dplyr)

obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public_anno.Rds")
obj <- JoinLayers(obj)

table(obj$orig.ident)
obj <- subset(obj, subset = orig.ident %in% c("JC", "LC"))
table(obj$celltype)
obj$celltype <- gsub("Mesophyll\\*", "Mesophyll", obj$celltype)


Idents(obj) <- "celltype"
# markers <- FindAllMarkers(obj, logfc.threshold = 1)
# markers <- markers[markers$p_val_adj <= 0.05, ]

ct <- unique(obj$celltype)
table(obj$orig.ident)

# a named list
vs <- lapply(ct, function(x) {
  subobj <- subset(obj, celltype == x)
  message(x)
  if (length(unique(subobj$orig.ident)) < 2) {
    return(NULL)
  }
  mm <- FindMarkers(
    subobj,
    group.by = "orig.ident",
    ident.1 = "JC",
    ident.2 = "LC"
  )
  # submarkers <- markers[markers$cluster == x, ]
  # mm[mm$p_val_adj <= 0.05 & rownames(mm) %in% submarkers$gene & abs(mm$avg_log2FC) > 1, ] %>%
  #   mutate(gene = rownames(.))

  mm %>% mutate(gene = rownames(.))
})

names(vs) <- ct

sapply(vs, nrow)

openxlsx::write.xlsx(vs, file = "./dataLib/athLeaf/JC_vs_LC_protolasted.xlsx")
