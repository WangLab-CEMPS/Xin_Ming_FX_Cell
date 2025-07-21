library(Seurat)

obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public_anno.Rds")
obj <- JoinLayers(obj)
Idents(obj) <- "celltype"

table(obj$orig.ident)

obj_jc <- subset(obj, subset = orig.ident == "JC")
obj_lc <- subset(obj, subset = orig.ident == "LC")
obj_lpz <- subset(obj, subset = orig.ident == "LPZ")


table(obj_jc$celltype)
table(obj_lc$celltype)
table(obj_lpz$celltype)


deg_jc <- FindAllMarkers(obj_jc, logfc.threshold = 1)
deg_lc <- FindAllMarkers(obj_lc, logfc.threshold = 1)
deg_lpz <- FindAllMarkers(obj_lpz, logfc.threshold = 1)

deg_jc <- deg_jc[deg_jc$p_val_adj < 0.05 & deg_jc$avg_log2FC >= 1, ]
dim(deg_jc)
deg_lc <- deg_lc[deg_lc$p_val_adj < 0.05 & deg_lc$avg_log2FC >= 1, ]
dim(deg_lc)
deg_lpz <- deg_lpz[deg_lpz$p_val_adj < 0.05 & deg_lpz$avg_log2FC >= 1, ]
dim(deg_lpz)


table(unique(deg_jc$gene) %in% unique(deg_lc$gene))
gs <- setdiff(deg_tq$gene, deg_fx$gene)
length(gs)

go <- read.csv("./dataLib/athLeaf/GO_response_to_wounding.tsv", sep = "\t", header = FALSE)
go <- go[[1]]
length(go)
table(gs %in% go)

gg <- gs[gs %in% go]


write.csv(deg_fx, "./dataLib/athRoot/FX_TQ_FindAllMarkers_celltype_FX.csv")
write.csv(deg_tq, "./dataLib/athRoot/FX_TQ_FindAllMarkers_celltype_TQ.csv")

write.csv(gg, "./dataLib/athRoot/wounding_genes.csv")

deg_tq_wound <- deg_tq[deg_tq$gene %in% gg, ]
rownames(deg_tq_wound) <- NULL
dim(deg_tq_wound)
write.csv(deg_tq_wound, "./dataLib/athRoot/TQ_wound.csv")
