library(Seurat)

obj <- readRDS("./dataLib/athRoot/FX_TQ_v5.Rds")


Idents(obj)
deg <- FindAllMarkers(obj, logfc.threshold = 1)

write.csv(deg, "./dataLib/athRoot/FX_TQ_FindAllMarkers_cluster.csv")

Idents(obj) <- "celltype"
deg2 <- FindAllMarkers(obj, logfc.threshold = 1)

write.csv(deg2, "./dataLib/athRoot/FX_TQ_FindAllMarkers_celltype.csv")


table(obj$tech, obj$celltype)
obj_fx <- subset(obj, subset = tech == "atroot")
obj_tq <- subset(obj, subset = tech == "tq")

table(obj_fx$celltype)
table(obj_tq$celltype)

obj_tq <- subset(obj_tq, subset = ! celltype %in% c("Unknown", "Epidermis"))

deg_fx <- FindAllMarkers(obj_fx, logfc.threshold = 1)
deg_tq <- FindAllMarkers(obj_tq, logfc.threshold = 1)

deg_fx <- deg_fx[deg_fx$p_val_adj < 0.05 & deg_fx$avg_log2FC >= 1, ]
dim(deg_fx)
deg_tq <- deg_tq[deg_tq$p_val_adj < 0.05 & deg_tq$avg_log2FC >= 1, ]
dim(deg_tq)

table(unique(deg_tq$gene) %in% unique(deg_fx$gene))
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
