library(Seurat)
library(dplyr)

maize <- readRDS("./dataLib/maize/maize.rds")
pub <- readRDS("./dataLib/maize/science_maize.rds")

colnames(pub@meta.data)

maize_pub  <- readRDS("./dataLib/maize/maize_pub.Rds")

subobj <- maize_pub[, maize_pub$celltype == "Root Cap"]
DimPlot(subobj, group.by = "orig.ident")
