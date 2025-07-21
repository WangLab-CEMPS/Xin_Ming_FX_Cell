library(Seurat)
library(ggplot2)
library(pigRutils)


obj <- Read10X("./rawdata/maize/root_pub") %>%
  CreateSeuratObject(
    min.cells = 3,
    min.features = 200,
    project = "pub"
  )

(p <- qc_plot(obj, min_nFeature = 800, max_nFeature = 8000, max_nCount = 60000, min_nCount = 3000, color = ""))

ggsave("./Plots/maize_pub_qc.pdf", p, width = 7, height = 6)


min(obj$nFeature_RNA)
max(obj$nFeature_RNA)
min(obj$nCount_RNA)
max(obj$nCount_RNA)

table(obj$nFeature_RNA < 8000 & obj$nFeature_RNA > 800 & obj$nCount_RNA > 1500 & obj$nCount_RNA < 60000)


(obj <- subset(
  obj,
  subset = nFeature_RNA < 8000 &
    nFeature_RNA > 800 &
    nCount_RNA > 1500 &
    nCount_RNA < 60000
))

saveRDS(obj, "./dataLib/maize/science_maize.rds")
