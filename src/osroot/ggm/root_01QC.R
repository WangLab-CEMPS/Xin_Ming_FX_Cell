library(Seurat)
library(magrittr)
library(pigRutils)

# 创建Seurat对象
sample <- c(
  # ggm1 = "./rawdata/osRoot/ggm1"
  ggm2 = "./rawdata/osRoot/ggm2"
)

(obj <- Read10X(data.dir = sample) %>%
  CreateSeuratObject(.,
    min.cells = 3,
    min.features = 200,
    project = names(sample)
  ))

# 添加细胞器基因比例
Mt <- read.table("~/wkdir/reference/Oryza_sativa/annotation/Mt_gene.txt")[[1]]
Pt <- read.table("~/wkdir/reference/Oryza_sativa/annotation/Pt_gene.txt")[[1]]

obj[["percent.mito"]] <- 100 * Matrix::colSums(
  GetAssayData(obj, layer = "counts")[Mt[Mt %in% rownames(obj)], ]
) / Matrix::colSums(GetAssayData(obj, layer = "counts"))

obj[["percent.chlo"]] <- 100 * Matrix::colSums(
  GetAssayData(obj, layer = "counts")[Pt[ Pt %in% rownames(obj)], ]
) / Matrix::colSums(GetAssayData(obj, layer = "counts"))


(p <- qc_plot(obj, max_nFeature = 6000, min_nFeature = 500))


ggplot2::ggsave(paste0("./Plots/", names(sample), "_qc.pdf"), p, width = 10, height = 8)

min(obj$nFeature_RNA)
max(obj$nFeature_RNA)
min(obj$nCount_RNA)
max(obj$nCount_RNA)

# 过滤数据
(obj <- subset(obj,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 6000 &
    nCount_RNA < 25000 &
    percent.mito < 5 & percent.chlo < 2
))

# 保存
saveRDS(obj, paste0("./dataLib/osRoot/", names(sample), ".rds"))
