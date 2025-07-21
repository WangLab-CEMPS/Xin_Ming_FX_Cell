library(Seurat)
library(magrittr)
library(pigRutils)

# 创建Seurat对象
sample <- c(
  # ggm1 = "./osRoot/osr1"
  ggm2 = "./osRoot/osr2"
)

(obj <- Read10X(data.dir = sample) %>%
  CreateSeuratObject(.,
    min.cells = 3,
    min.features = 200,
    project = names(sample)
  ))

# 添加细胞器基因比例
obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = "^ATMG")
obj[["percent.chlo"]] <- PercentageFeatureSet(obj, pattern = "^ATCG")

(p <- qc_plot(obj, max_nFeature = 8000, min_nFeature = 800, max_nCount = 100000))


ggplot2::ggsave(paste0("./Plots/", names(sample), "_qc.pdf"), p, width = 10, height = 8)

min(obj$nFeature_RNA)
max(obj$nFeature_RNA)
min(obj$nCount_RNA)
max(obj$nCount_RNA)

# 过滤数据
(
  obj <- subset(obj,
    subset = nFeature_RNA > 500 &
      nFeature_RNA < 8000 &
      nCount_RNA < 25000 &
      percent.mito < 5 & percent.chlo < 20
  )
)

# 保存
saveRDS(obj, paste0("./dataLib/athLeaf/", names(sample), ".rds"))
