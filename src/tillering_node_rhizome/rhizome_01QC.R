library(Seurat)
library(magrittr)
library(pigRutils)

# 创建Seurat对象
sample <- c(
  rhizome_S0 = "./rawdata/rhizome/stage0"
)

obj <- Read10X(data.dir = sample) %>%
  CreateSeuratObject(.,
    min.cells = 3,
    min.features = 200,
    project = names(sample)
  )

(p <- qc_plot(obj, max_nFeature = 4000, min_nFeature = 500, max_nCount = 10000, color = ""))


ggplot2::ggsave(paste0("./Plots/qc/", names(sample), "_qc.pdf"), p, width = 10, height = 8)

min(obj$nFeature_RNA)
max(obj$nFeature_RNA)
min(obj$nCount_RNA)
max(obj$nCount_RNA)

table(
  obj$nFeature_RNA > 499 &
    obj$nFeature_RNA < 4001 &
    obj$nCount_RNA > 998 &
    obj$nCount_RNA < 11001
)

# 过滤数据
obj <- subset(obj,
  subset = nFeature_RNA > 499 &
    nFeature_RNA < 4001 &
    nCount_RNA > 998 &
    nCount_RNA < 11001
)

# 保存
saveRDS(obj, paste0("./dataLib/tillering_node_rhizome/", names(sample), ".rds"))
