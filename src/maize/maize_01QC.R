library(Seurat)
library(magrittr)
library(pigRutils)

sample <- c(maize = "./rawdata/maize/root")

obj <- Read10X(data.dir = sample) %>%
  CreateSeuratObject(.,
    min.cells = 3,
    min.features = 200,
    project = "maize"
  )


(p <- qc_plot(obj, max_nFeature = 8000, min_nFeature = 800, max_nCount = 60000, min_nCount = 3000, color = ""))

ggplot2::ggsave(paste0("./Plots/maize_qc.pdf"), p, width = 7, height = 6)

# 打印统计信息
min(obj$nFeature_RNA)
max(obj$nFeature_RNA)
min(obj$nCount_RNA)
max(obj$nCount_RNA)

table(obj$nFeature_RNA < 8000 & obj$nFeature_RNA > 800 & obj$nCount_RNA > 1500 & obj$nCount_RNA < 60000)


# 过滤数据
(obj <- subset(obj,
  subset = nFeature_RNA > 800 &
    nFeature_RNA < 8000 &
    nCount_RNA < 60000 &
    nCount_RNA > 1500
))

# 保存处理后的对象
saveRDS(obj, paste0("./dataLib/maize/maize.rds"))
