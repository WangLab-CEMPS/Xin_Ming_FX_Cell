library(Seurat)
library(magrittr)
library(pigRutils)

samples <- c(
  refmethod1 = "./rawdata/ref_method/refmethod1",
  refmethod2 = "./rawdata/ref_method/refmethod2",
  snRNA1 = "./rawdata/ref_method/snRNA1",
  snRNA2 = "./rawdata/ref_method/snRNA2"
)


# 读取线粒体和叶绿体基因列表
Mt <- read.table("~/wkdir/reference/Oryza_sativa/annotation/Mt_gene.txt")[[1]]
Pt <- read.table("~/wkdir/reference/Oryza_sativa/annotation/Pt_gene.txt")[[1]]

# 处理每个样本
for (i in seq_along(samples)) {
  sample_name <- names(samples)[i]
  sample_path <- samples[i]
  
  obj <- Read10X(data.dir = sample_path) %>%
    CreateSeuratObject(.,
                       min.cells = 3,
                       min.features = 200,
                       project = sample_name
    )
  
  # 计算细胞器基因比例
  obj[["percent.mito"]] <- 100 * Matrix::colSums(
    GetAssayData(obj, layer = "counts")[Mt[Mt %in% rownames(obj)], ]
  ) / Matrix::colSums(GetAssayData(obj, layer = "counts"))
  
  obj[["percent.chlo"]] <- 100 * Matrix::colSums(
    GetAssayData(obj, layer = "counts")[Pt[Pt %in% rownames(obj)], ]
  ) / Matrix::colSums(GetAssayData(obj, layer = "counts"))
  
  # 绘制QC图并保存
  p <- qc_plot(obj, max_nFeature = 6000, min_nFeature = 500, max_nCount = 30000)
  ggplot2::ggsave(paste0("./Plots/", sample_name, "_qc.pdf"), p, width = 10, height = 8)
  
  # 打印统计信息
  cat(sprintf("\nSample: %s\n", sample_name))
  cat(sprintf("min nFeature_RNA: %d\n", min(obj$nFeature_RNA)))
  cat(sprintf("max nFeature_RNA: %d\n", max(obj$nFeature_RNA)))
  cat(sprintf("min nCount_RNA: %d\n", min(obj$nCount_RNA)))
  cat(sprintf("max nCount_RNA: %d\n", max(obj$nCount_RNA)))
  
  # 过滤数据
  obj <- subset(obj,
                subset = nFeature_RNA > 500 &
                  nFeature_RNA < 6000 &
                  nCount_RNA < 25000 &
                  percent.mito < 5 & 
                  percent.chlo < 2
  )
  
  # 保存处理后的对象
  saveRDS(obj, paste0("./dataLib/ref_method/", sample_name, ".rds"))
}
