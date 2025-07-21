library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(pigRutils)


obj <- readRDS("./dataLib/athLeaf/scsnRNA_harmony.Rds")
table(obj$orig.ident)


obj$orig.ident[obj$orig.ident == "SeuratProject"] <- "Public_data"
obj$orig.ident[obj$orig.ident == "LC"] <- "scLC"
obj$orig.ident[obj$orig.ident == "LPZ"] <- "scLPZ"
table(obj$orig.ident)


obj$tech <- obj$orig.ident
obj$tech[obj$tech %in% c("snLC", "snLPZ")] <- "snRNA"
obj$tech[obj$tech %in% c("scLC", "scLPZ")] <- "FX"
obj$tech[obj$tech == "Public_data"] <- "scRNA"
table(obj$tech)


obj$group <- obj$orig.ident
obj$group[obj$group %in% c("scLC", "snLC")] <- "LC"
obj$group[obj$group %in% c("scLPZ", "snLPZ")] <- "LPZ"
obj$group[obj$group == "Public_data"] <- "Common"
table(obj$group)


obj$tech <- factor(obj$tech, levels = c("snRNA", "FX", "scRNA"))
obj$group <- factor(obj$group, levels = c("LC", "LPZ", "Common"))



# add annotation ----------------
phloem_parenchyma <- c("AT3G48740", "AT5G23660", "AT3G11930") # "AT4G13770",
companion_cell <- c("AT1G22710", "AT5G57350") # "AT1G79430", "AT5G06850",
ab_PC <- c("AT2G26580")
ad_PC <- c("AT5G13930", "AT1G02205")
PC <- c("AT1G04040")
myrosinase_cell <- c("AT5G26000", "AT1G52342")
SE <- c("AT3G01680", "AT3G03270")
Spongy <- c("AT2G45190")
Palisade <- c("AT5G17790")
bundle_sheet <- c("AT5G41920", "AT1G77990", "AT2G37860")
GC <- c("AT2G46070", "AT1G08810", "AT1G22690") # "AT3G18040",
epidermis <- c("AT4G21750", "AT1G51500", "AT1G68530") # "AT1G27950"
xylem <- c("AT1G16410", "AT5G61480") # "AT5G12140",
hydathode <- c("AT3G54420", "AT3G16670", "AT1G62510") # "AT1G28230",
mesophyll <- c("AT3G01500", "AT4G12970", "AT1G29910") #  "AT2G05100"
cell_cycle <- c("AT1G18370", "AT3G25980", "AT5G25090")
trichome <- c("AT4G01060", "AT2G30420", "AT2G30424")

names(phloem_parenchyma) <- rep("Ph", length(phloem_parenchyma))
names(companion_cell) <- rep("Cc", length(companion_cell))
names(myrosinase_cell) <- rep("Mi", length(myrosinase_cell))
names(SE) <- rep("Se", length(SE))
names(bundle_sheet) <- rep("BSC", length(bundle_sheet))
names(GC) <- rep("GC", length(GC))
names(epidermis) <- rep("Ep", length(epidermis))
names(hydathode) <- rep("Hy", length(hydathode))
names(mesophyll) <- rep("Mp", length(mesophyll))
names(cell_cycle) <- rep("Dc", length(cell_cycle))
names(xylem) <- rep("Xy", length(xylem))


ff <- c(
  cell_cycle,
  companion_cell, SE, phloem_parenchyma, xylem,
  bundle_sheet, myrosinase_cell, GC,
  epidermis, hydathode, mesophyll
)


p2 <- DotPlot(obj, features = ff, cols = c("#eeeeee", "#ff4800"), col.min = 0, col.max = 3) +
  theme_bw(base_size = 11) +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "cm"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 80, vjust = 0.5, hjust = 0.5, size = 10)
  )

lev <- c(
  "12", "10", "19", "11",
  "6", "5", "18", "14",
  "15", "9", "1",
  "16", "17", "4", "0", "8",
  "2", "3",
  "7", "13"
)

p2$data$id <- factor(p2$data$id, levels = lev)
p2

ggsave("./figures/leaf/snRNA/anno_dot.pdf", p2, height = 5, width = 7.2)

# ------------------------
DimPlot(obj, group.by = "seurat_clusters", label = TRUE, label.box = TRUE)

obj[["celltype"]] <- obj$seurat_clusters
obj[["celltype"]] <- as.character(obj$celltype)

obj$celltype[obj$celltype %in% c("12")] <- "dividing cell"
obj$celltype[obj$celltype %in% c("10")] <- "companion cell"
obj$celltype[obj$celltype %in% c("19")] <- "sieve element cell"
obj$celltype[obj$celltype %in% c("11")] <- "phloem"
obj$celltype[obj$celltype %in% c("6")] <- "xylem"
obj$celltype[obj$celltype %in% c("5")] <- "bundle sheath cell"
obj$celltype[obj$celltype %in% c("18")] <- "myrosin idioblast" # myrosinase_cell
obj$celltype[obj$celltype %in% c("14")] <- "guard cell"
obj$celltype[obj$celltype %in% c("15", "9", "1")] <- "epidermis"
obj$celltype[obj$celltype %in% c("16", "17")] <- "hydathode"
obj$celltype[obj$celltype %in% c("13", "7", "3", "2")] <- "mesophyll(1)"
obj$celltype[obj$celltype %in% c("4", "8")] <- "mesophyll(2)"
obj$celltype[obj$celltype %in% c("0")] <- "mesophyll(3)"

col <- pigRutils::select_colors("col_33")
col <- col[c(1, 7, 5, 4, 3, 6, 2, 14, 9, 10, 8, 12, 13)]

p_umap <- DimPlot(
  obj,
  group.by = "celltype",
  reduction = "harmony_umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 2,
  cols = col,
  raster = TRUE
) &
  tidydr::theme_dr() &
  theme(panel.grid = element_blank(), title = element_blank()) &
  guides(
    color = guide_legend(
      keywidth = 0.2,
      override.aes = list(size = 2, stroke = 2)
    )
  )

ggsave(filename = "./figures/leaf/snRNA/anno_umap.pdf", p_umap, width = 8.2, height = 6)


col <- pigRutils::select_colors("col_33")
col <- col[c(1, 7, 5, 4, 3, 6, 2, 14, 9, 10, 8, 12, 13)]

pp <- p <- embSCdim(
  obj, "harmony_umap",
  group_by = "celltype",
  colors = col,
  save = TRUE,
  width = 9.2, height = 5.5,
  save_path = "./figures/leaf/snRNA/scsnRNA_harmony_umap_with_celltype.pdf"
)


# ----------
colnames(obj@meta.data)

p <- embSCdim(
  obj, "harmony_umap",
  group_by = "harmony_clusters",
  split_by = "group",
  colors = pigRutils::select_colors("col_33"),
  save = TRUE,
  width = 17.2, height = 5.5,
  save_path = "./figures/leaf/snRNA/scsnRNA_harmony_umap_with_cluster_split.pdf"
)

p <- embSCdim(
  obj, "harmony_umap",
  group_by = "celltype",
  split_by = "group",
  colors = col,
  save = TRUE,
  width = 17.2, height = 5.5,
  save_path = "./figures/leaf/snRNA/scsnRNA_harmony_umap_with_celltype_split.pdf"
)


p <- embSCdim(
  obj, "harmony_umap",
  group_by = "tech",
  split_by = "group",
  colors = pigRutils::select_colors("col_33"),
  save = TRUE,
  width = 17.2, height = 5.5,
  save_path = "./figures/leaf/snRNA/scsnRNA_harmony_umap_with_tech_split.pdf"
)


saveRDS(obj, "./dataLib/athLeaf/scsnRNA_harmony.Rds")