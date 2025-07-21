library(Seurat)
library(dplyr)
library(ggplot2)

obj <- readRDS("./dataLib/athRoot/FX_TQ_v5.rds")
obj <- FindClusters(obj, resolution = 0.45)
DimPlot(obj, label = TRUE)

subobj <- subset(obj, subset = seurat_clusters %in% c(21, 8))
DimPlot(subobj, label = TRUE)
subobj <- FindClusters(subobj, resolution = 0.1)
DimPlot(subobj, label = TRUE)
cc <- Cells(subobj)[subobj$seurat_clusters == 3]

obj[["tech"]] <- obj$orig.ident %>% gsub("[1-2]$", "", .)
table(obj$tech)
DimPlot(obj, label = TRUE, split.by = "tech")


# 1, 20
Root_cap <- c("AT1G33280", "AT1G79580")
# 9, 8 # "AT3G23580", "AT5G13840"
Dividing_Cells <- c("AT2G22490", "AT2G26760", "AT3G25980")
# 10
QC <- c("AT3G20840", "AT1G51190", "AT2G04025")
# 6
Root_hair <- c("AT5G49270", "AT1G16440", "AT2G03720")
# 14, 15, 12 # "AT5G07990", "AT3G21670", "AT1G09750"
Cortex <- c("AT1G29025", "AT1G62510", "AT5G64620")
# 18, 19, 4 # "AT1G34670", "AT5G06200", "AT3G11550", "AT2G36100"
Endodermis <- c("AT4G11290", "AT5G57620", "AT4G17280", "AT4G20140")
# 3 # "AT4G14020", "AT1G05760"
Phloem <- c("AT4G19840", "AT2G15310", "AT1G79430")
# 11
Xylem <- c("AT1G68810", "AT1G20850", "AT4G35350")
# 5
Lateral_root <- c("AT1G77690", "AT4G14650")
# 22
Epidermis <- c("AT2G42840", "AT4G21750")
# 21
photosynthetic_cell <- c("AT3G54890")
# 17, 2
Stem_cell_niche <- c("AT2G28790", "AT1G73590")
# 13
Atrichoblast <- c("AT1G79840", "AT5G40330", "AT2G37260")
# 0, 7, 16 # "AT1G68740",
Stele <- c("AT1G32450", "AT5G48070", "AT2G31083", "AT3G25710")
# 0, 16
procambium <- c("AT1G20700", "AT1G61660")
# Columella AT1G17400
# 0, 16
pericycle <- c("AT3G23430", "AT3G45700", "AT5G43180", "AT4G30450")

# -------

obj[["celltype"]] <- obj$seurat_clusters %>% as.character()
obj$celltype[Cells(obj) %in% cc] <- "22"

Idents(obj) <- obj$celltype

FeaturePlot(obj, features = Stele, order = TRUE, label = TRUE)
DotPlot(obj, features = Stele)


obj$celltype[obj$celltype %in% c("6")] <- "Root-Hair"
obj$celltype[obj$celltype %in% c("1", "20")] <- "Root-Cap"
obj$celltype[obj$celltype %in% c("13")] <- "Atrichoblast"
obj$celltype[obj$celltype %in% c("21")] <- "Unknown"
obj$celltype[obj$celltype %in% c("22")] <- "Epidermis"
obj$celltype[obj$celltype %in% c("5")] <- "Lateral-Root"
obj$celltype[obj$celltype %in% c("11")] <- "Xylem"
obj$celltype[obj$celltype %in% c("3")] <- "Phloem"
obj$celltype[obj$celltype %in% c("18", "19", "4")] <- "Endodermis"
obj$celltype[obj$celltype %in% c("14", "15", "12")] <- "Cortex"
obj$celltype[obj$celltype %in% c("0", "16", "7")] <- "Stele"
obj$celltype[obj$celltype %in% c("10")] <- "QC"
obj$celltype[obj$celltype %in% c("8", "9")] <- "Dividing-Cells"
obj$celltype[obj$celltype %in% c("17", "2")] <- "Stem-Cell-Niche"

DimPlot(obj, group.by = "celltype")

saveRDS(obj, "./dataLib/athRoot/FX_TQ_v5.rds")

# dotplot --------------
obj$seurat_clusters <- as.character(obj$seurat_clusters)
obj$seurat_clusters[obj$celltype == "Epidermis"] <- "22"


names(Epidermis) <- rep("Epi", length(Epidermis))
names(Dividing_Cells) <- rep("Div", length(Dividing_Cells))
names(QC) <- rep("QC", length(QC))
names(Root_hair) <- rep("RH", length(Root_hair))
names(Cortex) <- rep("Cor", length(Cortex))
names(Endodermis) <- rep("Endo", length(Endodermis))
names(Phloem) <- rep("Phlo", length(Phloem))
names(Xylem) <- rep("Xyle", length(Xylem))
names(Lateral_root) <- rep("LR", length(Lateral_root))
names(Root_cap) <- rep("RC", length(Root_cap))
names(Stem_cell_niche) <- rep("SCN", length(Stem_cell_niche))
names(Atrichoblast) <- rep("Atri", length(Atrichoblast))
names(Stele) <- rep("Stele", length(Stele))
names(photosynthetic_cell) <- rep("Un", length(photosynthetic_cell))

ff <- c(
  Dividing_Cells,
  QC, Root_cap, Root_hair, Atrichoblast, Cortex,
  Endodermis, Phloem, Xylem,
  Lateral_root,
  Stem_cell_niche, Stele, Epidermis, photosynthetic_cell
)

p2 <- DotPlot(obj, features = ff, cols = c("#eeeeee", "#ff4800"), col.min = 0, col.max = 3) +
  theme_test() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "cm"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 80, vjust = 0.5, hjust = 0.5, size = 10)
  )

lev <- c(
  "8", "9", "10", "1", "20",
  "6", "13", "15", "14", "12",
  "18", "19", "4",
  "3", "11", "5",
  "17", "2",
  "0", "7", "16",
  "22", "21"
)

p2$data$id <- factor(p2$data$id, levels = lev)
p2

ggsave("./figures/athroot/dotplot.pdf", p2, height = 4.6, width = 10)
