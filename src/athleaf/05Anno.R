library(Seurat)
library(ggplot2)

obj <- readRDS("./dataLib/athLeaf/LC_LPC_public.Rds")

phloem_parenchyma <- c("AT3G48740", "AT5G23660", "AT4G13770", "AT3G11930")
companion_cell <- c("AT1G22710", "AT5G06850", "AT1G79430", "AT5G57350")
ab_PC <- c("AT2G26580")
ad_PC <- c("AT5G13930", "AT1G02205")
PC <- c("AT1G04040")
myrosinase_cell <- c("AT5G26000", "AT1G52342")
SE <- c("AT3G01680", "AT3G03270")
Spongy <- c("AT2G45190")
Palisade <- c("AT5G17790")
bundle_sheet <- c("AT5G41920", "AT1G77990", "AT2G37860")
GC <- c("AT2G46070", "AT3G18040", "AT1G08810", "AT1G22690")
epidermis <- c("AT4G21750", "AT1G51500", "AT1G68530", "AT1G27950")
vasculature <- c("AT5G12140", "AT1G16410", "AT5G61480")
hydathode <- c("AT3G54420", "AT1G28230", "AT3G16670", "AT1G62510")
mesophyll <- c("AT1G29910", "AT3G01500", "AT4G12970") #  "AT2G05100"
cell_cycle <- c("AT1G18370", "AT3G25980", "AT5G25090")
trichome <- c("AT4G01060", "AT2G30420", "AT2G30424")

names(phloem_parenchyma) <- rep("PhloemParenchyma", length(phloem_parenchyma))
names(companion_cell) <- rep("CompanionCell", length(companion_cell))
names(ab_PC) <- rep("ab_pc", length(ab_PC))
names(ad_PC) <- rep("ad_pc", length(ad_PC))
names(PC) <- rep("PC", length(PC))
names(myrosinase_cell) <- rep("MI", length(myrosinase_cell))
names(SE) <- rep("SE", length(SE))
names(Spongy) <- rep("Spongy", length(Spongy))
names(Palisade) <- rep("Palisade", length(Palisade))
names(bundle_sheet) <- rep("bundle sheet", length(bundle_sheet))
names(GC) <- rep("GC", length(GC))
names(epidermis) <- rep("Epidermis", length(epidermis))
names(hydathode) <- rep("Hydathode", length(hydathode))
names(mesophyll) <- rep("Mesophyll", length(mesophyll))
names(cell_cycle) <- rep("CellCycle", length(cell_cycle))
names(vasculature) <- rep("Vasculature", length(vasculature))
names(trichome) <- rep("Trichome", length(trichome))


ff <- c(
  SE,
  companion_cell, phloem_parenchyma, vasculature, bundle_sheet, myrosinase_cell, GC,
  trichome, epidermis, ad_PC, ab_PC, PC, hydathode, mesophyll, Spongy, Palisade,
  cell_cycle
)


p2 <- DotPlot(obj, features = ff, cols = c("#eeeeee", "#ff4800"), col.min = 0, col.max = 3) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "cm"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 80, vjust = 0.5, hjust = 0.5, size = 10)
  )

lev <- c(
  "12", "11", "13", "4",
  "18", "14",
  "8", "0", "17", "16", "1", "10",
  "15", "6",
  "7", "2", "3", "5",
  "9"
)

p2$data$id <- factor(p2$data$id, levels = lev)
p2

ggsave("./Plots/LC_LPZ_public_anno_dot.pdf", p2, height = 5, width = 15)

# ------------------------
# obj[["percent.chlo"]] <- PercentageFeatureSet(obj, pattern = "^ATCG")

obj[["celltype"]] <- obj$seurat_clusters
obj[["celltype"]] <- as.character(obj$celltype)

obj$celltype[obj$celltype %in% c("9")] <- "Cell-Cycle"
obj$celltype[obj$celltype %in% c("12")] <- "Companion-Cell"
obj$celltype[obj$celltype %in% c("11")] <- "Phloem"
obj$celltype[obj$celltype %in% c("13")] <- "Xylem"
obj$celltype[obj$celltype %in% c("4")] <- "Bundle-Sheet"
obj$celltype[obj$celltype %in% c("18")] <- "Myrosin-Idioblast" # myrosinase_cell
obj$celltype[obj$celltype %in% c("14")] <- "Guard-Cell"
obj$celltype[obj$celltype %in% c("8")] <- "Epidermis(Trichome)"
obj$celltype[obj$celltype %in% c("0")] <- "Epidermis"
obj$celltype[obj$celltype %in% c("16", "17")] <- "Hydathode"
obj$celltype[obj$celltype %in% c("3", "2", "5", "7", "10")] <- "Mesophyll"
obj$celltype[obj$celltype %in% c("1", "6", "15")] <- "Mesophyll*"

col <- pigRutils::select_colors("col_33")
col <- col[c(1, 7, 5, 4, 3, 6, 2, 14, 9, 10, 8, 12)]

p_umap <- DimPlot(
  obj,
  group.by = "celltype",
  reduction = "harmony_umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4,
  cols = col
) &
  tidydr::theme_dr() &
  theme(panel.grid = element_blank(), title = element_blank()) &
  guides(
    color = guide_legend(
      keywidth = 0.2,
      override.aes = list(size = 2, stroke = 2)
    )
  )

p_umap


ggsave(filename = "./Plots/LC_LPZ_public_umap_Anno.pdf", p_umap, width = 8.2, height = 6)

saveRDS(obj, "./dataLib/athLeaf/LC_LPC_public.Rds")
