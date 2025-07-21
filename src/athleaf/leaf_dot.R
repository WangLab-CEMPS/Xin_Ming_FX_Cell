library(Seurat)
library(ggplot2)

obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public.Rds")

bundle_sheet <- c("AT5G41920", "AT1G77990", "AT2G37860")
cell_cycle <- c("AT1G18370", "AT3G25980", "AT5G25090")
epidermis <- c("AT4G21750", "AT1G51500", "AT1G68530") # "AT1G27950"
trichome <- c("AT4G01060", "AT2G30420") # "AT2G30424"
GC <- c("AT2G46070", "AT1G08810", "AT1G22690") # "AT3G18040"
hydathode <- c("AT3G54420", "AT3G16670", "AT1G62510") # "AT1G28230",
mesophyll <- c("AT1G29910", "AT3G01500", "AT4G12970")
myrosinase_cell <- c("AT5G26000", "AT1G52342")
phloem <- c("AT3G48740", "AT5G23660", "AT3G11930") # "AT4G13770"
SE <- c("AT1G22710", "AT1G79430", "AT5G57350") # "AT5G06850"
xylem <- c("AT1G16410", "AT5G61480") # "AT5G12140"

names(bundle_sheet) <- rep("BS", length(bundle_sheet))
names(cell_cycle) <- rep("Cyc", length(cell_cycle))
names(epidermis) <- rep("Epi", length(epidermis))
names(trichome) <- rep("Tri", length(trichome))
names(GC) <- rep("GC", length(GC))
names(hydathode) <- rep("Hyd", length(hydathode))
names(mesophyll) <- rep("Mes", length(mesophyll))
names(myrosinase_cell) <- rep("MI", length(myrosinase_cell))
names(phloem) <- rep("Phloem", length(phloem))
names(SE) <- rep("SE", length(SE))
names(xylem) <- rep("Xylem", length(xylem))

ff <- c(
  bundle_sheet, cell_cycle, epidermis,
  trichome, GC, hydathode, myrosinase_cell,
  phloem, SE, xylem,
  mesophyll
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
  "4", "9", "0", "8", "14",
  "17", "16", "18", "11", "12", "13",
  "1", "10",
  "15", "6",
  "7", "2", "3", "5"
)

p2$data$id <- factor(p2$data$id, levels = lev)
p2

ggsave("./figures/leaf/LC_LPZ_public_anno_dot.pdf", p2, height = 5, width = 9)
