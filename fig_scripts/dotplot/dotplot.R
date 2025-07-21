library(Seurat)
library(ggplot2)

# dotplot for leaf ---------------------
obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public.Rds")

bundle_sheet <- c("AT5G41920", "AT1G77990") # "AT2G37860"
cell_cycle <- c("AT1G18370", "AT3G25980", "AT5G25090")
epidermis <- c("AT4G21750", "AT1G51500", "AT1G68530") # "AT1G27950"
GC <- c("AT2G46070", "AT1G08810", "AT1G22690") # "AT3G18040"
hydathode <- c("AT3G54420", "AT3G16670", "AT1G62510") # "AT1G28230",
mesophyll <- c("AT1G29910", "AT3G01500", "AT4G12970")
myrosinase_cell <- c("AT5G26000", "AT1G52342")
phloem <- c("AT3G48740", "AT5G23660", "AT3G11930") # "AT4G13770"
CC <- c("AT1G22710", "AT1G79430", "AT5G57350") # "AT5G06850"
xylem <- c("AT1G16410", "AT5G61480") # "AT5G12140"

names(bundle_sheet) <- rep("BS", length(bundle_sheet))
names(cell_cycle) <- rep("Cyc", length(cell_cycle))
names(epidermis) <- rep("Epi", length(epidermis))
names(GC) <- rep("GC", length(GC))
names(hydathode) <- rep("Hyd", length(hydathode))
names(mesophyll) <- rep("Mes", length(mesophyll))
names(myrosinase_cell) <- rep("MI", length(myrosinase_cell))
names(phloem) <- rep("Phloem", length(phloem))
names(CC) <- rep("Companion", length(CC))
names(xylem) <- rep("Xylem", length(xylem))

ff <- c(
  bundle_sheet, cell_cycle, epidermis,
  GC, hydathode, myrosinase_cell,
  phloem, CC, xylem,
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

ggsave("./figures/leaf/LC_LPZ_public_anno_dot.pdf", p2, height = 4.6, width = 8)

# dotplot for tillering ---------------------
obj <- readRDS("./dataLib/tillering_node_rhizome/tillering.rds")

cyc <- c("Os03g0279200", "Os08g0490900", "Os05g0438700")
ep <- c("Os02g0776900", "Os02g0813600", "Os06g0653000")
xylem <- c("Os10g0467800", "Os09g0422500", "Os10g0532000", "Os02g0237100")
phloem <- c("Os06g0717200", "Os01g0971000", "Os10g0543800", "Os04g0105200")
co <- c("Os02g0662000", "Os01g0914300")
en <- c("Os12g0122000", "Os02g0745100")
ex <- c("Os04g0452700", "Os03g0107300")
pe <- c("Os03g0216700", "Os10g0524300", "Os02g0809800", "Os04g0556000")
rc <- c("Os04g0552000", "Os01g0248000")
cam <- c("Os04g0649400")
mer <- c("Os04g0653600", "Os06g0657500")

names(cyc) <- rep("Dc", length(cyc))
names(ep) <- rep("Ep", length(ep))
names(xylem) <- rep("Xy", length(xylem))
names(phloem) <- rep("Ph", length(phloem))
names(en) <- rep("En/Co", length(en))
names(co) <- rep("En/Co", length(co))
names(ex) <- rep("Ex", length(ex))
names(pe) <- rep("Pe", length(pe))
names(rc) <- rep("Rc", length(rc))
names(cam) <- rep("Ca", length(cam))
names(mer) <- rep("Me", length(mer))

# dotplot
ff <- c(cyc, mer, xylem, cam, phloem, pe, ep, en, co, ex, rc)
p <- DotPlot(obj, features = ff, cols = c("#eeeeee", "#ff4800"), col.min = 0, col.max = 3) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "cm"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 80, vjust = 0.5, hjust = 0.5, size = 10)
  )

lev <- c(
  "2", "0", "10", "7", "13", "8", "9",
  "4", "5", "1", "3", "12", "11", "6"
)

p$data$id <- factor(p$data$id, levels = lev)
p

ggsave("./figures/tillering_node_rhizome/tillering_dot.pdf", p, width = 10, height = 5)

# dotplot for rhizome ---------------------
obj <- readRDS("./dataLib/tillering_node_rhizome/rhizome_S0.rds")

cellcycle <- c("Bochr03G100050", "Bochr08G266560", "Bochr05G179090")
xylem <- c("Bochr10G305000", "Bochr09G282600")
protoxylem <- c("Bochr08G263530")
TrachearyElement <- c("Bochr01G048000", "Bochr05G161730")
phloem_cc <- c("Bochr06G210850", "Bochr06G214280", "Bochr01G047960")
phloem_pole_pericycle <- c("Bochr01G012110")
cortex <- c("Bochr01G032380", "Bochr02G076060")
endodermis <- c("Bochr02G082030", "Bochr04G156440")
exodermis <- c("Bochr04G144330", "Bochr04G130210")
pericycle <- c("Bochr09G286170", "Bochr01G013200")
cambium <- c("Bochr07G233300")
epidermis <- c("Bochr03G117850", "Bochr01G001790")
QC <- c("Bochr01G042780", "Bochr02G072640")
rc <- c("Bochr06G206520")
smc <- c("Bochr03G119130", "Bochr07G220090")


cellcycle <- cellcycle
xylem <- c(xylem, protoxylem)
TrachearyElement <- TrachearyElement
phloem_cc <- c(phloem_cc, phloem_pole_pericycle)
cortex <- cortex
en_ex_pe <- c(endodermis, exodermis, pericycle)
cambium <- cambium
lrc <- c(lateral_root_primordia)
qc_like_root_mer <- c(QC, rc)
smc <- smc

names(cellcycle) <- rep("Dc", length(cellcycle))
names(xylem) <- rep("Xy", length(xylem))
names(TrachearyElement) <- rep("Te", length(TrachearyElement))
names(phloem_cc) <- rep("Ph&Cc", length(phloem_cc))
names(cortex) <- rep("Co", length(cortex))
names(en_ex_pe) <- rep("En/Ex/Pe", length(en_ex_pe))
names(cambium) <- rep("Ca", length(cambium))
names(epidermis) <- rep("Epi", length(epidermis))
names(qc_like_root_mer) <- rep("QC", length(qc_like_root_mer))
names(smc) <- rep("SMC", length(smc))

# dotplot
ff <- c(
  cellcycle, xylem, TrachearyElement,
  phloem_cc, cortex, en_ex_pe, cambium,
  epidermis, qc_like_root_mer, smc
)

p <- DotPlot(obj, features = ff, cols = c("#eeeeee", "#ff4800"), col.min = 0, col.max = 3) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "cm"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 80, vjust = 0.5, hjust = 0.5, size = 10)
  )

lev <- c(
  "14", "18", "13", "11", "23", "5", "20", "16", "21", "19",
  "9", "15", "17", "8", "6", "4", "7", "22", "1", "2", "0", "10", "3", "12"
)

p$data$id <- factor(p$data$id, levels = lev)

ggsave("./figures/tillering_node_rhizome_fig4s4/rhizome_S0_dot.pdf", p, width = 8, height = 5)
