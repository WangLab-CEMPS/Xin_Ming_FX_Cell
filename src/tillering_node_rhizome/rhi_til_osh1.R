library(patchwork)
library(ggplot2)
library(ggrastr)
library(Seurat)
library(dplyr)

plt2 <- "Os06g0657500"
plt3 <- "Os02g0614300"
mer <- c("Os04g0653600", "Os03g0764900")
osh1 <- c("Os03g0727000")
osh15 <- c("Os07g0129700")

smc <- c("Bochr03G119130", "Bochr03G119150", "Bochr07G220090")

objrhi <- readRDS("./dataLib/tillering_node_rhizome/rhizome_S0.rds")
objtil <- readRDS("./dataLib/tillering_node_rhizome/tillering.rds")


ptil <- FeaturePlot(objtil, c(plt2, plt3, mer, osh1, osh15), order = TRUE)

ptil_list <- lapply(ptil, function(x) {
  pd <- x$data
  t <- x$labels$title

  ggplot(pd, aes(x = umap_1, y = umap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
    ggtitle(t)
})


pptil <- wrap_plots(ptil_list, ncol = 3) &
  scale_color_gradientn(colours = c("#c8c6c3", "#f6ee6c", "#f9a432", "#eb212f", "#88181d")) &
  theme_linedraw() &
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )


pdf("./figures/tillering_node_rhizome_fig4s4/tillering_mc_marker_umap.pdf", width = 16, height = 11.1)
pptil
dev.off()

prhi <- FeaturePlot(objrhi, smc, order = TRUE)

prhi_list <- lapply(prhi, function(x) {
  pd <- x$data
  t <- x$labels$title

  ggplot(pd, aes(x = umap_1, y = umap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
    ggtitle(t)
})


pprhi <- wrap_plots(prhi_list, ncol = 3) &
  scale_color_gradientn(colours = c("#c8c6c3", "#f6ee6c", "#f9a432", "#eb212f", "#88181d")) &
  theme_linedraw() &
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

pdf("./figures/tillering_node_rhizome_fig4s4/rhizome_smc_marker_umap.pdf", width = 14, height = 5)
pprhi
dev.off()
