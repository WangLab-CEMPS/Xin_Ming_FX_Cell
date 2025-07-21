library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)

# read data from ./dataLib/mRNA-QC.xlsx
pd <- read_excel("./dataLib/mRNA-QC.xlsx")
pd$groups <- rep(c("LS", "LSC", "LSF", "JB", "JBC", "JBF", "OSR", "OSRC", "OSRF"), each = 4)
pd$groups2 <- rep(c("LS", "JB", "OSR"), each = 12)
pd <- pd[c(2:6, 8, 10:12, 14:17, 19:20, 22:24, 26:30, 32:34, 36), ]

write.table(pd, file = "./dataLib/mRNA-QC.tsv", sep = "\t", row.names = FALSE)

# t.test group by groups2 and groups
rin_pvals <- pd %>%
  group_by(groups2) %>%
  rstatix::t_test(
    RIN ~ groups,
    p.adjust.method = "BH",
    var.equal = TRUE
  ) %>%
  rstatix::add_x_position(x = "groups", dodge = 0.9) %>% # dodge must match points
  mutate(label = p.adj.signif)

write.table(rin_pvals, file = "./dataLib/mRNA-QC-t-test.tsv", sep = "\t", row.names = FALSE)

# add sd and mean by groups2 and groups
pd <- pd %>%
  group_by(groups2, groups) %>%
  summarise(sd = sd(RIN), mean = mean(RIN)) %>%
  ungroup()

# Most basic error bar
p <- ggplot(pd) +
  geom_bar(aes(x = groups, y = mean, fill = groups2), stat = "identity", alpha = 0.7) +
  geom_errorbar(
    aes(x = groups, ymin = mean - sd, ymax = mean + sd),
    width = 0.4, colour = "#635e5e", alpha = 0.9, linewidth = 0.6
  ) +
  ggprism::theme_prism(
    base_line_size = 0.2
  ) +
  theme(
    axis.text.x = element_text(size = 11)
  ) +
  scale_fill_manual(values = c("#1dc1ab", "#56acfc", "#fc877f")) +
  scale_y_continuous(expand = c(0, 0.1))

pdf("./figures/mRNA-QC-RIN.pdf", width = 7, height = 4)
p
dev.off()
