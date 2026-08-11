###############################################################
# Beluga Project — Figure 3
# Fig3_PanelA: Heatmap of functional-category abundance (row-scaled),
#              ordered by Plastic:Water dispersion ratio, with Habitat
#              (top strip) and Significant/p_fdr<0.05 (left strip)
#              annotations -- built in native ggplot2.
# Fig3_PanelB: Paired station-level plot, aerobic corrin ring synthesis,
#              showing Plastic > Water at all 7 paired stations
###############################################################

## -----------------------------
## Required input files (working directory):
## -----------------------------
## 1) func_header_mat.rds  -- from beluga_MAG_functional_analysis.R (Section 6)
## 2) disp_by_header.rds   -- from beluga_MAG_functional_analysis.R (Section 7)
## 3) meta_func.rds        -- from the same script's Section 9 (save outputs)
## 4) corrin_df.rds        -- from the same script's Section 9
###############################################################

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(cowplot)

## -----------------------------
## 0) Load pre-computed objects
## -----------------------------
func_header_mat <- readRDS("func_header_mat.rds")
disp_by_header  <- readRDS("disp_by_header.rds")
meta_func       <- readRDS("meta_func.rds")
corrin_df       <- readRDS("corrin_df.rds")

habitat_labels <- c(ST = "Water", T = "Plastic")
habitat_colors <- c(Water = "#0072B2", Plastic = "#E69F00")
sig_colors     <- c("p_fdr < 0.05" = "grey30", "n.s." = "grey90")

## =============================================================
## Fig3_PanelA -- build three aligned ggplot layers (top strip,
## left strip, main heatmap) and compose them in a 2x2 grid.
## =============================================================

## Row order: Plastic:Water dispersion ratio, descending. This is a
## habitat-specific variability ranking, not an abundance ranking --
## the corrinoid-related categories are mid-to-low in raw and per-KO
## abundance, but stand out because they are differentially dispersed
## and directionally enriched between habitats.
header_order <- disp_by_header %>%
  arrange(desc(ratio_Plastic_Water)) %>%
  pull(header)

## Column order: Water block first, then Plastic block; within each,
## by River_Km_m ascending. Clean Habitat_Station labels.
col_meta <- meta_func %>%
  mutate(
    Habitat  = recode(Habitat, !!!habitat_labels),
    Habitat  = factor(Habitat, levels = c("Water", "Plastic")),
    ColLabel = paste0(Habitat, "_", Station)
  ) %>%
  arrange(Habitat, River_Km_m)

## Row-scaled (z-score) long-format data for the main heatmap
heat_long <- as.data.frame(func_header_mat) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "header", values_to = "abundance") %>%
  left_join(meta_func %>% select(SampleID, Habitat, Station), by = "SampleID") %>%
  mutate(Habitat = recode(Habitat, !!!habitat_labels)) %>%
  left_join(col_meta %>% select(SampleID, ColLabel), by = "SampleID") %>%
  group_by(header) %>%
  mutate(z = as.numeric(scale(abundance))) %>%
  ungroup() %>%
  mutate(
    ColLabel = factor(ColLabel, levels = col_meta$ColLabel),
    header   = factor(header, levels = rev(header_order))
  )

## Row annotation data (Significant / n.s.), aligned to same header levels
row_strip_df <- disp_by_header %>%
  mutate(Significant = ifelse(p_fdr < 0.05, "p_fdr < 0.05", "n.s.")) %>%
  filter(header %in% header_order) %>%
  mutate(header = factor(header, levels = rev(header_order)))

## Top annotation data (Habitat), aligned to same column levels
top_strip_df <- col_meta %>%
  mutate(ColLabel = factor(ColLabel, levels = col_meta$ColLabel))

## -- Main heatmap tile plot --
p_main <- ggplot(heat_long, aes(x = ColLabel, y = header, fill = z)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick3",
                       midpoint = 0, name = "Row z-score\n(abundance)") +
  scale_y_discrete(position = "left") +
  guides(fill = guide_colorbar(barheight = unit(3, "cm"), barwidth = unit(0.5, "cm"))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    panel.grid  = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 5, r = 5, b = 25, l = 5)   # extra bottom clearance for rotated labels
  )

## -- Top strip: Habitat --
p_top <- ggplot(top_strip_df, aes(x = ColLabel, y = 1, fill = Habitat)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = habitat_colors, name = "Habitat") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

## -- Left strip: Significant (p_fdr < 0.05) --
p_left <- ggplot(row_strip_df, aes(x = 1, y = header, fill = Significant)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = sig_colors, name = "Significant") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

## -- Extract legends, then strip them from the individual panels so
## each tile plot in the 2x2 grid is legend-free (legends assembled
## separately below, to avoid duplicated/misaligned legends).
legend_fill <- get_legend(p_main)
legend_hab  <- get_legend(p_top)
legend_sig  <- get_legend(p_left)

p_main_nolegend <- p_main + theme(legend.position = "none")
p_top_nolegend  <- p_top  + theme(legend.position = "none")
p_left_nolegend <- p_left + theme(legend.position = "none")

corner <- ggplot() + theme_void()

## -- 2x2 grid: [corner, top strip] / [left strip, main heatmap] --
## Row1 and Row2 share identical relative column widths so the strips
## align exactly with the heatmap columns/rows beneath/beside them.
grid_2x2 <- plot_grid(
  corner,          p_top_nolegend,
  p_left_nolegend, p_main_nolegend,
  ncol = 2,
  rel_widths  = c(0.06, 1),
  rel_heights = c(0.06, 1),
  align = "hv", axis = "tblr"
)

## -- Arrange the three legends in a single horizontal row, placed below
## the grid (frees horizontal width for the enlarged row labels).
legends_row <- plot_grid(legend_hab, legend_sig, legend_fill, nrow = 1)

Fig3_PanelA <- plot_grid(
  grid_2x2, legends_row,
  ncol = 1,
  rel_heights = c(1, 0.30)
)

## =============================================================
## Fig3_PanelB -- Paired station-level plot, aerobic corrin ring
## synthesis (Plastic > Water at all 7 stations, Wilcoxon W=3,
## p=0.004079)
## =============================================================

corrin_paired <- corrin_df %>%
  mutate(
    Habitat = recode(Habitat, !!!habitat_labels),
    Habitat = factor(Habitat, levels = c("Water", "Plastic"))
  )

Fig3_PanelB <- ggplot(corrin_paired, aes(x = Habitat, y = total_abundance, group = Station)) +
  geom_line(color = "grey60", linewidth = 0.5) +
  geom_point(aes(color = Habitat), size = 3) +
  scale_color_manual(values = habitat_colors) +
  labs(
    x = NULL,
    y = "Aerobic corrin ring synthesis\n(community-weighted abundance)",
    color = NULL
  ) +
  annotate(
    "text", x = 1.5, y = max(corrin_paired$total_abundance) * 1.05,
    label = "Wilcoxon exact test\nW = 3, p = 0.004079",
    size = 4.5, hjust = 0.5
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text  = element_text(size = 13),
    axis.title = element_text(size = 13),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

## =============================================================
## Save individual panels + composed figure
## =============================================================
ggsave("Figure3_PanelA_functional_heatmap.pdf", Fig3_PanelA, width = 8, height = 8, units = "in")
ggsave("Figure3_PanelA_functional_heatmap.png", Fig3_PanelA, width = 8, height = 8, units = "in", dpi = 300)

ggsave("Figure3_PanelB_corrinoid_paired.pdf", Fig3_PanelB, width = 4.5, height = 6, units = "in")
ggsave("Figure3_PanelB_corrinoid_paired.png", Fig3_PanelB, width = 4.5, height = 6, units = "in", dpi = 300)

fig3_composed <- plot_grid(
  Fig3_PanelA, Fig3_PanelB,
  labels = c("A)", "B)"),
  label_size = 14,
  ncol = 2,
  rel_widths = c(1.4, 1)
)

ggsave("Figure3_functional_heterogeneity_corrinoid.pdf", fig3_composed, width = 12.5, height = 8, units = "in")
ggsave("Figure3_functional_heterogeneity_corrinoid.png", fig3_composed, width = 12.5, height = 8, units = "in", dpi = 300)

cat("\nDone. Saved:\n")
cat(" - Figure3_PanelA_functional_heatmap.pdf/.png\n")
cat(" - Figure3_PanelB_corrinoid_paired.pdf/.png\n")
cat(" - Figure3_functional_heterogeneity_corrinoid.pdf/.png (composed)\n")
