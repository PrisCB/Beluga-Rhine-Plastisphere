###############################################
# Beluga Project — 16S (Prokaryotic) ASV Analysis
# Taxonomic community structure, nutrient envfit, alpha/beta diversity:
# Plastisphere (MT/L) vs Water (ST), Rhine transect
###############################################

## -----------------------------
## Required input files (working directory):
## -----------------------------
## 1) asv_table.txt
##      16S ASV count table (Cascabel/DADA2). First column "#OTU_ID",
##      remaining columns are NIOZ349.XXXX.YYYY sample IDs.
##
## 2) NIOZ349_mapping_file_correct_GP_3D.txt
##      Sample metadata for the 16S run. Must contain "#SampleID" and
##      "Description" columns.
##
## 3) representative_seq_set_tax_assignments.txt (16S version)
##      Taxonomy assignment table: Kingdom, Phylum, Class, Order, Family,
##      Genus, Species.
##
## 4) 201903-04_BelugaRhineNutrients.xlsx
##      Station-level nutrient concentrations (NO3NO2, Si, PO4)
###############################################

library(tidyverse)
library(readr)
library(readxl)
library(stringr)
library(tibble)
library(vegan)

## -----------------------------
## 1) Load ASV counts
## -----------------------------
asv_raw <- read_tsv("asv_table.txt", col_names = TRUE, show_col_types = FALSE)
cat("Raw ASV table dim:", dim(asv_raw), "\n")

names(asv_raw)[1] <- "ASV_ID"

X <- asv_raw %>%
  column_to_rownames("ASV_ID") %>%
  as.matrix() %>%
  t() %>%
  as.data.frame()

cat("Transposed count matrix dim (samples x ASVs):", dim(X), "\n")

## -----------------------------
## 2) Load and clean sample metadata
## -----------------------------
meta_asv <- read_tsv("NIOZ349_mapping_file_correct_GP_3D.txt", col_names = TRUE, show_col_types = FALSE)

meta2 <- meta_asv %>%
  mutate(
    SampleID = as.character(`#SampleID`),
    Description = str_trim(as.character(Description)),
    Description = str_replace_all(Description, ",+$", "")
  ) %>%
  filter(!str_detect(Description, regex("^Neg_PCR|^Mock", ignore_case = TRUE))) %>%
  filter(SampleID %in% rownames(X)) %>%
  mutate(
    Station_keepLetters = case_when(
      str_detect(Description, "GP_\\d{4}-\\d{2}[-_](MT|L|ST)") ~ {
        m <- str_match(Description, "(GP_\\d{4}-\\d{2})[-_](MT|L|ST)")
        paste0(m[,2], "_", m[,3])
      },
      str_detect(Description, "^GP0*\\d{2}ST_") ~ {
        m <- str_match(Description, "^GP0*(\\d{2})ST_")
        paste0("GP_2019-", m[,2], "_ST")
      },
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Station_keepLetters))

cat("\nUnique Station_keepLetters values:\n")
print(sort(unique(meta2$Station_keepLetters)))
cat("N unique station-habitat combinations:", length(unique(meta2$Station_keepLetters)), "\n")

## -----------------------------
## 3) Collapse technical replicates/subsamples -> station-habitat level
## -----------------------------
X_sub <- X[meta2$SampleID, , drop = FALSE]
X_station <- rowsum(as.matrix(X_sub), group = meta2$Station_keepLetters)
cat("\nStation-level count matrix dim:", dim(X_station), "\n")

## -----------------------------
## 4) Full NMDS across all 23 station-habitat combinations
## -----------------------------
X_hel <- decostand(X_station, method = "hellinger")
set.seed(999)
nmds <- metaMDS(X_hel, distance = "bray", k = 2, trymax = 200)
cat("\nFull (23-station) NMDS stress:", signif(nmds$stress, 3), "\n")

scores <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column("Station_keepLetters") %>%
  mutate(
    Station = str_extract(Station_keepLetters, "(?<=GP_2019-)\\d{2}") %>% factor(),
    Code    = str_extract(Station_keepLetters, "(?<=_)(MT|L|ST)"),
    Habitat = ifelse(Code == "ST", "Water", "Plastic") %>% factor(levels = c("Water", "Plastic")),
    Label   = paste0(Station, "_", Code)
  )

cat("Full NMDS scores dim:", dim(scores), "\n")

## -----------------------------
## 5) Subset + re-run NMDS on the 7 matched stations (Figure 2 taxonomic panel)
## -----------------------------
target_stations <- c("01", "06", "07", "08", "09", "13", "14")

meta2_matched <- meta2 %>%
  filter(
    str_extract(Station_keepLetters, "(?<=GP_2019-)\\d{2}") %in% target_stations,
    str_extract(Station_keepLetters, "(?<=_)(MT|ST)") %in% c("MT", "ST")
  )

X_station_matched <- X_station[unique(meta2_matched$Station_keepLetters), , drop = FALSE]
cat("\nMatched-station count matrix dim:", dim(X_station_matched), "\n")

X_hel_matched <- decostand(X_station_matched, method = "hellinger")
set.seed(999)
nmds_taxa <- metaMDS(X_hel_matched, distance = "bray", k = 2, trymax = 200)

scores_taxa <- as.data.frame(scores(nmds_taxa, display = "sites")) %>%
  rownames_to_column("Station_keepLetters") %>%
  mutate(
    Station  = str_extract(Station_keepLetters, "(?<=GP_2019-)\\d{2}") %>% factor(),
    Code     = str_extract(Station_keepLetters, "(?<=_)(MT|ST)"),
    Habitat  = ifelse(Code == "ST", "Water", "Plastic") %>% factor(levels = c("Water", "Plastic")),
    SampleID = Station_keepLetters
  )

cat("Matched-station NMDS scores dim:", dim(scores_taxa), "\n")
cat("\nStation x Habitat table (matched subset):\n")
print(table(scores_taxa$Station, scores_taxa$Habitat))
cat("\nMatched-station (taxonomic) NMDS stress:", signif(nmds_taxa$stress, 3), "\n")

## -----------------------------
## 6) Load and clean nutrient data
## -----------------------------
nut <- read_excel("201903-04_BelugaRhineNutrients.xlsx")

nut_clean <- nut %>%
  rename(
    SampleID = `Sample ID`,
    PO4      = `PO4 (µmol/L)`,
    Si       = `Si (µmol/L)`,
    NO3NO2   = `NO3+NO2 (µmol/L)`
  ) %>%
  mutate(
    across(c(NO3NO2, Si, PO4), ~ suppressWarnings(as.numeric(.x))),
    Station = str_pad(as.character(as.integer(
      str_extract(SampleID, "(?<=-)(\\d+)$"))), 2, pad = "0")
  ) %>%
  group_by(Station) %>%
  summarise(
    NO3NO2 = mean(NO3NO2, na.rm = TRUE),
    Si     = mean(Si,     na.rm = TRUE),
    PO4    = mean(PO4,    na.rm = TRUE),
    .groups = "drop"
  )

## -----------------------------
## 7) Match nutrients to the matched-subset NMDS scores and run envfit
## -----------------------------
scores_env <- scores_taxa %>%
  left_join(nut_clean, by = "Station")

stopifnot(sum(is.na(scores_env$NO3NO2)) == 0)

# --- envfit: plastisphere (MT) only ---
mt_idx <- scores_env$Code == "MT"
nmds_mt_scores <- scores(nmds_taxa, display = "sites")[mt_idx, ]
env_mt <- scores_env[mt_idx, c("NO3NO2", "Si", "PO4")]

set.seed(999)
ef_mt <- envfit(nmds_mt_scores, env_mt, permutations = 9999, na.rm = TRUE)
cat("\n=== PLASTISPHERE (MT) envfit ===\n")
print(ef_mt)

# --- envfit: water (ST) only ---
st_idx <- scores_env$Code == "ST"
nmds_st_scores <- scores(nmds_taxa, display = "sites")[st_idx, ]
env_st <- scores_env[st_idx, c("NO3NO2", "Si", "PO4")]

set.seed(999)
ef_st <- envfit(nmds_st_scores, env_st, permutations = 9999, na.rm = TRUE)
cat("\n=== WATER (ST) envfit ===\n")
print(ef_st)

## -----------------------------
## 8) Bootstrap (5,000 resamples) comparing envfit R² between habitats
## -----------------------------
set.seed(999)
n_boot  <- 5000
target  <- target_stations
r2_boot <- matrix(NA, nrow = n_boot, ncol = 2, dimnames = list(NULL, c("MT", "ST")))

for (i in seq_len(n_boot)) {
  stn_i <- sample(target, replace = TRUE)

  boot_mt <- scores_env[scores_env$Code == "MT" & scores_env$Station %in% stn_i, ]
  boot_st <- scores_env[scores_env$Code == "ST" & scores_env$Station %in% stn_i, ]

  if (length(unique(boot_mt$Station)) < 3) next
  if (length(unique(boot_st$Station)) < 3) next

  sc_mt <- scores(nmds_taxa, display = "sites")[
    rownames(scores(nmds_taxa, display = "sites")) %in% boot_mt$Station_keepLetters, , drop = FALSE]
  sc_st <- scores(nmds_taxa, display = "sites")[
    rownames(scores(nmds_taxa, display = "sites")) %in% boot_st$Station_keepLetters, , drop = FALSE]

  env_mt_i <- boot_mt[match(rownames(sc_mt), boot_mt$Station_keepLetters), c("NO3NO2", "Si", "PO4")]
  env_st_i <- boot_st[match(rownames(sc_st), boot_st$Station_keepLetters), c("NO3NO2", "Si", "PO4")]

  ef_mt_i <- tryCatch(envfit(sc_mt, env_mt_i, permutations = 0, na.rm = TRUE), error = function(e) NULL)
  ef_st_i <- tryCatch(envfit(sc_st, env_st_i, permutations = 0, na.rm = TRUE), error = function(e) NULL)

  if (!is.null(ef_mt_i) && !is.null(ef_st_i)) {
    r2_boot[i, "MT"] <- ef_mt_i$vectors$r["NO3NO2"]
    r2_boot[i, "ST"] <- ef_st_i$vectors$r["NO3NO2"]
  }
}

r2_boot <- na.omit(r2_boot)
cat("\n=== Bootstrap results (NO3NO2 R²) ===\n")
cat("MT median R²:", round(median(r2_boot[, "MT"]), 3), "\n")
cat("ST median R²:", round(median(r2_boot[, "ST"]), 3), "\n")
cat("ST > MT in", round(mean(r2_boot[, "ST"] > r2_boot[, "MT"]) * 100, 1), "% of resamples\n")

## -----------------------------
## 9) Habitat distinctiveness, alpha/beta diversity, phylum composition
## -----------------------------
matched_stations <- scores %>%
  group_by(Station) %>%
  summarise(habitats = paste(sort(unique(Habitat)), collapse = "+"), .groups = "drop") %>%
  filter(habitats == "Water+Plastic") %>%
  pull(Station) %>%
  as.character()

cat("\nStations with both Water and Plastic:", matched_stations, "\n")
cat("N matched stations:", length(matched_stations), "\n")

matched_ids <- scores %>%
  filter(Station %in% matched_stations) %>%
  pull(Station_keepLetters)

X_matched <- X_station[matched_ids, ]

# --- Alpha diversity ---
alpha_div <- data.frame(
  Station_keepLetters = rownames(X_matched),
  Shannon  = diversity(X_matched, index = "shannon"),
  Richness = specnumber(X_matched)
) %>%
  left_join(scores[, c("Station_keepLetters", "Habitat", "Station")], by = "Station_keepLetters")

cat("\n=== Alpha diversity by habitat ===\n")
print(alpha_div %>%
        group_by(Habitat) %>%
        summarise(
          Shannon_mean  = round(mean(Shannon), 2),
          Shannon_sd    = round(sd(Shannon), 2),
          Richness_mean = round(mean(Richness), 0),
          Richness_sd   = round(sd(Richness), 0),
          .groups = "drop"
        ))

print(wilcox.test(Shannon ~ Habitat, data = alpha_div))
print(wilcox.test(Richness ~ Habitat, data = alpha_div))

# --- Beta diversity dispersion (taxonomic ASV) ---
dist_taxa <- vegdist(X_hel_matched, method = "bray")
bd_taxa <- betadisper(dist_taxa, group = scores_taxa$Habitat)
cat("\n=== Betadisper (taxonomic) ANOVA ===\n")
print(anova(bd_taxa))

# --- PERMANOVA: formal test of habitat separation ---
cat("\n=== PERMANOVA (taxonomic, habitat) ===\n")
print(adonis2(dist_taxa ~ Habitat, data = scores_taxa, permutations = 9999))

## -----------------------------
## 10) Load taxonomy and compute top phyla by habitat
## -----------------------------
tax_raw <- read_tsv("representative_seq_set_tax_assignments.txt", col_names = TRUE, show_col_types = FALSE)
cat("\nTaxonomy table dim:", dim(tax_raw), "\n")

names(tax_raw)[1] <- "ASV_ID"
tax <- tax_raw

asv_long <- as.data.frame(t(X_station_matched)) %>%
  rownames_to_column("ASV_ID") %>%
  left_join(tax[, c("ASV_ID", "Phylum")], by = "ASV_ID") %>%
  pivot_longer(-c(ASV_ID, Phylum), names_to = "Station_keepLetters", values_to = "count") %>%
  left_join(scores_taxa[, c("Station_keepLetters", "Habitat")], by = "Station_keepLetters")

phylum_summary <- asv_long %>%
  group_by(Habitat, Phylum) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  group_by(Habitat) %>%
  mutate(rel_abund = total / sum(total) * 100) %>%
  arrange(Habitat, desc(rel_abund)) %>%
  slice_head(n = 5)

cat("\n=== Top 5 phyla by habitat ===\n")
print(phylum_summary)

## -----------------------------
## 11) Save key objects for downstream use / figure-building
## -----------------------------
saveRDS(scores, "scores_16S_full.rds")
saveRDS(scores_taxa, "scores_16S_matched.rds")
saveRDS(nmds_taxa, "nmds_taxa_16S.rds")
saveRDS(alpha_div, "alpha_div_16S.rds")
saveRDS(phylum_summary, "phylum_summary_16S.rds")
saveRDS(r2_boot, "r2_boot_16S.rds")

cat("\nDone. Key objects saved: scores_16S_full.rds, scores_16S_matched.rds, nmds_taxa_16S.rds, alpha_div_16S.rds, phylum_summary_16S.rds, r2_boot_16S.rds\n")

combos <- sort(unique(meta2$Station_keepLetters))

station_habitat_list <- data.frame(
  n = seq_along(combos),
  Station_Habitat = combos
)

print(station_habitat_list, row.names = FALSE)
cat("\nTotal unique station-habitat combinations:", length(combos), "\n")

## -----------------------------
## Station-to-station spatial heterogeneity in dissolved nutrients across
## all 14 stations
## -----------------------------
nutrient_range <- nut_clean %>%
  summarise(across(
    c(NO3NO2, Si, PO4),
    list(
      mean     = ~mean(.x, na.rm = TRUE),
      sd       = ~sd(.x, na.rm = TRUE),
      min      = ~min(.x, na.rm = TRUE),
      max      = ~max(.x, na.rm = TRUE),
      cv_pct   = ~sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE) * 100,
      fold_change = ~max(.x, na.rm = TRUE) / min(.x, na.rm = TRUE)
    )
  ))

print(t(nutrient_range))
