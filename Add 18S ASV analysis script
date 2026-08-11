###############################################
# Beluga Eukaryotic (18S) Analysis
# Diatom / ciliate / fungi / Dreissena enrichment + alpha diversity
# Plastic vs Water
###############################################

## -----------------------------
## Required input files (working directory):
## -----------------------------
## 1) euk_asv_table.txt
##      18S ASV count table (Cascabel, 951R-based "eukaryotic" primer set).
##      12,635 ASVs x 34 samples expected.
##
## 2) euk_representative_seq_set_tax_assignments.txt
##      Taxonomy assignment table, 9 ranks: Domain, Supergroup, Division,
##      Subdivision, Class, Order, Family, Genus, Species.
##
## 3) NIOZ349_mapping_file_corrected_GP_EUK.txt
##      Sample metadata for the 18S ("eukaryotic") run. Must contain
##      "#SampleID" and "Description" columns.
##
## -----------------------------
## Notes on data completeness:
## -----------------------------
## - Water samples at stations 13 and 14 have very low/zero read depth in
##   this 18S dataset:
##     * Station 13 (NIOZ349.9193.9094): zero reads, absent from the
##       counts table entirely.
##     * Station 14 (NIOZ349.9195.9096): 1 raw read total, present but
##       removed by the depth<100 filter in Section 3.
##   The equivalent 16S water samples at these two stations are present
##   with strong read depth, so this is a marker-specific outcome rather
##   than a missing sample.
###############################################

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(tibble)

## -----------------------------
## 0) Helper: fix "V1, V2, V3..." generic column names
## -----------------------------
## If a file was imported with the header row not recognized as column
## names, R auto-labels columns V1, V2, V3... and the true header ends up
## sitting in row 1 as data. This detects that pattern, promotes row 1 to
## column names, and drops it from the data.
fix_generic_header <- function(df) {
  if (all(grepl("^V[0-9]+$", names(df)))) {
    new_names <- as.character(unlist(df[1, ], use.names = FALSE))
    df <- df[-1, ]
    names(df) <- new_names
  }
  df
}

## -----------------------------
## 1) Load counts and taxonomy
## -----------------------------
counts_raw <- read_tsv("euk_asv_table.txt", col_names = TRUE, show_col_types = FALSE)
counts_raw <- fix_generic_header(counts_raw)
names(counts_raw)[1] <- "ASV_ID"

tax_raw <- read_tsv("euk_representative_seq_set_tax_assignments.txt", col_names = TRUE, show_col_types = FALSE)
tax_raw <- fix_generic_header(tax_raw)
names(tax_raw)[1] <- "ASV_ID"

## Numeric columns come back as character after fix_generic_header() promotes
## row 1 -> convert count columns to numeric
count_cols <- setdiff(names(counts_raw), "ASV_ID")
counts_raw <- counts_raw %>%
  mutate(across(all_of(count_cols), as.numeric))

cat("Counts dim:", dim(counts_raw), "\n")
cat("Tax dim:", dim(tax_raw), "\n")

## -----------------------------
## 2) Build sample metadata (Station/Habitat) from NIOZ349 mapping file
## -----------------------------
mapping_raw <- read_tsv("NIOZ349_mapping_file_corrected_GP_EUK.txt", col_names = TRUE, show_col_types = FALSE)
mapping_raw <- fix_generic_header(mapping_raw)
names(mapping_raw)[1] <- "SampleID"
names(mapping_raw)[ncol(mapping_raw)] <- "Description"

sample_ids <- names(counts_raw)[-1]  # drop ASV_ID column

meta <- mapping_raw %>%
  filter(SampleID %in% sample_ids) %>%
  select(SampleID, Description) %>%
  mutate(
    Habitat = case_when(
      str_detect(Description, "^Neg_PCR|^Mock")   ~ "Control",
      str_detect(Description, "MT|_L-")           ~ "Plastic",
      str_detect(Description, "_ST_")             ~ "Water",
      TRUE                                        ~ NA_character_
    ),
    Station = case_when(
      Habitat == "Plastic" ~ str_extract(Description, "(?<=GP_2019-)\\d{2}"),
      Habitat == "Water"   ~ str_extract(Description, "(?<=ST_)\\d{2}(?=_ST_a)"),
      TRUE ~ NA_character_
    )
  )

cat("\nHabitat counts:\n")
print(table(meta$Habitat, useNA = "ifany"))

## -----------------------------
## 3) Exclude control sample + failed/near-zero-depth sample(s)
## -----------------------------
counts_mat <- counts_raw %>%
  column_to_rownames("ASV_ID") %>%
  as.matrix()

depth <- colSums(counts_mat)
cat("\nRead depth per sample (sorted):\n")
print(sort(depth))

failed_samples <- names(depth)[depth < 100]
cat("\nFailed/near-zero-depth samples (excluded):", failed_samples, "\n")

exclude_samples <- unique(c(
  meta$SampleID[meta$Habitat == "Control"],
  failed_samples
))
cat("All excluded samples:", exclude_samples, "\n")

keep_samples <- setdiff(sample_ids, exclude_samples)
counts_mat <- counts_mat[, keep_samples]
meta <- meta %>% filter(SampleID %in% keep_samples) %>% arrange(match(SampleID, keep_samples))
stopifnot(identical(colnames(counts_mat), meta$SampleID))

## -----------------------------
## 4) Filter to true Eukaryota (drop Bacteria contamination + nucleomorph/plastid reads)
## -----------------------------
## Use %in% everywhere below, never == , when the column can contain NA.
## `NA == "x"` returns NA (not FALSE), and subsetting a matrix with a logical
## index containing NA inserts phantom all-NA rows that silently poison any
## downstream colSums()/rowSums(). `%in%` is NA-safe (always TRUE/FALSE).
euk_asv_ids <- tax_raw$ASV_ID[tax_raw$Domain %in% "Eukaryota"]
counts_euk <- counts_mat[rownames(counts_mat) %in% euk_asv_ids, , drop = FALSE]

cat("\nASVs retained (true Eukaryota):", nrow(counts_euk), "of", nrow(counts_mat), "\n")

tax_euk <- tax_raw %>% filter(ASV_ID %in% rownames(counts_euk)) %>%
  column_to_rownames("ASV_ID")
tax_euk <- tax_euk[rownames(counts_euk), ]  # align order

## -----------------------------
## 5) Relative abundance (proportion of eukaryotic reads per sample)
## -----------------------------
rel <- sweep(counts_euk, 2, colSums(counts_euk), FUN = "/")

division_abund <- rowsum(rel, group = tax_euk$Division)
division_mean_by_habitat <- t(apply(division_abund, 1, function(x) tapply(x, meta$Habitat, mean)))
cat("\n=== Division-level mean relative abundance (Plastic vs Water) ===\n")
print(round(division_mean_by_habitat[order(-division_mean_by_habitat[, "Plastic"]), ], 4))

## -----------------------------
## 6) Diatom classes specifically: Bacillariophyceae, Coscinodiscophyceae, Mediophyceae
## -----------------------------
gyrista_idx <- tax_euk$Subdivision %in% "Gyrista"
class_abund <- rowsum(rel[gyrista_idx, , drop = FALSE], group = tax_euk$Class[gyrista_idx])
class_mean_by_habitat <- t(apply(class_abund, 1, function(x) tapply(x, meta$Habitat, mean)))
cat("\n=== Class-level (Gyrista) mean relative abundance ===\n")
print(round(class_mean_by_habitat, 4))

## Total "diatoms" (combined 3 classes) -- shown for comparison only.
## This combined metric is confounded by planktonic Stephanodiscus
## dominance in water and is not the recommended enrichment metric (see
## pennate-specific test below).
diatom_classes <- c("Bacillariophyceae", "Coscinodiscophyceae", "Mediophyceae")
diatom_idx <- tax_euk$Class %in% diatom_classes
diatom_rel_per_sample <- colSums(rel[diatom_idx, , drop = FALSE])
cat("\n=== Combined diatom relative abundance per sample (for reference only) ===\n")
print(tapply(diatom_rel_per_sample, meta$Habitat, function(x) c(mean = mean(x), sd = sd(x), n = length(x))))

## -----------------------------
## 7) Pennate diatoms (Bacillariophyceae) -- the ecologically meaningful,
##    biofilm-forming group -- genus breakdown + enrichment test
## -----------------------------
bacillario_idx <- tax_euk$Class %in% "Bacillariophyceae"
genus_abund_bacillario <- rowsum(rel[bacillario_idx, , drop = FALSE], group = tax_euk$Genus[bacillario_idx])
genus_mean_bacillario <- t(apply(genus_abund_bacillario, 1, function(x) tapply(x, meta$Habitat, mean)))

cat("\n=== Genus breakdown within Bacillariophyceae (pennate diatoms) ===\n")
print(round(genus_mean_bacillario[order(-genus_mean_bacillario[, "Plastic"]), ], 4))

pennate_rel_per_sample <- colSums(rel[bacillario_idx, , drop = FALSE])

## Centric diatoms (Mediophyceae) -- planktonic/water-column specialists
medio_idx <- tax_euk$Class %in% "Mediophyceae"
genus_abund_medio <- rowsum(rel[medio_idx, , drop = FALSE], group = tax_euk$Genus[medio_idx])
genus_mean_medio <- t(apply(genus_abund_medio, 1, function(x) tapply(x, meta$Habitat, mean)))
cat("\n=== Genus breakdown within Mediophyceae (centric diatoms) ===\n")
print(round(genus_mean_medio[order(-genus_mean_medio[, "Water"]), ], 4))

## -----------------------------
## 8) Matched-station paired comparison (Bacillariophyceae / pennate diatoms)
## -----------------------------
plastic_stations <- unique(meta$Station[meta$Habitat == "Plastic"])
water_stations   <- unique(meta$Station[meta$Habitat == "Water"])
matched_stations_euk <- sort(intersect(plastic_stations, water_stations))
cat("\nMatched stations (eukaryotic, both habitats with usable data):", matched_stations_euk, "\n")

station_habitat_df <- data.frame(
  SampleID = meta$SampleID,
  Station  = meta$Station,
  Habitat  = meta$Habitat,
  Pennate  = as.numeric(pennate_rel_per_sample),  # positional -- rel/meta share row order
  stringsAsFactors = FALSE
)

station_means <- station_habitat_df %>%
  filter(Station %in% matched_stations_euk) %>%
  group_by(Station, Habitat) %>%
  summarise(mean_pennate = mean(Pennate, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Habitat, values_from = mean_pennate)

cat("\n=== Station-level mean pennate-diatom (Bacillariophyceae) relative abundance ===\n")
print(station_means)

wilcox_paired <- wilcox.test(station_means$Plastic, station_means$Water, paired = TRUE)
cat("\n=== Paired Wilcoxon signed-rank test (Plastic vs Water, matched stations, pennate diatoms) ===\n")
print(wilcox_paired)

wilcox_unpaired <- wilcox.test(
  pennate_rel_per_sample[meta$Habitat == "Plastic"],
  pennate_rel_per_sample[meta$Habitat == "Water"]
)
cat("\n=== Unpaired Mann-Whitney (all samples, pennate diatoms) ===\n")
print(wilcox_unpaired)

## -----------------------------
## 9) Ciliates: Peritrichia / Sessilida (Epistylis, Zoothamnium, etc.)
## -----------------------------
cil_idx <- tax_euk$Subdivision %in% "Ciliophora"
order_abund_cil <- rowsum(rel[cil_idx, , drop = FALSE], group = tax_euk$Order[cil_idx])
order_mean_cil <- t(apply(order_abund_cil, 1, function(x) tapply(x, meta$Habitat, mean)))
cat("\n=== Order-level breakdown within Ciliophora ===\n")
print(round(order_mean_cil[order(-order_mean_cil[, "Plastic"]), ], 4))

## Genus-level within the dominant sessile/attached order. In this taxonomy
## string, "Peritrichia_2" is the Order-level label for sessile peritrichs
## ("Sessilida" sits at Family rank instead, kept here as a fallback in
## case the local taxonomy table differs).
peritrichia_idx <- tax_euk$Order %in% c("Peritrichia_2", "Sessilida")
if (sum(peritrichia_idx) > 0) {
  genus_abund_perit <- rowsum(rel[peritrichia_idx, , drop = FALSE], group = tax_euk$Genus[peritrichia_idx])
  genus_mean_perit <- t(apply(genus_abund_perit, 1, function(x) tapply(x, meta$Habitat, mean)))
  cat("\n=== Genus breakdown within sessile peritrich ciliates ===\n")
  print(round(genus_mean_perit[order(-genus_mean_perit[, "Plastic"]), ], 4))
  
  perit_rel_per_sample <- colSums(rel[peritrichia_idx, , drop = FALSE])
  cat("\nSessile peritrich ciliates -- Mann-Whitney test:\n")
  print(wilcox.test(perit_rel_per_sample ~ meta$Habitat))
} else {
  cat("\nNo ASVs matched Peritrichia_2/Sessilida -- check unique(tax_euk$Order):\n")
  print(unique(tax_euk$Order))
}

## -----------------------------
## 10) Dreissena rostriformis (zebra mussel)
## -----------------------------
dre_idx <- tax_euk$Species %in% "Dreissena_rostriformis"
cat("\nN ASVs assigned to Dreissena rostriformis:", sum(dre_idx), "\n")
if (sum(dre_idx) > 0) {
  dre_rel_per_sample <- colSums(rel[dre_idx, , drop = FALSE])
  cat("Dreissena rostriformis relative abundance by habitat:\n")
  print(tapply(dre_rel_per_sample, meta$Habitat, function(x) c(mean = mean(x), max = max(x), n = length(x))))
  
  ## Means can look near-identical while the distributions differ sharply --
  ## Mann-Whitney tests the whole distribution, not just the mean. Always
  ## check median + zero-detection rate alongside the mean/SD.
  cat("\nDreissena median by habitat:\n")
  print(tapply(dre_rel_per_sample, meta$Habitat, median))
  cat("N zero-detection samples by habitat:\n")
  print(tapply(dre_rel_per_sample, meta$Habitat, function(x) sum(x == 0)))
  print(wilcox.test(dre_rel_per_sample ~ meta$Habitat))
  
  ## Station-level breakdown -- signal is concentrated at station 13
  ## (Strasbourg, major inland Rhine port) with a secondary elevation at
  ## station 07; where paired water samples exist (07, 08), Water > Plastic,
  ## suggesting local population density rather than substrate preference.
  ## No water sample exists at station 13 for comparison.
  dre_station_df <- data.frame(
    SampleID = meta$SampleID,
    Station  = meta$Station,
    Habitat  = meta$Habitat,
    Dreissena = round(as.numeric(dre_rel_per_sample), 4),
    stringsAsFactors = FALSE
  )
  cat("\n=== Dreissena rostriformis relative abundance by Station and Habitat ===\n")
  print(dre_station_df %>% arrange(Habitat, Station))
}

## -----------------------------
## 11a) Full community screen -- rank all classes by combined abundance,
##      done before drafting any composition text, rather than limiting
##      the screen to taxa named in advance.
## -----------------------------
class_abund_all <- rowsum(rel, group = tax_euk$Class)
class_mean_all <- t(apply(class_abund_all, 1, function(x) tapply(x, meta$Habitat, mean)))
class_sd_all <- t(apply(class_abund_all, 1, function(x) tapply(x, meta$Habitat, sd)))

## class_mean_all[, "Plastic"] is a named vector (inherits class names,
## including a literal NA, from the matrix rownames). Passing a named
## vector into data.frame() makes R try to use those names as row.names --
## and row.names cannot contain NA, which crashes with "row names contain
## missing values". Fix: unname() every vector going into the data.frame.
class_summary <- data.frame(
  Class        = rownames(class_mean_all),
  Plastic_mean = unname(class_mean_all[, "Plastic"]),
  Plastic_sd   = unname(class_sd_all[, "Plastic"]),
  Water_mean   = unname(class_mean_all[, "Water"]),
  Water_sd     = unname(class_sd_all[, "Water"]),
  row.names    = NULL,
  stringsAsFactors = FALSE
)
class_summary <- class_summary[order(-class_summary$Plastic_mean), ]
cat("\n=== FULL class-level ranking (top 20 by Plastic mean) ===\n")
print(round(head(class_summary, 20)[, -1], 4))
cat("\nSum of Plastic means across ALL classes (sanity check, should be ~1):", sum(class_mean_all[, "Plastic"], na.rm = TRUE), "\n")
cat("Sum of Water means across ALL classes (sanity check, should be ~1):", sum(class_mean_all[, "Water"], na.rm = TRUE), "\n")

## Flag: any class where SD >= mean in either habitat is likely
## outlier-driven (patchy detection in a few samples) rather than a real,
## consistent enrichment signal -- report for completeness but do not
## claim enrichment for these.
class_summary$flag_outlier_driven <- with(class_summary,
                                          (Plastic_sd >= Plastic_mean & Plastic_mean > 0) | (Water_sd >= Water_mean & Water_mean > 0)
)
cat("\n=== Classes flagged as likely outlier-driven (SD >= mean) ===\n")
print(class_summary[class_summary$flag_outlier_driven, c("Class", "Plastic_mean", "Plastic_sd", "Water_mean", "Water_sd")])

## -----------------------------
## 11a-ii) Wilcoxon tests for every remaining group reported in the
##         supplementary table.
##
##         Note on exact vs asymptotic p-values: wilcox.test() computes the
##         exact permutation p-value only when there are no ties in the
##         pooled data; with ties present it falls back to the asymptotic
##         normal approximation (with continuity correction) and prints a
##         warning. This matters at small n with complete separation
##         between groups, where exact and asymptotic p-values can differ
##         by orders of magnitude. Always check whether the tie warning
##         appeared to know which value is the correct one to report.
test_group <- function(idx_mask) {
  sub <- colSums(rel[idx_mask, , drop = FALSE])
  wilcox.test(sub ~ meta$Habitat)
}

cat("\n=== Mollusca (incl. Dreissena) -- Wilcoxon ===\n")
print(test_group(tax_euk$Class %in% "Mollusca"))

cat("\n=== Other ciliates (Spirotrichea) -- Wilcoxon ===\n")
print(test_group(tax_euk$Class %in% "Spirotrichea"))

cat("\n=== Chrysophyceae -- Wilcoxon ===\n")
print(test_group(tax_euk$Class %in% "Chrysophyceae"))

cat("\n=== Centric diatoms (Mediophyceae) -- Wilcoxon ===\n")
print(test_group(tax_euk$Class %in% "Mediophyceae"))

cat("\n=== Cryptophyceae -- Wilcoxon ===\n")
print(test_group(tax_euk$Class %in% "Cryptophyceae"))

## -----------------------------
## 11b) Fungi (Ascomycota + Basidiomycota)
## -----------------------------
fungi_idx <- tax_euk$Class %in% c("Ascomycota", "Basidiomycota")
cat("\nN ASVs assigned to Fungi (Ascomycota+Basidiomycota):", sum(fungi_idx), "\n")
if (sum(fungi_idx) > 0) {
  genus_abund_fungi <- rowsum(rel[fungi_idx, , drop = FALSE], group = tax_euk$Genus[fungi_idx])
  genus_mean_fungi <- t(apply(genus_abund_fungi, 1, function(x) tapply(x, meta$Habitat, mean)))
  cat("\n=== Genus breakdown within Fungi (expect Trechispora on top) ===\n")
  print(round(genus_mean_fungi[order(-genus_mean_fungi[, "Plastic"]), ], 4))
  
  fungi_rel_per_sample <- colSums(rel[fungi_idx, , drop = FALSE])
  cat("\nFungi relative abundance by habitat:\n")
  print(tapply(fungi_rel_per_sample, meta$Habitat, function(x) c(mean = mean(x), sd = sd(x), n = length(x))))
  print(wilcox.test(fungi_rel_per_sample ~ meta$Habitat))
}

## -----------------------------
## 11c) High-mean-but-non-significant groups (Embryophyceae, Annelida,
##      Nematoda) -- included for completeness, but flagged as
##      outlier-driven (SD >= mean), not enrichment.
## -----------------------------
patchy_classes <- c("Embryophyceae", "Annelida", "Nematoda")
for (cl in patchy_classes) {
  idx <- tax_euk$Class %in% cl
  if (sum(idx) > 0) {
    sub_rel <- colSums(rel[idx, , drop = FALSE])
    cat("\n===", cl, "(flagged outlier-driven -- see class_summary$flag_outlier_driven) ===\n")
    print(tapply(sub_rel, meta$Habitat, function(x) c(mean = mean(x), sd = sd(x), n = length(x))))
    print(wilcox.test(sub_rel ~ meta$Habitat))
  }
}

## -----------------------------
## 12) Alpha diversity: richness & Shannon (Eukaryota-filtered), Plastic vs Water
## -----------------------------
richness <- colSums(counts_euk > 0)
shannon_fn <- function(col) {
  p <- col[col > 0] / sum(col)
  -sum(p * log(p))
}
shannon_div <- apply(counts_euk, 2, shannon_fn)

alpha_df <- data.frame(
  SampleID = names(richness),
  Richness = as.numeric(richness),
  Shannon  = as.numeric(shannon_div),
  Habitat  = meta$Habitat[match(names(richness), meta$SampleID)]
)

cat("\n=== Alpha diversity by habitat (raw counts) ===\n")
print(aggregate(cbind(Richness, Shannon) ~ Habitat, data = alpha_df,
                FUN = function(x) c(mean = mean(x), sd = sd(x))))

wilcox_richness <- wilcox.test(Richness ~ Habitat, data = alpha_df)
wilcox_shannon  <- wilcox.test(Shannon  ~ Habitat, data = alpha_df)
cat("\n=== Richness Wilcoxon ===\n"); print(wilcox_richness)
cat("\n=== Shannon Wilcoxon ===\n");  print(wilcox_shannon)

## -----------------------------
## 13) Rarefaction check -- rule out sequencing depth as a confound
## -----------------------------
depth_per_sample <- colSums(counts_euk)
depth_habitat <- meta$Habitat[match(names(depth_per_sample), meta$SampleID)]
cat("\n=== Eukaryota-filtered read depth by habitat ===\n")
print(tapply(depth_per_sample, depth_habitat,
             function(x) c(mean = mean(x), sd = sd(x), min = min(x), max = max(x))))
print(wilcox.test(depth_per_sample ~ depth_habitat))

min_depth <- min(depth_per_sample)
cat("\nRarefying to min depth:", min_depth, "\n")

set.seed(42)

rarefy_one <- function(col, depth, n_iter = 50) {
  counts_int <- as.integer(col)
  asv_pool <- rep(seq_along(counts_int), counts_int)
  rich <- numeric(n_iter); shan <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    samp <- sample(asv_pool, size = depth, replace = FALSE)
    tab <- table(samp)
    rich[i] <- length(tab)
    p <- tab / sum(tab)
    shan[i] <- -sum(p * log(p))
  }
  c(rare_richness = mean(rich), rare_shannon = mean(shan))
}

rare_mat <- sapply(colnames(counts_euk), function(s) rarefy_one(counts_euk[, s], min_depth))
rare_df <- as.data.frame(t(rare_mat))
rare_df$Habitat <- meta$Habitat[match(rownames(rare_df), meta$SampleID)]

cat("\n=== Rarefied alpha diversity by habitat ===\n")
print(aggregate(cbind(rare_richness, rare_shannon) ~ Habitat, data = rare_df,
                FUN = function(x) c(mean = mean(x), sd = sd(x))))
print(wilcox.test(rare_richness ~ Habitat, data = rare_df))
print(wilcox.test(rare_shannon  ~ Habitat, data = rare_df))

## -----------------------------
## 14) Save key objects for downstream use / figure-building
## -----------------------------
saveRDS(meta, "scores_euk.rds")
saveRDS(rel, "rel_euk.rds")
saveRDS(tax_euk, "tax_euk.rds")
saveRDS(station_means, "pennate_diatom_station_means.rds")
saveRDS(class_summary, "class_summary_full_euk.rds")
saveRDS(alpha_df, "alpha_df_euk.rds")
saveRDS(rare_df, "rare_df_euk.rds")

cat("\nDone. Key objects saved: scores_euk.rds, rel_euk.rds, tax_euk.rds, pennate_diatom_station_means.rds, class_summary_full_euk.rds, alpha_df_euk.rds, rare_df_euk.rds\n")

###############################################
## Appendix: station 13/14 water sample diagnostic
##
## Documents why stations 13 and 14 have no/near-zero eukaryotic water
## data. Both failures occur at Cascabel demultiplexing, before DADA2.
## Station 13 yields exactly zero reads; station 14 yields exactly one.
## The equivalent 16S water samples at both stations are unaffected.
###############################################

## Confirm station 13's water sample exists in the raw mapping file
## (i.e. was collected) but never reaches counts_raw
mapping_raw %>%
  filter(str_detect(Description, "ST_13_ST_a")) %>%
  select(SampleID, Description)

cat("\nWas NIOZ349.9193.9094 (station 13 water) ever a column in the raw ASV table?",
    "NIOZ349.9193.9094" %in% names(counts_raw), "\n")

## Cross-check against Cascabel's DADA2 filter_summary.out for the 18S run
## (expects a "reads.in"/"filtered" style column; adjust names if the local
## file differs). Confirms station 13 water absent from this file entirely;
## station 14 water present with reads.in = 1.
if (file.exists("filter_summary_euk.out")) {
  filter_summary_euk <- read_tsv("filter_summary_euk.out", show_col_types = FALSE)
  cat("\nStation 13 water in filter_summary_euk.out:\n")
  print(filter_summary_euk %>% filter(str_detect(.[[1]], "9193")))
  cat("\nStation 14 water in filter_summary_euk.out:\n")
  print(filter_summary_euk %>% filter(str_detect(.[[1]], "9195")))
}
