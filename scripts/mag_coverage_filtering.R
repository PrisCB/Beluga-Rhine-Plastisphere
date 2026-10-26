# Reconstruct the 120-MAG relative-abundance matrix from archived inputs.
# Run from the repository root: Rscript scripts/mag_coverage_filtering.R
# Requires base R only. This does NOT validate upstream coverage assignments.
options(stringsAsFactors = FALSE)
coverage_file <- "scicore/mag_coverage_byMAG_v2.tsv"
qc_file <- "data/metadata/mag_qc.csv"
selected_file <- "data/metadata/mag_metadata_120.csv"
samples_file <- "data/metadata/sample_metadata.csv"
out_dir <- "results/mag_processing"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
coverage <- read.delim(coverage_file, check.names = FALSE)
qc <- read.csv(qc_file, check.names = FALSE)
selected <- read.csv(selected_file, check.names = FALSE)
samples <- read.csv(samples_file, check.names = FALSE)
stopifnot("MAG_ID" %in% names(coverage), all(c("MAG", "Completeness", "Contamination") %in% names(qc)),
          "MAG" %in% names(selected), "SampleID_cov" %in% names(samples))
stopifnot(nrow(coverage) == 302L, nrow(qc) == 302L, nrow(selected) == 120L,
          !anyDuplicated(coverage$MAG_ID), !anyDuplicated(qc$MAG),
          !anyDuplicated(selected$MAG), setequal(coverage$MAG_ID, qc$MAG))
sample_ids <- names(coverage)[-1]
stopifnot(length(sample_ids) == 14L, setequal(sample_ids, samples$SampleID_cov))
coverage_mat <- as.matrix(coverage[, -1, drop = FALSE])
rownames(coverage_mat) <- coverage$MAG_ID
storage.mode(coverage_mat) <- "double"
stopifnot(all(is.finite(coverage_mat)), all(coverage_mat >= 0), all(colSums(coverage_mat) > 0))
# Important: normalise across all 302 MAGs BEFORE selecting 120.
rel_mat <- sweep(coverage_mat, 2, colSums(coverage_mat), "/")
qc_aligned <- qc[match(rownames(rel_mat), qc$MAG), , drop = FALSE]
max_rel <- apply(rel_mat, 1, max)
eligible <- qc_aligned$Completeness >= 70 & qc_aligned$Contamination <= 10 &
            max_rel >= 0.001 & rowSums(coverage_mat > 0) >= 1
eligible_ids <- rownames(rel_mat)[eligible]
# Retain the archived 120-MAG set and order; independently verify eligibility.
stopifnot(all(selected$MAG %in% rownames(rel_mat)),
          all(selected$Completeness >= 70), all(selected$Contamination <= 10),
          all(selected$MAG %in% eligible_ids))
if (!setequal(eligible_ids, selected$MAG)) {
  warning("The archived 120 MAGs do not equal the full set passing all thresholds. Inspect eligibility_report.csv.")
}
rel_mat_filt <- rel_mat[selected$MAG, , drop = FALSE]
write.csv(rel_mat, file.path(out_dir, "relative_abundance_302_MAGs.csv"))
write.csv(rel_mat_filt, file.path(out_dir, "relative_abundance_120_MAGs.csv"))
report <- data.frame(MAG = rownames(rel_mat), completeness = qc_aligned$Completeness,
  contamination = qc_aligned$Contamination, max_relative_abundance = max_rel,
  detected_samples = rowSums(coverage_mat > 0), eligible = eligible,
  in_archived_120 = rownames(rel_mat) %in% selected$MAG)
write.csv(report, file.path(out_dir, "eligibility_report.csv"), row.names = FALSE)
cat("MAGs:", nrow(rel_mat), "| eligible:", sum(eligible),
    "| archived selection:", nrow(rel_mat_filt), "\n")
cat("Outputs written to", out_dir, "\n")
