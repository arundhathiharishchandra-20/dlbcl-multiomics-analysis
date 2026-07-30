BASE_DIR <- "/mnt/d/lymphoma_project"   # <-- change this to your own path
library(methylKit)

# ─── STEP 1: GET LIST OF ALL COV.GZ FILES ────────────────────────

cov_files <- list.files(
    paste0(BASE_DIR, "/WGBS/bismark_output/methylation/"),
    pattern = "*.bismark.cov.gz",
    full.names = TRUE
)

cat("Number of coverage files found:", length(cov_files), "\n")

# Get sample names from filenames
sample_ids <- gsub("_bismark_bt2.deduplicated.bismark.cov.gz", "", 
                   basename(cov_files))

cat("Sample IDs:\n")
print(sample_ids)

# ─── STEP 2: READ INTO METHYLKIT ─────────────────────────────────

# Since we have no treatment groups, set all as 0
treatment <- rep(0, length(cov_files))

myobj <- methRead(
    as.list(cov_files),
    sample.id = as.list(sample_ids),
    assembly = "hg38",
    pipeline = "bismarkCoverage",
    treatment = treatment,
    mincov = 10
)

cat("methylKit object created!\n")
cat("Number of samples:", length(myobj), "\n")

# ─── STEP 3: QUALITY CHECK ───────────────────────────────────────

# Coverage statistics per sample
png(paste0(BASE_DIR, "/WGBS/methylkit_coverage_stats.png"),
    width = 2000, height = 1000, res = 150)
getCoverageStats(myobj[[1]], plot = TRUE, both.strands = FALSE)
dev.off()

# Methylation statistics
png(paste0(BASE_DIR, "/WGBS/methylkit_methylation_stats.png"),
    width = 2000, height = 1000, res = 150)
getMethylationStats(myobj[[1]], plot = TRUE, both.strands = FALSE)
dev.off()

cat("QC plots saved!\n")

# ─── STEP 4: FILTER BY COVERAGE ──────────────────────────────────

filtered_obj <- filterByCoverage(
    myobj,
    lo.count = 10,
    lo.perc = NULL,
    hi.count = NULL,
    hi.perc = 99.9
)

cat("Coverage filtering done!\n")

# ─── STEP 5: MERGE ALL SAMPLES ───────────────────────────────────

meth <- unite(filtered_obj, destrand = FALSE)

cat("Samples merged!\n")
cat("CpG sites covered in all samples:", nrow(meth), "\n")

# ─── STEP 6: PCA OF METHYLATION DATA ─────────────────────────────

png(paste0(BASE_DIR, "/WGBS/methylkit_PCA.png"),
    width = 1200, height = 1000, res = 150)
PCASamples(meth, screeplot = FALSE)
dev.off()

cat("PCA plot saved!\n")

# ─── STEP 7: CLUSTERING ──────────────────────────────────────────

png(paste0(BASE_DIR, "/WGBS/methylkit_clustering.png"),
    width = 1200, height = 1000, res = 150)
clusterSamples(meth, dist = "correlation", method = "ward",
               plot = TRUE)
dev.off()

cat("Clustering plot saved!\n")

# Save workspace
save(myobj, filtered_obj, meth,
     file = paste0(BASE_DIR, "/WGBS/methylkit_workspace.RData"))

cat("ALL DONE! Workspace saved!\n")
