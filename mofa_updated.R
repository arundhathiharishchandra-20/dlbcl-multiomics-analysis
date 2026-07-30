# ═══════════════════════════════════════════════════════════════════
# MOFA2 Multi-Omics Integration: RNA-seq + WGBS for DLBCL
# ═══════════════════════════════════════════════════════════════════

if(!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
if(!requireNamespace("MOFA2", quietly=TRUE))
    BiocManager::install("MOFA2")

library(MOFA2)
library(data.table)
library(ggplot2)

# ── Paths ──────────────────────────────────────────────────────────
COUNTS_DIR <- "/mnt/d/lymphoma_project/RNA_seq/counts"
OUTDIR     <- "/mnt/d/lymphoma_project/mofa_output"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ── Load data ──────────────────────────────────────────────────────
cat("Loading data...\n")
rna    <- fread(file.path(COUNTS_DIR, "Norm_rnaseq.csv"), data.table=FALSE)
meth   <- fread(file.path(COUNTS_DIR, "Gene_Methylation.csv"), data.table=FALSE)
paired <- read.csv("/mnt/d/lymphoma_project/WGBS/mofa_paired_samples.csv")

rownames(rna)  <- rna[,1];  rna  <- rna[,-1]
rownames(meth) <- meth[,1]; meth <- meth[,-1]

cat("RNA-seq:", nrow(rna), "genes x", ncol(rna), "samples\n")
cat("Methylation:", nrow(meth), "genes x", ncol(meth), "samples\n")

# Convert methylation gene symbols to ENSEMBL IDs
library(org.Hs.eg.db)
symbols <- rownames(meth)
mapping <- select(org.Hs.eg.db,
    keys = symbols,
    columns = "ENSEMBL",
    keytype = "SYMBOL")
mapping <- mapping[!is.na(mapping$ENSEMBL), ]
mapping <- mapping[!duplicated(mapping$SYMBOL), ]
mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
cat("Mapped genes:", nrow(mapping), "\n")
meth <- meth[mapping$SYMBOL, ]
rownames(meth) <- mapping$ENSEMBL
cat("Methylation after ID conversion:", nrow(meth), "genes\n")# ── Match samples ──────────────────────────────────────────────────
cat("\nMatching samples...\n")
paired_avail <- paired[
    paired$WGBS_GSM %in% colnames(meth) &
    paired$TRIM_ID  %in% colnames(rna), ]
cat("Matched paired samples:", nrow(paired_avail), "\n")

rna_matched  <- rna[,  paired_avail$TRIM_ID,  drop=FALSE]
meth_matched <- meth[, paired_avail$WGBS_GSM, drop=FALSE]

# Rename methylation to TRIM IDs
colnames(meth_matched) <- paired_avail$TRIM_ID

# ── Find common genes ──────────────────────────────────────────────
commonGenes <- intersect(rownames(rna_matched), rownames(meth_matched))
cat("Common genes:", length(commonGenes), "\n")

rna_matched  <- rna_matched[commonGenes,  ]
meth_matched <- meth_matched[commonGenes, ]

# ── Remove genes with too many NAs in methylation ──────────────────
na_frac <- rowMeans(is.na(meth_matched))
keep    <- na_frac < 0.5
rna_matched  <- rna_matched[keep,  ]
meth_matched <- meth_matched[keep, ]
cat("Genes after NA filter:", nrow(rna_matched), "\n")

# Fill remaining NAs with row mean
for(i in 1:nrow(meth_matched)){
    row_mean <- mean(as.numeric(meth_matched[i,]), na.rm=TRUE)
    meth_matched[i, is.na(meth_matched[i,])] <- row_mean
}

# ── Create MOFA data list ──────────────────────────────────────────
cat("\nCreating MOFA object...\n")
data_list <- list(
    RNA         = as.matrix(rna_matched),
    Methylation = as.matrix(meth_matched)
)

cat("RNA matrix:", nrow(data_list$RNA), "x", ncol(data_list$RNA), "\n")
cat("Methylation matrix:", nrow(data_list$Methylation), "x", ncol(data_list$Methylation), "\n")

# ── Create and train MOFA model ────────────────────────────────────
MOFAobject <- create_mofa(data_list)
print(MOFAobject)

data_opts  <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

train_opts <- get_default_training_options(MOFAobject)
train_opts$maxiter          <- 1000
train_opts$convergence_mode <- "medium"
train_opts$seed             <- 42

MOFAobject <- prepare_mofa(MOFAobject,
    data_options     = data_opts,
    model_options    = model_opts,
    training_options = train_opts)

outfile <- file.path(OUTDIR, "DLBCL_MOFA2.hdf5")
cat("\nTraining MOFA model...\n")
MOFAobject <- run_mofa(MOFAobject, outfile, use_basilisk=TRUE)
cat("Training complete!\n")

saveRDS(MOFAobject, file.path(OUTDIR, "MOFA_model.rds"))

# ── Load and plot results ──────────────────────────────────────────
model <- load_model(outfile)

# Variance explained
pdf(file.path(OUTDIR, "VarianceExplained.pdf"), width=10, height=6)
print(plot_variance_explained(model))
dev.off()
cat("Variance explained plot saved!\n")

# Factor scatter
pdf(file.path(OUTDIR, "Factor_Scatter.pdf"), width=8, height=6)
print(plot_factors(model, factors=1:min(5, model@dimensions$K)))
dev.off()

# Factor correlation
pdf(file.path(OUTDIR, "Factor_Correlation.pdf"), width=8, height=7)
print(plot_factor_cor(model))
dev.off()

# ── Extract weights ────────────────────────────────────────────────
cat("\nExtracting weights...\n")
rna_weights  <- get_weights(model, views="RNA",         factors="all")[[1]]
meth_weights <- get_weights(model, views="Methylation", factors="all")[[1]]

write.csv(rna_weights,  file.path(OUTDIR, "RNA_Weights.csv"))
write.csv(meth_weights, file.path(OUTDIR, "Methylation_Weights.csv"))

# Top RNA genes for Factor 1
cat("\nTop RNA genes for Factor 1:\n")
top_rna <- rna_weights[order(abs(rna_weights[,1]), decreasing=TRUE), , drop=FALSE]
print(head(top_rna, 20))

# Top methylation genes for Factor 1
cat("\nTop methylation genes for Factor 1:\n")
top_meth <- meth_weights[order(abs(meth_weights[,1]), decreasing=TRUE), , drop=FALSE]
print(head(top_meth, 20))

# Weight plots
pdf(file.path(OUTDIR, "Weights_RNA_Factor1.pdf"), width=8, height=6)
print(plot_weights(model, view="RNA", factor=1, nfeatures=20))
dev.off()

pdf(file.path(OUTDIR, "Weights_Methylation_Factor1.pdf"), width=8, height=6)
print(plot_weights(model, view="Methylation", factor=1, nfeatures=20))
dev.off()

# Heatmaps
pdf(file.path(OUTDIR, "Heatmap_RNA_Factor1.pdf"), width=10, height=8)
print(plot_data_heatmap(model, view="RNA", factor=1, features=25,
    show_rownames=TRUE, show_colnames=FALSE))
dev.off()

pdf(file.path(OUTDIR, "Heatmap_Methylation_Factor1.pdf"), width=10, height=8)
print(plot_data_heatmap(model, view="Methylation", factor=1, features=25,
    show_rownames=TRUE, show_colnames=FALSE))
dev.off()

# Save factor values
factors <- get_factors(model, factors="all")[[1]]
write.csv(factors, file.path(OUTDIR, "MOFA_Factors.csv"))

cat("\n=== MOFA2 Analysis Complete! ===\n")
cat("All outputs saved to:", OUTDIR, "\n")
cat("\nFiles generated:\n")
cat("  DLBCL_MOFA2.hdf5              - Trained model\n")
cat("  MOFA_model.rds                - R model object\n")
cat("  MOFA_Factors.csv              - Factor values\n")
cat("  RNA_Weights.csv               - RNA feature weights\n")
cat("  Methylation_Weights.csv       - Methylation weights\n")
cat("  VarianceExplained.pdf         - Variance per factor\n")
cat("  Factor_Scatter.pdf            - Factor scatter plot\n")
cat("  Factor_Correlation.pdf        - Factor correlations\n")
cat("  Weights_RNA_Factor1.pdf       - Top RNA features\n")
cat("  Weights_Methylation_Factor1.pdf\n")
cat("  Heatmap_RNA_Factor1.pdf\n")
cat("  Heatmap_Methylation_Factor1.pdf\n")
