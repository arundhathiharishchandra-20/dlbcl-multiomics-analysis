# ═══════════════════════════════════════════════════════════════════
# Multi-Omics Analysis: RNA-seq + WGBS Methylation for DLBCL
# ═══════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
library(pheatmap)
library(igraph)

# ── Paths ──────────────────────────────────────────────────────────
COUNTS_DIR <- "/mnt/d/lymphoma_project/RNA_seq/counts"
OUTDIR     <- "/mnt/d/lymphoma_project/multiomics_output"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ── Load data ──────────────────────────────────────────────────────
cat("Loading data...\n")
expr    <- read.csv(file.path(COUNTS_DIR, "Norm_rnaseq.csv"),
                    row.names=1, check.names=FALSE)
meth    <- read.csv(file.path(COUNTS_DIR, "Gene_Methylation.csv"),
                    row.names=1, check.names=FALSE)
modules <- read.csv(file.path(COUNTS_DIR, "WGCNA_ModuleGenes.csv"))
paired  <- read.csv("/mnt/d/lymphoma_project/WGBS/mofa_paired_samples.csv")

cat("RNA-seq:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("Methylation:", nrow(meth), "genes x", ncol(meth), "samples\n")
cat("Paired samples available:", nrow(paired), "\n")

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
cat("Methylation after ID conversion:", nrow(meth), "genes\n")# ── Match samples between RNA-seq and WGBS ─────────────────────────
cat("\nMatching samples...\n")

# Filter paired to only samples present in both datasets
paired_avail <- paired[
    paired$WGBS_GSM %in% colnames(meth) &
    paired$TRIM_ID  %in% colnames(expr), ]

cat("Matched paired samples:", nrow(paired_avail), "\n")

# Subset and rename columns to match
expr_matched <- expr[, paired_avail$TRIM_ID, drop=FALSE]
meth_matched <- meth[, paired_avail$WGBS_GSM, drop=FALSE]

# Rename methylation columns to TRIM IDs for consistency
colnames(meth_matched) <- paired_avail$TRIM_ID

cat("Expression matrix:", nrow(expr_matched), "x", ncol(expr_matched), "\n")
cat("Methylation matrix:", nrow(meth_matched), "x", ncol(meth_matched), "\n")

# ── Find common genes ──────────────────────────────────────────────
commonGenes <- intersect(rownames(expr_matched), rownames(meth_matched))
cat("Common genes:", length(commonGenes), "\n")

expr_matched <- expr_matched[commonGenes, ]
meth_matched <- meth_matched[commonGenes, ]

# ── Correlation analysis ───────────────────────────────────────────
cat("\nRunning correlation analysis...\n")
cor_results <- data.frame()

for(gene in commonGenes){
    expr_vals <- as.numeric(expr_matched[gene, ])
    meth_vals <- as.numeric(meth_matched[gene, ])

    # Skip if too many NAs
    valid <- !is.na(expr_vals) & !is.na(meth_vals)
    if(sum(valid) < 3) next

    r <- tryCatch(
        cor.test(expr_vals[valid], meth_vals[valid], method="pearson"),
        error = function(e) NULL
    )
    if(is.null(r)) next

    cor_results <- rbind(cor_results, data.frame(
        Gene        = gene,
        Correlation = r$estimate,
        Pvalue      = r$p.value
    ))
}

cat("Genes tested:", nrow(cor_results), "\n")

# ── FDR correction ─────────────────────────────────────────────────
cor_results$FDR <- p.adjust(cor_results$Pvalue, method="BH")

# ── Significant genes ──────────────────────────────────────────────
sigGenes <- cor_results %>%
    filter(Pvalue < 0.05 & abs(Correlation) > 0.3)
cat("Significant genes (FDR<0.05, |r|>0.3):", nrow(sigGenes), "\n")

# Save results
write.csv(cor_results, file.path(OUTDIR, "All_Correlation_Results.csv"), row.names=FALSE)
write.csv(sigGenes,    file.path(OUTDIR, "Significant_Correlated_Genes.csv"), row.names=FALSE)

# ── Volcano plot ───────────────────────────────────────────────────
cat("Plotting volcano plot...\n")
cor_results$Significant <- cor_results$FDR < 0.05

p1 <- ggplot(cor_results,
    aes(Correlation, -log10(FDR), color=Significant)) +
    geom_point(alpha=0.5, size=1) +
    scale_color_manual(values=c("grey", "red")) +
    theme_bw() +
    labs(title="Methylation-Expression Correlation",
         x="Pearson Correlation", y="-log10(FDR)") +
    geom_vline(xintercept=c(-0.3, 0.3), linetype="dashed")

ggsave(file.path(OUTDIR, "Volcano_Plot.pdf"), p1, width=8, height=6)
cat("Volcano plot saved!\n")

# ── Merge with WGCNA modules ───────────────────────────────────────
cat("\nMerging with WGCNA modules...\n")
merged <- merge(sigGenes, modules, by="Gene")
cat("Significant genes with module assignment:", nrow(merged), "\n")
print(table(merged$Module))
write.csv(merged, file.path(OUTDIR, "SigGenes_with_Modules.csv"), row.names=FALSE)

# ── Heatmap of top correlated genes ───────────────────────────────
if(nrow(sigGenes) >= 10){
    cat("Plotting heatmap...\n")
    topGenes <- head(sigGenes$Gene[order(abs(sigGenes$Correlation), decreasing=TRUE)], 50)
    topGenes <- topGenes[topGenes %in% rownames(expr_matched)]

    pdf(file.path(OUTDIR, "Heatmap_Top_Genes.pdf"), width=12, height=10)
    pheatmap(expr_matched[topGenes, ],
        scale="row",
        show_rownames=TRUE,
        show_colnames=TRUE,
        main="Top Correlated Genes - Expression")
    dev.off()
    cat("Heatmap saved!\n")
}

# ── GO enrichment ──────────────────────────────────────────────────
if(nrow(merged) >= 10){
    cat("\nRunning GO enrichment...\n")
    tryCatch({
        ego <- enrichGO(
            gene          = merged$Gene,
            OrgDb         = org.Hs.eg.db,
            keyType       = "ENSEMBL",
            ont           = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff  = 0.05
        )
        if(!is.null(ego) && nrow(ego@result) > 0){
            write.csv(as.data.frame(ego),
                file.path(OUTDIR, "GO_Results.csv"))
            cat("GO enrichment done:", nrow(ego@result), "terms\n")

            pdf(file.path(OUTDIR, "GO_Dotplot.pdf"), width=10, height=8)
            print(dotplot(ego, showCategory=20))
            dev.off()
        }
    }, error = function(e) cat("GO enrichment failed:", e$message, "\n"))
}

# ── STRING network ─────────────────────────────────────────────────
if(nrow(merged) >= 5){
    cat("\nBuilding STRING network...\n")
    tryCatch({
        string_db <- STRINGdb$new(
            version         = "12",
            species         = 9606,
            score_threshold = 400
        )
        mapped <- string_db$map(merged, "Gene", removeUnmappedRows=TRUE)
        if(nrow(mapped) > 0){
            ppi <- string_db$get_interactions(mapped$STRING_id)
            write.csv(ppi, file.path(OUTDIR, "STRING_Network.csv"), row.names=FALSE)
            cat("STRING network saved:", nrow(ppi), "interactions\n")

            pdf(file.path(OUTDIR, "STRING_Network.pdf"), width=12, height=10)
            string_db$plot_network(mapped$STRING_id)
            dev.off()
        }
    }, error = function(e) cat("STRING network failed:", e$message, "\n"))
}

cat("\n=== Multi-Omics Analysis Complete! ===\n")
cat("All outputs saved to:", OUTDIR, "\n")
cat("\nFiles generated:\n")
cat("  All_Correlation_Results.csv\n")
cat("  Significant_Correlated_Genes.csv\n")
cat("  SigGenes_with_Modules.csv\n")
cat("  Volcano_Plot.pdf\n")
cat("  Heatmap_Top_Genes.pdf\n")
cat("  GO_Results.csv\n")
cat("  GO_Dotplot.pdf\n")
cat("  STRING_Network.csv\n")
cat("  STRING_Network.pdf\n")

