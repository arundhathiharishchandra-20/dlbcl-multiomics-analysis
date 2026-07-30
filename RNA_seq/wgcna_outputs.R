BASE_DIR <- "/mnt/d/lymphoma_project"   # <-- change this to your own path
library(WGCNA)
library(gprofiler2)

# ─── MODULE SUMMARY TABLE ────────────────────────────────────────
module_summary <- data.frame()

for(mod in 1:11){
    genes <- names(net$colors[net$colors == mod])
    n_genes <- length(genes)
    
    top_pathway <- "No significant pathway"
    if(!is.null(all_results[[mod]])){
        res <- all_results[[mod]]$result
        res <- res[order(res$p_value), ]
        top_pathway <- res$term_name[1]
    }
    
    module_summary <- rbind(module_summary, data.frame(
        Module = mod,
        Genes = n_genes,
        Top_Pathway = top_pathway
    ))
}

write.csv(module_summary, paste0(BASE_DIR, "/RNA_seq/counts/module_summary.csv"), row.names = FALSE)
cat("Module summary saved!\n")

# ─── TOP 20 HUB GENES ────────────────────────────────────────────
MEs <- moduleEigengenes(datExpr_hvg, net$colors)$eigengenes
kME <- cor(datExpr_hvg, MEs, use = "p")

hub_genes_all <- data.frame()
for(mod in 1:11){
    mod_genes <- names(net$colors[net$colors == mod])
    mod_kme <- kME[mod_genes, paste0("ME", mod)]
    top20 <- names(sort(mod_kme, decreasing = TRUE))[1:min(20, length(mod_genes))]
    top20_kme <- sort(mod_kme, decreasing = TRUE)[1:min(20, length(mod_genes))]
    hub_genes_all <- rbind(hub_genes_all, data.frame(
        Module = mod,
        Gene = top20,
        kME = round(top20_kme, 3)
    ))
}

write.csv(hub_genes_all, paste0(BASE_DIR, "/RNA_seq/counts/top20_hub_genes_all_modules.csv"), row.names = FALSE)
cat("Top 20 hub genes saved!\n")

# ─── MODULE EIGENGENE CLUSTERING PLOT ────────────────────────────
MEs_plot <- moduleEigengenes(datExpr_hvg, net$colors)$eigengenes

png(paste0(BASE_DIR, "/RNA_seq/counts/module_eigengene_clustering.png"), width = 1200, height = 800, res = 150)
plotEigengeneNetworks(MEs_plot, "Module Eigengene Clustering", marDendro = c(0,4,2,0), marHeatmap = c(3,4,2,2), plotDendrograms = TRUE, xLabelsAngle = 90)
dev.off()
cat("Module eigengene clustering saved!\n")

cat("ALL OUTPUTS SAVED!\n")
