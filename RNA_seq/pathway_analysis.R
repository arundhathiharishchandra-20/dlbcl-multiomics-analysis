BASE_DIR <- "/mnt/d/lymphoma_project"   # <-- change this to your own path
library(gprofiler2)
library(WGCNA)

module_genes <- list()
for(mod in 1:11){
    module_genes[[mod]] <- names(net$colors[net$colors == mod])
}

all_results <- list()

for(mod in 1:11){
    cat("Running Module", mod, "(", length(module_genes[[mod]]), "genes)\n")
    genes <- module_genes[[mod]]
    if(length(genes) > 10){
        gost_result <- gost(query = genes, organism = "hsapiens", sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG"), significant = TRUE, user_threshold = 0.05)
        all_results[[mod]] <- gost_result
        if(!is.null(gost_result)){
            cat("Module", mod, "- Pathways:", nrow(gost_result$result), "\n")
        } else {
            cat("Module", mod, "- No significant pathways\n")
        }
    }
}

cat("All done!\n")
save.image(paste0(BASE_DIR, "/RNA_seq/counts/workspace.RData"))
