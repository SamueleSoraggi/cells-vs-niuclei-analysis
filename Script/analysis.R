# Function to convert ENTREZ IDs in core_enrichment to gene symbols

convert_core_enrichment <- function(enrich_result, conversion_table) {
    symbolVector <- sapply(
        enrich_result$core_enrichment,
        function(x) {
            ids <- unlist(strsplit(x, "/"))
            #cat(ids, "\n")
            symbols <- conversion_table$SYMBOL[match(ids, conversion_table$ENTREZID)]
            #cat(symbols, "\n")
            paste(symbols, collapse = "/")
        }
    )
    return(symbolVector)
}


# Function to create GSEA plots for first N or last N GO terms
plot_gsea_terms <- function(gsea_result, table, n_terms = 5, orderby="rank", decreasing = FALSE, use_first = TRUE, ncol = 2, save_folder = FALSE) {
    #object name
    objname <- deparse(substitute(gsea_result))

    # Get the indices of terms to plot
    if (use_first) {
        indices <- 1:min(n_terms, nrow(table))
        plot_title <- paste0("GSEA_plots_top_", n_terms, "_terms_", orderby,"_", objname)
    } else {
        indices <- max(1, nrow(table) - n_terms + 1):nrow(table)
        plot_title <- paste0("GSEA_plots_bottom_", n_terms, "_terms_", orderby,"_", objname)
    }
    
    # Create list to store plots
    plot_list <- list()

    # order the table by ranking. Then find the indices again which were defined above 
    table = table[order(table[, orderby], decreasing = decreasing), ]
    IDS = table$ID[indices]
    #cat(IDS)
    new_indices = which(gsea_result$ID %in% IDS)
    #cat(new_indices)

    # Generate plots for each term
    for (i in new_indices) {
        p1 <- gseaplot(gsea_result, geneSetID = i, by = "runningScore", 
                      title = gsea_result$Description[i])
        p2 <- gseaplot(gsea_result, geneSetID = i, by = "preranked", 
                      title = gsea_result$Description[i])
        plot_list <- c(plot_list, list(p1, p2))
    }
    
    # Create grid of plots
    n_plots <- length(plot_list)
    cowplot::plot_grid(plotlist = plot_list, ncol = ncol)

    # save plot
    if(save_folder!=FALSE){
        ggsave(filename = paste0(save_folder, "/", plot_title, ".pdf"), 
               width = 20, height = 8 * ceiling(n_plots / ncol))
        ggsave(filename = paste0(save_folder, "/", plot_title, ".png"), 
               width = 20, height = 8 * ceiling(n_plots / ncol))
    

    cowplot::plot_grid(plotlist = plot_list, ncol = ncol)

    }
}