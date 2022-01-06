#' Compute enrichment
#'
#' @param dbs
#' @param sobj
#' @param res_search
#' @param output_path
#' @import enrichR
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @return 
#' @export
compute_enrichment <- function(dbs, sobj, res_search, output_path){
    enrich_list <- list()
    enrich_list2 <- list()
    markers_list <- list()
    # TODO: Update res_search parameter for loop. Clanky!
    for (res in res_search){
        Idents(sobj) <- paste0('res.', res)
        df <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
        markers_list[[paste0('res.', res)]] <- df
        for (c in unique(df$cluster)){
            c_df <- df %>% filter(cluster %in% c)
            enriched <- enrichr(c_df$gene, dbs)
            enrich_list2[[paste0('res.', res, '_', c)]] <- enriched
            enriched <- enriched[[dbs]]
            enriched <- enriched[,c('Term','Combined.Score')]
            colnames(enriched) <- c('gene.set', paste0('res.', res, '_', c))
            enrich_list[[paste0('res.', res, '_', c)]] <- enriched
        }
    }
    top_genes_df <- select_geneset_markers(markers_list)
    write.csv(top_genes_df, paste0(output_path, 'top_genes.csv'))
    
    enriched_combined <- Reduce(function(...) merge(..., by='gene.set',all=T), enrich_list)
    enriched_combined[is.na(enriched_combined)] <- 0
    enriched_combined <- gather(enriched_combined, 'res', 'combined.score', -gene.set)
    colnames(enriched_combined) <- c('gene.set', 'res_cluster', 'combined.score')
    write.csv(enriched_combined, paste0(output_path,'enrichment.csv'))
}

#' Select geneset markers 
#' Default is the top 5 genes
#'
#' @param markers_list
#' @import enrichR
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @return 
#' @export
select_geneset_markers <- function(markers_list){
    for (res in names(markers_list)){
        markers_list[[res]][,'cluster'] <- paste0(res, '_', markers_list[[res]][,'cluster'])
    }
    markers_df <- do.call("rbind", markers_list)
    top_genes_df <- data.frame(res=c(), top_genes=c())
    for (c in unique(markers_df$cluster)){
        c_df <- markers_df %>% filter(cluster %in% c)
        top_genes <- paste(c_df$gene[1:5], collapse=',')
        row <- data.frame(res=c(c), top_genes = c(top_genes))
        top_genes_df <- rbind(top_genes_df, row)
    }
    colnames(top_genes_df) <- c('res_cluster', 'top_genes')
    return(top_genes_df)
}





                            