
compute_enrichment <- function(){
    dbs <- "GO_Biological_Process_2018"
    enrich_list <- list()
    enrich_list2 <- list()
    markers_list <- list() 
    for (res in seq(0.1, 0.5, by=0.1)){
        Idents(cl.prep$sobj) <- paste0('res.', res)
        df <- FindAllMarkers(cl.prep$sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
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
}