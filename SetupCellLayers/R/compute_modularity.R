#' Compute modularity
#'
#' @param 
#' @return
#' @export
compute_modularity <- function(sobj,res_search){
    # Resolution grid search
    modularity_list <- list()
    for (i in seq_along(res_search)) {
        out <- capture.output(sobj <- FindClusters(sobj, resolution=res_search[[i]]), type =  "output")
        sobj@meta.data$'seurat_clusters' <- NULL
        res_name <- paste('res.',res_search[[i]],sep='')
        modularity_list[[res_name]] <- as.numeric(unlist(strsplit(x = out[7], split = ":"))[2])
        names(sobj@meta.data)[length(names(sobj@meta.data))]<-res_name
    }
    
    mod_df <- stack(modularity_list)
    colnames(mod_df) <- c('modularity', 'resolution')    
    mod_df <- mod_df[,c('resolution', 'modularity')]
    write.csv(mod_df, '../Data/PBMC/pbmc_modularity.csv')
    sil_list <- compute_silhouette_score(as.character(mod_df$'resolution'),sobj)
    return(list('sobj'=sobj, 'mod'=mod_df, 'sil'=sil_list))
}