#' Compute multi-resolution search
#'
#' Compute clustering for all resolution parameters.
#' 
#' @param sobj A Seurat object
#' @param res_search  A vector of type double that specifies the resolution parameters
#' @param output_path A character denoting the output path
#' @return
#' @export
compute_multi_res_search <- function(sobj, res_search, output_path){
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
    write.csv(mod_df, output_path)
    cl.setup <- list('sobj'=sobj, 'mod'=mod_df)
    return(cl.setup)
}