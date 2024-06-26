#' Compute silhouette scores
#' 
#' @description
#' Compute silhouette scores using the PCA embedding space (default 30 dimensions) for each resolution's communities.
#'
#' @param sobj A Seurat object
#' @param res_search A vector of type double that specifies the resolution parameters
#' @param output_path A character denoting the output path                                 
#' @return sil_df silhouette score summary statistics for each resolution's communities
#' @export
compute_silhouette_scores <- function(sobj, res_search, output_path){
    sil_list <- list()
    for (j in 1:length(res_search)) {
        Idents(sobj) <- res_search[j]
        clusters <- Idents(sobj)
        dist.matrix <- dist(x = Embeddings(object = sobj[['pca']])[, 1:30])
        sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
        sil_summary <- aggregate(sil[,3], list(sil[,1]), mean)
        sil_summary$Group.1 <- sil_summary$Group.1 - 1
        colnames(sil_summary) <- c(res_search[j], 'sil')
        sil_list[[res_search[j]]] <- sil_summary
    }
    sil_df <- wrangle_silhouette_scores(sil_list)
    write.csv(sil_df, output_path)
    return(sil_df)
}

#' Wrangle silhouette scores
#'
#' @param sil_list A list of silhouette score summary statistics for each resolution's communities
#' @return sil_df A dataframe containing the silhouette score summary statistics for each resolution's communities
#' @export
wrangle_silhouette_scores <- function(sil_list){
    for(i in seq_along(sil_list)){
        sil_list[i][[1]][,1] <- paste(colnames(sil_list[i][[1]])[1], sil_list[i][[1]][,1], sep='_')
        colnames(sil_list[i][[1]])[1] <- "res"
        sil_list[i][[1]] <- sil_list[i][[1]]
    }
    sil_df <- Reduce(function(...) full_join(..., by = c('sil', 'res')), sil_list)
    colnames(sil_df) <- c('res_cluster', 'sil')
    return(sil_df)
}