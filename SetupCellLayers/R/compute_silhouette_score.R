#' Compute silhouette scores
#'
#' @param 
#' @return
#' @export
compute_silhouette_score <- function(res,sobj){
    sil_list <- list()
    for (j in 1:length(res)) {
        Idents(sobj) <- res[j]
        clusters <- Idents(sobj)
        dist.matrix <- dist(x = Embeddings(object = sobj[['pca']])[, 1:30])
        sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
        sil_summary <- aggregate(sil[,3], list(sil[,1]), mean)
        sil_summary$Group.1 <- sil_summary$Group.1 - 1
        colnames(sil_summary) <- c(res[j], 'sil')
        sil_list[[res[j]]] <- sil_summary
    }
    return(sil_list)
}