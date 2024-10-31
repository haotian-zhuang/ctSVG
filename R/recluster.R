#' Reassign cells to clusters
#'
#' @param d The output from \code{\link{runseurat}}
#' @param n.max Number of reclustering times
#' @param resolution The parameter used in \code{\link{runseurat}} (1.2 by default)
#' @import Seurat
#' @return A list including the matrix of Jaccard indices, and the list of reassigned clusters
#' @export
#'
recluster <- function(d, n.max = 1e3, resolution = 1.2) {
  sc <- Idents(d)
  sccell <- split(names(sc), sc)
  
  clu.list <- vector(mode = "list", length = length(sccell))
  names(clu.list) <- names(sccell)
  
  jac.mat <- matrix(nrow = length(sccell), ncol = n.max)
  rownames(jac.mat) <- names(sccell)
  
  for (i in seq_len(n.max)) {
    d.tmp <- FindClusters(d, verbose = FALSE, resolution = resolution, random.seed = i)
    sc.tmp <- Idents(d.tmp)
    sccell.tmp <- split(names(sc.tmp), sc.tmp)
    jac.tmp <- PairWiseJaccardSets(sc, sc.tmp)
    
    jac.mat[, i] <- apply(jac.tmp, 1, max)
    match.tmp <- apply(jac.tmp, 1, which.max)
    for (j in names(match.tmp)) {
      clu.list[[j]] <- c(clu.list[[j]], list(sccell.tmp[[match.tmp[j]]]))
    }
  }
  return(list(jac.mat = jac.mat, clu.list = clu.list))
}

JaccardSets <- function(set1, set2){
  length(intersect(set1, set2)) / length(unique(c(set1, set2)))
}

PairWiseJaccardSets <- function(ident1, ident2){
  ident1.list <- split(names(ident1), ident1)
  ident2.list <- split(names(ident2), ident2)
  res <- matrix(nrow = length(ident1.list), ncol = length(ident2.list),
                dimnames = list(names(ident1.list), names(ident2.list)))
  for (i in seq_along(ident1.list)){
    res[i, ] <- purrr::map_dbl(ident2.list, ~JaccardSets(ident1.list[[i]], .x))
  }
  return(res)
}
