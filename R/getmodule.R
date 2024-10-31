#' Identify gene modules within each cell cluster
#'
#' @param d The output from \code{\link{runseurat}}
#' @param clutestRes The output from \code{\link{ctsvg_test}}
#' @param clufitRes The output from \code{\link{ctsvg_fit}}
#' @param seed Random seed for k-means clustering (2024 by default)
#' @param clunum Number of modules (\code{NULL} by default, automatically selected)
#' @import Seurat
#' @return A list with each element representing a cell cluster and indicating which module each gene belongs to
#' @export
#'
ctsvg_module <- function(d, clutestRes, clufitRes, seed = 2024, clunum = NULL) {
  sc <- Idents(d)
  sccell <- split(names(sc), sc)
  ctsvg <- clutestRes[clutestRes$fdr < 0.05, ]
  
  set.seed(seed = seed)
  module <- sapply(unique(ctsvg$cluster), function(k) {
    gene.tmp <- ctsvg[ctsvg$cluster == k, "gene"]
    pl <- clufitRes[gene.tmp, intersect(colnames(clufitRes), sccell[[k]])]
    ml <- Matrix::rowMeans(pl, na.rm = TRUE)
    sl <- matrixStats::rowSds(pl, na.rm = TRUE)
    mat <- (pl - ml) / sl
    mykmeans(mat = mat, clunum = clunum)
  }, simplify = FALSE)
  return(module)
}

mykmeans <- function(mat, clunum = NULL) {
  if(is.null(clunum)) {
    per <- sapply(1:30, function(clunum) {
      tmp <- kmeans(mat, centers = clunum, iter.max = 1e3)
      tmp$tot.withinss / tmp$totss
    })
    clunum <- findPC::findPC(sort(per, decreasing = TRUE), number = 30, figure = FALSE)[1, 1]
  }
  
  kl <- kmeans(mat, centers = clunum, iter.max = 1e3)
  kl$cluster
}
