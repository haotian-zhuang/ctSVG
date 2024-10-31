#' Generate single-cell gene expression profiles
#'
#' @param cellassignRes The output from \code{\link{cellassign}}
#' @param visiumpath The path to the Visium HD output folder
#' @import Seurat Matrix
#' @return The single-cell gene expression count matrix
#' @export
#'
getexpr <- function(cellassignRes, visiumpath) {
  mat <- Read10X_h5(paste0(visiumpath, "/filtered_feature_bc_matrix.h5"))
  mat <- mat[, names(cellassignRes$assign)]
  mat <- t(mat)
  g <- cut(1:ncol(mat), ncol(mat) / 2000)
  t(do.call(cbind, sapply(unique(g), function(i) {
    as(SparseArray::rowsum(mat[, g == i], cellassignRes$assign, reorder = FALSE), Class = "dgCMatrix")
  }, simplify = FALSE)))
}
