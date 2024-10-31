#' Run Seurat pipeline
#'
#' @param count The output from \code{\link{getexpr}}
#' @param min.genes The minimum number of genes that are expressed in a cell for it to be retained (300 by default)
#' @param min.cells The minimum proportion of cells a gene is expressed in to be retained (0.01 by default)
#' @param npcs Number of PCs to use (10 by default)
#' @param resolution The parameter in the \code{FindClusters} function (1.2 by default)
#' @import Seurat
#' @return A Seurat object with cell clustering results
#' @export
#'
runseurat <- function(count, min.genes = 300, min.cells = 0.01, npcs = 10, resolution = 1.2) {
  count <- count[, colSums(count > 0) >= min.genes]
  count <- count[rowMeans(count > 0) > min.cells, ]
  d <- CreateSeuratObject(count)
  d <- NormalizeData(d, verbose = FALSE)
  d <- FindVariableFeatures(d, verbose = FALSE)
  d <- ScaleData(d, verbose = FALSE)
  d <- RunPCA(d, verbose = FALSE, npcs = npcs)
  d <- FindNeighbors(d, verbose = FALSE, reduction = "pca", dims = 1:npcs)
  d <- FindClusters(d, verbose = FALSE, resolution = resolution, random.seed = 0)
  d <- RunUMAP(d, verbose = FALSE, reduction = "pca", dims = 1:npcs)
  d
}
