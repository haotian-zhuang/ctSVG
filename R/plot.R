#' Visualize spatial patterns of a single gene
#'
#' @param d The output from \code{\link{runseurat}}
#' @param clufitRes The output from \code{\link{ctsvg_fit}} (\code{NULL} by default, original gene expression matrix)
#' @param coord The matrix of spatial locations from \code{\link{cellassign}}
#' @param gene Which gene to plot
#' @param cluster A character value indicating which cell cluster to plot (\code{NULL} by default, all cells)
#' @param background Whether to plot all cells (\code{FALSE} by default)
#' @param ... The parameters in the \code{geom_point(...)} function
#' @import ggplot2
#' @return A ggplot2 plot
#' @export
#'
plotGene <- function(d, clufitRes = NULL, coord, gene, cluster = NULL, background = FALSE, ...) {
  sc <- Idents(d)
  sccell <- split(names(sc), sc)
  
  if(is.null(clufitRes)) {
    expr <- d@assays$RNA$data
  } else {
    expr <- clufitRes
  }
  
  if(is.null(cluster)) {
    i <- rownames(coord) 
  } else {
    i <- sccell[[cluster]]
  }
  i <- intersect(intersect(i, rownames(coord)), colnames(expr))
  df <- data.frame(expr = expr[gene, i], row = coord[i, "row"], col = coord[i, "col"])
  l <- df$expr
  l[l > quantile(l, 0.99, na.rm = TRUE)] <- quantile(l, 0.99, na.rm = TRUE)
  l[l < quantile(l, 0.01, na.rm = TRUE)] <- quantile(l, 0.01, na.rm = TRUE)
  df$expr <- l
  
  p <- ggplot(df, aes(x = row, y = col, color = expr))
  if(background) {
    p <- p + geom_point(data = data.frame(row = coord[, "row"], col = coord[, "col"]),
                        color = "gray", size = 0.01, alpha = 0.1, stroke = 0.2)
  }
  p <- p +
    geom_point(size = 0.01, alpha = 0.8, ...) +
    coord_fixed() +
    scale_color_gradientn(colors = SpatialColors(100), name = gene) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size = 6),
          legend.margin = margin(0, 0, 0, 0, "cm"),
          legend.title = element_text(size = 8, face = "italic"),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  p
}

#' Visualize spatial patterns of metagenes within a cell cluster
#'
#' @param d The output from \code{\link{runseurat}}
#' @param clufitRes The output from \code{\link{ctsvg_fit}}
#' @param moduleRes The output from \code{\link{ctsvg_module}}
#' @param coord The matrix of spatial locations from \code{\link{cellassign}}
#' @param cluster A character value indicating which cell cluster to plot
#' @param background Whether to plot all cells (\code{FALSE} by default)
#' @param ... The parameters in the \code{geom_point(...)} function
#' @import ggplot2
#' @return A ggplot2 plot
#' @export
#'
plotMetagene <- function(d, clufitRes, moduleRes, coord, cluster, background = FALSE, ...) {
  sc <- Idents(d)
  sccell <- split(names(sc), sc)
  
  gene.agg <- sapply(names(moduleRes), function(k) {
    clu.tmp <- moduleRes[[k]]
    pl <- clufitRes[names(clu.tmp), intersect(colnames(clufitRes), sccell[[k]])]
    ml <- Matrix::rowMeans(pl, na.rm = TRUE)
    sl <- matrixStats::rowSds(pl, na.rm = TRUE)
    mat <- (pl - ml) / sl
    
    clu.list <- split(names(clu.tmp), clu.tmp)
    do.call(rbind, sapply(names(clu.list), function(l) {
      expr.ave <- Matrix::colMeans(mat[clu.list[[l]], ], na.rm = TRUE)
      data.frame(module = paste("Module", l), cell = names(expr.ave), expr = expr.ave)
    }, simplify = FALSE))
  }, simplify = FALSE)
  
  df <- gene.agg[[cluster]]
  df$module <- factor(df$module, levels = unique(df$module))
  l <- df$expr
  l[l > quantile(l, 0.99, na.rm = TRUE)] <- quantile(l, 0.99, na.rm = TRUE)
  l[l < quantile(l, 0.01, na.rm = TRUE)] <- quantile(l, 0.01, na.rm = TRUE)
  df$expr <- l
  df[, "row"] <- coord[df$cell, "row"]
  df[, "col"] <- coord[df$cell, "col"]
  
  p <- ggplot(df, aes(x = row, y = col, color = expr))
  if(background) {
    p <- p + geom_point(data = data.frame(row = coord[, "row"], col = coord[, "col"]),
                        color = "gray", size = 0.01, alpha = 0.1, stroke = 0.2)
  }
  p <- p +
    geom_point(size = 0.01, alpha = 0.8, ...) +
    facet_wrap(~module, ncol = 2) +
    coord_fixed() +
    scale_color_gradientn(colors = SpatialColors(100), name = NULL) +
    theme_minimal() +
    theme(legend.title = element_text(face = "italic"),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  p
}

SpatialColors <- colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
