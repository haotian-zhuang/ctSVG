#' Identify spatially variable genes within each cell cluster
#'
#' @param d The output from \code{\link{runseurat}}
#' @param coord The matrix of spatial locations from \code{\link{cellassign}}
#' @param cell.filter Whether to remove spatially isolated cells within each cell cluster (\code{TRUE} by default)
#' @param clu.list The list of reassigned clusters from \code{\link{recluster}}
#' @param jac.mat The matrix of Jaccard indices from \code{\link{recluster}}
#' @param min.cells The minimum proportion of cells within each cell cluster a gene is expressed in to be retained (0.01 by default)
#' @param seed Random seed for permutation (2024 by default)
#' @param verbose Whether to print the cell cluster in progress (\code{FALSE} by default)
#' @param n.permute Number of permutation times (reassigned clusters with the highest Jaccard indices, 100 by default)
#' @import Seurat
#' @return A data frame with p-values and FDRs for each gene in each cell cluster
#' @export
#'
ctsvg_test <- function(d, coord, cell.filter = TRUE, clu.list = NULL, jac.mat = NULL, min.cells = 0.01, seed = 2024, verbose = FALSE, n.permute = 100) {
  expr <- d@assays$RNA$data
  if(!identical(colnames(expr), rownames(coord))) { coord <- coord[intersect(colnames(expr), rownames(coord)), ] }
  
  sc <- Idents(d)
  sccell <- split(names(sc), sc)
  if(cell.filter) { sccell <- lapply(sccell, function(k) { k[!detectOutlier(coord = coord[k, ])] }) }
  
  if(is.null(clu.list)) {
    clutestRes <- do.call(rbind, sapply(names(sccell), function(k) {
      i <- sccell[[k]]
      expr_clu <- expr[, i]
      res <- spatialTest(expr = expr_clu[rowMeans(expr_clu > 0) > min.cells, ], coord = coord[i, ], coord_permute = NULL, knot = FALSE)
      data.frame(cluster = k, gene = rownames(res), res)
    }, simplify = FALSE))
  } else {
    set.seed(seed = seed)
    clutestRes <- do.call(rbind, sapply(names(sccell), function(k) {
      if(verbose) { print(k) }
      top.index <- order(jac.mat[k, ], decreasing = TRUE)[1:n.permute]
      cell.index <- clu.list[[k]][top.index]
      
      i <- sccell[[k]]
      coord_permute <- lapply(cell.index, function(j) {
        coord.tmp <- coord[j, ]
        rownames(coord.tmp) <- sample(rownames(coord.tmp), replace = FALSE)
        coord.tmp
      })
      
      expr_clu <- expr[, i]
      res <- spatialTest(expr = expr[rowMeans(expr_clu > 0) > min.cells, ], coord = coord[i, ], coord_permute = coord_permute, knot = FALSE)
      data.frame(cluster = k, gene = rownames(res), res)
    }, simplify = FALSE))
  }
  return(clutestRes)
}

#' Fit spatial gene expression patterns within each cell cluster
#'
#' @param d The output from \code{\link{runseurat}}
#' @param coord The matrix of spatial locations from \code{\link{cellassign}}
#' @param cell.filter Whether to remove spatially isolated cells within each cell cluster (\code{TRUE} by default)
#' @import Seurat
#' @return The fitted cell-cluster-specific gene expression matrix
#' @export
#'
ctsvg_fit <- function(d, coord, cell.filter = TRUE) {
  expr <- d@assays$RNA$data
  if(!identical(colnames(expr), rownames(coord))) { coord <- coord[intersect(colnames(expr), rownames(coord)), ] }
  
  sc <- Idents(d)
  sccell <- split(names(sc), sc)
  if(cell.filter) { sccell <- lapply(sccell, function(k) { k[!detectOutlier(coord = coord[k, ])] }) }
  
  clufitRes <- do.call(cbind, sapply(names(sccell), function(k) {
    i <- sccell[[k]]
    spatialFit(expr = expr[, i], coord = coord[i, ], knot = FALSE)
  }, simplify = FALSE))
  return(clufitRes)
}

detectOutlier <- function(coord, prop = 0.01, min.cells = 10, max.cells = 50, slope = 6) {
  distcell <- as.matrix(stats::dist(coord))
  score <- apply(distcell, 1, function(i) mean(head(sort(i), min(max(prop * ncol(distcell), min.cells), max.cells))))
  return(score > mean(score) + slope * sd(score))
}

spatialTest <- function(expr, coord, coord_permute = NULL, knot = F, maxknotallowed = 5) {
  
  if(is.null(coord_permute)) {
    return(spatialTestFixed(expr = expr, coord = coord,
                            knot = knot, maxknotallowed = maxknotallowed))
  } else {
    res.ori <- spatialTestFixed(expr = expr, coord = coord,
                                knot = knot, maxknotallowed = maxknotallowed)
    fstat.ori <- setNames(res.ori[, "fstat"], nm = rownames(res.ori))
    
    fstat.perm <- sapply(coord_permute, function(i) {
      res.perm <- spatialTestFixed(expr = expr, coord = i,
                                   knot = knot, maxknotallowed = maxknotallowed)
      return(setNames(res.perm[, "fstat"], nm = rownames(res.perm))[names(fstat.ori)])
    })
    
    pval.empirical <- (Matrix::rowSums(fstat.perm >= fstat.ori) + 1)/(ncol(fstat.perm) + 1)
    
    pval.parametric <- sapply(1:nrow(fstat.perm), function(i) {
      if(fstat.ori[i] == 0) { return(setNames(1, nm = names(fstat.ori)[i])) }
      
      if(!'try-error'%in%class(try(capture.output(fitdistrplus::fitdist(fstat.perm[i, ], 'gamma')), silent = T))) {
        fit.gamma <- fitdistrplus::fitdist(fstat.perm[i, ], 'gamma')
        return(stats::pgamma(fstat.ori[i], shape = fit.gamma$estimate[1], rate = fit.gamma$estimate[2], lower.tail = F))
      } else {
        fit.gamma <- fitdistrplus::fitdist(fstat.perm[i, ], 'gamma', method = "mme")
        return(stats::pgamma(fstat.ori[i], shape = fit.gamma$estimate[1], rate = fit.gamma$estimate[2], lower.tail = F))
      }
    })
    
    fdr.empirical <- stats::p.adjust(pval.empirical, method = 'fdr')
    fdr.parametric <- stats::p.adjust(pval.parametric, method = 'fdr')
    res <- data.frame(fdr = fdr.parametric, pval = pval.parametric,
                      fdr.empirical = fdr.empirical, pval.empirical = pval.empirical,
                      fstat.ori = fstat.ori)
    res <- res[order(res$pval, res$pval.empirical, -res$fstat.ori), ]
    return(res)
  }
}

spatialTestFixed <- function(expr, coord, knot = F, maxknotallowed = 5) {
  
  expr <- expr[, rownames(coord), drop = F]
  if(knot == F) {
    
    knotnum <- rep(0, nrow(expr))
    names(knotnum) <- rownames(expr)
    xrow <- cbind(1, splines::bs(coord[,'row'], intercept = F, df = 3))
    ycol <- cbind(1, splines::bs(coord[,'col'], intercept = F, df = 3))
    
    B <- xrow[,rep(1:ncol(xrow), ncol(xrow))] * ycol[,rep(1:ncol(ycol), each = ncol(ycol))]
    B <- B[, which(matrixStats::colSds(B)>0), drop = F]
    B <- cbind(1, B)
    rownames(B) <- colnames(expr)
    colnames(B) <- NULL
    tBB <- crossprod(B)
    rownames(tBB) <- colnames(tBB) <- NULL
    
    pred <- as.matrix(expr %*% B %*% tcrossprod(chol2inv(chol(tBB)), B))
    
    SSE <- Matrix::rowSums((expr - pred)^2)
    SST <- Matrix::rowSums(sweep(expr, 1, Matrix::rowMeans(expr), FUN = '-')^2)
    
    if (any(SST<SSE)) print(names(which(SST<SSE)))
    
    fstat <- ((SST - SSE)/(ncol(B) - 1))/(SSE/(nrow(B) - ncol(B)))
    fstat[which(Matrix::rowSums(expr) == 0)] <- 0
    
    pval <- stats::pf(q = fstat, df1 = ncol(B) - 1, df2 = nrow(B) - ncol(B), lower.tail = F)
    logpval <- stats::pf(q = fstat, df1 = ncol(B) - 1, df2 = nrow(B) - ncol(B), lower.tail = F, log.p = T)
    
    res <- data.frame(fstat = fstat, pval = pval, logpval = logpval)
    
  } else {
    
    knotnum0 <- 0:maxknotallowed
    names(knotnum0) <- knotnum0
    
    Blist <- lapply(knotnum0, function(numknot) {
      
      xrow <- cbind(1, splines::bs(coord[,'row'], intercept = F, df = numknot+3))
      ycol <- cbind(1, splines::bs(coord[,'col'], intercept = F, df = numknot+3))
      
      B <- xrow[,rep(1:ncol(xrow), ncol(xrow))] * ycol[,rep(1:ncol(ycol), each = ncol(ycol))]
      B <- B[, which(matrixStats::colSds(B)>0), drop = F]
      B <- cbind(1, B)
      rownames(B) <- colnames(expr)
      colnames(B) <- NULL
      tBB <- crossprod(B)
      rownames(tBB) <- colnames(tBB) <- NULL
      list(B = B, tBB = tBB)
    })
    names(Blist) <- as.character(knotnum0)
    
    testpos <- sapply(knotnum0, function(numknot) {
      
      tBB <- Blist[[as.character(numknot)]][['tBB']]
      !'try-error'%in%class(try(chol(tBB), silent = T))
      #matrixcalc::is.positive.definite(tBB)
    })
    
    if(mean(testpos) != 1) {
      maxknot <- which(testpos == F)[1] - 2
      knotnum0 <- 0:maxknot
      names(knotnum0) <- knotnum0 }
    
    expr <- Matrix::t(expr)
    
    bic <- sapply(knotnum0, Calbic, Blist = Blist, expr = expr)
    
    knotnum <- knotnum0[apply(bic, 1, which.min)]
    names(knotnum) <- rownames(bic)
    
    res <- lapply(unique(knotnum), function(k) {
      
      B <- Blist[[as.character(k)]][['B']]
      tBB <- Blist[[as.character(k)]][['tBB']]
      
      beta <- as.matrix(tcrossprod(chol2inv(chol(tBB)), B) %*% expr[, which(knotnum == k), drop = F])
      pred <- B %*% beta
      
      expr.sub <- expr[, colnames(pred), drop = F]
      
      SSE <- Matrix::colSums((expr.sub - pred)^2)
      SST <- Matrix::colSums(sweep(expr.sub, 2, Matrix::colMeans(expr.sub), FUN = '-')^2)
      
      if (any(SST<SSE)) print(names(which(SST<SSE)))
      
      fstat <- ((SST - SSE)/(ncol(B) - 1))/(SSE/(nrow(B) - ncol(B)))
      fstat[which(Matrix::colSums(expr.sub) == 0)] <- 0
      
      pval <- stats::pf(q = fstat, df1 = ncol(B) - 1, df2 = nrow(B) - ncol(B), lower.tail = F)
      logpval <- stats::pf(q = fstat, df1 = ncol(B) - 1, df2 = nrow(B) - ncol(B), lower.tail = F, log.p = T)
      
      data.frame(fstat = fstat, pval = pval, logpval = logpval)
    })
    
    res <- do.call(rbind, res)
    res <- res[colnames(expr), ]
  }
  
  res$fdr <- stats::p.adjust(res$pval, method = 'fdr')
  res$knotnum <- knotnum
  res <- res[, c("fdr", "logpval", "pval", "fstat", "knotnum")]
  res <- res[order(res$logpval, -res$fstat), ]
  return(res)
}

spatialFit <- function(expr, coord, knot = F, maxknotallowed = 5) {
  
  expr <- expr[, rownames(coord), drop = F]
  if(knot == F) {
    
    knotnum <- rep(0, nrow(expr))
    names(knotnum) <- rownames(expr)
    xrow <- cbind(1, splines::bs(coord[,'row'], intercept = F, df = 3))
    ycol <- cbind(1, splines::bs(coord[,'col'], intercept = F, df = 3))
    
    B <- xrow[,rep(1:ncol(xrow), ncol(xrow))] * ycol[,rep(1:ncol(ycol), each = ncol(ycol))]
    B <- B[, which(matrixStats::colSds(B)>0), drop = F]
    B <- cbind(1, B)
    rownames(B) <- colnames(expr)
    colnames(B) <- NULL
    tBB <- crossprod(B)
    rownames(tBB) <- colnames(tBB) <- NULL
    
    pred <- as.matrix(expr %*% B %*% tcrossprod(chol2inv(chol(tBB)), B))
    return(pred)
  } else {
    
    knotnum0 <- 0:maxknotallowed
    names(knotnum0) <- knotnum0
    
    Blist <- lapply(knotnum0, function(numknot) {
      
      xrow <- cbind(1, splines::bs(coord[,'row'], intercept = F, df = numknot+3))
      ycol <- cbind(1, splines::bs(coord[,'col'], intercept = F, df = numknot+3))
      
      B <- xrow[,rep(1:ncol(xrow), ncol(xrow))] * ycol[,rep(1:ncol(ycol), each = ncol(ycol))]
      B <- B[, which(matrixStats::colSds(B)>0), drop = F]
      B <- cbind(1, B)
      rownames(B) <- colnames(expr)
      colnames(B) <- NULL
      tBB <- crossprod(B)
      rownames(tBB) <- colnames(tBB) <- NULL
      list(B = B, tBB = tBB)
    })
    names(Blist) <- as.character(knotnum0)
    
    testpos <- sapply(knotnum0, function(numknot) {
      
      tBB <- Blist[[as.character(numknot)]][['tBB']]
      !'try-error'%in%class(try(chol(tBB), silent = T))
      #matrixcalc::is.positive.definite(tBB)
    })
    
    if(mean(testpos) != 1) {
      maxknot <- which(testpos == F)[1] - 2
      knotnum0 <- 0:maxknot
      names(knotnum0) <- knotnum0 }
    
    expr <- Matrix::t(expr)
    
    bic <- sapply(knotnum0, Calbic, Blist = Blist, expr = expr)
    
    knotnum <- knotnum0[apply(bic, 1, which.min)]
    names(knotnum) <- rownames(bic)
    print(table(knotnum))
    
    pred <- lapply(unique(knotnum), function(k) {
      
      B <- Blist[[as.character(k)]][['B']]
      tBB <- Blist[[as.character(k)]][['tBB']]
      
      beta <- as.matrix(tcrossprod(chol2inv(chol(tBB)), B) %*% expr[, which(knotnum == k), drop = F])
      pred <- B %*% beta
      pred
    })
    
    pred <- do.call(cbind, pred)
    pred <- pred[, colnames(expr)]
    pred <- Matrix::t(pred)
    return(pred)
  }
}

Calbic <- function(numknot, Blist, expr) {
  
  B <- Blist[[as.character(numknot)]][['B']]
  tBB <- Blist[[as.character(numknot)]][['tBB']]
  
  beta <- as.matrix(tcrossprod(chol2inv(chol(tBB)), B) %*% expr)
  pred <- B %*% beta
  mse <- Matrix::colMeans((expr - pred)^2)
  bic <- nrow(B)*(1+log(2*pi)+log(mse)) + log(nrow(B))*(ncol(B)+1) # stats::AIC(lm(), k = log(nrow(B)))
  return(bic)
}
