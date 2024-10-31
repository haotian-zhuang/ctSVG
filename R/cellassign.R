#' Assign each square to a cell
#'
#' @param visiumpath The path to the Visium HD output folder
#' @param segmentationpath The path to the StarDist segmentation file
#' @param cellnucratio The ratio of cell size and nucleus size (2 by default)
#' @import sf
#' @return A list including a vector indicating which cell each square belongs to, and the matrix of spatial locations
#' @export
#'
cellassign <- function(visiumpath,segmentationpath,cellnucratio=2) {
  sf::sf_use_s2(FALSE)
  # read in data
  segdf <- seg <- data.table::fread(segmentationpath,data.table=F)
  segdf[,1] <- seg[,1] <- paste0('Cell:',seg[,1])
  spotcenter <- as.data.frame(arrow::read_parquet(paste0(visiumpath,'/spatial/tissue_positions.parquet')))
  spotcenter <- spotcenter[spotcenter$in_tissue==1,]
  spotrowcol <- spotcenter[,c('barcode','array_row','array_col')]
  spotcenter <- spotcenter[,c('barcode','pxl_row_in_fullres','pxl_col_in_fullres')]
  rad <- rjson::fromJSON(file=paste0(visiumpath,'/spatial/scalefactors_json.json'))$spot_diameter_fullres/2
  
  # make cell segmentation into sf
  seglist <- lapply(split(seg, seg$Cell), function(df) {
    df <- rbind(df, df[1, ])
    st_polygon(list(cbind(df$X, df$Y)))
  })
  logarea <- log2(sapply(seglist,st_area))
  cut <- mean(logarea[logarea!= - Inf])+2*sd(logarea[logarea!= - Inf])
  seglist <- seglist[logarea < cut]
  seg <- st_sf(data.frame(cell=names(seglist)), geometry = do.call(st_sfc, seglist), crs = 4326)
  segdf <- segdf[segdf[,1]%in%seg$cell,]
  
  # make spot into sf
  constmat <- rbind(c(-rad,-rad),c(-rad,rad),c(rad,rad),c(rad,-rad),c(-rad,-rad)) 
  spotlist <- apply(spotcenter[,c('pxl_row_in_fullres','pxl_col_in_fullres')],1,function(i) {
    st_polygon(list(matrix(rep(i,each=5),ncol=2)+constmat))
  })
  names(spotlist) <- spotcenter[,1]
  spot <- st_sf(data.frame(spot=names(spotlist)), geometry = do.call(st_sfc, spotlist), crs = 4326)
  
  # check intersections between cell segmentation and spot
  intersections <- st_intersects(spot,seg)
  names(intersections) <- spot$spot
  len <- sapply(intersections,length)
  intersections <- intersections[len > 0]
  intersections <- sapply(intersections,function(i) seg$cell[i])
  
  # record uniquely assigned spots
  assign <- rep(NA,nrow(spotcenter))
  names(assign) <- spotcenter$barcode
  uniinter <- unlist(intersections[names(which(len==1))])
  assign[names(uniinter)] <- uniinter
  
  # process spots assigned to multiple cells by calculating area of intersection
  multiinter <- intersections[names(which(len > 1))]
  multiinterspot <- names(multiinter)
  multiintercell <- unique(unlist(multiinter))
  
  subseg <- seg[seg$cell%in%multiintercell,]
  subspotcenter <- spotcenter[spotcenter[,1]%in%multiinterspot,]
  gridid <- expand.grid(0:10,0:10)
  overlap <- pbapply::pbsapply(1:nrow(gridid),function(i) {
    spotgrid <- subspotcenter
    spotgrid[,'pxl_row_in_fullres']  <- spotgrid[,'pxl_row_in_fullres']-rad+rad/5*gridid[i,1]
    spotgrid[,'pxl_col_in_fullres']  <- spotgrid[,'pxl_col_in_fullres']-rad+rad/5*gridid[i,2]
    spotgrid <- st_as_sf(spotgrid, coords = c("pxl_row_in_fullres", "pxl_col_in_fullres"), crs = 4326)
    spotgridwithin <- as.data.frame(st_join(spotgrid, subseg, join = st_within))
    spotgridwithin <- spotgridwithin[!is.na(spotgridwithin$cell),]
    paste0(spotgridwithin[,1],':-_-:',spotgridwithin[,2])
  },simplify = F)
  overlap <- table(unlist(overlap))
  overlap <- sort(overlap,decreasing = T)
  overlap <- data.frame(spot=sub(':-_-:.*','',names(overlap)),cell=sub('.*:-_-:','',names(overlap)),count=as.vector(overlap))
  overlap <- overlap[!duplicated(overlap[,1]),]
  assign[overlap[,1]] <- overlap[,2]
  nucassign <- assign
  
  # expand the nuclei and redo assignment
  scale_factor <- sqrt(cellnucratio)
  segcentroid <- st_centroid(seg)
  cellname <- segcentroid$cell
  segcentroid <- do.call(rbind,segcentroid$geometry)
  segcentroid <- segcentroid[match(segdf$Cell,cellname),]
  expandseg <- segcentroid + (segdf[,2:3] - segcentroid) * scale_factor
  
  expandseg <- lapply(split(expandseg, segdf$Cell), function(df) {
    df <- rbind(df, df[1, ])
    st_polygon(list(cbind(df$X, df$Y)))
  })
  #expandarea <- sapply(expandseg,st_area) # check if the polygon expand correctly
  expandseg <- st_sf(data.frame(cell=names(expandseg)), geometry = do.call(st_sfc, expandseg), crs = 4326)
  
  remainspot <- names(assign)[which(is.na(assign))]
  remainspot <- st_sf(data.frame(spot=remainspot), geometry = do.call(st_sfc, spotlist[remainspot]), crs = 4326)
  intersections <- st_intersects(remainspot,expandseg)
  names(intersections) <- remainspot$spot
  len <- sapply(intersections,length)
  intersections <- intersections[len > 0]
  intersections <- sapply(intersections,function(i) expandseg$cell[i])
  
  uniinter <- unlist(intersections[names(which(len==1))])
  assign[names(uniinter)] <- uniinter
  
  multiinter <- intersections[names(which(len > 1))]
  multiinterspot <- names(multiinter)
  multiintercell <- unique(unlist(multiinter))
  
  
  subseg <- expandseg[expandseg$cell%in%multiintercell,]
  subspotcenter <- spotcenter[spotcenter[,1]%in%multiinterspot,]
  gridid <- expand.grid(0:10,0:10)
  overlap <- pbapply::pbsapply(1:nrow(gridid),function(i) {
    spotgrid <- subspotcenter
    spotgrid[,'pxl_row_in_fullres']  <- spotgrid[,'pxl_row_in_fullres']-rad+rad/5*gridid[i,1]
    spotgrid[,'pxl_col_in_fullres']  <- spotgrid[,'pxl_col_in_fullres']-rad+rad/5*gridid[i,2]
    spotgrid <- st_as_sf(spotgrid, coords = c("pxl_row_in_fullres", "pxl_col_in_fullres"), crs = 4326)
    spotgridwithin <- as.data.frame(st_join(spotgrid, subseg, join = st_within))
    spotgridwithin <- spotgridwithin[!is.na(spotgridwithin$cell),]
    paste0(spotgridwithin[,1],':-_-:',spotgridwithin[,2])
  },simplify = F)
  overlap <- table(unlist(overlap))
  overlap <- sort(overlap,decreasing = T)
  overlap <- data.frame(spot=sub(':-_-:.*','',names(overlap)),cell=sub('.*:-_-:','',names(overlap)),count=as.vector(overlap))
  overlap <- overlap[!duplicated(overlap[,1]),]
  assign[overlap[,1]] <- overlap[,2]
  
  ### refine results
  # only keep cells that have nuclei spots
  assign <- assign[!is.na(assign)]
  nucassign <- nucassign[!is.na(nucassign)]
  assign <- assign[assign %in% nucassign]
  
  # remove cells whose covered area is less than half of total area
  coverarea <- table(nucassign)*(rad*2)^2
  coverprop <- as.vector(coverarea)/2^logarea[names(coverarea)]
  rmcell <- names(which(coverprop < 0.5))
  assign <- assign[!assign%in%rmcell]
  nucassign <- nucassign[!nucassign%in%rmcell]
  
  # check validity of cell segmentation and remove disconnected spots
  spotrowcol <- spotrowcol[spotrowcol[,1]%in%names(assign),]
  spotrowcol$cell <- assign[spotrowcol[,1]]
  spotrowcollist <- split(spotrowcol,spotrowcol$cell)
  
  rmspot <- NULL
  for (i in names(spotrowcollist)) {
    i <- spotrowcollist[[i]]
    spotname <- i[,1]
    d <- as.matrix(dist(i[,2:3]))
    e <- which(d<=sqrt(2),arr.ind=T)
    e <- data.frame(spotname[e[,1]],spotname[e[,2]])
    g <- igraph::graph_from_data_frame(e,directed = F,vertices = spotname)
    comp <- igraph::components(g)
    mem <- comp$membership
    if (comp$no > 1) {
      nucprop <- tapply(names(mem)%in%names(nucassign),list(mem),mean)
      rmspot <- c(rmspot,names(mem)[mem%in%as.numeric(names(which(nucprop==0)))])
    }
  }
  
  assign[rmspot] <- NA
  assign <- assign[!is.na(assign)]
  
  centroid <- st_centroid(seg)
  coord <- do.call(rbind,centroid$geometry)
  coord <- data.frame(cell=centroid$cell,X=coord[,1],Y=coord[,2])
  coord <- coord[match(unique(assign),coord[,1]),]
  coord <- as.matrix(coord[,2:3])
  rownames(coord) <- unique(assign)
  colnames(coord) <- c('row','col')
  
  list(assign = assign, coord = coord)
}
