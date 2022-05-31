try(suppressMessages(library(dplyr)))
try(suppressMessages(library(dbplyr)))
try(suppressMessages(library(plyr)))
try(suppressMessages(library(scales)))
try(suppressMessages(library(ComplexHeatmap)))
try(suppressMessages(library(circlize)))
try(suppressMessages(library(stringr)))
try(suppressMessages(library(RColorBrewer)))
try(suppressMessages(library(car)))
try(suppressMessages(library(igraph)))
try(suppressMessages(library(MSigDB)))
try(suppressMessages(library(GenomicRanges)))
try(suppressMessages(library(biomaRt)))
try(suppressMessages(library(caTools)))
try(suppressMessages(library(Matrix)))
try(suppressMessages(library(AUCell)))
try(suppressMessages(library(GENIE3)))
try(suppressMessages(library(RcisTarget)))
try(suppressMessages(library(reshape2)))
try(suppressMessages(library(visNetwork)))
try(suppressMessages(library(shiny)))
try(suppressMessages(library(rsvd)))
try(suppressMessages(library(ggplot2)))
try(suppressMessages(library(rjags)))
try(suppressMessages(library(infercnv)))
try(suppressMessages(library(Seurat)))
try(suppressMessages(library(reticulate)))
try(suppressMessages(library(mHG)))
try(suppressMessages(library(irlba)))
try(suppressMessages(library(NMF)))
try(suppressMessages(library(VennDiagram)))
try(suppressMessages(library(tsne)))
try(suppressMessages(library(SingleR)))
try(suppressMessages(library(SCENIC)))
try(suppressMessages(library(parallel)))
try(suppressMessages(library(ggsignif)))
try(suppressMessages(library(tidyverse)))
try(suppressMessages(library(nichenetr)))
try(suppressMessages(library(cowplot)))
try(suppressMessages(library(ggpubr)))
try(suppressMessages(library(DiagrammeR)))
try(suppressMessages(library(SeuratDisk)))
try(suppressMessages(library(SeuratWrappers)))
try(suppressMessages(library(nnls)))
try(suppressMessages(library(colormap)))
try(suppressMessages(library(genesorteR)))
try(suppressMessages(library(dplyr)))
try(suppressMessages(library(survival)))
try(suppressMessages(library(survminer)))
try(suppressMessages(library(tidyr)))
try(suppressMessages(library(caret)))


#### Miscellaneous ####

# Source: Seurat utilities
'%||%' = function (lhs, rhs){
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

# To load an RData object and specify name
# Source: https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
loadRData = function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#' Plots a series of barplots and connects them 
#' Modified from https://stackoverflow.com/questions/22560850/barplot-with-connected-series
#' 
#' @param dat NxM matrix with N rows as features and M columns as samples
#' @param color Vector of N colors
#' @param space Space between barplots
#' @param alpha Alpha for area connecting barplots
#' 
#' @examples
#' dat <- matrix(rnorm(100),10,10)
#' dat <- abs(matrix(rnorm(100),10,10)) 
#' connectedBarplot(dat, color=rainbow(nrow(dat)))
#'
connectedBarplot <- function(dat, color=rainbow(nrow(dat)), space=1, alpha=0.5, ...) {  
  b <- barplot(dat, col=color, space = space, ...)                     
  
  for (i in seq_len(ncol(dat) - 1)) {     
    lines(c(b[i]+0.5, b[i+1]-0.5), c(0, 0)) ## bottom line       
    
    for (j in seq_len(nrow(dat))) {     
      if (j == 1) {                   
        lines(c(b[i]+0.5, b[i+1]-0.5), c(dat[j,i], dat[j,i+1]))                       
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(0, dat[j,i], dat[j,i+1], 0),               
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
      if (j == 2) {                   
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))                      
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(dat[1,i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], dat[1,i+1]),                       
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
      if (j > 2) {                    
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))                      
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(colSums(dat[1:(j-1),])[i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], colSums(dat[1:(j-1),])[i+1]),              
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
    }          
  }              
}   

# To summarize multiple matrices of same dimensions element-wise
SummarizeMatrices = function(
  list, 
  fun,
  ...){
  a = nrow(list[[1]])
  b = ncol(list[[2]])
  mat = t(sapply(1:a, function(i){
    sapply(1:b, function(j){
      x = sapply(list, '[', i, j)
      return(fun(x, ...))
    })
  }))
  rownames(mat) = rownames(list[[1]])
  colnames(mat) = colnames(list[[2]])
  return(mat)
}

# To identify clusters using a flexible cutoff (cut the tree at several places, then use correlation and length thresholds)
# Modified from Neftel et al. identification of differentially expressed genes
FindSets = function(
  data,
  stat = median,
  is.cor = FALSE,
  k.range = 2:50,
  c.thresh_low = 0,
  c.thresh_high = 1,
  l.thresh_low = 0,
  l.thresh_high = 1,
  do.reduce.overlap = TRUE,
  do.sort = TRUE,
  clustering_distance_columns = 'pearson', 
  clustering_distance_rows = 'pearson',
  ...
){
  if (!is.cor){
    dend = hclust(as.dist(1-cor(data)))
  }
  if (is.cor){
    dend = hclust(data)
  }
  
  if (l.thresh_high <= 1){
    l.thresh_high = l.thresh_high*dim(data)[2]
  }
  if (l.thresh_low < 1){
    l.thresh_low = l.thresh_low*dim(data)[2]
  }
  
  k.range = k.range[k.range <= dim(data)[2]]
  
  sets = unlist(lapply(c(k.range), function(k){
    groupings = as.factor(cutree(as.hclust(dend), k = k))
    sets_k = lapply(levels(groupings), function(x){
      names(groupings)[groupings == x]
    })
    sets_k = sets_k[sapply(sets_k, function(x){
      l = length(x)
      if (l == 1){
        return(l >= l.thresh_low)
      } else {
        if (is.cor){
          c = data[x,x]
        } else {
          c = cor(data[,x], method = 'pearson')
        }
        diag(c) = NA
        m = stat(c, na.rm = TRUE)
        return(m >= c.thresh_low & m <= c.thresh_high & l >= l.thresh_low & l <= l.thresh_high)
      }
    })]
    return(sets_k)
  }), recursive = FALSE)
  
  if (do.reduce.overlap){
    i = 1
    while (i < length(sets)){
      set1 = sets[[i]]
      o = sapply(sets, function(set2){
        length(intersect(set1, set2)) > 0
        #length(intersect(set1, set2))/sqrt(length(set1)*length(set2)) > 0.1
        #phyper(length(intersect(set1, set2)), length(set1), length(set1) - length(intersect(set1, set2)), length(set2), lower.tail = FALSE) < 0.01
      })
      sets[[i]] = Reduce(union, sets[o])
      o[i] = FALSE
      sets = sets[!o]
      i = i+1
    }
    names(sets) = sapply(sets, '[', 1)
  }
  
  if (do.sort){
    sets = sets[order(sapply(sets, function(x){
      if (is.cor){
        c = data[x,x]
      } else {
        c = cor(data[,x], method = 'pearson')
      }
      diag(c) = NA
      m = stat(c, na.rm = TRUE)
    }), decreasing = TRUE)]
  }
  
  return(sets)
}

# Modify ComplexHeatmap defaults
Heatmap = function(
  matrix,
  colors = NULL,
  breaks = NULL,
  col = NULL,
  clustering_method_rows = 'ward.D2', 
  clustering_method_columns = 'ward.D2',
  clustering_distance_rows = 'pearson',
  clustering_distance_columns = 'pearson',
  show_row_names = FALSE,
  show_column_names = FALSE,
  is.symmetric = FALSE,
  ...
){
  if (is.null(breaks)){
    mi = min(matrix, na.rm = TRUE)
    ma = max(matrix, na.rm = TRUE)
    if (mi < 0 & ma > 0){
      breaks = c(mi, 0, ma)
    } else {
      breaks = c(mi, (mi + ma)/2, ma)
    }
  }
  if (is.null(colors)){
    colors = c('blue','white','red')
  }
  if (is.null(col)){
    col = colorRamp2(breaks = breaks, colors = colors)
  }
  if (is.symmetric == TRUE){
    h = ComplexHeatmap::Heatmap(matrix, 
                                col = col,
                                clustering_method_rows = clustering_method_rows, 
                                clustering_distance_rows = clustering_distance_rows,
                                cluster_columns = FALSE,
                                show_row_names = show_row_names,
                                show_column_names = show_column_names,
                                ...)
    o = unlist(row_order(h))
    return(ComplexHeatmap::Heatmap(matrix, 
                                   col = col,
                                   clustering_method_rows = clustering_method_rows, 
                                   clustering_distance_rows = clustering_distance_rows,
                                   cluster_columns = FALSE,
                                   column_order = o,
                                   show_row_names = show_row_names,
                                   show_column_names = show_column_names, 
                                   ...))
  } else {
    return(ComplexHeatmap::Heatmap(matrix, 
                                   col = col,
                                   clustering_method_rows = clustering_method_rows, 
                                   clustering_method_columns = clustering_method_columns,
                                   clustering_distance_rows = clustering_distance_rows,
                                   clustering_distance_columns = clustering_distance_columns,
                                   show_row_names = show_row_names,
                                   show_column_names = show_column_names,
                                   ...))
  }
  
}

PCAPlot3d = function(
  srt, 
  group.by = 'ident',
  color.scheme = NULL,
  dims = c(1,2,3),
  ...
){
  embeddings = Embeddings(srt, 'pca')
  PC1 = embeddings[,dims[1]]
  PC2 = embeddings[,dims[2]]
  PC3 = embeddings[,dims[3]]
  groups = factor(FetchData(srt, group.by)[,1])
  if (is.null(color.scheme)){
    color.scheme = hue_pal()(nlevels(groups))
  }
  point.col = color.scheme[groups]
  print(unique(point.col))
  scatter3d(PC1, PC2, PC3, surface = FALSE, point.col = point.col, ...)
}

#### SCENIC modifications

# Modified from SCENIC
runCorrelation <- function(
  exprMat_filtered,
  scenicOptions
){
  corrMat <- cor(t(exprMat_filtered), method="spearman")
  saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))
}

# Modified from SCENIC
runSCENIC_1_coexNetwork2modules <- function (scenicOptions){
  linkList <- loadInt(scenicOptions, "genie3ll")
  if (!all(colnames(linkList)[1:3] == c("TF", "Target", "weight"))) 
    stop("The link list colnames should be \"TF\", \"Target\", \"weight\"")
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tCreating TF modules")
  if (getSettings(scenicOptions, "verbose")) 
    message(msg)
  quantile(linkList$weight, probs = c(0.75, 0.9))
  SCENIC:::.openDev(fileName = getIntName(scenicOptions, "genie3weighPlot"), 
                    devType = getSettings(scenicOptions, "devType"))
  plot(linkList$weight[1:1e+06], type = "l", ylim = c(0, max(linkList$weight)), 
       main = "Weight of the links", ylab = "Weight", xlab = "Links sorted decreasingly")
  abline(h = 0.001, col = "blue")
  dev.off()
  linkList_001 <- linkList[which(linkList[, "weight"] > getSettings(scenicOptions, 
                                                                    "modules/weightThreshold")), ]
  if (getSettings(scenicOptions, "verbose")) 
    message("Number of links between TFs and targets: ", 
            nrow(linkList_001))
  tfModules <- list()
  linkList_001$TF <- as.character(linkList_001$TF)
  linkList_001$Target <- as.character(linkList_001$Target)
  tfModules[["w001"]] <- split(linkList_001$Target, factor(linkList_001$TF))
  llminW <- linkList_001[which(linkList_001[, "weight"] > 0.005), 
                         ]
  tfModules[["w005"]] <- split(llminW$Target, factor(llminW$TF))
  tfModules[["top50"]] <- lapply(tfModules[["w001"]], function(x) x[1:(min(length(x), 
                                                                           50))])
  linkList_001_byTarget <- split(linkList_001, factor(linkList_001$Target))
  nTopTfs <- c(5, 10, 50)
  nTopTfs <- setNames(nTopTfs, paste("top", nTopTfs, "perTarget", 
                                     sep = ""))
  topTFsperTarget <- lapply(linkList_001_byTarget, function(llt) {
    nTFs <- nTopTfs[which(nTopTfs <= nrow(llt))]
    reshape2::melt(lapply(nTFs, function(x) llt[1:x, "TF"]))
  })
  topTFsperTarget <- topTFsperTarget[which(!sapply(sapply(topTFsperTarget, 
                                                          nrow), is.null))]
  topTFsperTarget.asDf <- data.frame(data.table::rbindlist(topTFsperTarget, 
                                                           idcol = TRUE))
  colnames(topTFsperTarget.asDf) <- c("Target", "TF", "method")
  head(topTFsperTarget.asDf)
  tfModules.melted <- reshape2::melt(tfModules)
  colnames(tfModules.melted) <- c("Target", "TF", "method")
  tfModules <- rbind(tfModules.melted, topTFsperTarget.asDf)
  rm(tfModules.melted)
  if (getSettings(scenicOptions, "verbose")) 
    print(rbind(nTFs = length(unique(tfModules$TF)), nTargets = length(unique(tfModules$Target)), 
                nGeneSets = nrow(unique(tfModules[, c("TF", "method")])), 
                nLinks = nrow(tfModules)))
  corrMat <- loadInt(scenicOptions, "corrMat")
  tfs <- as.character(unique(tfModules$TF))
  missingTFs <- tfs[which(!tfs %in% rownames(corrMat))]
  if (length(missingTFs) > 0) {
    warning("The following TFs are missing from the correlation matrix: ", 
            paste(missingTFs, collapse = ", "))
    tfs <- tfs[which(tfs %in% rownames(corrMat))]
    corrMat <- corrMat[tfs, ]
  }
  tfModules_byTF <- split(tfModules, factor(tfModules$TF))
  tfModules_withCorr_byTF <- lapply(tfModules_byTF[tfs], function(tfGeneSets) {
    tf <- as.character(unique(tfGeneSets$TF))
    targets <- as.character(tfGeneSets$Target)
    targets = targets[which(targets %in% rownames(corrMat))]
    cbind(tfGeneSets[match(targets, tfGeneSets$Target), ], 
          corr = c(as.numeric(corrMat[tf, targets] > 0.03) - as.numeric(corrMat[tf, targets] < -0.03)))
  })
  tfModules_withCorr <- data.frame(data.table::rbindlist(tfModules_withCorr_byTF))
  if (length(missingTFs) > 0) {
    tfModules_withCorr <- rbind(tfModules_withCorr, data.frame(tfModules[tfModules$TF %in% 
                                                                           missingTFs, ], corr = NA))
  }
  saveRDS(tfModules_withCorr, file = getIntName(scenicOptions, 
                                                "tfModules_asDF"))
}


#### Seurat modifications ####

# Modified from Seurat
AddModuleScore = function (object, features, pool = NULL, nbin = 24, ctrl = 100, 
                           k = FALSE, assay = NULL, name = "Cluster", seed = 1) 
{
  set.seed(seed = seed)
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetData(object)
  features.old <- features
  if (k) {
    .NotYetUsed(arg = "k")
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == 
                                         i))
    }
    cluster.length <- length(x = features)
  }
  else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
      return(intersect(x = x, y = rownames(x = object)))
    })
    cluster.length <- length(x = features)
  }
  if (!all(Seurat:::LengthCheck(values = features))) {
    warning(paste("Could not find enough features in the object from the following feature lists:", 
                  paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))), 
                  "Attempting to match case..."))
    features <- lapply(X = features.old, FUN = CaseMatch, 
                       match = rownames(x = object))
  }
  if (!all(Seurat:::LengthCheck(values = features))) {
    stop(paste("The following feature lists do not have enough features present in the object:", 
               paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))), 
               "exiting..."))
  }
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                         n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                                                                              data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
                        ncol = ncol(x = object))
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, 
                                                        ])
  }
  features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
                            ncol = ncol(x = object))
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- name
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  Seurat:::CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}

# Modified from Seurat
merge2 = function (x = NULL, y = NULL, add.cell.ids = NULL, merge.data = TRUE, 
                   merge.dr = NULL, project = "SeuratProject", na.rm = FALSE, ...){
  objects <- c(x, y)
  if (!is.null(x = add.cell.ids)) {
    if (length(x = add.cell.ids) != length(x = objects)) {
      stop("Please provide a cell identifier for each object provided to merge")
    }
    for (i in 1:length(x = objects)) {
      objects[[i]] <- RenameCells(object = objects[[i]], 
                                  add.cell.id = add.cell.ids[i])
    }
  }
  objects <- Seurat:::CheckDuplicateCellNames(object.list = objects)
  assays <- lapply(X = objects, FUN = Seurat:::FilterObjects, classes.keep = "Assay")
  fake.feature <- Seurat:::RandomName(length = 17)
  assays <- unique(x = unlist(x = assays, use.names = FALSE))
  combined.assays <- vector(mode = "list", length = length(x = assays))
  names(x = combined.assays) <- assays
  for (assay in assays) {
    assays.merge <- lapply(X = objects, FUN = function(object) {
      return(tryCatch(expr = object[[assay]], error = function(e) {
        return(CreateAssayObject(counts = Matrix(data = 0, 
                                                 ncol = ncol(x = object), dimnames = list(fake.feature, 
                                                                                          colnames(x = object)), sparse = TRUE)))
      }))
    })
    if (all(Seurat:::IsSCT(assay = assays.merge))) {
      scaled.features <- unique(x = unlist(x = lapply(X = assays.merge, 
                                                      FUN = function(x) rownames(x = GetAssayData(object = x, 
                                                                                                  slot = "scale.data")))))
      for (ob in 1:length(x = objects)) {
        if (assay %in% Seurat:::FilterObjects(object = objects[[ob]], 
                                              classes.keep = "Assay")) {
          objects[[ob]] <- suppressWarnings(GetResidual(object = objects[[ob]], 
                                                        features = scaled.features, assay = assay, 
                                                        verbose = FALSE))
          assays.merge[[ob]] <- objects[[ob]][[assay]]
        }
      }
      scaled.features <- names(x = which(x = table(x = unlist(x = lapply(X = assays.merge, 
                                                                         FUN = function(x) rownames(x = GetAssayData(object = x, 
                                                                                                                     slot = "scale.data"))))) == length(x = assays.merge)))
      if (length(x = scaled.features) > 0) {
        for (a in 1:length(x = assays.merge)) {
          assays.merge[[a]] <- SetAssayData(object = assays.merge[[a]], 
                                            slot = "scale.data", new.data = GetAssayData(object = assays.merge[[a]], 
                                                                                         slot = "scale.data")[scaled.features, ])
        }
      }
      else {
        for (a in 1:length(x = assays.merge)) {
          assays.merge[[a]] <- SetAssayData(object = assays.merge[[a]], 
                                            slot = "scale.data", new.data = new(Class = "matrix"))
        }
      }
    }
    merged.assay <- merge2.Assay(x = assays.merge[[1]], y = assays.merge[2:length(x = assays.merge)], 
                                 merge.data = merge.data, na.rm = na.rm, ...)
    merged.assay <- subset(x = merged.assay, features = rownames(x = merged.assay)[rownames(x = merged.assay) != 
                                                                                     fake.feature])
    if (length(x = Key(object = merged.assay)) == 0) {
      Key(object = merged.assay) <- paste0(assay, "_")
    }
    combined.assays[[assay]] <- merged.assay
  }
  combined.meta.data <- data.frame(row.names = colnames(x = combined.assays[[1]]))
  new.idents <- c()
  for (object in objects) {
    old.meta.data <- object[[]]
    if (any(!colnames(x = old.meta.data) %in% colnames(x = combined.meta.data))) {
      cols.to.add <- colnames(x = old.meta.data)[!colnames(x = old.meta.data) %in% 
                                                   colnames(x = combined.meta.data)]
      combined.meta.data[, cols.to.add] <- NA
    }
    i <- sapply(X = old.meta.data, FUN = is.factor)
    old.meta.data[i] <- lapply(X = old.meta.data[i], FUN = as.vector)
    combined.meta.data[rownames(x = old.meta.data), colnames(x = old.meta.data)] <- old.meta.data
    new.idents <- c(new.idents, as.vector(Idents(object = object)))
  }
  names(x = new.idents) <- rownames(x = combined.meta.data)
  new.idents <- factor(x = new.idents)
  if (DefaultAssay(object = x) %in% assays) {
    new.default.assay <- DefaultAssay(object = x)
  }
  else if (DefaultAssay(object = y) %in% assays) {
    new.default.assay <- DefaultAssay(object = y)
  }
  else {
    new.default.assay <- assays[1]
  }
  combined.images <- vector(mode = "list", length = length(x = unlist(x = lapply(X = objects, 
                                                                                 FUN = Images))))
  index <- 1L
  for (i in 1:length(x = objects)) {
    object <- objects[[i]]
    for (image in Images(object = object)) {
      image.obj <- object[[image]]
      if (image %in% names(x = combined.images)) {
        image <- if (is.null(x = add.cell.ids)) {
          make.unique(names = c(na.omit(object = names(x = combined.images)), 
                                image))[index]
        }
        else {
          paste(image, add.cell.ids[i], sep = "_")
        }
      }
      combined.images[[index]] <- image.obj
      names(x = combined.images)[index] <- image
      index <- index + 1L
    }
  }
  combined.reductions <- list()
  if (!is.null(x = merge.dr)) {
    for (dr in merge.dr) {
      drs.to.merge <- list()
      for (i in 1:length(x = objects)) {
        if (!dr %in% Reductions(object = objects[[i]])) {
          warning("The DimReduc ", dr, " is not present in all objects being ", 
                  "merged. Skipping and continuing.", call. = FALSE, 
                  immediate. = TRUE)
          break
        }
        drs.to.merge[[i]] <- objects[[i]][[dr]]
      }
      if (length(x = drs.to.merge) == length(x = objects)) {
        combined.reductions[[dr]] <- merge(x = drs.to.merge[[1]], 
                                           y = drs.to.merge[2:length(x = drs.to.merge)])
      }
    }
  }
  merged.object <- new(Class = "Seurat", assays = combined.assays, 
                       reductions = combined.reductions, images = combined.images, 
                       meta.data = combined.meta.data, active.assay = new.default.assay, 
                       active.ident = new.idents, project.name = project, version = packageVersion(pkg = "SeuratObject"))
  VariableFeatures(merged.object) = union(VariableFeatures(x), VariableFeatures(y))
  merged.object@meta.data[] = lapply(merged.object@meta.data, function(i){
    if (is.character(i)){ i = as.factor(i) }
    return(i)
  })
  return(merged.object)
}

# Modified from Seurat
merge2.Assay = function (x = NULL, y = NULL, add.cell.ids = NULL, merge.data = TRUE, ...) 
{
  assays <- c(x, y)
  if (!is.null(x = add.cell.ids)) {
    for (i in seq_along(along.with = assays)) {
      assays[[i]] <- RenameCells(object = assays[[i]], 
                                 new.names = add.cell.ids[i])
    }
  }
  counts.mats <- lapply(X = assays, FUN = Seurat:::ValidateDataForMerge, 
                        slot = "counts")
  keys <- sapply(X = assays, FUN = Key)
  merged.counts <- SeuratObject:::RowMergeSparseMatrices(mat1 = counts.mats[[1]], 
                                                         mat2 = counts.mats[2:length(x = counts.mats)])
  combined.assay <- CreateAssayObject(counts = merged.counts, 
                                      min.cells = -1, min.features = -1)
  if (length(x = unique(x = keys)) == 1) {
    Key(object = combined.assay) <- keys[1]
  }
  if (merge.data) {
    data.mats <- lapply(X = assays, FUN = Seurat:::ValidateDataForMerge, 
                        slot = "data")
    merged.data <- SeuratObject:::RowMergeSparseMatrices(mat1 = data.mats[[1]], 
                                                         mat2 = data.mats[2:length(x = data.mats)])
    if (!all.equal(target = colnames(x = combined.assay), 
                   current = colnames(x = merged.data))) {
      merged.data <- merged.data[, colnames(x = combined.assay)]
    }
    combined.assay <- SetAssayData(object = combined.assay, 
                                   slot = "data", new.data = merged.data)
  }
  if (all(Seurat:::IsSCT(assay = assays))) {
    vst.set.new <- list()
    idx <- 1
    umi.assay.new <- list()
    for (i in 1:length(x = assays)) {
      vst.set.old <- Misc(object = assays[[i]], slot = "vst.set")
      umi.assay.old <- Misc(object = assays[[i]], slot = "umi.assay")
      if (!is.null(x = vst.set.old) && length(x = vst.set.old) > 
          1) {
        for (j in 1:length(x = vst.set.old)) {
          vst.set.new[[idx]] <- vst.set.old[[j]]
          umi.assay.new[[idx]] <- umi.assay.old[[j]]
          idx <- idx + 1
        }
      }
      else if (!is.null(x = Misc(object = assays[[i]], 
                                 slot = "vst.out"))) {
        vst.set.new[[idx]] <- Misc(object = assays[[i]], 
                                   slot = "vst.out")
        umi.assay.new[[idx]] <- Misc(object = assays[[i]], 
                                     slot = "umi.assay")
        idx <- idx + 1
      }
    }
    Misc(object = combined.assay, slot = "vst.set") <- vst.set.new
    Misc(object = combined.assay, slot = "umi.assay") <- umi.assay.new
    scale.data <- do.call(what = cbind, args = lapply(X = assays, 
                                                      FUN = GetAssayData, slot = "scale.data"))
    combined.assay <- SetAssayData(object = combined.assay, 
                                   slot = "scale.data", new.data = scale.data)
  }
  return(combined.assay)
}

# Modified from Seurat
GetResidual = function (object, features, assay = NULL, umi.assay = NULL, clip.range = NULL, 
                        replace.value = FALSE, na.rm = TRUE, verbose = TRUE) 
{
  assay <- assay %||% DefaultAssay(object = object)
  #if (Seurat:::IsSCT(assay = object[[assay]])) {
  #  object[[assay]] <- as(object[[assay]], "SCTAssay")
  #}
  #if (!inherits(x = object[[assay]], what = "SCTAssay")) {
  #  stop(assay, " assay was not generated by SCTransform")
  #}
  sct.models <- levels(x = object[[assay]])
  if (length(x = sct.models) == 0) {
    warning("SCT model not present in assay", call. = FALSE, 
            immediate. = TRUE)
    return(object)
  }
  possible.features <- unique(x = unlist(x = lapply(X = sct.models, 
                                                    FUN = function(x) {
                                                      rownames(x = SCTResults(object = object[[assay]], 
                                                                              slot = "feature.attributes", model = x))
                                                    })))
  bad.features <- setdiff(x = features, y = possible.features)
  if (length(x = bad.features) > 0) {
    warning("The following requested features are not present in any models: ", 
            paste(bad.features, collapse = ", "), call. = FALSE)
    features <- intersect(x = features, y = possible.features)
  }
  features.orig <- features
  if (na.rm) {
    features <- names(x = which(x = table(unlist(x = lapply(X = sct.models, 
                                                            FUN = function(x) {
                                                              rownames(x = SCTResults(object = object[[assay]], 
                                                                                      slot = "feature.attributes", model = x))
                                                            }))) == length(x = sct.models)))
    if (length(x = features) == 0) {
      return(object)
    }
  }
  features <- intersect(x = features.orig, y = features)
  if (length(x = sct.models) > 1 & verbose) {
    message("This SCTAssay contains multiple SCT models. Computing residuals for cells using")
  }
  new.residuals <- lapply(X = sct.models, FUN = function(x) {
    Seurat:::GetResidualSCTModel(object = object, assay = assay, SCTModel = x, 
                                 new_features = features, replace.value = replace.value, 
                                 clip.range = clip.range, verbose = verbose)
  })
  existing.data <- GetAssayData(object = object, slot = "scale.data", 
                                assay = assay)
  all.features <- union(x = rownames(x = existing.data), y = features)
  new.scale <- matrix(data = NA, nrow = length(x = all.features), 
                      ncol = ncol(x = object), dimnames = list(all.features, 
                                                               Cells(x = object)))
  if (nrow(x = existing.data) > 0) {
    new.scale[1:nrow(x = existing.data), ] <- existing.data
  }
  if (length(x = new.residuals) == 1 & is.list(x = new.residuals)) {
    new.residuals <- new.residuals[[1]]
  }
  else {
    new.residuals <- Reduce(cbind, new.residuals)
  }
  new.scale[rownames(x = new.residuals), colnames(x = new.residuals)] <- new.residuals
  if (na.rm) {
    new.scale <- new.scale[!rowAnyNAs(x = new.scale), ]
  }
  object <- SetAssayData(object = object, assay = assay, slot = "scale.data", 
                         new.data = new.scale)
  if (any(!features.orig %in% rownames(x = new.scale))) {
    bad.features <- features.orig[which(!features.orig %in% 
                                          rownames(x = new.scale))]
    warning("Residuals not computed for the following requested features: ", 
            paste(bad.features, collapse = ", "), call. = FALSE)
  }
  return(object)
}

# Modified from Seurat
FindAllMarkers2 = function(
  srt,
  assay = NULL,
  group.by = NULL,
  test.use = 'wilcox',
  latent.vars = NULL, 
  features = VariableFeatures(srt), 
  logfc.threshold = 0.25, 
  min.pct = 0.1, 
  return.thresh = 0.01,
  return.thresh.adj = 1,
  only.pos = TRUE,
  do.enrichment = TRUE,
  enrichment.matrix = NULL,
  enrichment.type = 'GO',
  do.plot = TRUE,
  cols = NULL,
  do.print = TRUE,
  print.bar = FALSE,
  ...
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  if (!is.null(group.by)){
    Idents(srt) = group.by
  }
  if (!all(srt$species == 'H|Hs')){
    do.enrichment = FALSE
  }
  genes_all = FindAllMarkers(srt, assay = assay, latent.vars = latent.vars, features = features, test.use = test.use, logfc.threshold = logfc.threshold, 
                             min.pct = min.pct, return.thresh = return.thresh, return.thresh.adj = return.thresh.adj, only.pos = only.pos, print.bar = print.bar, ...)
  genes_100 = top_n(group_by(top_n(group_by(genes_all, cluster), -100, p_val_adj), cluster), 100, avg_log2FC)
  genes_10 = top_n(group_by(top_n(group_by(genes_all, cluster), -10, p_val_adj), cluster), 10, avg_log2FC)
  genes_1 = top_n(group_by(top_n(group_by(genes_all, cluster), -1, p_val_adj), cluster), 1, avg_log2FC)
  if (do.print){
    print(genes_1)
  }
  if (do.plot){
    #h = FeaturePlot(srt, features = genes_1$gene)
    #print(h)
    h = DoHeatmap(srt, features = genes_10$gene, group.by = group.by, assay = assay, cols = cols, label = FALSE, ...)
    print(h)
  }
  if (do.enrichment){
    if (is.null(enrichment.matrix)){
      enrichment.matrix = BuildEnrichmentMatrix(rownames(srt), type = enrichment.type)
      #enrichment.matrix = BuildEnrichmentMatrix(genes.use, type = enrichment.type)
    }
    go_all = sapply(split(genes_all, genes_all$cluster), function(tbl){
      Enrichment(geneset = tbl$gene, enrichment.matrix = enrichment.matrix)
    })
    names(go_all) = levels(Idents(srt))
    go_100 = lapply(colnames(go_all), function(cluster){
      sort(go_all[,cluster])[1:100]
    })
    names(go_100) = levels(Idents(srt))
    go_10 = lapply(go_100, '[', 1:10)
    names(go_10) = levels(Idents(srt))
    go_1 = lapply(go_10, '[', 1)
    names(go_1) = levels(Idents(srt))
  } else {
    go_all = NA
    go_100 = NA
    go_10 = NA
    go_1 = NA
  }
  if (do.print){
    print(go_10)
  }
  return(list('genes_all' = genes_all, 'genes_100' = genes_100, 'genes_10' = genes_10, 'genes_1' = genes_1, 
              'go_all' = go_all, 'go_100' = go_100, 'go_10' = go_10, 'go_1' = go_1))
}

# Modified from Seurat
FindAllConservedMarkers2 = function(
  srt,
  assay = NULL,
  group.by = NULL,
  grouping.var = NULL,
  features = VariableFeatures(srt), 
  logfc.threshold = 0.25, 
  min.pct = 0.1, 
  return.thresh = 0.1,
  return.thresh.adj = 0.1,
  only.pos = TRUE,
  do.enrichment = TRUE,
  enrichment.matrix = NULL,
  enrichment.type = 'GO',
  do.plot = TRUE,
  cols = NULL,
  do.print = TRUE,
  print.bar = FALSE,
  ...
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  if (!is.null(group.by)){
    Idents(srt) = group.by
  }
  if (!all(srt$species == 'H|Hs')){
    do.enrichment = FALSE
  }
  genes.list = lapply(levels(Idents(srt)), function(l){
    g = FindConservedMarkers(srt, assay = assay, ident.1 = l, grouping.var = grouping.var, features = features, logfc.threshold = logfc.threshold, 
                             min.pct = min.pct, return.thresh = return.thresh, return.thresh.adj = return.thresh.adj, only.pos = only.pos, print.bar = print.bar, ...)
    g$cluster = l
    g$gene = rownames(g)
    return(g)
  })
  genes_all = Reduce(rbind, genes.list)
  genes_100 = top_n(group_by(genes_all, cluster), -100, max_pval)
  genes_10 = top_n(group_by(genes_all, cluster), -10, max_pval)
  genes_1 = top_n(group_by(genes_all, cluster), -1, max_pval)
  if (do.print){
    print(genes_1)
  }
  if (do.plot){
    #h = FeaturePlot(srt, features = genes_1$gene)
    #print(h)
    h = DoHeatmap(srt, features = genes_10$gene, group.by = group.by, assay = assay, cols = cols, ...)
    print(h)
  }
  if (do.enrichment){
    if (is.null(enrichment.matrix)){
      enrichment.matrix = BuildEnrichmentMatrix(rownames(srt), type = enrichment.type)
      #enrichment.matrix = BuildEnrichmentMatrix(genes.use, type = enrichment.type)
    }
    go_all = sapply(split(genes_all, genes_all$cluster), function(tbl){
      Enrichment(geneset = tbl$gene, enrichment.matrix = enrichment.matrix)
    })
    names(go_all) = levels(Idents(srt))
    go_100 = lapply(colnames(go_all), function(cluster){
      sort(go_all[,cluster])[1:100]
    })
    names(go_100) = levels(Idents(srt))
    go_10 = lapply(go_100, '[', 1:10)
    names(go_10) = levels(Idents(srt))
    go_1 = lapply(go_10, '[', 1)
    names(go_1) = levels(Idents(srt))
  } else {
    go_all = NA
    go_100 = NA
    go_10 = NA
    go_1 = NA
  }
  if (do.print){
    print(go_10)
  }
  return(list('genes_all' = genes_all, 'genes_100' = genes_100, 'genes_10' = genes_10, 'genes_1' = genes_1, 
              'go_all' = go_all, 'go_100' = go_100, 'go_10' = go_10, 'go_1' = go_1))
}

# Modified from Seurat
DoHeatmap = function (object, features = NULL, cells = NULL, group.by = "ident", 
                      group.bar = TRUE, disp.min = -2.5, disp.max = NULL, slot = "scale.data", 
                      assay = NULL, label = TRUE, size = 5.5, hjust = 0, angle = 45, 
                      raster = TRUE, draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02, 
                      combine = TRUE, cols = NULL) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[group.by]][cells, , drop = FALSE]
  plots <- vector(mode = "list", length = ncol(x = groups.use))
  for (i in 1:ncol(x = groups.use)) {
    data.group <- data
    group.use <- groups.use[, i, drop = TRUE]
    if (!is.factor(x = group.use)) {
      group.use <- factor(x = group.use)
    }
    names(x = group.use) <- cells
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use)) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- rep(x = levels(x = group.use), 
                                times = lines.width)
      group.levels <- levels(x = group.use)
      names(x = placeholder.groups) <- placeholder.cells
      group.use <- as.vector(x = group.use)
      names(x = group.use) <- cells
      group.use <- factor(x = c(group.use, placeholder.groups), 
                          levels = group.levels)
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = names(x = sort(x = group.use)), group.by = group.use)
    if (group.bar) {
      group.use2 <- sort(x = group.use)
      if (draw.lines) {
        na.group <- Seurat:::RandomName(length = 20)
        levels(x = group.use2) <- c(levels(x = group.use2), 
                                    na.group)
        group.use2[placeholder.cells] <- na.group
        if (is.null(cols)){
          cols <- c(hue_pal()(length(x = levels(x = group.use))), 
                    "#FFFFFF")
        }
      }
      else {
        if (is.null(cols)){
          cols <- c(hue_pal()(length(x = levels(x = group.use))))
        }
      }
      pbuild <- ggplot_build(plot = plot)
      if (length(cols) < nlevels(group.use2)){
        cols = cols[1:nlevels(group.use2)]
      }
      if (length(cols) > nlevels(group.use2)){
        cols = rep(cols, nlevels(group.use2))[1:nlevels(group.use2)]
      }
      names(x = cols) <- levels(x = group.use2)
      y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
      y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + 
        y.range * 0.015
      y.max <- y.pos + group.bar.height * y.range
      plot <- plot + annotation_raster(raster = t(x = cols[group.use2]), 
                                       , xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
        coord_cartesian(ylim = c(0, y.max), clip = "off")
      if (label) {
        x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
        x.divs <- pbuild$layout$panel_params[[1]]$x.major
        x <- data.frame(group = sort(x = group.use), 
                        x = x.divs)
        label.x.pos <- tapply(X = x$x, INDEX = x$group, 
                              FUN = median) * x.max
        label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                  label.x.pos)
        plot <- plot + geom_text(stat = "identity", data = label.x.pos, 
                                 aes_string(label = "group", x = "label.x.pos"), 
                                 y = y.max + y.max * 0.03 * 0.5, angle = angle, 
                                 hjust = hjust, size = size)
        plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                 y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) * 
                                                                   size), clip = "off"))
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

#### Transcriptome ####

# To get full expression matrix
GetData = function(
  srt,
  genes = NULL,
  slot = 'scale.data',
  assay = NULL
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  if (is.null(genes)){
    if ('RNA' %in% names(srt@assays)){
      genes = rownames(GetAssayData(srt, assay = 'RNA', slot = 'counts'))
    } else if  ('Spatial' %in% names(srt@assays)){
      genes = rownames(GetAssayData(srt, assay = 'Spatial', slot = 'counts'))
    } else {
      genes = rownames(GetAssayData(srt, assay = assay, slot = 'counts'))
    }
  }
  data = GetAssayData(srt, assay = assay, slot = slot)
  missing = setdiff(genes, rownames(data))
  add = matrix(0, nrow = length(missing), ncol = ncol(data))
  rownames(add) = missing
  data = rbind(data, add)
  data = data[genes, ]
  return(data)
}

FilterCells = function(
  srt,
  ngene.min = 0,
  numi.min = 0,
  mt.min = 0,
  rb.min = 0,
  ngene.max = Inf,
  numi.max = Inf,
  mt.max = 0.3,
  rb.max = 0.3,
  do.plot = FALSE,
  assay = NULL,
  ...
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  
  print(paste('Starting from', ncol(srt), 'cells'))
  
  genes = rownames(srt)
  counts = GetAssayData(srt, assay = assay, slot = 'counts')
  
  if (!'percent.mito' %in% colnames(srt@meta.data)){
    mt.pattern = "^MT-|^MTRNR2L"
    mito.genes = grep(pattern = mt.pattern, x = genes, value = TRUE, ignore.case = TRUE)
    srt$percent.mito = Matrix::colSums(counts[mito.genes, ])/Matrix::colSums(counts)
  }
  
  if (!'percent.ribo' %in% colnames(srt@meta.data)){
    rb.pattern = "RPL|RPS|RPP"
    ribo.genes = grep(pattern = rb.pattern, x = genes, value = TRUE, ignore.case = TRUE)
    srt$percent.ribo = Matrix::colSums(counts[ribo.genes, ])/Matrix::colSums(counts)
  }
  
  if (do.plot){
    par(mfrow = c(2,2))
    hist(srt@meta.data[,paste0('nFeature_',assay)], 100, xlab = 'nGene', main = '')
    abline(v = ngene.min, col = 'darkgreen')
    abline(v = ngene.max, col = 'darkred')
    hist(srt@meta.data[,paste0('nCount_',assay)], 100, xlab = 'nUMI', main = '')
    abline(v = numi.min, col = 'darkgreen')
    abline(v = numi.max, col = 'darkred')
    hist(srt@meta.data$percent.mito, 100, xlab = 'percent.mito', xlim = c(0,1), main = '')
    abline(v = mt.min, col = 'darkgreen')
    abline(v = mt.max, col = 'darkred')
    hist(srt@meta.data$percent.ribo, 100, xlab = 'percent.ribo', xlim = c(0,1), main = '')
    abline(v = rb.min, col = 'darkgreen')
    abline(v = rb.max, col = 'darkred')
    par(mfrow = c(1,1))
  }
  
  srt = srt[, srt@meta.data[, paste('nFeature', assay, sep = '_')] >= ngene.min &
              srt@meta.data[, paste('nFeature', assay, sep = '_')] <= ngene.max &
              srt@meta.data[, paste('nCount', assay, sep = '_')] >= numi.min &
              srt@meta.data[, paste('nCount', assay, sep = '_')] <= numi.max &
              srt@meta.data$percent.mito >= mt.min &
              srt@meta.data$percent.mito <= mt.max &
              srt@meta.data$percent.ribo >= rb.min &
              srt@meta.data$percent.ribo <= rb.max]
  
  print(paste('Filtering to', ncol(srt), 'cells'))
  
  return(srt)
}

FilterGenes = function(
  srt, 
  genes.exclude = c(),
  assay = NULL
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  genes = rownames(srt)
  mt.pattern = "^MT-|^MTRNR2L"
  mito.genes = grep(pattern = mt.pattern, x = genes, value = TRUE, ignore.case = TRUE)
  rb.pattern = "RPL|RPS|RPP"
  ribo.genes = grep(pattern = rb.pattern, x = genes, value = TRUE, ignore.case = TRUE)
  genes = !duplicated(genes) & !genes%in%mito.genes & !genes%in%ribo.genes & !genes %in% genes.exclude
  srt = srt[genes, ]
  return(srt)
}

FindArchetypes = function(
  srt,
  assay = NULL,
  slot = 'data',
  arch.k = 4,
  reduction = 'pca',
  dims = 1:3,
  neighborhood.k = NULL,
  do.search = FALSE,
  do.plot = TRUE,
  ...
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  if (is.null(neighborhood.k)){
    neighborhood.k = ncol(srt)/(2*arch.k)
  }
  
  data.pca = Embeddings(srt, reduction = reduction, assay = assay)[,dims]
  
  if (do.search){
    k.max = max(arch.k, max(dims)) + 5
    steps = stepArchetypes(data.pca, k = 1:k.max, verbose = FALSE, nrep = 10, maxIterations = 100)
    if (do.plot){
      ElbowPlot(srt)
      screeplot(steps)
    }
    a = bestModel(steps)[[arch.k]]
    if (is.null(a$archetypes)){a = archetypes(data.pca, k = arch.k, maxIterations = 100)}
  } else {
    steps = stepArchetypes(data.pca, k = arch.k, verbose = FALSE, nrep = 50, maxIterations = 100)
    a = bestModel(steps)
    if (is.null(a$archetypes)){a = archetypes(data.pca, k = arch.k, maxIterations = 100)}
  }
  
  archetypes = apply(a$alphas, 1, order, decreasing = TRUE)
  top_archetype = archetypes[1,]
  srt$archetype = top_archetype
  srt$archetype = factor(srt$archetype)
  
  # srt$neighborhood = apply(a$alphas, 1, function(cell){
  #   s = sort(cell, decreasing = TRUE)
  #   if (s[1] > 2*s[2]){
  #     return(which.max(cell))
  #   } else {
  #     return(0)
  #   }
  # })
  # srt$neighborhood = factor(srt$neighborhood)
  
  neighborhoods = apply(a$alphas, 2, order, decreasing = TRUE)
  top_neighborhood = neighborhoods[1:neighborhood.k, ]
  srt$neighborhood = 0
  srt$neighborhood[top_neighborhood] = top_archetype[top_neighborhood]
  srt$neighborhood = factor(srt$neighborhood)
  
  if (do.plot){
    xyplot(a, data.pca, chull = chull(data.pca), atypes.col = hue_pal()(arch.k))
    text(a$archetypes*1.05, pch = 20, col = hue_pal()(arch.k))
    for (i in seq(1, length(dims)-1, by = 2)){
      plot(data.pca[, i:(i+1)], pch = 20, col = hue_pal()(arch.k)[srt$archetype])
      points(a$archetypes[, i:(i+1)], pch = 23, bg = hue_pal()(arch.k), col = 'black')
      text(a$archetypes[, i:(i+1)], pch = 20, col = hue_pal()(arch.k), pos = 4)
      plot(data.pca[, i:(i+1)], pch = 20, col = c('darkgrey', hue_pal()(arch.k))[as.numeric(srt$neighborhood)])
      points(a$archetypes[, i:(i+1)], pch = 23, bg = hue_pal()(arch.k), col = 'black')
      text(a$archetypes[, i:(i+1)], pch = 20, col = hue_pal()(arch.k), pos = 4)
    }
  }
  
  Idents(srt) = 'archetype'
  
  return(srt)
}

CompareGroupings = function(
  srt,
  group1 = 'orig.ident',
  group2 = 'orig.ident',
  ...
){
  try({
    h = DimPlot(srt, group.by = group1)
    print(h)
    h = DimPlot(srt, group.by = group2)
    print(h)
  })
  
  id1 = unlist(FetchData(srt, group1))
  id2 = unlist(FetchData(srt, group2))
  
  if (!any(is.na(as.integer(id1))) | is.factor(id1)
      & (!any(is.na(as.integer(id2))) | is.factor(id2))){
    tab = table(id1, id2)
    names(dimnames(tab)) = c(group1, group2)
    barplot(height = tab, legend = FALSE, col = hue_pal()(dim(tab)[1]), args.legend = list(title = group1), xlab = group2, ylab = 'Frequency',...)
    barplot(height = t(t(tab)/colSums(tab)), legend = FALSE, col = hue_pal()(dim(tab)[1]), args.legend = list(title = group1), xlab = group2, ylab = 'Proportion',...)
    barplot(height = t(tab), legend = FALSE, col = hue_pal()(dim(tab)[2]), args.legend = list(title = group2), xlab = group1, ylab = 'Frequency',...)
    barplot(height = t(tab/rowSums(tab)), legend = FALSE, col = hue_pal()(dim(tab)[2]), args.legend = list(title = group2), xlab = group1, ylab = 'Proportion',...)
    print(tab)
    print(chisq.test(tab))
    return(tab)
  } else {
    if (!any(is.na(as.integer(id1))) | is.factor(id1)){
      VlnPlot(srt, group2, group.by = group1)
      boxplot(id2 ~ as.factor(id1))
      print(anova(lm(id2 ~ as.factor(id1))))
      return()
    }
    if (!any(is.na(as.integer(id2))) | is.factor(id2)){
      VlnPlot(srt, group1, group.by = group2)
      boxplot(id1 ~ as.factor(id2))
      print(anova(lm(id1 ~ as.factor(id2))))
      return()
    }
  }
}

#### Gene modules ####

FindMSigDB = function(
  type
){
  if (type == 'GO'){
    db = MSigDB$C5_GENE_ONTOLOGY
  } else if (type == 'HALLMARK'){
    db = MSigDB$HALLMARK
  } else if (type == 'MOTIF'){
    db = MSigDB$C3_MOTIF[-grep('UNKNOWN|MIR',names(MSigDB$C3_MOTIF))]
  } else if (type == 'PATHWAYS'){
    db = MSigDB$C2_CURATED[grep('BIOCARTA|REACTOME|KEGG',names(MSigDB$C2_CURATED))]
  } else if (type == 'BIOCARTA'){
    db = MSigDB$C2_CURATED[grep('BIOCARTA',names(MSigDB$C2_CURATED))]
  } else if (type == 'KEGG'){
    db = MSigDB$C2_CURATED[grep('KEGG',names(MSigDB$C2_CURATED))]
  } else if (type == 'REACTOME'){
    db = MSigDB$C2_CURATED[grep('REACTOME',names(MSigDB$C2_CURATED))]
  } else {
    db = MSigDB$C5_GENE_ONTOLOGY[grep(type, names(MSigDB$C5_GENE_ONTOLOGY), ignore.case = TRUE, value = TRUE)]
  }
  return(db)
}

# From GO-PCA
MEnrichment = function(
  loadings,
  sign = 'positive',
  genes = NULL,
  srt = NULL,
  type = 'GO',
  db = NULL,
  do.sort = FALSE,
  pval.thresh = 1,
  pval.n = Inf,
  do.plot = FALSE,
  do.print = FALSE,
  ...
){
  if (is.null(genes) & is.null(db)){
    genes = rownames(srt)
    db = FindMSigDB(type)
    db = lapply(db, intersect, genes)
  } 
  if (!is.null(genes) & is.null(db)){
    db = FindMSigDB(type)
    db = lapply(db, intersect, genes)
  }
  
  if (sign == 'positive'){
    geneset.rank = names(sort(loadings, decreasing = TRUE))
  } else if (sign == 'negative'){
    geneset.rank = names(sort(loadings, decreasing = FALSE))
  }
  pvals = sapply(db, function(set){
    lambdas = geneset.rank %in% set
    t = mHG.test(lambdas)
    return(t$p.value)
  })
  
  if (do.sort){
    pvals = sort(pvals)
  }  
  pvals = pvals[pvals <= pval.thresh]
  if (length(pvals) >= pval.n){
    pvals = sort(pvals)
    pvals = pvals[1:pval.n]
  }
  if (do.print){
    print(pvals)
  }
  if (do.plot){
    PlotEnrichment(pvals, do.sort = TRUE, ...)
  }
  return(pvals)
  
}

# From GO-PCA
ReduceDB = function(
  srt,
  db,
  assay = NULL,
  dims = 1:5,
  l_max = 1 / 8,
  do.print = TRUE,
  do.plot = FALSE
){
  
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  } 
  DefaultAssay(srt) = assay
  srt = RunPCA(srt, features = rownames(srt), verbose = FALSE)
  srt = ProjectDim(srt, overwrite = TRUE, verbose = FALSE)
  
  signatures = lapply(db, function(set){
    
    all = t(sapply(dims, function(pc){
      gene.loadings = Loadings(srt, reduction = 'pca', projected = TRUE)[,pc]
      sapply(c(-1,1), function(pc.sign){
        genes.rank = names(sort(gene.loadings*pc.sign, decreasing = TRUE))
        lambdas = genes.rank %in% set
        stat = mHG.statistic.calc(lambdas = lambdas, n_max = length(lambdas) * l_max)
        pval = mHG.pval.calc(p = stat@mHG, N = length(lambdas), B = length(set), n_max = length(lambdas) * l_max)
        return(pval)
      })
    }))
    pval = min(all)
    a = which(all == pval, arr.ind = TRUE)[1,]
    pc = dims[a[1]]
    pc.sign = c(-1,1)[a[2]]
    gene.loadings = Loadings(srt, reduction = 'pca', projected = TRUE)[,pc]
    genes.rank = names(sort(gene.loadings*pc.sign, decreasing = TRUE))
    lambdas = genes.rank %in% set
    stat = mHG.statistic.calc(lambdas = lambdas, n_max = length(lambdas) * l_max)
    n = stat@n
    db = intersect(genes.rank[1:n], set)
    
    return(list('db' = db, 'stats' = c('pval' = pval, 'n' = n, 'length' = length(db), 'fraction' = length(db)/length(set), 'pc' = pc, 'pc.sign' = pc.sign)))
  })
  
  dbs = lapply(signatures, '[[', 'db')
  stats = as.data.frame(t(sapply(signatures, '[[', 'stats')))
  stats$pval.log = -log10(stats$pval)
  stats$pval.adjust = p.adjust(stats$pval, method = 'BH')
  stats$pval.adjust.log = -log10(stats$pval.adjust)
  
  return(list('dbs' = dbs, 'stats' = stats))
  
}

BuildEnrichmentMatrix = function(
  genes,
  type = 'GO',
  db = NULL
){
  if (is.null(db)){
    db = FindMSigDB(type)
  }
  terms = names(db)
  enrichment.matrix = sapply(db, function(term){
    genes %in% term
  })
  rownames(enrichment.matrix) = genes
  #enrichment.matrix = enrichment.matrix[, colSums(enrichment.matrix) > 0]
  enrichment.matrix[, colSums(enrichment.matrix) == 0] = NA
  return(enrichment.matrix)
}

Enrichment = function(
  geneset,
  genes = NULL,
  enrichment.matrix = NULL,
  srt = NULL,
  type = 'GO',
  do.sort = FALSE,
  pval.thresh = 1,
  pval.n = Inf,
  do.plot = FALSE,
  do.print = FALSE,
  ...){
  
  if (is.null(genes) & is.null(enrichment.matrix)){
    genes = rownames(srt)
    enrichment.matrix = BuildEnrichmentMatrix(genes, type = type)
  } 
  if (!is.null(genes) & is.null(enrichment.matrix)){
    enrichment.matrix = BuildEnrichmentMatrix(genes, type = type)
  }
  if (is.null(genes) & !is.null(enrichment.matrix)){
    genes = rownames(enrichment.matrix)
  } 
  
  geneset.log = (genes %in% geneset)
  present = geneset.log %*% enrichment.matrix
  pvals = phyper(present, colSums(enrichment.matrix), colSums(1 - enrichment.matrix), sum(geneset.log), lower.tail = FALSE, log.p = FALSE)
  names(pvals) = colnames(enrichment.matrix)
  
  if (do.sort){
    pvals = sort(pvals)
  }  
  pvals = pvals[pvals <= pval.thresh]
  if (length(pvals) >= pval.n){
    pvals = sort(pvals)
    pvals = pvals[1:pval.n]
  }
  if (do.print){
    print(pvals)
  }
  if (do.plot){
    PlotEnrichment(pvals, do.sort = TRUE, ...)
  }
  return(pvals)
}

PlotEnrichment = function(
  pvals,
  do.sort = FALSE,
  pval.thresh = 1,
  pval.n = Inf,
  do.plot = TRUE,
  do.print = FALSE,
  ...
){
  pvals = pvals[pvals < pval.thresh]
  if (do.sort){
    pvals = sort(pvals)
  }  
  pvals = pvals[pvals < pval.thresh]
  if (length(pvals) >= pval.n){
    pvals = pvals[1:pval.n]
  }
  if (do.print){
    print(pvals)
  }
  z = (pvals == 0)
  pvals[z] = 10^-10
  names(pvals)[z] = paste0(names(pvals)[z], '*')
  if (do.plot & (length(pvals) > 0)){ 
    barplot(-log10(pvals), horiz = TRUE, xlab = '-log10(p-value)', space = 0, names.arg = '', ...)
    text(0, (1:length(pvals)) - 0.5, labels = names(pvals), pos = 4, cex = 0.7)
  }
}

# Make equivalent random modules
MakeRand = function(
  srt,
  db,
  assay = NULL,
  nrand = 3,
  nbin = 25
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  data = GetData(srt, slot = 'data')
  db = lapply(db, intersect, rownames(data))
  data.avg = sort(rowMeans(x = data))
  data.cut = cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                        n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) = names(x = data.avg)
  binned = split(names(data.cut), data.cut)
  db_rand = lapply(names(db), function(m){
    lapply(1:10^nrand, function(i){
      used = vector()
      unused = binned
      for (g in db[[m]]){
        pool = data.cut[g]
        new = sample(unused[[pool]], 1)
        used = c(used, new)
        unused[[pool]] = setdiff(unused[[pool]], new)
      }
      return(used)
    })
  })
  names(db_rand) = names(db)
  return(db_rand)
}

# Modules to cells
GeneToEnrichment = function(
  srt,
  type = 'GO',
  db = NULL,
  method = 'rand',
  genes = NULL,
  assay = NULL,
  do.rescale = FALSE,
  min.cells = 0,
  min.genes = 0,
  min.var = 0,
  min.var.rescaled = 0,
  auc_percentile = 0.05,
  db_rand = NULL,
  nrand = 4,
  nbin = 25,
  ...
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  if (is.null(db)){
    db = FindMSigDB(type)
  }
  
  counts = as.matrix(GetData(srt, assay = assay, slot = 'counts'))
  genes = rownames(counts)
  genes.expr = rownames(counts)[rowSums(counts) > min.cells]
  
  if (method == 'metagene'){
    
    data = as.matrix(GetAssayData(srt, assay = assay, slot = 'scale.data'))
    
    db = lapply(db, intersect, genes.expr)
    
    enrichment.profile = t(sapply(names(db), function(m){
      colMeans(data[db[[m]], ], na.rm = TRUE)
    }))
    
    enrichment.profile = enrichment.profile[sapply(names(db), function(x){
      v = var(enrichment.profile[x, ])
      l = length(db[[x]])
      return(l > min.genes
             && v > min.var
             && v*l^2 > min.var.rescaled)
    }), ]
    
    if (do.rescale){
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn) / sqrt(v)
    }
    
    srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  
  if (method == 'auc'){
    
    data = as.matrix(GetData(srt, assay = assay, slot = 'data'))
    
    cells_rankings = AUCell_buildRankings(data)
    cells_AUC = AUCell_calcAUC(db, cells_rankings, aucMaxRank=nrow(cells_rankings)*auc_percentile)
    enrichment.profile = getAUC(cells_AUC)
    
    if (do.rescale){
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn) / sqrt(v)
    }
    
    srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  
  if (method == 'score'){
    
    temp = AddModuleScore(srt, features = db, assay = assay, name = names(db), nbin = nbin, ...)
    
    enrichment.profile = t(temp@meta.data[, names(db)])
    
    if (do.rescale){
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn) / sqrt(v)
    }
    
    srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  
  if (method == 'rand'){
    
    data = as.matrix(GetData(srt, assay = assay, slot = 'scale.data'))
    
    db = lapply(db, intersect, genes)
    
    if (is.null(db_rand)){
      db_rand = MakeRand(srt, db, nrand = nrand, nbin = nbin)
    } else {
      nrand = log10(length(db_rand[[1]]))
    }
    
    enrichment.profile = t(sapply(names(db), function(m){
      ra = sapply(db_rand[[m]], function(i){
        colMeans(data[i, ], na.rm = TRUE)
      })
      re = colMeans(data[db[[m]], ], na.rm = TRUE)
      p = rowMeans(ra >= re)
      p = -log10(p)
      return(p)
    }))
    enrichment.profile[is.infinite(enrichment.profile)] = nrand
    enrichment.profile = enrichment.profile/nrand
    
    srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  
  return(srt)
}

# Modules to dataset statistics
ScoreModule = function(
  srt,
  db = NULL,
  type = 'GO',
  method = 'rand',
  genes = NULL,
  assay = NULL,
  do.rescale = FALSE,
  min.cells = 0,
  min.genes = 0,
  min.var = 0,
  min.var.rescaled = 0,
  auc_percentile = 0.05,
  db_rand = NULL,
  nrand = 4,
  nbin = 25,
  ...
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  if (is.null(db)){
    db = FindMSigDB(type)
  }
  
  counts = as.matrix(GetData(srt, assay = assay, slot = 'counts'))
  genes = rownames(counts)
  genes.expr = rownames(counts)[rowSums(counts) > min.cells]
  
  if (method == 'rand'){
    
    data = as.matrix(GetData(srt, assay = assay, slot = 'scale.data'))
    
    db = lapply(db, intersect, genes)
    
    if (is.null(db_rand)){
      db_rand = MakeRand(srt, db, nrand = nrand, nbin = nbin)
    } else {
      nrand = log10(length(db_rand[[1]]))
    }
    
    v = sapply(names(db), function(m){
      ra = sapply(db_rand[[m]], function(i){
        res = c(NA,NA)
        tryCatch(expr = {
          srt_temp = RunPCA(srt, features = i, verbose = FALSE, npcs = 2)
          res = srt_temp@reductions$pca@stdev
        }, error = function(e){})
        return(res)
      })
      re = c(NA,NA)
      tryCatch(expr = {
        srt_temp = RunPCA(srt, features = modules[[m]], verbose = FALSE, npcs = 2)
        re = srt_temp@reductions$pca@stdev
      }, error = function(e){})
      presence = -log10(mean(ra[1,] >= re[1]))
      presence[is.infinite(presence)] = nrand
      coordination = -log10(mean(ra[1,]/ra[2,] >= re[1]/re[2]))
      coordination[is.infinite(coordination)] = nrand
      return(c('presence' = presence, 'coordination' = coordination))
    })
    
  }
  
  return(v)
  
}

NMFToModules = function(
  res,
  gmin = 5
){
  
  scores = basis(res)
  coefs = coefficients(res)
  
  # Remove if fewer than gmin genes
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  l = sapply(modules, length)
  keep = (l >= gmin)
  scores = scores[, keep]
  coefs = coefs[keep, ]
  
  # Find modules
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  
  names(modules) = sapply(modules, '[', 1)
  names(modules) = paste('m', names(modules), sep = '_')
  names(modules) = gsub('-','_',names(modules))
  
  return(modules)
}

BasisToModules = function(
  basis,
  gmin = 5
){
  
  scores = basis
  
  # Remove if fewer than gmin genes
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  l = sapply(modules, length)
  keep = (l >= gmin)
  scores = scores[, keep]
  
  # Find modules
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  
  names(modules) = sapply(modules, '[', 1)
  
  return(modules)
}


#### Genome ####

FindGranges = function(
  srt = NULL,
  genes = NULL,
  ...
){
  if (is.null(genes)){
    genes = rownames(srt)
  }
  species = unique(srt@meta.data$species)
  if (species == 'H|Hs'){
    mart = useMart("ensembl", dataset = 'hsapiens_gene_ensembl')
  } else if (species == 'M|Mm'){
    mart = useMart("ensembl", dataset = 'mmusculus_gene_ensembl')
  } else if (species == 'D|Dr'){
    mart = useMart("ensembl", dataset = 'drerio_gene_ensembl')
  }
  granges = getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"), 
                  filters = "external_gene_name", values = genes, mart = mart, uniqueRows = TRUE)
  granges = dplyr:::filter(granges, 
                           external_gene_name %in% genes & !duplicated(external_gene_name) & !is.na(as.numeric(chromosome_name)))
  granges = granges[order(as.numeric(granges$chromosome_name)),]
  granges = GRanges(seqnames = granges$chromosome_name, ranges = IRanges(start = granges$start_position, end = granges$end_position), name = granges$external_gene_name)
  granges = sort(granges)
  return(granges)
}

ExportGrangesForCNV = function(
  granges
){
  positions = data.frame(granges, stringsAsFactors = FALSE)
  rownames(positions) = positions$name
  positions = positions[, c('seqnames','start','end')]
  colnames(positions) = c('C_CHR', 'C_START', 'C_STOP')
  write.table(positions, file = 'gene_order.tsv', quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)
}

ConvertSpecies = function(
  srt1 = NULL,
  genes1 = NULL,
  species1 = NULL,
  species2 = 'Hs',
  map = NULL,
  assay = NULL,  
  ...
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  if (is.null(species1)){
    species1 = unique(srt1@meta.data$species)
  }
  # tryCatch(
  #   load(paste0('~/Documents/R/map_', species1, '_', species2, '.RData')),
  #   error = function(e){c()})
  if (is.null(map)){
    if (is.null(genes1)){
      genes1 = rownames(srt1)
    }
    if (species1 == 'H|Hs'){
      mart = useMart("ensembl", dataset = 'hsapiens_gene_ensembl')
    } else if (species1 == 'M|Mm'){
      mart = useMart("ensembl", dataset = 'mmusculus_gene_ensembl')
    } else if (species1 == 'D|Dr'){
      mart = useMart("ensembl", dataset = 'drerio_gene_ensembl')
    }
    if (species2 == 'H|Hs'){
      attributes = c("external_gene_name", grep('hsapiens', listAttributes(mart)$name, value = TRUE))
    } else if (species2 == 'M|Mm'){
      attributes = c("external_gene_name", grep('mmusculus', listAttributes(mart)$name, value = TRUE))
    } else if (species2 == 'D|Dr'){
      attributes = c("external_gene_name", grep('drerio', listAttributes(mart)$name, value = TRUE))
    }
    map = getBM(attributes = attributes, 
                filters = "external_gene_name", values = genes1, mart = mart, uniqueRows = TRUE)
    map = dplyr:::filter(map, external_gene_name %in% genes1 & !duplicated(external_gene_name))
    genes1 = map[, grep('gene_name', colnames(map), value = TRUE)[1]]
    genes2 = map[, grep('gene_name', colnames(map), value = TRUE)[2]]
    map = map[!(genes2 == '' | duplicated(genes2)), ]
    save(map, file = paste0('~/Documents/R/map_', species1, '_', species2, '.RData'))
  }
  if (!is.null(srt1)){
    genes1 = map[, grep('gene_name', colnames(map), value = TRUE)[1]]
    genes2 = map[, grep('gene_name', colnames(map), value = TRUE)[2]]
    data1 = GetAssayData(srt1, assay = assay, slot = 'counts')[genes1, ]
    data2 = data1
    rownames(data2) = genes2
    srt2 = CreateSeuratObject(counts = data2)
    srt2 = FilterGenes(srt2)
    srt2@meta.data = srt1@meta.data
    srt2$species = species2
    return(srt2)
  } else {
    save(map, file = paste0('~/Documents/R/map_', species1, '_', species2, '.RData'))
  }
}

#### Spatial transcriptomics ####

FindSTNeighbors = function(
  st, 
  d_max, 
  d_min = 0
){
  coord = st@images$slice1@coordinates[,c('imagerow','imagecol')]
  distances = as.matrix(dist(coord))
  distances = distances/sort(distances[distances > 0], decreasing = FALSE)[10]
  #distances = round(distances, digits = 1)
  neighbors = apply(distances, 1, function(d){
    names(d)[d >= d_min & d <= d_max]
  })
  names(neighbors) = rownames(coord)
  return(neighbors)
}

MakeSTRand = function(
  st
){
  data = GetData(st, assay = 'SCT', slot = 'data')
  data = t(apply(data, 1, function(row){
    sample(row, size = length(row), replace = FALSE)
  }))
  colnames(data) = colnames(st)
  st@assays$SCT@data = data
  st = ScaleData(st, assay = 'SCT', do.center = TRUE, do.scale = FALSE)
}
