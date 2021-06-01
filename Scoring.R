#!/usr/bin/env Rscript
source('~/Documents/R/seurat_functions.R')
setwd('~/Documents/Analysis/Tumors/')

load('parameters.RData')
modules = loadRData('modules.RData')
colors.module = loadRData('colors.module.RData')

## All

mclapply(files.all, function(f){
  setwd(f)
  srt = loadRData('srt.malignant.RData')

  modules_rand = MakeRand(srt, db = modules, nrand = nrand, nbin = nbin)

  ini = matrix(0,nrow = ncol(srt), ncol = length(modules))
  rownames(ini) = colnames(srt)
  colnames(ini) = names(modules)
  srt@meta.data[,names(modules)] = as.data.frame(ini)

  for (m in names(modules)){
    tryCatch(expr = {
      srt = GeneToEnrichment(srt, db = modules[m], method = 'rand', db_rand = modules_rand[m])
    }, error = function(e){c()})
  }

  scores = srt@meta.data[,names(modules)]
  frequency = colMeans(scores > 0.5, na.rm = TRUE)
  srt$state = apply(srt@meta.data[,names(modules)], 1, function(x){
    top = which.max(x)
    return(names(x)[top])
  })
  srt$state = factor(srt$state, levels = names(modules))
  save(srt, file = 'srt.scored.RData')

  pdf('scoring.pdf')
  plot.list = lapply(names(modules), function(m){
    h = FeaturePlot(srt, features = m, pt.size = 1, cols = c('grey',colors.module[m])) +
      NoAxes() + NoLegend() +
      theme(plot.title = element_text(size = 10))
  })
  print(CombinePlots(plot.list, ncol = 4))
  h = DimPlot(srt, group.by = 'state', cols = colors.module)
  print(h)
  # Co-occurence
  occ = colMeans(scores > 0.5)
  coocc = sapply(names(modules), function(m1){
    sapply(names(modules), function(m2){
      pval = phyper(
        sum((scores[, m1] > 0.5) * (scores[, m2] > 0.5)),
        sum(scores[, m2] > 0.5),
        dim(scores)[1] - sum(scores[, m2] > 0.5),
        sum(scores[, m1] > 0.5),
        lower.tail = FALSE)
      if (pval <= 0.5){
        return(-log10(pval))
      }
      if (pval > 0.5){
        return(log10(1-pval))
      }
    })
  })
  coocc[is.infinite(coocc) & coocc > 0] = max(coocc[is.finite(coocc)])
  coocc[is.infinite(coocc) & coocc < 0] = min(coocc[is.finite(coocc)])
  coocc[which(occ <= 0.01),] = 0
  coocc[, which(occ <= 0.01)] = 0
  coocc[coocc > pval_thresh & coocc < - pval_thresh] = 0
  coocc[coocc > 10] = 10
  coocc[coocc < -10] = -10
  g = graph_from_adjacency_matrix(coocc, diag = FALSE, weighted = TRUE, mode = 'undirected')
  coords = layout.circle(g)
  E(g)$color = colorRamp2(breaks = seq(-10, 10, length = 9), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 9))(E(g)$weight)
  E(g)$width = abs(E(g)$weight)/5
  V(g)$color = colors.module[names(V(g))]
  V(g)$size = occ*200
  plot(g,
       layout = coords,
       vertex.shape = 'circle',
       vertex.label.size = 1, vertex.label.color = 'black',
       frame=FALSE)
  dev.off()

  stats = ScoreModule(srt, db = modules, method = 'rand', db_rand = modules_rand)
  stats = rbind(stats, frequency)
  save(stats, file = 'stats.RData')

  setwd('~/Documents/Analysis/Tumors/')
  return()
}, mc.cores = ncores)


## Paired

mclapply(files.paired.list, function(f){
  
  setwd(paste0(f[['normal']],'Paired'))
  
  srt.normal = loadRData(paste0(f['normal'],'srt.malignant.RData'))
  srt.tumor = loadRData(paste0(f['tumor'],'srt.malignant.RData'))
  
  stats.normal = ScoreModule(srt.normal, db = modules, method = 'rand', nrand = nrand, nbin = nbin)
  rownames(stats.normal) = paste0(rownames(stats.normal),'.normal')
  stats.tumor = ScoreModule(srt.tumor, db = modules, method = 'rand', nrand = nrand, nbin = nbin)
  rownames(stats.tumor) = paste0(rownames(stats.tumor),'.tumor')

  srt = merge(srt.normal, srt.tumor)
  srt$cat = 'tumor'
  srt$cat[1:ncol(srt.normal)] = 'normal'
  
  srt = SCTransform(srt, return.only.var.genes = FALSE)
  srt = RunPCA(srt)
  srt = RunUMAP(srt, dims = 1:10)
  
  modules_rand = MakeRand(srt, db = modules, nrand = nrand, nbin = nbin)
  
  ini = matrix(0,nrow = ncol(srt), ncol = length(modules))
  rownames(ini) = colnames(srt)
  colnames(ini) = names(modules)
  srt@meta.data[,names(modules)] = as.data.frame(ini)
  
  for (m in names(modules)){
    tryCatch(expr = {
      srt = GeneToEnrichment(srt, db = modules[m], method = 'rand', db_rand = modules_rand[m])
    }, error = function(e){c()})
  }
  
  scores = srt@meta.data[,names(modules)]
  frequency.normal = colMeans(scores[srt$cat == 'normal',] > 0.5, na.rm = TRUE)
  frequency.tumor = colMeans(scores[srt$cat == 'tumor',] > 0.5, na.rm = TRUE)
  srt$state = apply(srt@meta.data[,names(modules)], 1, function(x){
    top = which.max(x)
    return(names(x)[top])
  })
  srt$state = factor(srt$state, levels = names(modules))
  save(srt, file = 'srt.scored.RData')
  
  pdf('scoring.pdf')
  plot.list = lapply(names(modules), function(m){
    h = FeaturePlot(srt, features = m, pt.size = 1, cols = c('grey',colors.module[m]), split.by = 'cat') +
      NoAxes() + NoLegend() +
      theme(plot.title = element_text(size = 10))
  })
  print(CombinePlots(plot.list, ncol = 4))
  h = DimPlot(srt, group.by = 'state', cols = colors.module, split.by = 'cat')
  print(h)
  h = DimPlot(srt, group.by = 'cat', cols = colors.cat)
  print(h)
  # Co-occurence
  for (ca in c('normal','tumor')){
    occ = colMeans(scores[srt$cat == ca,] > 0.5)
    coocc = sapply(names(modules), function(m1){
      sapply(names(modules), function(m2){
        pval = phyper(
          sum((scores[srt$cat == ca, m1] > 0.5) * (scores[srt$cat == ca, m2] > 0.5)), 
          sum(scores[srt$cat == ca, m2] > 0.5), 
          dim(scores[srt$cat == ca,])[1] - sum(scores[srt$cat == ca, m2] > 0.5),
          sum(scores[srt$cat == ca, m1] > 0.5), 
          lower.tail = FALSE)
        if (pval <= 0.5){
          return(-log10(pval))
        }
        if (pval > 0.5){
          return(log10(1-pval))
        }
      })
    })
    coocc[is.infinite(coocc) & coocc > 0] = max(coocc[is.finite(coocc)])
    coocc[is.infinite(coocc) & coocc < 0] = min(coocc[is.finite(coocc)])
    coocc[which(occ <= 0.01),] = 0
    coocc[, which(occ <= 0.01)] = 0
    coocc[coocc > pval_thresh & coocc < - pval_thresh] = 0
    coocc[coocc > 10] = 10
    coocc[coocc < -10] = -10
    g = graph_from_adjacency_matrix(coocc, diag = FALSE, weighted = TRUE, mode = 'undirected')
    coords = layout.circle(g)
    E(g)$color = colorRamp2(breaks = seq(-10, 10, length = 9), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 9))(E(g)$weight)
    E(g)$width = abs(E(g)$weight)/5
    V(g)$color = colors.module[names(V(g))]
    V(g)$size = occ*200
    plot(g, 
         layout = coords,
         vertex.shape = 'circle', 
         vertex.label.size = 1, vertex.label.color = 'black',
         frame=FALSE) 
  }
  dev.off()
  
  stats = rbind(stats.normal, stats.tumor, frequency.normal, frequency.tumor)
  save(stats, file = 'stats.RData')
  
  setwd('~/Documents/Analysis/Tumors/')
  return()
}, mc.cores = ncores)

