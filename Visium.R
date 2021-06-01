#!/usr/bin/env Rscript
source('~/Documents/R/seurat_functions.R')
setwd('~/Documents/Analysis/Visium/')

load('~/Documents/Analysis/Tumors/parameters.RData')
w = which(names(colors) == 'LIHCHs1')
l = colors[['LIHCHs1']]
colors = c(colors, 'LIHCHs1' = l, 'LIHCHs1A' = l,  'LIHCHs1B' = l)
colors = colors[-w]
nrand = 2
nbin = 10
cats = c('Malignant','Both','Normal')
colors.cat = brewer_pal(palette = 'RdYlGn')(8)[c(2,5,7)]
names(colors.cat) = cats
method = 'score'
pval_thresh = 1

for (prox in c('inverse')){

modules = loadRData('../Tumors/modules.RData')
colors.module = loadRData('../Tumors/colors.module.RData')
types.module = loadRData('../Tumors/types.module.RData')
modules = modules[!names(modules) %in% c('AC','OPC','NPC')]
colors.module = colors.module[names(modules)]
types.module = types.module[names(modules)]
ma = names(modules)

colors.pop = loadRData('../Tumors/colors.pop.ours.RData')

colors.bin = rep(Seurat:::SpatialColors(n = 3)[c(1,3)],2)
names(colors.bin) = c('FALSE','TRUE','0','1')

srt.tme = loadRData('~/Documents/Analysis/Tumors/srt.ours.tme.RData')
exclude = c('Embryonic_stem_cells','iPS_cells','Tissue_stem_cells','Smooth_muscle_cells','Chondrocytes','Epithelial_cells','Hepatocytes','Keratinocytes','Astrocyte','Neurons','Neuroepithelial_cell')
srt.tme = srt.tme[, !srt.tme$pop %in% exclude]
rm = (srt.tme$pop == 'B_cell' & !colSums(GetData(srt.tme, slot = 'data')[c('CD4','GZMB'),]) == 0) # rm errors in B cell annotations
srt.tme = srt.tme[, !rm]
# keep = sample(1:ncol(srt.tme), size = 5000, replace = FALSE)
# keep = unlist(lapply(levels(srt.tme$pop), function(ct){
#   x = sample(colnames(srt.tme)[srt.tme$pop == ct])
#   if (length(x) > 1000){
#     x = sample(x, size = 1000, replace = FALSE)
#   }
#   return(x)
# }))
# srt.tme = srt.tme[,keep]

pdf('macrophages.pdf')
#srt.mac = srt.tme
srt.mac = loadRData('~/Documents/Tumors/OVCAHs/OVCAHs1/srt.normal.RData')
srt.mac = srt.mac[, srt.mac$pop == 'Macrophage']
srt.mac = SCTransform(srt.mac, return.only.var.genes = FALSE)
srt.mac = RunPCA(srt.mac)
srt.mac = RunUMAP(srt.mac, dims = 1:10)
srt.mac = FindNeighbors(srt.mac)
srt.mac = FindClusters(srt.mac)
srt.mac$mac = factor(srt.mac$seurat_clusters, labels = c('M1','M2'))
h = DimPlot(srt.mac, group.by = 'mac', cols = c('darkorange','darkred'))
print(h)
markers = FindAllMarkers2(srt.mac, group.by = 'mac', do.print = FALSE, do.enrichment = FALSE, do.plot = TRUE)
db_mac = split(markers$genes_100$gene, markers$genes_100$cluster)
print(db_mac)
data = GetData(srt.mac, slot = 'scale.data')[unlist(db_mac),]
top_ann = HeatmapAnnotation(df = data.frame('mac' = srt.mac$mac),
                            col = list('mac' = c('M1' = 'darkorange','M2' = 'darkred')),
                            which = 'column')
h = Heatmap(data, top_annotation = top_ann,
        cluster_rows = FALSE, cluster_columns = FALSE, 
        breaks = c(-2,0,2), colors = c('blue','white','red'), 
        column_split = srt.mac$mac,
        row_split = unlist(lapply(names(db_mac), function(x){rep(x,length(db_mac[[x]]))})))
labels = c('IFIT1','TNF','IL1B','C1QA','CD14','CD163','APOE','HLA-DRA','MS4A4A')
r = rowAnnotation(link = anno_mark(at = which(rownames(data) %in% labels), labels = labels))
print(h+r)
dev.off()

files = c(
  'OVCAHs/OVCAHs1',
  #'OVCAHs/OVCAHs2',
  'OVCAHs/OVCAHs3',
  'BRCAHs/BRCAHs0',
  'BRCAHs/BRCAHs1',
  'BRCAHs/BRCAHs2',
  'UCECHs/UCECHs3',
  'LIHCHs/LIHCHs1',
  #'LIHCHs/LIHCHs1A',
  #'LIHCHs/LIHCHs1B',
  'PDACHs/PDACHs1',
  #'PDACHs/PDACHs2',
  'GISTHs/GISTHs1',
  'GISTHs/GISTHs2')
names(files) = sapply(strsplit(files, '/', fixed = TRUE), '[', 2)

meta.list = mclapply(names(files), function(fa){
  
  pdf(paste0(fa,'.',prox,'.pdf'), height = 10, width = 15)
  f = files[[fa]]
  meta = c()
  
  #### Loading data ####
  
  # ST
  
  st = loadRData(paste0('~/Documents/Visium/',f,'/st.RData'))
  
  res = loadRData(paste0('~/Documents/Visium/',f,'/res.RData'))
  modules_st = NMFToModules(res, gmin = 5)
  nmf = names(modules_st)
  
  st = AddMetaData(st, metadata = st@images$slice1@coordinates[,c('row','col')], col.name = c('row','col'))
  st$axis = st$row
  if (fa == 'OVCAHs1'){
    st$axis = st$row
  }
  if (fa == 'OVCAHs2'){
    st$axis = st$col
  }
  if (fa == 'OVCAHs3'){
    st$axis = sqrt((st$row - 0)^2 + (st$col - 70)^2)
  }
  if (fa == 'BRCAHs0'){
    st$axis = st$row - 1/3*st$col
  }
  if (fa == 'BRCAHs1'){
    st$axis = st$row
  }
  if (fa == 'BRCAHs2'){
    st$axis = - sqrt((st$col - 80)^2 + ((st$row - 25)*1.5)^2)
  }
  if (fa == 'UCECHs3'){
    st$axis = st$col - 1.2*st$row
  }
  if (fa == 'LIHCHs1'){
    st$axis = st$row
  }
  if (fa == 'PDACHs1'){
    st$axis = st$col - st$row
  }
  if (fa == 'PDACHs2'){
    st$axis = st$col
  }
  if (fa == 'GISTHs1'){
    st$axis = - st$row - 3/4*st$col
  }
  if (fa == 'GISTHs2'){
    st$axis = st$row + 1/3*st$col
  }
  st$axis = (st$axis - min(st$axis))/(max(st$axis) - min(st$axis))
  h = SpatialFeaturePlot(st, 'axis')
  print(h)
  
  h = SpatialFeaturePlot(st, c('nCount_Spatial','nCount_SCT','nFeature_Spatial','nFeature_SCT'))
  print(h)
  
  colors.clone = brewer_pal(palette = 'Accent')(nlevels(st$clone))
  names(colors.clone) = levels(st$clone)
  
  h = SpatialDimPlot(st, group.by = 'clone', cols = colors.clone)
  print(h)
  
  # Single cell
  
  if (fa %in% c('LIHCHs1A','LIHCHs1B')){
    srt.malignant = loadRData(paste0('~/Documents/Tumors/LIHCHs/LIHCHs1/srt.malignant.RData'))
  } else {
    srt.malignant = loadRData(paste0('~/Documents/Tumors/',f,'/srt.malignant.RData'))
  }
  if (!all(sapply(ma, function(x){x %in% names(srt.malignant@meta.data)}))){
    srt.malignant = GeneToEnrichment(srt.malignant, db = modules, method = 'rand', nrand = nrand, nbin = nbin)
  }
  ma_bin = paste0(ma, '_bin')
  scores = srt.malignant@meta.data[,ma]
  scores_bin = scores
  scores_bin[] = as.numeric(scores_bin > 0.5)
  srt.malignant = AddMetaData(srt.malignant, metadata = scores_bin, col.name = ma_bin)
  
  if (fa %in% c('LIHCHs1A','LIHCHs1B')){
    srt.normal = loadRData(paste0('~/Documents/Tumors/LIHCHs/LIHCHs1/srt.normal.RData'))
  } else {
    srt.normal = loadRData(paste0('~/Documents/Tumors/',f,'/srt.normal.RData'))
  }
  if (f %in% c('GISTHs/GISTHs1')){
    exclude = c('Embryonic_stem_cells','iPS_cells','Tissue_stem_cells','Smooth_muscle_cells','Chondrocytes','Hepatocytes','Keratinocytes','Astrocyte','Neurons','Neuroepithelial_cell')
  } else {
    exclude = c('Embryonic_stem_cells','iPS_cells','Tissue_stem_cells','Smooth_muscle_cells','Chondrocytes','Epithelial_cells','Hepatocytes','Keratinocytes','Astrocyte','Neurons','Neuroepithelial_cell')
  }
  srt.normal = srt.normal[, !srt.normal$pop %in% exclude]
  ta = table(srt.normal$pop)
  spec = names(ta)[ta > 20]
  gen = setdiff(levels(srt.tme$pop), spec)
  srt.normal = merge2(srt.normal[,srt.normal$pop %in% spec], srt.tme[,srt.tme$pop %in% gen])
  srt.normal = SCTransform(srt.normal, return.only.var.genes = FALSE)
  
  srt = merge2(srt.malignant, srt.normal)
  srt = SCTransform(srt, return.only.var.genes = FALSE)
  srt$pop = as.character(srt$pop)
  srt$pop = factor(srt$pop, levels = intersect(names(colors.pop), srt$pop))
  srt$pop = relevel(srt$pop, ref = 'Malignant')
  Idents(srt) = 'pop'
  srt = SCTransform(srt, return.only.var.genes = FALSE)
  srt = RunPCA(srt)
  srt = RunUMAP(srt, dims = 1:10)
  srt = RunTSNE(srt, dims = 1:10)
  h = DimPlot(srt, reduction = 'umap', group.by = 'pop', cols = colors.pop)
  print(h)
  h = DimPlot(srt, reduction = 'tsne', group.by = 'pop', cols = colors.pop)
  print(h)
  
  markers = FindAllMarkers2(srt, group.by = 'pop')
  genes.srt = markers$genes_100$gene
  
  #### NMF ####
  
  plot.list = lapply(nmf, function(x){
    h = SpatialFeaturePlot(st, x)
  })
  print(CombinePlots(plot.list))
  
  #### Anchoring ####
  
  # Find pure cell types from integration labels
  
  options(future.globals.maxSize = 5000*1024^2)
  object.list = list('SC' = srt, 'ST' = st)
  genes.use = intersect(SpatiallyVariableFeatures(st), VariableFeatures(srt))
  #genes.use = intersect(SpatiallyVariableFeatures(st), genes.srt)
  object.list = PrepSCTIntegration(object.list = object.list, 
                                   anchor.features = genes.use, 
                                   verbose = FALSE)
  anchors = FindTransferAnchors(reference = object.list$SC, query = object.list$ST, 
                                normalization.method = 'SCT',
                                features = genes.use,
                                verbose = FALSE)
  predictions = TransferData(anchorset = anchors, refdata = srt$pop)
  predictions = predictions[,!colnames(predictions) %in% c('predicted.id','prediction.score.max')]
  colnames(predictions) = gsub('prediction.score','pred',colnames(predictions))
  pred = colnames(predictions)
  st = AddMetaData(st, metadata = predictions, col.name = pred)
  plot.list = lapply(pred, function(x){
    h = SpatialFeaturePlot(st, features = x)
  })
  print(CombinePlots(plot.list))
  
  pred = colnames(predictions)[apply(predictions, 2, var) > 0]
  pred = pred[order(pred)]
  pred = c('pred.Malignant', setdiff(pred, 'pred.Malignant'))
  
  # Binarize
  
  predictions_bin = predictions
  predictions_bin[] = as.numeric(predictions_bin > 0.9)
  pred_bin = paste0(pred, '_bin')
  st = AddMetaData(st, predictions_bin, col.name = pred_bin)
  plot.list = lapply(pred_bin, function(x){
    h = SpatialDimPlot(st, x, cols = colors.bin)
  })
  print(CombinePlots(plot.list))
  
  #### NNLS from single-cell ####
  
  # Data
  
  genes.use = intersect(SpatiallyVariableFeatures(st), VariableFeatures(srt))
  #genes.use = intersect(SpatiallyVariableFeatures(st), genes.srt)
  prof = AverageExpression(srt, assay = 'SCT', slot = 'data')$SCT[genes.use,]
  data = as.matrix(GetAssayData(st, assay = 'SCT', slot = 'data'))[genes.use,]
  
  # Regression
  
  coef = t(apply(data, 2, function(y){
    coef(nnls(as.matrix(prof), y))
  }))
  colnames(coef) = colnames(prof)
  nnls = colnames(coef)[colSums(coef > 0) >= 20]
  prof = prof[,nnls]
  coef = t(apply(data, 2, function(y){
    coef(nnls(as.matrix(prof), y))
  }))
  colnames(coef) = nnls
  st = AddMetaData(st, coef, col.name = nnls)
  plot.list = lapply(nnls, function(x){
    h = SpatialFeaturePlot(st, x)
  })
  print(CombinePlots(plot.list))
  plot.list = lapply(nnls, function(x){
    h = SpatialFeaturePlot(st, x)
  })
  print(CombinePlots(plot.list))
  
  # Scaling from null from pure cell types
  
  coef_scaled = sapply(nnls, function(x){
    vec = coef[,x]
    y = paste0('pred.',x)
    spots.use = colnames(st)[st@meta.data[,y] == 0]
    #spots.use = colnames(st)[st@meta.data[,y] == 0 & rowSums(predictions_bin) == 1]
    if (length(spots.use) < 5){
      spots.use = colnames(st)[order(st@meta.data[,y])[1:5]]
    }
    data_rand = Reduce(cbind, lapply(1:10^nrand, function(i){
      t(apply(data[,spots.use], 1, function(expr){
        sample(expr, length(expr), replace = FALSE)
      }))
    }))
    coef_rand = t(apply(data_rand, 2, function(y){
      coef(nnls(as.matrix(prof), y))
    }))
    colnames(coef_rand) = nnls
    if (sd(coef_rand[,x]) == 0){
      return((coef[,x] - mean(coef_rand[,x]))/min(coef[coef[,x] > 0,x]))
    } else {
      return((coef[,x] - mean(coef_rand[,x]))/sd(coef_rand[,x]))
    }
  })
  nnls_scaled = paste0(nnls, '_scaled')
  st = AddMetaData(st, coef_scaled, col.name = nnls_scaled)
  plot.list = lapply(nnls_scaled, function(x){
    h = SpatialFeaturePlot(st, x)
  })
  print(CombinePlots(plot.list))
  
  # Binarizing
  
  coef_bin = coef_scaled
  coef_bin[] = (coef_bin > 2)
  # for (x in rownames(coef_bin)){
  #   if (sum(coef_bin[x,]) == 0){
  #     y = coef_scaled[x,]
  #     y = sort(y, decreasing = TRUE)
  #     if (y[1] > 1 & y[1] - y[2] > 1){
  #       coef_bin[x,names(y)[1]] = 1
  #     }
  #   }
  # }
  
  nnls_bin = paste0(nnls, '_bin')
  st = AddMetaData(st, coef_bin, col.name = nnls_bin)
  plot.list = lapply(nnls_bin, function(x){
    h = SpatialDimPlot(st, x, cols = colors.bin)
  })
  print(CombinePlots(plot.list))
  
  nnls = setdiff(nnls, 'Epithelial_cells')
  nnls_scaled = paste0(nnls, '_scaled')
  nnls_bin = paste0(nnls, '_bin')
  
  # M1 vs M2
  
  st.mac = st[, st$Macrophage_bin == 1]
  st.mac = SCTransform(st.mac, return.only.var.genes = FALSE, assay = 'Spatial')
  st.mac = GeneToEnrichment(st.mac, db = db_mac, method = 'score', nbin = nbin)
  st.mac$M1M2 = st.mac$M1 - st.mac$M2
  st.mac$M1_bin = as.numeric(st.mac$M1 > st.mac$M2)
  st.mac$M2_bin = as.numeric(st.mac$M2 > st.mac$M1)
  st.mac$M1M2_bin = c('M1','M2')[1 + (st.mac$M1 > st.mac$M2)]
  
  h = SpatialFeaturePlot(st.mac, c('M1','M2','M1M2'))
  print(h)
  
  srt = GeneToEnrichment(srt, db = db_mac, method = 'score', nbin = nbin)
  srt$M1M2 = srt$M1 - srt$M2
  h = VlnPlot(srt, 'M1M2', group.by = 'pop', cols = colors.pop)
  print(h)
  
  st@meta.data[,c('M1','M2','M1M2','M1M2_bin')] = NA
  st@meta.data[,c('M1_bin','M2_bin')] = 0
  st@meta.data[colnames(st.mac),c('M1','M2','M1M2','M1M2_bin','M1_bin','M2_bin')] = st.mac@meta.data[,c('M1','M2','M1M2','M1M2_bin','M1_bin','M2_bin')]
  
  nnls = c(nnls, 'M1', 'M2')
  nnls_scaled = paste0(nnls, '_scaled')
  nnls_bin = paste0(nnls, '_bin')
  
  nnls1 = setdiff(nnls, c('M1', 'M2'))
  nnls1_scaled = paste0(nnls1, '_scaled')
  nnls1_bin = paste0(nnls1, '_bin')
  
  nnls2 = setdiff(nnls, 'Macrophage')
  nnls2_scaled = paste0(nnls2, '_scaled')
  nnls2_bin = paste0(nnls2, '_bin')
  
  # Average neighborhood coefficient
  
  nei = FindSTNeighbors(st, d_min = 0, d_max = 1.5)
  coef_nei = sapply(nnls_bin, function(x){
    sapply(nei, function(spot){
      y = st@meta.data[spot,x]
      y[y < 0] = 0
      mean(y, na.rm = TRUE)
    })
  })
  colnames(coef_nei) = nnls
  nnls_nei = paste0(nnls, '_nei')
  nnls1_nei = paste0(nnls1, '_nei')
  nnls2_nei = paste0(nnls2, '_nei')
  st = AddMetaData(st, coef_nei, col.name = nnls_nei)
  h = SpatialFeaturePlot(st, nnls_nei)
  print(h)
  
  # M1 M2
  
  M1M2_nei = sapply(nei, function(spot){
    y = st@meta.data[spot,'M1M2']
    mean(y, na.rm = TRUE)
  })
  M1M2_nei[is.nan(M1M2_nei)] = NA
  st = AddMetaData(st, M1M2_nei, col.name = 'M1M2_nei')
  h = SpatialFeaturePlot(st, c('M1','M2','M1M2','M1_nei','M2_nei','M1M2_nei'))
  print(h)
  
  #### Distances ####
  
  # Add distance to cell types
  
  coord = st@images$slice1@coordinates[,c('imagerow','imagecol')]
  distances = as.matrix(dist(coord))
  distances = distances/min(distances[distances > 0])
  #distances = round(distances, digits = 1)
  coef_dist = sapply(nnls_bin, function(x){
    w = colnames(st)[as.logical(st@meta.data[,x])]
    w = w[!is.na(w)]
    if (length(w) == 1){
      mi = as.numeric(distances[,w])
    } else {
      mi = apply(distances[,w], 1, function(x){min(as.numeric(x), na.rm = TRUE)})
    }
    return(mi)
  })
  if (prox == 'inverse'){
    coef_dist = 1/(1+coef_dist)
  }
  if (prox == 'opposite'){
    coef_dist = -coef_dist
  }
  colnames(coef_dist) = nnls
  nnls_dist = paste0(nnls, '_dist')
  nnls1_dist = paste0(nnls1, '_dist')
  nnls2_dist = paste0(nnls2, '_dist')
  st = AddMetaData(st, coef_dist, col.name = nnls_dist)
  plot.list = lapply(nnls_dist, function(x){
    h = SpatialFeaturePlot(st, x)
  })
  print(CombinePlots(plot.list))
  
  # Malignant
  
  depth = apply(coef_dist[,!colnames(coef_dist)=='Malignant'], 1, min)
  if (prox == 'inverse'){
    depth = 1/(1+depth)
  }
  if (prox == 'opposite'){
    depth = -depth
  }  
  st = AddMetaData(st, depth, col.name = 'Depth')
  h = SpatialFeaturePlot(st, 'Depth')
  print(h)
  
  # M1 M2
  
  st$M1M2_dist = st$M1_dist - st$M2_dist
  h = SpatialFeaturePlot(st, c('M1_dist','M2_dist','M1M2_dist'))
  print(h)
  
  #### Splitting ####
  
  # Categories
  
  st$cat = apply(st@meta.data[,nnls_bin], 1, function(x){
    if (x['Malignant_bin']){
      if (sum(x) == 1){
        return('Malignant')
      } else {
        return('Both')
      }
    } else {
      if (sum(x) == 0){
        return(NA)
      } else {
        return('Normal')
      }
    }
  })
  st$cat = factor(st$cat, levels = cats)
  h = SpatialDimPlot(st, 'cat', pt.size.factor = 1, stroke = 0)
  print(h + scale_fill_manual(values = colors.cat))
  print(table(st$cat))
  
  # Compare all metrics
  
  h = VlnPlot(st, nmf, group.by = 'cat', pt.size = 0)
  print(h)
  
  h = VlnPlot(st, pred, group.by = 'cat', pt.size = 0)
  print(h)
  
  h = VlnPlot(st, nnls1_scaled, group.by = 'cat', pt.size = 0)
  print(h)
  
  plot.list = lapply(nmf, function(x){
    h = FeatureScatter(st, 'axis', x, group.by = 'cat', cols = colors.cat, span = 0.3, plot.cor = FALSE)
  })
  print(CombinePlots(plot.list))
  
  plot.list = lapply(pred, function(x){
    h = FeatureScatter(st, 'axis', x, group.by = 'cat', cols = colors.cat, span = 0.3, plot.cor = FALSE)
  })
  print(CombinePlots(plot.list))
  
  plot.list = lapply(nnls1_scaled, function(x){
    h = FeatureScatter(st, 'axis', x, group.by = 'cat', cols = colors.cat, span = 0.3, plot.cor = FALSE)
  })
  print(CombinePlots(plot.list))
  
  tryCatch(expr = {
    co = sapply(nmf, function(x){
      sapply(pred, function(y){
        cor(st@meta.data[,x], st@meta.data[,y])
      })
    })
    h = Heatmap(co, name = 'Correlation',
                row_dend_reorder = 1:nrow(co), column_dend_reorder = 1:ncol(co),
                show_row_names = TRUE, show_column_names = TRUE,
                breaks = seq(-0.5, 0.5, length = 11), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 11))
    print(h)
  }, error = function(e){c()})
  
  tryCatch(expr = {
    co = sapply(pred, function(x){
      sapply(nnls1_scaled, function(y){
        cor(st@meta.data[,x], st@meta.data[,y])
      })
    })
    h = Heatmap(co, name = 'Correlation',
                row_dend_reorder = 1:nrow(co), column_dend_reorder = 1:ncol(co),
                show_row_names = TRUE, show_column_names = TRUE,
                breaks = seq(-1, 1, length = 11), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 11))
    print(h)
  }, error = function(e){c()})
  
  tryCatch(expr = {
    co = sapply(nmf, function(x){
      sapply(nnls1_scaled, function(y){
        cor(st@meta.data[,x], st@meta.data[,y])
      })
    })
    h = Heatmap(co, name = 'Correlation',
                row_dend_reorder = 1:nrow(co), column_dend_reorder = 1:ncol(co),
                show_row_names = TRUE, show_column_names = TRUE,
                breaks = seq(-0.5, 0.5, length = 11), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 11))
    print(h)
  }, error = function(e){c()})
  
  # Entropy
  
  ent = apply(st@meta.data[,nnls_bin], 1, function(x){
    x = x/sum(x)
    x = x[!x==0]
    return(sum(x*log(x)))
  })
  st = AddMetaData(st, ent, col.name = 'ent')
  h = SpatialFeaturePlot(st, 'ent', min.cutoff = -3, max.cutoff = 0)
  print(h)
  
  # Split
  
  st.list = SplitObject(st[,!is.na(st$cat)], 'cat')
  st.list = st.list[cats]
  st.malignant = st.list$Malignant
  st.malignant = SCTransform(st.malignant, return.only.var.genes = FALSE, assay = 'Spatial')  
  # st.malignant = FindSpatiallyVariableFeatures(st.malignant, selection.method = 'markvariogram')
  st.normal = st.list$Normal
  st.normal = SCTransform(st.normal, return.only.var.genes = FALSE, assay = 'Spatial')  
  # st.normal = FindSpatiallyVariableFeatures(st.normal, selection.method = 'markvariogram')
  
  h = SpatialFeaturePlot(st.malignant, nnls_nei)
  print(h)
  
  h = SpatialFeaturePlot(st.malignant, nnls_dist)
  print(h)
  
  
  #### Modules ####
  
  # Add modules
  
  st.malignant = GeneToEnrichment(st.malignant, db = modules, method = method, nrand = nrand, nbin = nbin)
  st.normal = GeneToEnrichment(st.normal, db = modules, method = method, nrand = nrand, nbin = nbin)
  
  plot.list = lapply(ma, function(m){
    h = SpatialFeaturePlot(st.malignant, m)
  })
  print(CombinePlots(plot.list))
  plot.list = lapply(ma, function(x){
    h = FeatureScatter(st.malignant, 'axis', x, group.by = 'orig.ident', cols = colors, span = 0.3, plot.cor = FALSE) + NoLegend()
  })
  print(CombinePlots(plot.list))
  
  h = Heatmap(cor(st.malignant@meta.data[,ma]), 
              show_row_names = TRUE, is.symmetric = TRUE,
              breaks = seq(-1,1, length = 7), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 7))
  print(h)
  h = Heatmap(cor(st.malignant@meta.data[,ma]), 
              show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE,
              breaks = seq(-1,1, length = 7), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 7))
  print(h)
  
  ma_bin = paste0(ma, '_bin')
  scores_bin = st.malignant@meta.data[,ma] > 0
  scores_bin[] = as.numeric(scores_bin)
  st.malignant = AddMetaData(st.malignant, scores_bin, col.name = ma_bin)
  plot.list = lapply(ma_bin, function(m){
    h = SpatialDimPlot(st.malignant, m, cols = colors.bin)
  })
  print(CombinePlots(plot.list))
  plot.list = lapply(ma_bin, function(x){
    h = FeatureScatter(st.malignant, 'axis', x, group.by = 'orig.ident', cols = colors, span = 0.3, plot.cor = FALSE) + NoLegend()
  })
  print(CombinePlots(plot.list))
  
  st@meta.data[,ma] = NA
  st@meta.data[colnames(st.malignant),ma] = st.malignant@meta.data[,ma]
  st@meta.data[,ma_bin] = 0
  st@meta.data[colnames(st.malignant),ma_bin] = st.malignant@meta.data[,ma_bin]
  
  # Module neighborhood
  
  nei = FindSTNeighbors(st, d_min = 0, d_max = 1.5)
  coef_nei = sapply(ma, function(x){
    sapply(nei, function(spot){
      y = st@meta.data[spot,x]
      mean(y, na.rm = TRUE)
    })
  })
  colnames(coef_nei) = ma
  ma_nei = paste0(ma, '_nei')
  st = AddMetaData(st, coef_nei, col.name = ma_nei)
  
  su = st@meta.data[,c(ma_nei, nnls1_nei)]
  su = su[apply(su, 1, function(x){!any(is.na(x))}),]
  co = cor(su)
  h = Heatmap(co, is.symmetric = TRUE,
              show_row_names = TRUE, show_column_names = TRUE,
              breaks = seq(-0.3,0.3,length = 7), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 7))
  print(h)
  co = cor(su)[ma_nei, nnls1_nei]
  h = Heatmap(co, 
              show_row_names = TRUE, show_column_names = TRUE,
              breaks = seq(-0.3,0.3,length = 7), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 7))
  print(h)
  
  su = st.malignant@meta.data[,c(ma, nnls1_nei)]
  su = su[apply(su, 1, function(x){!any(is.na(x))}),]
  co = cor(su)
  h = Heatmap(co, name = 'Correlation', is.symmetric = TRUE,
              show_row_names = TRUE, show_column_names = TRUE,
              breaks = seq(-0.3,0.3,length = 7), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 7))
  groups = cutree(as.hclust(row_dend(h)), k = 5)
  groups = groups[row_order(h)]
  groups = factor(groups, levels = unique(groups))
  h = Heatmap(co, name = 'Correlation', row_order = row_order(h), is.symmetric = TRUE,
              show_row_names = TRUE, show_column_names = TRUE,
              row_names_gp = gpar(col = c(colors.module, rep('black', length(nnls1_nei)))),
              breaks = seq(-0.3,0.3,length = 7), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 7))
  print(h)
  decorate_heatmap_body('Correlation', {
    l = c(0,cumsum(table(groups))[1:length(table(groups))])
    for (k in 1:length(l)){
      i = unit(l[k:(k+1)]/ncol(co), 'npc')
      j = unit(1-l[k:(k+1)]/nrow(co), 'npc')
      grid.lines(i, j[1], gp = gpar(lwd = 1, col = 'black'))
      grid.lines(i, j[2], gp = gpar(lwd = 1, col = 'black'))
      grid.lines(i[1], j, gp = gpar(lwd = 1, col = 'black'))
      grid.lines(i[2], j, gp = gpar(lwd = 1, col = 'black'))
    }
  })
  
  co = cor(su)[ma, nnls1_nei]
  h = Heatmap(co, 
              show_row_names = TRUE, show_column_names = TRUE,
              breaks = seq(-0.2,0.2,length = 7), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 7))
  print(h)
  
  # Module distance
  
  coord = st@images$slice1@coordinates[,c('imagerow','imagecol')]
  distances = as.matrix(dist(coord))
  distances = distances/min(distances[distances > 0])
  #distances = round(distances, digits = 1)
  coef_dist = sapply(ma_bin, function(x){
    w = colnames(st)[as.logical(st@meta.data[,x])]
    w = w[!is.na(w)]
    if (length(w) == 1){
      mi = as.numeric(distances[,w])
    } else {
      mi = apply(distances[,w], 1, function(x){min(as.numeric(x), na.rm = TRUE)})
    }
    return(mi)
  })
  coef_dist = -coef_dist
  colnames(coef_dist) = ma
  ma_dist = paste0(ma, '_dist')
  st = AddMetaData(st, coef_dist, col.name = ma_dist)
  plot.list = lapply(ma_dist, function(x){
    h = SpatialFeaturePlot(st, x)
  })
  print(CombinePlots(plot.list))
  
  # Plot by axis
  
  plot.list = lapply(ma, function(m){
    sub = st.malignant@meta.data
    model = lm(sub[,m] ~ sub[,'axis'])
    est = round(summary(model)$coefficients[2,1], digits = 1)
    pval = round(-log10(summary(model)$coefficients[2,4]), digits = 1)
    if (pval > pval_thresh){
      p = ggplot(sub, aes_string(x = 'axis', y = m, color = 'orig.ident')) +
        geom_jitter(show.legend = FALSE, height = 0.05, width = 0) +
        #geom_point(show.legend = FALSE) +
        geom_smooth(method='lm', show.legend = FALSE,
                    fullrange = TRUE) + 
        annotate('text', label = paste('pval', pval), x = 0.5, y = 0.5) + 
        scale_color_manual(values = colors.module[[m]], breaks = unique(st.malignant$orig.ident)) +
        ggtitle(m) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 0),
          axis.text.y = element_text(angle = 90),
          axis.title.y = element_text(size = 0),
          axis.title.x = element_text(size = 0))
    } else {
      p = ggplot(sub, aes_string(x = 'axis', y = m, color = 'orig.ident')) +
        geom_jitter(show.legend = FALSE, height = 0.05, width = 0) +
        #geom_point(show.legend = FALSE) +
        annotate('text', label = 'NS', x = 0.5, y = 0.5) + 
        scale_color_manual(values = colors.module[[m]], breaks = unique(st.malignant$orig.ident)) +
        ggtitle(m) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 0),
          axis.text.y = element_text(angle = 90),
          axis.title.y = element_text(size = 0),
          axis.title.x = element_text(size = 0))
    }
    return(p)
  })
  print(CombinePlots(plot.list, ncol = 4))
  
  #### Integrating ####
  
  tryCatch(expr = {
    
    # Integrating
    
    options(future.globals.maxSize = 5000*1024^2)
    object.list = list('SC' = srt, 'ST' = st)
    genes.use = intersect(VariableFeatures(srt), SpatiallyVariableFeatures(st))
    object.list = PrepSCTIntegration(object.list = object.list,
                                     anchor.features = genes.use,
                                     verbose = FALSE)
    anchors = FindIntegrationAnchors(object.list = object.list, reference = 1,
                                     normalization.method =  'SCT',
                                     anchor.features = genes.use,
                                     verbose = FALSE)
    integrated = IntegrateData(anchorset = anchors,
                               normalization.method = "SCT",
                               verbose = FALSE)
    integrated = RunPCA(integrated)
    integrated = RunUMAP(integrated, dims = 1:10)
    integrated = RunTSNE(integrated, dims = 1:10)
    
    # Plotting
    
    h = DimPlot(integrated, reduction = 'umap', group.by = c('pop'), label = TRUE, repel = TRUE, cols = colors.pop, pt.size = 2)
    print(h)
    h = DimPlot(integrated, reduction = 'umap', group.by = c('cat'), label = TRUE, repel = TRUE, cols = colors.cat, pt.size = 2)
    print(h)
    h = FeaturePlot(integrated, reduction = 'umap', features = 'axis', cols = brewer_pal(palette = 'Purples', direction = 1)(n = 9), pt.size = 2)
    print(h)
    plot.list = lapply(nnls, function(x){
      h = FeaturePlot(integrated, reduction = 'umap', x)
    })
    print(CombinePlots(plot.list))
    plot.list = lapply(nnls_bin, function(x){
      h = DimPlot(integrated, reduction = 'umap', group.by = x, cols = colors.bin)
    })
    print(CombinePlots(plot.list))
    plot.list = lapply(pred, function(x){
      h = FeaturePlot(integrated, reduction = 'umap', features = x)
    })
    print(CombinePlots(plot.list))
    plot.list = lapply(pred_bin, function(x){
      h = DimPlot(integrated, reduction = 'umap', group.by = x, cols = colors.bin)
    })
    print(CombinePlots(plot.list))
    plot.list = lapply(nmf, function(x){
      h = FeaturePlot(integrated, reduction = 'umap', x)
    })
    print(CombinePlots(plot.list))
    for (m in ma){
      h = FeaturePlot(integrated, reduction = 'umap', features = m, cols = c('lightgrey',colors.module[m]), split.by = 'technique', pt.size = 2)
      print(h)
    }
    
    h = DimPlot(integrated, reduction = 'tsne', group.by = c('pop'), label = TRUE, repel = TRUE, cols = colors.pop, pt.size = 2)
    print(h)
    h = DimPlot(integrated, reduction = 'tsne', group.by = c('cat'), label = TRUE, repel = TRUE, cols = colors.cat, pt.size = 2)
    print(h)
    h = FeaturePlot(integrated, reduction = 'tsne', features = 'axis', cols = brewer_pal(palette = 'Purples', direction = 1)(n = 9), pt.size = 2)
    print(h)
    plot.list = lapply(nnls, function(x){
      h = FeaturePlot(integrated, reduction = 'tsne', x)
    })
    print(CombinePlots(plot.list))
    plot.list = lapply(nnls_bin, function(x){
      h = DimPlot(integrated, reduction = 'tsne', group.by = x, cols = colors.bin)
    })
    print(CombinePlots(plot.list))
    plot.list = lapply(pred, function(x){
      h = FeaturePlot(integrated, reduction = 'tsne', features = x)
    })
    print(CombinePlots(plot.list))
    plot.list = lapply(pred_bin, function(x){
      h = DimPlot(integrated, reduction = 'tsne', group.by = x, cols = colors.bin)
    })
    print(CombinePlots(plot.list))
    plot.list = lapply(nmf, function(x){
      h = FeaturePlot(integrated, reduction = 'tsne', x)
    })
    print(CombinePlots(plot.list))
    for (m in ma){
      h = FeaturePlot(integrated, reduction = 'tsne', features = m, cols = c('lightgrey',colors.module[m]), split.by = 'technique', pt.size = 2)
      print(h)
    }
    
  }, error = function(e){warning(e)})
  
  dev.off()
  print(fa)
  meta = st@meta.data
  return(meta)
  
}, mc.cores = length(files))

names(meta.list) = names(files)
save(meta.list, file = paste0('meta.list.',prox,'.RData'))

}

for (prox in c('inverse')){

pdf(paste0('Visium.',prox,'.pdf'), height = 7, width = 7)

load('~/Documents/Analysis/Tumors/parameters.RData')
w = which(names(colors) == 'LIHCHs1')
l = colors[['LIHCHs1']]
colors = c(colors, 'LIHCHs1' = l, 'LIHCHs1A' = l,  'LIHCHs1B' = l)
colors = colors[-w]
nrand = 2
nbin = 10
cats = c('Malignant','Both','Normal')
colors.cat = brewer_pal(palette = 'RdYlGn')(8)[c(2,5,7)]
names(colors.cat) = cats
method = 'score'
pval_thresh = 1

modules = loadRData('../Tumors/modules.RData')
colors.module = loadRData('../Tumors/colors.module.RData')
types.module = loadRData('../Tumors/types.module.RData')
modules = modules[!names(modules) %in% c('AC','OPC','NPC')]
colors.module = colors.module[names(modules)]
types.module = types.module[names(modules)]
ma = names(modules)
ma_bin = paste0(ma, '_bin')
ma_nei = paste0(ma, '_nei')
ma_dist = paste0(ma, '_dist')

load(paste0('meta.list.',prox,'.RData'))
meta.list = meta.list[sapply(meta.list, length) > 0]
meta.list = meta.list[!names(meta.list) %in% c('LIHCHs1A','LIHCHs1B')]

# Make meta, adding NAs
keep = Reduce(union, sapply(meta.list, colnames))
meta.list = lapply(meta.list, function(meta){
  if (all(keep %in% colnames(meta))){
    return(meta)
  } else {
    add = matrix(NA, nrow = nrow(meta), ncol = length(setdiff(keep, colnames(meta))))
    colnames(add) = setdiff(keep, colnames(meta))
    rownames(add) = rownames(meta)
    meta = cbind(meta, add)
    return(meta)
  }
})
meta = Reduce(rbind, lapply(meta.list, function(meta){meta[, keep]}))
meta = meta[!is.na(meta$cat),]
meta$orig.ident = factor(meta$orig.ident, levels = intersect(names(colors), unique(meta$orig.ident)))
meta$cat = factor(meta$cat, levels = cats)
nnls_dist = setdiff(grep('_dist', colnames(meta), value = TRUE), ma_dist)
nnls_nei = grep('_nei', colnames(meta), value = TRUE)
nnls = gsub('_dist','',nnls_dist)

nnls1 = setdiff(nnls, c('M1', 'M2', 'M1M2'))
nnls1_dist = paste0(nnls1, '_dist')
nnls1_nei = paste0(nnls1, '_nei')
nnls1_bin = paste0(nnls1, '_bin')

nnls2 = setdiff(nnls, 'Macrophage')
nnls2_dist = paste0(nnls2, '_dist')
nnls2_nei = paste0(nnls2, '_nei')
nnls2_bin = paste0(nnls2, '_bin')

# Sample composition

for (o in levels(meta$orig.ident)){
  b = meta[meta$orig.ident == o,]
  sp = t(sapply(c('Normal','Both'), function(s){
    a = b[b$cat == s, nnls1_bin]
    a[a < 0] = 0
    a = colMeans(a, na.rm = TRUE)
    a = a[!is.nan(a)]
    a = a[!names(a) == 'Malignant_bin']
    a = a/sum(a)
  }))
  colnames(sp) = gsub('_bin', '', colnames(sp))
  sp = sp[,rev(colnames(sp))]
  connectedBarplot(t(rbind(sp['Normal',], sp['Both',])), space = 0.2,
                   main = o, names.arg = c('',''), axes = FALSE,
                   las = 2, beside = FALSE, col = colors.pop[colnames(sp)])
}


comp.list = lapply(levels(meta$orig.ident), function(o){
  b = meta[meta$orig.ident == o,]
  sp = t(sapply(c('Normal','Both'), function(s){
    a = b[b$cat == s, nnls1_bin]
    a[a < 0] = 0
    a = colMeans(a, na.rm = TRUE)
    a = a[!is.nan(a)]
    a = a[!names(a) == 'Malignant_bin']
    a = a/sum(a)
  }))
  colnames(sp) = gsub('_bin', '', colnames(sp))
  return(sp)
})
names(comp.list) = levels(meta$orig.ident)
comp = melt(comp.list)
names(comp) = c('cat','ct','value','orig.ident')
p.list = lapply(levels(comp$ct), function(ct){
  sub = comp[comp$ct == ct,]
  sub = sub[order(sub$cat), ]
  p = ggpaired(sub, x = 'cat', y = 'value',
               color = 'cat', line.color = 'grey', line.size = 0.4) +
    stat_compare_means(paired = TRUE) +
    scale_color_manual(values = colors.cat) +
    ggtitle(ct) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      legend.position = 'none')
})
print(CombinePlots(p.list), ncol = 4)

comp.list = lapply(levels(meta$orig.ident), function(o){
  b = meta[meta$orig.ident == o,]
  sp = t(sapply(c('Normal','Both'), function(s){
    a = b[b$cat == s, nnls1_bin]
    a[a < 0] = 0
    a = colMeans(a, na.rm = TRUE)
    a = a[!is.nan(a)]
    a = a[!names(a) == 'Malignant_bin']
    a = a/sum(a)
  }))
  colnames(sp) = gsub('_bin', '', colnames(sp))
  log2(sp['Both',]/sp['Normal',])
})
keep = Reduce(union, sapply(comp.list, names))
comp.list = lapply(comp.list, function(comp){
  if (all(keep %in% names(comp))){
    return(comp)
  } else {
    add = matrix(NA, nrow = 1, ncol = length(setdiff(keep, names(comp))))
    names(add) = setdiff(keep, names(comp))
    comp = c(comp, add)
    return(comp)
  }
})
comp = Reduce(rbind, lapply(comp.list, function(comp){comp[keep]}))
comp[is.infinite(comp)] = 5
rownames(comp) = levels(meta$orig.ident)
barplot(comp, beside = TRUE, col = colors[rownames(comp)], las = 2)

# Modules by distance

val = Reduce(rbind,lapply(ma, function(m){
  Reduce(rbind,lapply(nnls1, function(ct){
    res = sapply(levels(meta$orig.ident), function(o){
      val = NA
      tryCatch(expr = {
        ct_dist = paste0(ct, '_dist')
        sub = meta[meta$orig.ident == o & meta$cat == 'Malignant',]
        model = lm(sub[,m] ~ sub[,ct_dist])
        est = summary(model)$coefficients[2,1]
        pval = -log10(summary(model)$coefficients[2,4])
        val = pval*sign(est)
        return(val)
      }, error = function(e){return(NA)})
    })
    #res[abs(res)<pval_thresh] = 0
    #res = sign(res)
    #med = sum(res, na.rm = TRUE)
    #amp = sum(res == sign(med), na.rm = TRUE)
    med = median(res, na.rm = TRUE)
    amp = mean(sign(res) == sign(med), na.rm = TRUE)
    return(data.frame('module' = m, 'ct' = ct, 'med' = med, 'amp' = amp))
    return(data.frame('module' = m, 'ct' = ct, 'med' = med, 'amp' = amp))
  }))
}))
val$module = factor(val$module, levels = rev(ma))
val$ct = factor(val$ct, levels = nnls)
val$highlight = (val$amp >= 0.5) & (abs(val$med) >= 0.75)
val$med[val$med > 4] = 4
val$med[val$med < -1] = -1
p = ggplot(val, aes(x=ct,y=module)) +
  geom_point(shape=21,aes(size=amp,fill=med,colour=highlight)) +
  scale_fill_gradient2(low = brewer_pal(palette = 'RdBu', direction = -1)(n = 3)[1],
                       mid = brewer_pal(palette = 'RdBu', direction = -1)(n = 3)[2], 
                       high = brewer_pal(palette = 'RdBu', direction = -1)(n = 3)[3]) +
  scale_colour_manual(breaks=c('FALSE','TRUE'),values=c('white','black')) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size=14, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.text.x = element_text(size=12, colour = "black", angle = 90, hjust = 1),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
print(p)

for (ct in nnls){
  val = Reduce(rbind,lapply(ma, function(m){
    res = sapply(levels(meta$orig.ident), function(o){
      val = NA
      tryCatch(expr = {
        ct_dist = paste0(ct, '_dist')
        sub = meta[meta$orig.ident == o & meta$cat == 'Malignant',]
        model = lm(sub[,m] ~ sub[,ct_dist])
        est = summary(model)$coefficients[2,1]
        pval = -log10(summary(model)$coefficients[2,4])
        pval[pval > 10] = 10
        val = pval*sign(est)
        return(val)
      }, error = function(e){return(NA)})
    })
    return(data.frame('orig.ident' = levels(meta$orig.ident), 'module' = m, 'res' = res))
  }))
  val$module = factor(val$module, levels = ma)
  val$orig.ident = factor(val$orig.ident)
  p = ggplot(val, aes(x=module,y=res)) +
    geom_boxplot(coef = 0, outlier.shape = NA) +
    #geom_violin() + 
    geom_jitter(aes(colour = orig.ident, size = 2), width = 0.25) +
    geom_hline(yintercept = log10(0.05), linetype = 2) + 
    geom_hline(yintercept = log10(0.05), linetype = 2) + 
    scale_colour_manual(values = colors, breaks = names(colors)) +
    ggtitle(ct) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size=14, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.text.x = element_text(size=12, colour = "black", angle = 90, hjust = 1),
          axis.title=element_blank()) + NoLegend()
  print(p)
}

for (m in ma){
  val = Reduce(rbind,lapply(nnls1, function(ct){
    res = sapply(levels(meta$orig.ident), function(o){
      val = NA
      tryCatch(expr = {
        ct_dist = paste0(ct, '_dist')
        sub = meta[meta$orig.ident == o & meta$cat == 'Malignant',]
        model = lm(sub[,m] ~ sub[,ct_dist])
        est = summary(model)$coefficients[2,1]
        pval = -log10(summary(model)$coefficients[2,4])
        pval[pval > 10] = 10
        val = pval*sign(est)
        return(val)
      }, error = function(e){return(NA)})
    })
    return(data.frame('orig.ident' = levels(meta$orig.ident), 'ct' = ct, 'res' = res))
  }))
  val$ct = factor(val$ct, levels = nnls)
  val$orig.ident = factor(val$orig.ident)
  p = ggplot(val, aes(x=ct,y=res)) +
    geom_boxplot(coef = 0, outlier.shape = NA) +
    #geom_violin() + 
    geom_jitter(aes(colour = orig.ident, size = 2), width = 0.25) +
    geom_hline(yintercept = log10(0.05), linetype = 2) + 
    geom_hline(yintercept = -log10(0.05), linetype = 2) + 
    scale_colour_manual(values = colors, breaks = names(colors)) +
    ggtitle(m) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size=14, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.text.x = element_text(size=12, colour = "black", angle = 90, hjust = 1),
          axis.title=element_blank()) + NoLegend()
  print(p)
}

for (ct in nnls){
  for (m in ma){
    
    ct_nei = paste0(ct, '_nei')
    ct_dist = paste0(ct, '_dist')
    
    plot.list = lapply(levels(meta$orig.ident), function(o){
      p = NULL
      tryCatch(expr = {
        sub = meta[meta$orig.ident == o & meta$cat == 'Malignant',]
        model = lm(sub[,m] ~ sub[,ct_dist])
        est = round(summary(model)$coefficients[2,1], digits = 1)
        pval = round(-log10(summary(model)$coefficients[2,4]), digits = 1)
        if (pval > pval_thresh){
          p = ggplot(sub, aes_string(x = ct_dist, y = m, color = 'orig.ident')) +
            geom_jitter(show.legend = FALSE, height = 0.05, width = 0) +
            #geom_point(show.legend = FALSE) +
            geom_smooth(method='lm', show.legend = FALSE,
                        fullrange = TRUE) + 
            annotate('text', label = paste('pval', pval), x = 0.05, y = 0.05) + 
            scale_color_manual(values = colors) +
            ggtitle(o) +
            theme_classic() +
            theme(
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size = 0),
              axis.text.y = element_text(angle = 90),
              axis.title.y = element_text(),
              axis.title.x = element_text())
        } else {
          p = ggplot(sub, aes_string(x = ct_dist, y = m, color = 'orig.ident')) +
            geom_jitter(show.legend = FALSE, height = 0.05, width = 0) +
            #geom_point(show.legend = FALSE) +
            annotate('text', label = 'NS', x = 0.05, y = 0.05) + 
            scale_color_manual(values = colors) +
            ggtitle(o) +
            theme_classic() +
            theme(
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size = 0),
              axis.text.y = element_text(angle = 90),
              axis.title.y = element_text(),
              axis.title.x = element_text())
        }
      }, error = function(e){c()})
      return(p)
    })
    print(CombinePlots(plot.list, ncol = 3))
  }
}

# Modules by neighborhood

val = Reduce(rbind,lapply(ma, function(m){
  Reduce(rbind,lapply(nnls1, function(ct){
    res = sapply(levels(meta$orig.ident), function(o){
      val = NA
      tryCatch(expr = {
        ct_nei = paste0(ct, '_nei')
        sub = meta[meta$orig.ident == o & meta$cat == 'Malignant',]
        model = lm(sub[,m] ~ sub[,ct_nei])
        est = summary(model)$coefficients[2,1]
        pval = -log10(summary(model)$coefficients[2,4])
        val = pval*sign(est)
        return(val)
      }, error = function(e){return(NA)})
    })
    #res[abs(res)<pval_thresh] = 0
    #res = sign(res)
    #med = sum(res, na.rm = TRUE)
    #amp = sum(res == sign(med), na.rm = TRUE)
    med = median(res, na.rm = TRUE)
    amp = mean(sign(res) == sign(med), na.rm = TRUE)
    return(data.frame('module' = m, 'ct' = ct, 'med' = med, 'amp' = amp))
    return(data.frame('module' = m, 'ct' = ct, 'med' = med, 'amp' = amp))
  }))
}))
val$module = factor(val$module, levels = rev(ma))
val$ct = factor(val$ct, levels = nnls)
val$highlight = (val$amp >= 0.5) & (abs(val$med) >= 0.75)
val$med[val$med > 4] = 4
val$med[val$med < -1] = -1
p = ggplot(val, aes(x=ct,y=module)) +
  geom_point(shape=21,aes(size=amp,fill=med,colour=highlight)) +
  scale_fill_gradient2(low = brewer_pal(palette = 'RdBu', direction = -1)(n = 3)[1],
                         mid = brewer_pal(palette = 'RdBu', direction = -1)(n = 3)[2], 
                         high = brewer_pal(palette = 'RdBu', direction = -1)(n = 3)[3]) +
  scale_colour_manual(breaks=c('FALSE','TRUE'),values=c('white','black')) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size=14, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.text.x = element_text(size=12, colour = "black", angle = 90, hjust = 1),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
print(p)


for (ct in nnls){
  val = Reduce(rbind,lapply(ma, function(m){
    res = sapply(levels(meta$orig.ident), function(o){
      val = NA
      tryCatch(expr = {
        ct_nei = paste0(ct, '_nei')
        sub = meta[meta$orig.ident == o & meta$cat == 'Malignant',]
        model = lm(sub[,m] ~ sub[,ct_nei])
        est = summary(model)$coefficients[2,1]
        pval = -log10(summary(model)$coefficients[2,4])
        pval[pval > 10] = 10
        val = pval*sign(est)
        return(val)
      }, error = function(e){return(NA)})
    })
    return(data.frame('orig.ident' = levels(meta$orig.ident), 'module' = m, 'res' = res))
  }))
  val$module = factor(val$module, levels = ma)
  val$orig.ident = factor(val$orig.ident)
  p = ggplot(val, aes(x=module,y=res)) +
    geom_boxplot(coef = 0, outlier.shape = NA) +
    #geom_violin() + 
    geom_jitter(aes(colour = orig.ident, size = 2), width = 0.25) +
    geom_hline(yintercept = log10(0.05), linetype = 2) + 
    geom_hline(yintercept = -log10(0.05), linetype = 2) + 
    scale_colour_manual(values = colors, breaks = names(colors)) +
    ggtitle(ct) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size=14, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.text.x = element_text(size=12, colour = "black", angle = 90, hjust = 1),
          axis.title=element_blank()) + NoLegend()
  print(p)
}

for (m in ma){
  val = Reduce(rbind,lapply(nnls1, function(ct){
    res = sapply(levels(meta$orig.ident), function(o){
      val = NA
      tryCatch(expr = {
        ct_nei = paste0(ct, '_nei')
        sub = meta[meta$orig.ident == o & meta$cat == 'Malignant',]
        model = lm(sub[,m] ~ sub[,ct_nei])
        est = summary(model)$coefficients[2,1]
        pval = -log10(summary(model)$coefficients[2,4])
        pval[pval > 10] = 10
        val = pval*sign(est)
        return(val)
      }, error = function(e){return(NA)})
    })
    return(data.frame('orig.ident' = levels(meta$orig.ident), 'ct' = ct, 'res' = res))
  }))
  val$ct = factor(val$ct, levels = nnls)
  val$orig.ident = factor(val$orig.ident)
  p = ggplot(val, aes(x=ct,y=res)) +
    geom_boxplot(coef = 0, outlier.shape = NA) +
    #geom_violin() + 
    geom_jitter(aes(colour = orig.ident, size = 2), width = 0.25) +
    geom_hline(yintercept = log10(0.05), linetype = 2) + 
    geom_hline(yintercept = log10(0.05), linetype = 2) + 
    scale_colour_manual(values = colors, breaks = names(colors)) +
    ggtitle(m) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size=14, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.text.x = element_text(size=12, colour = "black", angle = 90, hjust = 1),
          axis.title=element_blank()) + NoLegend()
  print(p)
}


for (ct in nnls){
  for (m in ma){
    
    ct_nei = paste0(ct, '_nei')
    ct_dist = paste0(ct, '_dist')
    
    plot.list = lapply(levels(meta$orig.ident), function(o){
      tryCatch(expr = {
        p = NULL
        sub = meta[meta$orig.ident == o & meta$cat == 'Malignant',]
        model = lm(sub[,m] ~ sub[,ct_nei])
        est = round(summary(model)$coefficients[2,1], digits = 1)
        pval = round(-log10(summary(model)$coefficients[2,4]), digits = 1)
        if (pval > pval_thresh){
          p = ggplot(sub, aes_string(x = ct_nei, y = m, color = 'orig.ident')) +
            geom_jitter(show.legend = FALSE, height = 0.05, width = 0) +
            #geom_point(show.legend = FALSE) +
            geom_smooth(method='lm', show.legend = FALSE,
                        fullrange = TRUE) + 
            annotate('text', label = paste('pval', pval), x = 0.05, y = 0.05) + 
            scale_color_manual(values = colors) +
            ggtitle(o) +
            theme_classic() +
            theme(
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size = 0),
              axis.text.y = element_text(angle = 90),
              axis.title.y = element_text(),
              axis.title.x = element_text())
        } else {
          p = ggplot(sub, aes_string(x = ct_nei, y = m, color = 'orig.ident')) +
            geom_jitter(show.legend = FALSE, height = 0.05, width = 0) +
            #geom_point(show.legend = FALSE) +
            annotate('text', label = 'NS', x = 0.05, y = 0.05) + 
            scale_color_manual(values = colors) +
            ggtitle(o) +
            theme_classic() +
            theme(
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size = 0),
              axis.text.y = element_text(angle = 90),
              axis.title.y = element_text(),
              axis.title.x = element_text())
        }
      }, error = function(e){c()})
      return(p)
    })
    print(CombinePlots(plot.list, ncol = 3))
  }
}

# M1 vs M2 distances

plot.list = lapply(levels(meta$orig.ident), function(o){
  sub = meta[meta$orig.ident == o,]
  p = ggplot(sub, aes_string(x = 'M1_dist', y = 'M2_dist', fill = 'cat', col = 'cat')) + 
    geom_point() +
    geom_smooth(method = 'lm') + 
    scale_fill_manual(values = colors.cat, breaks = cats) +
    scale_color_manual(values = colors.cat, breaks = cats) +
    ggtitle(o) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(angle = 90),
      axis.title.y = element_text(),
      axis.title.x = element_text())
})
print(CombinePlots(plot.list, ncol = 3))

sub = meta[,c('orig.ident','cat','M1M2_dist')]
sub = melt(sub, id.vars = c('orig.ident','cat'), measure.vars = 'M1M2_dist')
p.values <- sapply(split(sub, sub$orig.ident), function(x){wilcox.test(value~cat, x, subset = cat %in% c('Malignant','Normal'))$p.value})
labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*","n.s."))
y.values <- sapply(split(sub, sub$orig.ident), function(x){max(sapply(split(x, x$cat), function(xx){boxplot(xx$value, plot=F)$stats[4, ]}), na.rm = TRUE)}+0.5)
p = ggplot(sub, aes_string(x = 'orig.ident', y = 'value', fill = 'cat')) + 
  geom_boxplot(coef = 0, outlier.size = 0, outlier.shape = NA, varwidth = FALSE, show.legend = FALSE) +
  geom_signif(annotations = labels, xmin = (1:nlevels(sub$orig.ident))-0.2, xmax = (1:nlevels(sub$orig.ident))+0.2, y_position = y.values, tip_length = 0) + 
  scale_fill_manual(values = colors.cat, breaks = names(colors.cat)) +
  ggtitle('') +
  ylim(-2,4) + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(angle = 90),
    axis.text.x = element_text(angle = 90),
    axis.title.y = element_text(size = 0),
    axis.title.x = element_text(size = 0))
print(p)


# M1 vs M2 neighborhood

plot.list = lapply(levels(meta$orig.ident), function(o){
  sub = meta[meta$orig.ident == o,]
  p = ggplot(sub, aes_string(x = 'M1_nei', y = 'M2_nei', fill = 'cat', col = 'cat')) + 
    geom_point() +
    geom_smooth(method = 'lm') + 
    scale_fill_manual(values = colors.cat, breaks = names(colors.cat)) +
    scale_color_manual(values = colors.cat, breaks = names(colors.cat)) +
    ggtitle(o) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(angle = 90),
      axis.title.y = element_text(),
      axis.title.x = element_text())
})
print(CombinePlots(plot.list, ncol = 3))

sub = meta[,c('orig.ident','cat','M1M2_nei')]
sub = melt(sub, id.vars = c('orig.ident','cat'), measure.vars = 'M1M2_nei')
p.values <- sapply(split(sub, sub$orig.ident), function(x){wilcox.test(value~cat, x, subset = cat %in% c('Malignant','Normal'))$p.value})
labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*","n.s."))
y.values <- sapply(split(sub, sub$orig.ident), function(x){max(sapply(split(x, x$cat), function(xx){boxplot(xx$value, plot=F)$stats[4, ]}), na.rm = TRUE)}+0.05)
p = ggplot(sub, aes_string(x = 'orig.ident', y = 'value', fill = 'cat')) + 
  geom_boxplot(coef = 0, outlier.size = 0, outlier.shape = NA, varwidth = FALSE, show.legend = FALSE) +
  geom_signif(annotations = labels, xmin = (1:nlevels(sub$orig.ident))-0.2, xmax = (1:nlevels(sub$orig.ident))+0.2, y_position = y.values, tip_length = 0) + 
  scale_fill_manual(values = colors.cat, breaks = names(colors.cat)) +
  ggtitle('') +
  ylim(-0.2,0.4) + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(angle = 90),
    axis.text.x = element_text(angle = 90),
    axis.title.y = element_text(size = 0),
    axis.title.x = element_text(size = 0))
print(p)

# M1 vs M2

plot.list = lapply(levels(meta$orig.ident), function(o){
  sub = meta[meta$orig.ident == o,]
  p = ggplot(sub, aes_string(x = 'M1', y = 'M2', fill = 'cat', col = 'cat')) + 
    geom_point() +
    geom_smooth(method = 'lm') + 
    scale_fill_manual(values = colors.cat, breaks = names(colors.cat)) +
    scale_color_manual(values = colors.cat, breaks = names(colors.cat)) +
    ggtitle(o) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(angle = 90),
      axis.title.y = element_text(),
      axis.title.x = element_text())
})
print(CombinePlots(plot.list, ncol = 3))

sub = meta[,c('orig.ident','cat','M1M2')]
sub = sub[sub$cat %in% c('Both','Normal') & !is.na(sub$M1M2),]
sub = melt(sub, id.vars = c('orig.ident','cat'), measure.vars = 'M1M2')
sub$cat = factor(as.character(sub$cat), levels = c('Normal','Both'))
p.values <- sapply(split(sub, sub$orig.ident), function(x){wilcox.test(value~cat, x, subset = cat %in% c('Normal','Both'))$p.value})
labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*","n.s."))
y.values <- sapply(split(sub, sub$orig.ident), function(x){max(sapply(split(x, x$cat), function(xx){boxplot(xx$value, plot=F)$stats[4, ]}), na.rm = TRUE)}+0.05)
p = ggplot(sub, aes_string(x = 'orig.ident', y = 'value', fill = 'cat')) + 
  geom_boxplot(coef = 0, outlier.size = 0, outlier.shape = NA, varwidth = FALSE, show.legend = FALSE) +
  geom_signif(annotations = labels, xmin = (1:nlevels(sub$orig.ident))-0.2, xmax = (1:nlevels(sub$orig.ident))+0.2, y_position = y.values, tip_length = 0) + 
  scale_fill_manual(values = colors.cat, breaks = names(colors.cat)) +
  ggtitle('') +
  ylim(-0.2,0.4) + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(angle = 90),
    axis.text.x = element_text(angle = 90),
    axis.title.y = element_text(size = 0),
    axis.title.x = element_text(size = 0))
print(p)

dev.off()

}
