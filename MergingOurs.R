#!/usr/bin/env Rscript
source('~/Documents/R/seurat_functions.R')
setwd('~/Documents/Analysis/Tumors/')

load('parameters.RData')

#### Ours ####

pdf('MergingOurs.pdf', height = 10, width = 12)

modules = loadRData('modules.RData')
colors.module = loadRData('colors.module.RData')
types.module = loadRData('types.module.RData')
modules = modules[!names(modules) %in% c('AC','OPC','NPC')]
colors.module = colors.module[names(modules)]
types.module = types.module[names(modules)]

#### Malignant ####

## Srt

srt.list.ours = lapply(files.ours, function(f){
  loadRData(paste0(f,'srt.scored.RData'))
})
cancer = factor(apply(sapply(strsplit(names(files.ours), ''), '[', 1:4), 2, paste0, collapse = ''))
system = factor(corr[as.character(cancer),'system'])
germlayer = factor(corr[as.character(cancer),'germlayer'])

for (s in names(srt.list.ours)){
  srt = srt.list.ours[[s]]
  h = DimPlot(srt, group.by = 'orig.ident', cols = colors)
  print(h)
  plot.list = lapply(names(modules), function(m){
    h = FeaturePlot(srt, features = m, pt.size = 1, cols = c('grey',colors.module[m])) +
      NoAxes() + NoLegend() +
      theme(plot.title = element_text(size = 10))
  })
  print(CombinePlots(plot.list, ncol = 4))
  tryCatch(expr = {
    srt = RunTSNE(srt, dims = 1:10)
    plot.list = lapply(names(modules), function(m){
      h = FeaturePlot(srt, features = m, pt.size = 1, cols = c('grey',colors.module[m]), split.by = 'cat', reduction = 'tsne') +
        NoAxes() + NoLegend() +
        theme(plot.title = element_text(size = 10))
    })
    print(CombinePlots(plot.list, ncol = 4))
  }, error = function(e){c()})
  srt$state = apply(srt@meta.data[,names(modules)], 1, function(x){
    top = which.max(x)
    return(names(x)[top])
  })
  srt$state = factor(srt$state, levels = names(modules))
  h = DimPlot(srt, group.by = 'state', cols = colors.module)
  print(h)
  srt$num = rowSums(srt@meta.data[,names(modules)] > 0.5)
  h = FeaturePlot(srt, 'num', cols = c('grey','darkred'))
  print(h)
  srt.list.ours[[s]] = srt
}

srt.ours = Reduce(merge2, srt.list.ours)
srt.ours = RunPCA(srt.ours)
srt.ours = RunUMAP(srt.ours, dims = 1:10)

h = DimPlot(srt.ours, group.by = 'orig.ident', cols = colors, label = TRUE)
print(h)
h = DimPlot(srt.ours, group.by = 'orig.ident', cols = colors)
print(h)
plot.list = lapply(names(modules), function(m){
  h = FeaturePlot(srt.ours, features = m, pt.size = 1, cols = c('grey',colors.module[m])) +
    NoAxes() + NoLegend() +
    theme(plot.title = element_text(size = 10))
})
print(CombinePlots(plot.list, ncol = 4))
srt.ours$state = apply(srt.ours@meta.data[,names(modules)], 1, function(x){
  top = which.max(x)
  return(names(x)[top])
})
srt.ours$state = factor(srt.ours$state, levels = names(modules))
h = DimPlot(srt.ours, group.by = 'state', cols = colors.module)
print(h)

save(srt.ours, file = 'srt.ours.RData')

srt.ours = SCTransform(srt.ours, return.only.var.genes = FALSE)
srt.ours = RunPCA(srt.ours)
srt.ours = RunUMAP(srt.ours, dims = 1:10)
#srt.ours = GeneToEnrichment(srt.ours, db = modules, method = 'rand', nrand = nrand, nbin = nbin)

h = DimPlot(srt.ours, group.by = 'orig.ident', cols = colors, label = TRUE)
print(h)
h = DimPlot(srt.ours, group.by = 'orig.ident', cols = colors)
print(h)
plot.list = lapply(names(modules), function(m){
  h = FeaturePlot(srt.ours, features = m, pt.size = 0.5, cols = c('grey',colors.module[m])) +
    NoAxes() + NoLegend() +
    theme(plot.title = element_text(size = 10))
})
print(CombinePlots(plot.list, ncol = 4))
srt.ours$state = apply(srt.ours@meta.data[,names(modules)], 1, function(x){
  top = which.max(x)
  return(names(x)[top])
})
srt.ours$state = factor(srt.ours$state, levels = names(modules))
h = DimPlot(srt.ours, group.by = 'state', cols = colors.module)
print(h)

## TSNE

meta = Reduce(rbind, lapply(srt.list.ours, function(srt){
  srt@meta.data[, names(modules)]
}))
ts = tsne(meta, max_iter = max_iter, perplexity = perplexity)

ts = data.frame(ts)
colnames(ts) = c('x','y')
ts$state = unlist(lapply(srt.list.ours, function(srt){srt$state}))
ts$num = unlist(lapply(srt.list.ours, function(srt){srt$num}))
ts$orig.ident = unlist(lapply(srt.list.ours, function(srt){srt$orig.ident}))
for (m in names(modules)){
  ts[,m] = as.numeric(unlist(lapply(srt.list.ours, function(srt){srt@meta.data[,m]})))
}

distances = as.matrix(dist(t(dist(ts[,c('x','y')])), method = 'euclidean'))
ts$entropy = sapply(1:nrow(ts), function(c){
  nn = order(distances[,c])[1:20]
  ta = table(ts$orig.ident[nn])
  ta = ta/sum(ta)
  ta = ta[!ta == 0]
  return(sum(ta*log(ta)))
})
ts$mentropy = sapply(1:nrow(ts), function(c){
  nn = order(distances[,c])[1:20]
  ta = colSums(ts[nn,names(modules)] > 0.5)
  ta = ta/sum(ta)
  ta = ta[!ta == 0]
  return(sum(ta*log(ta)))
})

p = ggplot(ts, aes(x = x, y = y, fill = state, color = state)) +
  geom_point() +
  scale_fill_manual(values = colors.module) +
  scale_color_manual(values = colors.module) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
print(p)
p = ggplot(ts, aes(x = x, y = y, fill = orig.ident, color = orig.ident)) +
  geom_point() +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
print(p)
p = ggplot(ts, aes_string(x = 'x', y = 'y', color = 'entropy')) +
  geom_point() +
  scale_color_gradientn(colors = c("lightgrey", "slateblue4")) +
  ggtitle('entropy') +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
print(p)
p = ggplot(ts, aes_string(x = 'x', y = 'y', color = 'num')) +
  geom_point() +
  scale_color_gradientn(colors = c("lightgrey", "darkred")) +
  ggtitle('entropy') +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
print(p)
p = ggplot(ts, aes(x = x, y = y, fill = orig.ident, color = orig.ident)) +
  geom_point() +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  facet_wrap(~ orig.ident, ncol = 5) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
print(p)
p = ggplot(ts, aes(x = x, y = y, fill = state, color = state)) +
  geom_point() +
  scale_fill_manual(values = colors.module) +
  scale_color_manual(values = colors.module) +
  facet_wrap(~ orig.ident, ncol = 6) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
print(p)
plot.list = lapply(names(modules), function(m){
  p = ggplot(ts, aes_string(x = 'x', y = 'y', color = m)) +
    geom_point(size = 0.25) +
    scale_color_gradientn(colors = c("lightgrey", colors.module[m])) +
    ggtitle(m) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
})
print(CombinePlots(plot.list, ncol = 4))

## Stats

stats.list.ours = lapply(files.ours, function(f){
  stats = loadRData(paste0(f,'stats.RData'))
  stats = stats[,names(modules)]
  return(stats)
})

for (s in names(stats.list.ours)){
  stats = stats.list.ours[[s]]
  srt = srt.ours[, srt.ours$orig.ident == s]
  scores = srt@meta.data[,names(modules)]
  fraction = colMeans(scores > 0.5, na.rm = TRUE)
  data = GetData(srt, slot = 'scale.data')
  average = sapply(modules, function(mod){
    mean(data[intersect(mod,rownames(data)),])
  })
  num = mean(srt$num, na.rm = TRUE)
  stats = rbind(stats, fraction, average, num)
  stats.list.ours[[s]] = stats
}

stats.ours = melt(lapply(stats.list.ours, function(stats){
  melt(stats)
}))
names(stats.ours) = c('measure','module','variable','value','sample')
stats.ours$cancer = factor(apply(sapply(strsplit(stats.ours$sample, ''), '[', 1:4), 2, paste0, collapse = ''))
stats.ours$system = factor(corr[as.character(stats.ours$cancer),'system'])
stats.ours$germlayer = factor(corr[as.character(stats.ours$cancer),'germlayer'])
save(stats.ours, file = 'stats.ours.RData')

for (measure in levels(stats.ours$measure)){
  
  mat = stats.ours[stats.ours$measure == measure, c('value','module','sample')]
  mat = reshape(mat, timevar = 'sample', idvar = 'module', direction = 'wide')
  rownames(mat) = mat[,1]
  mat = mat[,-1]
  colnames(mat) = gsub('value.', '', colnames(mat), fixed = TRUE)
  
  df = data.frame('subsample' = colnames(mat), stringsAsFactors = FALSE)
  df$sample = sapply(strsplit(df$subsample, split = '.', fixed = TRUE), '[', 1)
  df$cancer = sapply(as.character(df$sample), function(x){
    y = sapply(strsplit(x, ''), '[', 1:4)
    paste0(y, collapse = '')
  })
  df$system = factor(corr[as.character(df$cancer),'system'])
  df$germlayer = factor(corr[as.character(df$cancer),'germlayer'])
  top_ann = HeatmapAnnotation(df = df[,c('sample','system','germlayer')],
                              col = list('sample' = colors, 'system' = colors.system, 'germlayer' = colors.germlayer),
                              which = 'column')
  h = Heatmap(mat, name = measure,
              top_annotation = top_ann,
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = TRUE, show_column_names = TRUE,
              row_names_gp = gpar(col = colors.module[rownames(mat)]), column_names_gp = gpar(col = colors[colnames(mat)]),
              breaks = seq(quantile(unlist(mat), seq(1:10)/10, na.rm = TRUE)['10%'], quantile(unlist(mat), seq(1:10)/10, na.rm = TRUE)['90%'], length = 11), colors = viridis_pal(begin = 0, end = 1, option = 'magma')(11))
  print(h)
  decorate_heatmap_body(measure, {
    #l = c(0,cumsum(rle(as.character(cancer))$lengths))
    l = 0:ncol(mat)
    for (k in 1:length(l)){
      i = unit(l[k]/max(l), 'npc')
      grid.lines(i, c(0,1), gp = gpar(lwd = 1, col = 'black'))
    }
    l = 0:nrow(mat)
    for (k in 1:length(l)){
      i = unit(l[k]/max(l), 'npc')
      grid.lines(c(0,1), i, gp = gpar(lwd = 1, col = 'black'))
    }
  })
  
}

## Co-occurence

par(mfrow = c(2,2))
coocc.list = list()
occ.list = list()
for (s in names(srt.list.ours)){
  
  srt = srt.list.ours[[s]]
  stats = stats.list.ours[[s]]
  meta = srt@meta.data[, names(modules)]
  
  occ = stats['frequency',]
  
  coocc = sapply(names(modules), function(m1){
    sapply(names(modules), function(m2){
      pval = phyper(
        sum((meta[, m1] > 0.5) * (meta[, m2] > 0.5)), 
        sum(meta[, m2] > 0.5), 
        dim(meta)[1] - sum(meta[, m2] > 0.5),
        sum(meta[, m1] > 0.5), 
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
  
  occ.list[[s]] = occ
  coocc.list[[s]] = coocc
  
  keep = names(which(stats['presence',] > 0.5 & stats['frequency',] > 0.05))
  occ = occ[keep]
  coocc = coocc[keep, keep]
  coocc[ coocc > -0.5 & coocc < 0.5] = 0
  coocc[coocc > 10] = 10
  coocc[coocc < -10] = -10
  
  tryCatch(expr = {
    g = graph_from_adjacency_matrix(coocc, diag = FALSE, weighted = TRUE, mode = 'undirected')
    V(g)$type = (types.module[names(V(g))] == 'Type')
    coords = layout.circle(g)
    E(g)$color = colorRamp2(breaks = seq(-10, 10, length = 9), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 9))(E(g)$weight)
    E(g)$width = 5*abs(E(g)$weight)/max(E(g)$weight)
    V(g)$color = colors.module[names(V(g))]
    V(g)$size = 50*occ/max(occ)
    # Plot graph
    plot(g, main = s,
         layout = coords,
         vertex.shape = 'circle', 
         vertex.label.size = 1, vertex.label.color = 'black',
         frame=FALSE)
  }, error = function(e){c()})
}

for (can in levels(cancer)){
  coocc = SummarizeMatrices(coocc.list[cancer == can], median)
  occ = sapply(names(modules), function(m){median(sapply(occ.list[cancer == can], '[', m))})
  keep = (names(which(occ >= 0.05)))
  occ = occ[keep]
  coocc = coocc[keep, keep]
  coocc[ coocc > -0.5 & coocc < 0.5] = 0
  coocc[coocc > 10] = 10
  coocc[coocc < -10] = -10
  tryCatch(expr = {
    g = graph_from_adjacency_matrix(coocc, diag = FALSE, weighted = TRUE, mode = 'undirected')
    V(g)$type = (types.module[names(V(g))] == 'Type')
    coords = layout.circle(g)
    E(g)$color = colorRamp2(breaks = seq(-10, 10, length = 9), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 9))(E(g)$weight)
    E(g)$width = abs(E(g)$weight)/5
    V(g)$color = colors.module[names(V(g))]
    V(g)$size = occ*200
    # Plot graph
    plot(g, main = can,
         layout = coords,
         vertex.shape = 'circle', 
         vertex.label.size = 1, vertex.label.color = 'black',
         frame=FALSE)
  }, error = function(e){c()})
}

for (sys in levels(system)){
  coocc = SummarizeMatrices(coocc.list[system == sys], median)
  occ = sapply(names(modules), function(m){median(sapply(occ.list[system == sys], '[', m))})
  keep = (names(which(occ >= 0.05)))
  occ = occ[keep]
  coocc = coocc[keep, keep]
  coocc[ coocc > -0.5 & coocc < 0.5] = 0
  coocc[coocc > 10] = 10
  coocc[coocc < -10] = -10
  tryCatch(expr = {
    g = graph_from_adjacency_matrix(coocc, diag = FALSE, weighted = TRUE, mode = 'undirected')
    V(g)$type = (types.module[names(V(g))] == 'Type')
    coords = layout.circle(g)
    E(g)$color = colorRamp2(breaks = seq(-10, 10, length = 9), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 9))(E(g)$weight)
    E(g)$width = abs(E(g)$weight)/5
    V(g)$color = colors.module[names(V(g))]
    V(g)$size = occ*200
    # Plot graph
    plot(g, main = sys,
         layout = coords,
         vertex.shape = 'circle', 
         vertex.label.size = 1, vertex.label.color = 'black',
         frame=FALSE)
  }, error = function(e){c()})
}

coocc = SummarizeMatrices(coocc.list, median)
occ = sapply(names(modules), function(m){median(sapply(occ.list, '[', m))})
keep = (names(which(occ >= 0.05)))
occ = occ[keep]
coocc = coocc[keep, keep]
coocc[coocc > -0.5 & coocc < 0.5] = 0
coocc[coocc > 10] = 10
coocc[coocc < -10] = -10
tryCatch(expr = {
  g = graph_from_adjacency_matrix(coocc, diag = FALSE, weighted = TRUE, mode = 'undirected')
  V(g)$type = (types.module[names(V(g))] == 'Type')
  coords = layout.circle(g)
  E(g)$color = colorRamp2(breaks = seq(-10, 10, length = 9), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 9))(E(g)$weight)
  E(g)$width = abs(E(g)$weight)/5
  V(g)$color = colors.module[names(V(g))]
  V(g)$size = occ*200
  # Plot graph
  plot(g, main = 'All',
       layout = coords,
       vertex.shape = 'circle', 
       vertex.label.size = 1, vertex.label.color = 'black',
       frame=FALSE)
}, error = function(e){c()})

par(mfrow = c(1,1))


#### All ####

srt.list.ours.all = lapply(files.ours, function(f){
  loadRData(paste0(f,'srt.all.RData'))
})

srt.ours.all = Reduce(merge2, srt.list.ours.all)
srt.ours.all = RunPCA(srt.ours.all)
srt.ours.all = RunUMAP(srt.ours.all, dims = 1:10)
save(srt.ours.all, file = 'srt.ours.all.RData')

colors.pop = hue_pal()(nlevels(srt.ours.all$pop))
names(colors.pop) = levels(srt.ours.all$pop)
colors.pop['Malignant'] = 'black'
save(colors.pop, file = 'colors.pop.ours.RData')

o = sample(1:ncol(srt.ours.all), size = ncol(srt.ours.all), replace = FALSE)
h = DimPlot(srt.ours.all, group.by = 'orig.ident', cols = colors, order = o, pt.size = 1, label = TRUE)
print(h)
h = DimPlot(srt.ours.all, group.by = 'orig.ident', cols = colors, order = o, pt.size = 1)
print(h)
h = DimPlot(srt.ours.all, group.by = 'type', cols = colors.type, order = o, pt.size = 1, label = TRUE)
print(h)
h = DimPlot(srt.ours.all, group.by = 'type', cols = colors.type, order = o, pt.size = 1)
print(h)
h = DimPlot(srt.ours.all, group.by = 'pop', cols = colors.pop, order = o, pt.size = 1, label = TRUE)
print(h)
h = DimPlot(srt.ours.all, group.by = 'pop', cols = colors.pop, order = o, pt.size = 1, label = FALSE)
print(h)

srt.ours.tme = srt.ours.all[, !srt.ours.all$pop == 'Malignant']
srt.ours.tme = RunPCA(srt.ours.tme)
srt.ours.tme = RunUMAP(srt.ours.tme, dims = 1:10)
save(srt.ours.tme, file = 'srt.ours.tme.RData')

# Frequency regression

freq = t(sapply(levels(srt.ours$orig.ident), function(or){
  colMeans(srt.ours@meta.data[srt.ours$orig.ident == or, names(modules)] > 0.5, na.rm = TRUE)
}))
tme = table(srt.ours.tme$orig.ident, srt.ours.tme$pop)
exclude = c('Malignant','Embryonic_stem_cells','iPS_cells','Tissue_stem_cells','Smooth_muscle_cells','Chondrocytes','Epithelial_cells','Hepatocytes','Keratinocytes','Astrocyte','Neurons','Neuroepithelial_cell',
            'CMP','Erythroblast','GMP','Gametocytes','GMP','HSC_CD34+','Pre-B_cell_CD34-','Pro-B_cell_CD34+','Pro-Myelocyte')
tme = tme[, setdiff(colnames(tme), exclude)]
tme = tme/rowSums(tme)

# Multiple regression
test = sapply(names(modules), function(m){
  fit = summary(lm(freq[,m] ~ tme))$coefficients
  pval = -log10(fit[,4])
  pval[pval < -log10(pval_thresh)] = 0
  names(pval) = gsub('tme','',names(pval))
  est = fit[,1]
  names(est) = gsub('tme','',names(est))
  val = pval*sign(est)
  return(val)
})
h = Heatmap(test, 
            show_row_names = TRUE, show_column_names = TRUE,
            cluster_rows = FALSE, cluster_columns = FALSE,
            breaks = seq(-3, 3, length = 11), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 11))
print(h)

# Plot hits
par(mfrow = c(4,4))
for (ct in colnames(tme)){
  for (m in names(modules)){
    fit = summary(lm(freq[,m] ~ tme[,ct]))$coefficients
    pval = -log10(fit[2,4])
    plot(freq[,m] ~ tme[,ct], col = colors[rownames(tme)], pch = 20, xlab = ct, ylab = m, main = round(pval, digits = 2))
    if (pval > 0.5){
      abline(a = fit[1,1], b = fit[2,1], lty = 2)
    }
  }
}
par(mfrow = c(1,1))

dev.off()

