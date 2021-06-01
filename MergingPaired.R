#!/usr/bin/env Rscript
source('~/Documents/R/seurat_functions.R')
setwd('~/Documents/Analysis/Tumors/')

load('parameters.RData')

#### Paired ####

pdf('MergingPaired.pdf', height = 10, width = 20)

modules = loadRData('modules.RData')
colors.module = loadRData('colors.module.RData')
types.module = loadRData('types.module.RData')
modules = modules[!names(modules) %in% c('AC','OPC','NPC')]
colors.module = colors.module[names(modules)]
types.module = types.module[names(modules)]

#### Malignant ####

## Srt

srt.list.paired = lapply(files.paired.list, function(f){
  loadRData(paste0(f[['normal']],'Paired/srt.scored.RData'))
})
cancer = factor(apply(sapply(strsplit(names(files.paired.list), ''), '[', 1:4), 2, paste0, collapse = ''))
system = factor(corr[as.character(cancer),'system'])
germlayer = factor(corr[as.character(cancer),'germlayer'])

for (s in names(srt.list.paired)){
  srt = srt.list.paired[[s]]
  h = DimPlot(srt, group.by = 'orig.ident', cols = colors)
  print(h)
  h = DimPlot(srt, group.by = 'cat', cols = colors.cat)
  print(h)
  plot.list = lapply(names(modules), function(m){
    h = FeaturePlot(srt, features = m, pt.size = 1, cols = c('grey',colors.module[m]), split.by = 'cat') +
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
  srt$num = rowSums(srt@meta.data[,names(modules)] > 0.5)
  h = FeaturePlot(srt, 'num', cols = c('grey','darkred'))
  print(h)
  h = VlnPlot(srt, 'num', group.by = 'cat', cols = colors.cat, pt.size = 0)
  print(h)
  h = DimPlot(srt, group.by = 'state', cols = colors.module)
  print(h)
  srt.list.paired[[s]] = srt
}

srt.paired = Reduce(merge2, srt.list.paired)
srt.paired = RunPCA(srt.paired)
srt.paired = RunUMAP(srt.paired, dims = 1:10)

h = DimPlot(srt.paired, group.by = 'orig.ident', cols = colors, label = TRUE)
print(h)
h = DimPlot(srt.paired, group.by = 'orig.ident', cols = colors)
print(h)
plot.list = lapply(names(modules), function(m){
  h = FeaturePlot(srt.paired, features = m, pt.size = 0.5, cols = c('grey',colors.module[m])) +
    NoAxes() + NoLegend() +
    theme(plot.title = element_text(size = 10))
})
print(CombinePlots(plot.list, ncol = 4))
srt.paired$state = apply(srt.paired@meta.data[,names(modules)], 1, function(x){
  top = which.max(x)
  return(names(x)[top])
})
srt.paired$state = factor(srt.paired$state, levels = names(modules))
h = DimPlot(srt.paired, group.by = 'state', cols = colors.module)
print(h)
h = FeaturePlot(srt.paired, 'num', cols = c('grey','darkred'))
print(h)
h = DimPlot(srt.paired, group.by = 'cat', cols = colors.cat)
print(h)

save(srt.paired, file = 'srt.paired.RData')

srt.paired = SCTransform(srt.paired, return.only.var.genes = FALSE)
srt.paired = RunPCA(srt.paired)
srt.paired = RunUMAP(srt.paired, dims = 1:10)
srt.paired = GeneToEnrichment(srt.paired, db = modules, method = 'rand', nrand = nrand, nbin = nbin)

h = DimPlot(srt.paired, group.by = 'orig.ident', cols = colors, label = TRUE)
print(h)
h = DimPlot(srt.paired, group.by = 'orig.ident', cols = colors)
print(h)
h = DimPlot(srt.paired, group.by = 'cat', cols = colors.cat)
print(h)
plot.list = lapply(names(modules), function(m){
  h = FeaturePlot(srt.paired, features = m, pt.size = 1, cols = c('grey',colors.module[m])) +
    NoAxes() + NoLegend() +
    theme(plot.title = element_text(size = 10))
})
print(CombinePlots(plot.list, ncol = 4))
srt.paired$state = apply(srt.paired@meta.data[,names(modules)], 1, function(x){
  top = which.max(x)
  return(names(x)[top])
})
srt.paired$state = factor(srt.paired$state, levels = names(modules))
h = DimPlot(srt.paired, group.by = 'state', cols = colors.module)
print(h)

# ## TSNE
# 
# meta = Reduce(rbind, lapply(srt.list.paired, function(srt){
#   srt@meta.data[, names(modules)]
# }))
# ts = tsne(meta, max_iter = max_iter, perplexity = perplexity)
# 
# ts = data.frame(ts)
# colnames(ts) = c('x','y')
# ts$state = unlist(lapply(srt.list.paired, function(srt){srt$state}))
# ts$num = unlist(lapply(srt.list.paired, function(srt){srt$num}))
# ts$orig.ident = unlist(lapply(srt.list.paired, function(srt){srt$orig.ident}))
# for (m in names(modules)){
#   ts[,m] = as.numeric(unlist(lapply(srt.list.paired, function(srt){srt@meta.data[,m]})))
# }
# 
# distances = as.matrix(dist(t(dist(ts[,c('x','y')])), method = 'euclidean'))
# ts$entropy = sapply(1:nrow(ts), function(c){
#   nn = order(distances[,c])[1:20]
#   ta = table(ts$orig.ident[nn])
#   ta = ta/sum(ta)
#   ta = ta[!ta == 0]
#   return(sum(ta*log(ta)))
# })
# ts$mentropy = sapply(1:nrow(ts), function(c){
#   nn = order(distances[,c])[1:20]
#   ta = colSums(ts[nn,names(modules)] > 0.5)
#   ta = ta/sum(ta)
#   ta = ta[!ta == 0]
#   return(sum(ta*log(ta)))
# })
# 
# p = ggplot(ts, aes(x = x, y = y, fill = state, color = state)) +
#   geom_point() +
#   scale_fill_manual(values = colors.module) +
#   scale_color_manual(values = colors.module) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(fill = NA),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank())
# print(p)
# p = ggplot(ts, aes(x = x, y = y, fill = orig.ident, color = orig.ident)) +
#   geom_point() +
#   scale_fill_manual(values = colors) +
#   scale_color_manual(values = colors) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(fill = NA),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank())
# print(p)
# p = ggplot(ts, aes_string(x = 'x', y = 'y', color = 'entropy')) +
#   geom_point() +
#   scale_color_gradientn(colors = c("lightgrey", "slateblue4")) +
#   ggtitle('entropy') +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(fill = NA),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank())
# print(p)
# p = ggplot(ts, aes_string(x = 'x', y = 'y', color = 'num')) +
#   geom_point() +
#   scale_color_gradientn(colors = c("lightgrey", "darkred")) +
#   ggtitle('entropy') +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(fill = NA),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank())
# print(p)
# p = ggplot(ts, aes(x = x, y = y, fill = orig.ident, color = orig.ident)) +
#   geom_point() +
#   scale_fill_manual(values = colors) +
#   scale_color_manual(values = colors) +
#   facet_wrap(~ orig.ident, ncol = 6) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(fill = NA),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank())
# print(p)
# p = ggplot(ts, aes(x = x, y = y, fill = state, color = state)) +
#   geom_point() +
#   scale_fill_manual(values = colors.module) +
#   scale_color_manual(values = colors.module) +
#   facet_wrap(~ orig.ident, ncol = 6) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(fill = NA),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank())
# print(p)
# plot.list = lapply(names(modules), function(m){
#   p = ggplot(ts, aes_string(x = 'x', y = 'y', color = m)) +
#     geom_point(size = 0.25) +
#     scale_color_gradientn(colors = c("lightgrey", colors.module[m])) +
#     ggtitle(m) +
#     theme_classic() +
#     theme(
#       panel.border = element_rect(fill = NA),
#       legend.position = "none",
#       plot.title = element_text(hjust = 0.5),
#       axis.line.x = element_blank(),
#       axis.line.y = element_blank(),
#       axis.ticks.x = element_blank(),
#       axis.ticks.y = element_blank(),
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title.x = element_blank(),
#       axis.title.y = element_blank())
# })
# print(CombinePlots(plot.list, ncol = 4))

## Stats

stats.list.paired = lapply(files.paired.list, function(f){
  stats = loadRData(paste0(f[['normal']],'Paired/stats.RData'))
  stats = stats[,names(modules)]
  return(stats)
})

for (s in names(stats.list.paired)){
  stats = stats.list.paired[[s]]
  srt = srt.paired[, srt.paired$orig.ident == s]
  scores = srt@meta.data[,names(modules)]
  fraction.normal = colMeans(scores[srt$cat == 'normal',] > 0.5, na.rm = TRUE)
  fraction.tumor = colMeans(scores[srt$cat == 'tumor',] > 0.5, na.rm = TRUE)
  data = GetData(srt, slot = 'scale.data')
  average.normal = sapply(modules, function(mod){
    mean(data[intersect(mod,rownames(data)),srt$cat == 'normal'])
  })
  average.tumor = sapply(modules, function(mod){
    mean(data[intersect(mod,rownames(data)),srt$cat == 'tumor'])
  })
  num.normal = mean(srt@meta.data[srt$cat == 'normal','num'], na.rm = TRUE)
  num.tumor = mean(srt@meta.data[srt$cat == 'tumor','num'], na.rm = TRUE)
  stats = rbind(stats, fraction.normal, fraction.tumor, average.normal, average.tumor, num.normal, num.tumor)
  stats.list.paired[[s]] = stats
}

print(stats.list.paired$PDACHs1[c('num.normal','num.tumor'),1])
print(stats.list.paired$PDACHs2[c('num.normal','num.tumor'),1])
print(stats.list.paired$PDACHs3[c('num.normal','num.tumor'),1])

stats.paired = melt(lapply(stats.list.paired, function(stats){
  melt(stats)
}))
names(stats.paired) = c('measure','module','variable','value','sample')
stats.paired$cat = factor(sapply(strsplit(as.character(stats.paired$measure), '.', fixed = TRUE), '[', 2), levels = names(colors.cat))
stats.paired$measure = factor(sapply(strsplit(as.character(stats.paired$measure), '.', fixed = TRUE), '[', 1))
stats.paired$cancer = factor(apply(sapply(strsplit(stats.paired$sample, ''), '[', 1:4), 2, paste0, collapse = ''))
stats.paired$system = factor(corr[as.character(stats.paired$cancer),'system'])
stats.paired$germlayer = factor(corr[as.character(stats.paired$cancer),'germlayer'])
stats.paired$sample.cat = paste(stats.paired$sample, stats.paired$cat, sep = '.')
save(stats.paired, file = 'stats.paired.RData')

for (measure in levels(stats.paired$measure)){
  
  mat = stats.paired[stats.paired$measure == measure, c('value','module','sample.cat')]
  mat = reshape(mat, timevar = 'sample.cat', idvar = 'module', direction = 'wide')
  rownames(mat) = mat[,1]
  mat = mat[,-1]
  colnames(mat) = gsub('value.', '', colnames(mat), fixed = TRUE)
  
  df = data.frame('subsample' = colnames(mat), stringsAsFactors = FALSE)
  df$sample = sapply(strsplit(df$subsample, split = '.', fixed = TRUE), '[', 1)
  df$cat = sapply(strsplit(df$subsample, split = '.', fixed = TRUE), '[', 2)
  df$cancer = sapply(as.character(df$sample), function(x){
    y = sapply(strsplit(x, ''), '[', 1:4)
    paste0(y, collapse = '')
  })
  df$system = factor(corr[as.character(df$cancer),'system'])
  df$germlayer = factor(corr[as.character(df$cancer),'germlayer'])
  top_ann = HeatmapAnnotation(df = df[,c('sample','cat')],
                              col = list('sample' = colors, 'cat' = colors.cat),
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
  
  p.list = lapply(c('Squamous','Glandular','pEMT','Interferon'), function(m){
    sub = stats.paired[stats.paired$measure == measure & stats.paired$module == m,]
    sub = sub[order(sub$cat), ]
    p = ggpaired(sub, x = 'cat', y = 'value',
               color = 'cat', line.color = 'grey', line.size = 0.4, facet.by = 'cancer') +
      stat_compare_means(paired = TRUE) +
      scale_color_manual(values = colors.cat) +
      ggtitle(m) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0),
        legend.position = 'none')
  })
  print(CombinePlots(p.list), ncol = 4)
  
  p.list = lapply(c('Cycle','Stress','Hypoxia',
                    'Oxphos','Metal','Mesenchymal',
                    'Alveolar','Basal','Ciliated'), function(m){
    sub = stats.paired[stats.paired$measure == measure & stats.paired$module == m,]
    sub = sub[order(sub$cat), ]
    p = ggpaired(sub, x = 'cat', y = 'value',
                 color = 'cat', line.color = 'grey', line.size = 0.4, facet.by = 'cancer') +
      stat_compare_means(paired = TRUE) +
      scale_color_manual(values = colors.cat) +
      ggtitle(m) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0),
        legend.position = 'none')
  })
  print(CombinePlots(p.list), ncol = 3)
  
}


# #### All ####
# 
# srt.list.paired.all = lapply(files.paired, function(f){
#   loadRData(paste0(f,'srt.all.RData'))
# })
# 
# srt.paired.all = Reduce(merge2, srt.list.paired.all)
# srt.paired.all = RunPCA(srt.paired.all)
# srt.paired.all = RunUMAP(srt.paired.all, dims = 1:10)
# save(srt.paired.all, file = 'srt.paired.all.RData')
# 
# colors.pop = hue_pal()(nlevels(srt.paired.all$pop))
# names(colors.pop) = levels(srt.paired.all$pop)
# colors.pop['Malignant'] = 'black'
# save(colors.pop, file = 'colors.pop.paired.RData')
# 
# o = sample(1:ncol(srt.paired.all), size = ncol(srt.paired.all), replace = FALSE)
# h = DimPlot(srt.paired.all, group.by = 'orig.ident', cols = colors, order = o, pt.size = 1, label = TRUE)
# print(h)
# h = DimPlot(srt.paired.all, group.by = 'orig.ident', cols = colors, order = o, pt.size = 1)
# print(h)
# h = DimPlot(srt.paired.all, group.by = 'type', cols = colors.type, order = o, pt.size = 1, label = TRUE)
# print(h)
# h = DimPlot(srt.paired.all, group.by = 'type', cols = colors.type, order = o, pt.size = 1)
# print(h)
# h = DimPlot(srt.paired.all, group.by = 'pop', cols = colors.pop, order = o, pt.size = 1, label = TRUE)
# print(h)
# h = DimPlot(srt.paired.all, group.by = 'pop', cols = colors.pop, order = o, pt.size = 1, label = FALSE)
# print(h)
# 
# srt.paired.tme = srt.paired.all[, !srt.paired.all$pop == 'Malignant']
# srt.paired.tme = RunPCA(srt.paired.tme)
# srt.paired.tme = RunUMAP(srt.paired.tme, dims = 1:10)
# save(srt.paired.tme, file = 'srt.paired.tme.RData')


dev.off()