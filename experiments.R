#!/usr/bin/env Rscript
source('~/Documents/R/seurat_functions.R')
setwd('~/Documents/Analysis/Experiments/')
pdf('Experiments.pdf', height = 20, width = 20)

nrand = 3
nbin = 25
range = 8:12
gmin = 5
pval_thresh = 0.1

dir.list = c('Adaptive/RagWT',
             'Adaptive/RagKO',
             'PancreasPeritoneum/Pancreas',
             'PancreasPeritoneum/Peritoneum',
             'LiverPeritoneum/Peritoneum',
             'LiverPeritoneum/Liver')

col_ec = c(brewer.pal(name = 'Purples', n = 3)[3], brewer.pal(name = 'Blues', n = 3)[3],brewer.pal(name = 'RdYlGn', n = 11)[c(1,4,7,10)])
names(col_ec) = c('Adaptive-RagWT','Adaptive-RagKO',
                  'PancreasPeritoneum-Pancreas','PancreasPeritoneum-Peritoneum',
                  'LiverPeritoneum-Peritoneum','LiverPeritoneum-Liver')

col_experiment = col_ec[c(1,3,5)]
names(col_experiment) = c('Adaptive','PancreasPeritoneum','LiverPeritoneum')

col_eci = c(rep(col_ec[1], 4), rep(col_ec[2], 4), rep(col_ec[3], 2), rep(col_ec[4], 2), rep(col_ec[5], 2), rep(col_ec[6], 2))
names(col_eci) = c("Adaptive-RagWT-F1", "Adaptive-RagWT-F2",
                   "Adaptive-RagWT-M1", "Adaptive-RagWT-M2",
                   "Adaptive-RagKO-F1", "Adaptive-RagKO-F2",
                   "Adaptive-RagKO-M1", "Adaptive-RagKO-M2",
                   "PancreasPeritoneum-Pancreas-1", "PancreasPeritoneum-Pancreas-2", 
                   "PancreasPeritoneum-Peritoneum-1", "PancreasPeritoneum-Peritoneum-2",
                   "LiverPeritoneum-Peritoneum-3", "LiverPeritoneum-Peritoneum-4",
                   "LiverPeritoneum-Liver-5", "LiverPeritoneum-Liver-6"
)

#### Part1 ####

srt.list = lapply(dir.list, function(dir){
  srt = loadRData(paste('~/Documents/Experiments', dir, 'srt.RData', sep = '/'))
  return(srt)
})
save(srt.list, file = 'srt.list.RData')
srt = Reduce(merge2, srt.list)
srt = srt[,!is.na(srt$condition)]
srt = srt[,!is.na(srt$experiment)]

srt$eci = factor(srt$eci, levels = c(
  "Adaptive-RagWT-F1", "Adaptive-RagWT-F2",
  "Adaptive-RagWT-M1", "Adaptive-RagWT-M2",
  "Adaptive-RagKO-F1", "Adaptive-RagKO-F2",
  "Adaptive-RagKO-M1", "Adaptive-RagKO-M2",
  "PancreasPeritoneum-Pancreas-1", "PancreasPeritoneum-Pancreas-2",
  "PancreasPeritoneum-Peritoneum-1", "PancreasPeritoneum-Peritoneum-2",
  "LiverPeritoneum-Peritoneum-3", "LiverPeritoneum-Peritoneum-4",
  "LiverPeritoneum-Liver-5", "LiverPeritoneum-Liver-6"))
srt$ec = factor(srt$ec, levels = c(
  'Adaptive-RagWT','Adaptive-RagKO',
  'PancreasPeritoneum-Pancreas','PancreasPeritoneum-Peritoneum',
  'LiverPeritoneum-Peritoneum','LiverPeritoneum-Liver'))
srt$experiment = factor(srt$experiment, levels = c(
  'Adaptive','PancreasPeritoneum','LiverPeritoneum'))

srt = RunPCA(srt)
srt = RunUMAP(srt, dims = 1:10)

srt = FindNeighbors(srt)
srt = FindClusters(srt)
srt$cluster = Idents(srt)
col_cluster = hue_pal()(nlevels(srt$cluster))
names(col_cluster) = levels(srt$cluster)
markers = FindAllMarkers2(srt, label = FALSE)
write.table(markers$genes_10, file = 'markers.csv', sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)

save(srt, file = 'srt.RData')

#### Part 2 ####

load('srt.RData')
col_cluster = hue_pal()(nlevels(srt$cluster))
names(col_cluster) = levels(srt$cluster)

# Cell types
srt$singler = 'Malignant'
srt$singler[srt$cluster %in% c(7,8)] = 'Macrophage'
srt$singler[srt$cluster %in% c(12)] = 'DC'
srt$singler[srt$cluster %in% c(15)] = 'TCell'
srt$singler[srt$cluster %in% c(17)] = 'Endothelial'
srt$singler[srt$cluster %in% c(18)] = 'Fibroblast'
srt$singler[srt$cluster %in% c(19)] = 'BCell'
srt$singler[srt$cluster %in% c(13)] = 'Neutrophil'
srt$singler = factor(srt$singler)
col_singler = c(brewer.pal(8, 'Dark2'), brewer.pal(8, 'Accent'))[1:nlevels(srt$singler)]
names(col_singler) = levels(srt$singler)
# Finding cells
cluster.malignant = unique(as.character(srt$cluster)[srt$singler == 'Malignant'])
cluster.normal = setdiff(levels(srt$cluster), cluster.malignant)
cells.malignant = WhichCells(srt, expression = cluster %in% cluster.malignant)
cells.normal = WhichCells(srt, expression = cluster %in% cluster.normal)
cells.normal = setdiff(cells.normal, cells.malignant)
# Filtering
test = c('Ccl5','Ptprc','Apoe','Tyrobp','Lyz2')
data = GetAssayData(srt, assay = 'RNA', slot = 'counts')
for (g in test){
  cells.malignant = cells.malignant[data[g, cells.malignant] < 1]
}
# Assigning
srt$type = 'undetermined'
srt@meta.data[cells.malignant, 'type'] = 'malignant'
srt@meta.data[cells.normal, 'type'] = 'normal'
srt$type = factor(srt$type)
print(table(srt$type))
print(table(srt$singler, srt$cluster, srt$type))
col_type = c(brewer_pal(palette = 'Blues')(5)[4], brewer_pal(palette = 'Greens')(5)[4], 'grey')[1:nlevels(srt$type)]
names(col_type) = levels(srt$type)
# Pop
srt$pop = as.character(srt$type)
srt$pop[srt$type == 'normal'] = as.character(srt$singler[srt$type == 'normal'])
srt$pop = gsub('malignant','Malignant',srt$pop)
srt$pop = factor(srt$pop)
srt$pop = relevel(srt$pop, ref = 'Malignant')
col_pop = col_singler[intersect(levels(srt$pop),names(col_singler))]
col_pop = c(col_pop, 'Malignant' = 'grey20')
markers = FindAllMarkers2(srt, group.by = 'pop', label = FALSE)
write.table(markers$genes_10, file = 'markers_singler.csv', sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)

# Plotting
o = sample(1:ncol(srt), size = ncol(srt), replace = FALSE)
h = DimPlot(srt, pt.size = 2, group.by = 'type', cols = col_type, order = o)
print(h)
h = DimPlot(srt, pt.size = 2, group.by = 'cluster', cols = col_cluster, label = TRUE, order = o)
print(h)
h = DimPlot(srt, pt.size = 2, group.by = 'singler', cols = col_singler, order = o)
print(h)
h = DimPlot(srt, pt.size = 2, group.by = 'pop', cols = col_pop, order = o)
print(h)
h = DimPlot(srt, pt.size = 2, group.by = 'ec', cols = col_ec, order = o)
print(h)
h = DimPlot(srt, pt.size = 2, group.by = 'experiment', cols = col_experiment, order = o)
print(h)

o = sample(1:ncol(srt), size = ncol(srt), replace = FALSE)
h = DimPlot(srt, pt.size = 1, group.by = 'type', cols = col_type, order = o, split.by = 'ec')
print(h)
h = DimPlot(srt, pt.size = 1, group.by = 'cluster', cols = col_cluster, label = TRUE, order = o, split.by = 'ec')
print(h)
h = DimPlot(srt, pt.size = 1, group.by = 'singler', cols = col_singler, order = o, split.by = 'ec')
print(h)
h = DimPlot(srt, pt.size = 1, group.by = 'pop', cols = col_pop, order = o, split.by = 'ec')
print(h)
h = DimPlot(srt, pt.size = 1, group.by = 'ec', cols = col_ec, order = o, split.by = 'ec')
print(h)
h = DimPlot(srt, pt.size = 1, group.by = 'experiment', cols = col_experiment, order = o, split.by = 'ec')
print(h)

tme = table(srt$eci, srt$pop)
tme = tme[c(
  "Adaptive-RagWT-F1", "Adaptive-RagWT-F2",
  "Adaptive-RagWT-M1", "Adaptive-RagWT-M2",
  "Adaptive-RagKO-F1", "Adaptive-RagKO-F2",
  "Adaptive-RagKO-M1", "Adaptive-RagKO-M2",
  "PancreasPeritoneum-Pancreas-1", "PancreasPeritoneum-Pancreas-2",
  "PancreasPeritoneum-Peritoneum-1", "PancreasPeritoneum-Peritoneum-2",
  "LiverPeritoneum-Peritoneum-3", "LiverPeritoneum-Peritoneum-4",
  "LiverPeritoneum-Liver-5", "LiverPeritoneum-Liver-6"),]
exclude = c('undetermined')
tme = tme[, setdiff(colnames(tme), exclude)]
tme = tme/rowSums(tme)

barplot(tme, beside = TRUE, col = col_eci[rownames(tme)])
barplot(t(tme), beside = FALSE, col = col_singler[colnames(tme)], las = 2)

# Splitting
save(srt, file = 'srt.RData')
srt.malignant = srt[, srt$type == 'malignant']
save(srt.malignant, file = 'srt.malignant.RData')
srt.normal = srt[, srt$type == 'normal']
save(srt.normal, file = 'srt.normal.RData')

#### NMF ####
# 
# srt.malignant = loadRData('srt.malignant.RData')
# srt.malignant = srt.malignant[, srt.malignant$condition == 'RagWT']
# srt.malignant = SCTransform(srt.malignant, return.only.var.genes = FALSE)
# srt.malignant = RunPCA(srt.malignant)
# srt.malignant = RunUMAP(srt.malignant, dims = 1:10)
# 
# data = as.matrix(GetAssayData(srt.malignant, slot = 'scale.data'))
# data = data[VariableFeatures(srt.malignant),]
# data[data < 0] = 0
# data = data[apply(data, 1, var) > 0, ]
# print(dim(data))
# 
# res.list = mclapply(range, function(r){
#   nmf(data, rank = r, nrun = 1, seed = 'ica', method = 'nsNMF')
# }, mc.cores = ncores)
# names(res.list) = range
# modules.list = lapply(res.list, NMFToModules, gmin = gmin)
# print(sapply(modules.list,length))
# r = names(which(sapply(modules.list, length) == as.numeric(names(modules.list))))
# r = r[length(r)]
# print(r)
# res = res.list[[r]]
# save(res.list, file = 'res.list.RData')
# save(res, file = 'res.RData')
# 
# modules = NMFToModules(res, gmin = gmin)
# print(modules)
# modules = c(modules, list('Stress' = sort(Reduce(union, modules[c('m_Egr1','m_Fos','m_Hspb1')]))))
# modules = modules[c('m_Hmgb2','Stress','m_H2_K1',
#                     'm_Mt1','m_Rbp4')]
# names(modules) = c('Cycle','Stress','Interferon',
#                    'Hypoxia','Glandular')
# ma = names(modules)
# save(modules, file = 'modules.Mm.RData')

#### Human mouse ####

srt.malignant = loadRData('srt.malignant.RData')
genes.mm = rownames(GetData(srt.malignant))

# mart = useMart("ensembl", dataset = 'mmusculus_gene_ensembl')
# attributes = c("external_gene_name", grep('hsapiens', listAttributes(mart)$name, value = TRUE))
# map = getBM(attributes = attributes,
#             filters = "external_gene_name", values = genes.mm, mart = mart, uniqueRows = TRUE)
# map = dplyr:::filter(map, external_gene_name %in% genes.mm & !duplicated(external_gene_name))
# conv = data.frame(map[, grep('gene_name', colnames(map), value = TRUE)[c(1,2)]], stringsAsFactors = FALSE)
# colnames(conv) = c('mm','hs')
# conv = conv[!(conv$hs == '' | duplicated(conv$hs)), ]
load('conv.RData')

types = c('GO','KEGG','REACTOME','HALLMARK')
enrichment.matrix.list = lapply(types, function(type){
  BuildEnrichmentMatrix(genes = conv$hs, type = type)
})
names(enrichment.matrix.list) = types

hs_hs = loadRData('../Tumors/modules.RData')
hs_hs = hs_hs[1:(length(hs_hs) - 4)]
hs_mm = lapply(hs_hs, function(mod){
  a = conv[match(mod, conv$hs),'mm']
  a = a[!is.na(a)]
})

colors.hs = loadRData('../Tumors/colors.module.RData')
colors.hs = colors.hs[names(hs_hs)]

mm_mm = loadRData('modules.Mm.RData')
ma = names(mm_mm)
mm_hs = lapply(mm_mm, function(mod){
  a = conv[match(mod, conv$mm),'hs']
  a = a[!is.na(a)]
})
save(mm_hs, file = 'modules.Mm_Hs.RData')

colors.mm = colors.hs[names(mm_mm)]
save(colors.mm, file = 'colors.module.Mm.RData')

ovlp = sapply(mm_hs, function(m1){
  sapply(hs_hs, function(m2){
    phyper(length(intersect(m1,m2)), length(m1), nrow(conv) - length(m1), length(m2), lower.tail = FALSE)
  })
})
ovlp = -log10(ovlp)
ovlp[is.infinite(ovlp)] = max(ovlp[is.finite(ovlp)])
df = data.frame('top' = apply(ovlp, 2, which.max))
rownames(df) = colnames(ovlp)

h = Heatmap(name = 'Module overlap pvalue', ovlp,
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = TRUE, show_column_names = TRUE,
            row_names_gp = gpar(col = colors.hs), column_names_gp = gpar(col = colors.mm),
            breaks = seq(0, 8, length = 9), colors = brewer_pal(palette = 'Oranges')(n = 9))
print(h)

h = Heatmap(name = 'Module overlap pvalue', ovlp[names(mm_hs),],
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = TRUE, show_column_names = TRUE,
            row_names_gp = gpar(col = colors.hs), column_names_gp = gpar(col = colors.mm),
            breaks = seq(0, 8, length = 9), colors = brewer_pal(palette = 'Oranges')(n = 9))
print(h)

# Venn diagrams
for (n in names(mm_hs)){
  if (n %in% names(hs_hs)){
    venn.diagram(x = list('mm' = mm_hs[[n]], 'hs' = hs_hs[[n]]), filename = n)
  }
}

##### Scoring cells #####

load('srt.malignant.RData')
modules = loadRData('modules.Mm.RData')
modules = c(modules, list('InterferonRed' = setdiff(modules$Interferon, c('B2m','H2-D1','H2-K1'))))
ma = names(modules)
colors.module = loadRData('colors.module.Mm.RData' )
colors.module['InterferonRed'] = colors.module['Interferon']
colors.module = colors.module[ma]

# Score cells by experiment
srt.malignant.list = SplitObject(srt.malignant, 'experiment')
srt.malignant.list = mclapply(srt.malignant.list, function(srt){
  srt = SCTransform(srt, return.only.var.genes = FALSE)
  srt = GeneToEnrichment(srt, db = modules, method = 'rand', nrand = nrand)
}, mc.cores = ncores)
for (srt in srt.malignant.list){
  srt = SCTransform(srt, return.only.var.genes = FALSE)
  srt = RunPCA(srt)
  srt = RunUMAP(srt, dims = 1:10)
  o = sample(colnames(srt), size = ncol(srt), replace = FALSE)
  h = DimPlot(srt, group.by = 'ec', cols = col_ec, order = o, pt.size = 2)
  print(h)
  h = DimPlot(srt, group.by = 'experiment', cols = col_experiment, order = o, pt.size = 2)
  print(h)
  plot.list = lapply(ma, function(n){
    FeaturePlot(srt, n, pt.size = 1, cols = c('grey',colors.module[n]))
  })
  print(CombinePlots(plot.list, ncol = 3))
  save(srt, file = paste0(unique(srt$experiment),'.malignant.RData'))
}
save(srt.malignant.list, file = 'srt.malignant.list.RData')
srt.malignant = Reduce(merge2, srt.malignant.list)

# Test modules by condition
srt.malignant.list = SplitObject(srt.malignant, 'ec')
vals = mclapply(srt.malignant.list, function(srt){
  v = NULL
  modules_rand = MakeRand(srt, db = modules, nrand = nrand, nbin = nbin)
  v = sapply(ma, function(m){
    ra = sapply(modules_rand[[m]], function(mod){
      res = c(NA,NA)
      tryCatch(expr = {
        srt_temp = RunPCA(srt, features = mod, verbose = FALSE, npcs = 2)
        res = srt_temp@reductions$pca@stdev
      }, error = function(e){})
      return(res)
    })
    re = c(NA,NA)
    tryCatch(expr = {
      srt_temp = RunPCA(srt, features = modules[[m]], verbose = FALSE, npcs = 2)
      re = srt_temp@reductions$pca@stdev
    }, error = function(e){})
    dispersion = -log10(mean(ra[1,] >= re[1]))
    coordination = -log10(mean(ra[1,]/ra[2,] >= re[1]/re[2]))
    return(c('dispersion' = dispersion, 'coordination' = coordination))
  })
}, mc.cores = ncores)
dispersion = sapply(vals, '[', 'dispersion', ma)
dispersion[is.infinite(dispersion)] = nrand
top_ann = HeatmapAnnotation(df = data.frame('ec' = names(srt.malignant.list)), 
                            col = list('ec' = col_ec), 
                            which = 'column')
h = Heatmap(dispersion, name = 'Dispersion',
            top_annotation = top_ann,
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = TRUE, show_column_names = TRUE,
            row_names_gp = gpar(col = colors.module[rownames(dispersion)]), column_names_gp = gpar(col = col_ec[colnames(dispersion)]),
            breaks = seq(0, nrand, length = 11), col = viridis_pal(begin = 0, end = 1, option = 'magma')(11))
print(h)
decorate_heatmap_body('Dispersion', {
  #l = c(0,cumsum(rle(as.character(cancer))$lengths))
  l = 0:ncol(dispersion)
  for (k in 1:length(l)){
    i = unit(l[k]/max(l), 'npc')
    grid.lines(i, c(0,1), gp = gpar(lwd = 1, col = 'black'))
  }
  l = 0:nrow(dispersion)
  for (k in 1:length(l)){
    i = unit(l[k]/max(l), 'npc')
    grid.lines(c(0,1), i, gp = gpar(lwd = 1, col = 'black'))
  }
})
cors = mclapply(srt.malignant.list, function(srt){
  data = GetData(srt, slot = 'data')
  v = sapply(ma, function(m){
    re = cor(t(as.matrix(data[modules[[m]],])))
    re[is.na(re)] = 0
    diag(re) = NA
    return(re)
  })
}, mc.cores = ncores)
cors = lapply(ma, function(m){
  a = sapply(cors, '[', m)
  names(a) = names(srt.malignant.list)
  return(a)
})
names(cors) = ma
for (m in ma){
  h.list = NULL
  for (s in names(srt.malignant.list)){
    h.list = h.list + Heatmap(cors[[m]][[s]], name = m,
                                cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE,
                              breaks = seq(-1, 1, length = 11), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 11))
  }
  print(h.list)
}
par(mfrow = c(4,4))
for (m in ma){
  a = cors[[m]]
  a = as.vector(a)
  boxplot(a, main = m, col = col_ec[names(cors[[m]])], pch = 20, cex = 0.1, las = 2, ylim = c(-0.5,1), names = rep('', length(names(cors[[m]]))))
  #abline(h = 0)
}
par(mfrow = c(1,1))
for (srt in srt.malignant.list){
  srt = SCTransform(srt, return.only.var.genes = FALSE)
  srt = RunPCA(srt)
  srt = RunUMAP(srt, dims = 1:10)
  o = sample(colnames(srt), size = ncol(srt), replace = FALSE)
  h = DimPlot(srt, group.by = 'ec', cols = col_ec, order = o, pt.size = 2)
  print(h)
  h = DimPlot(srt, group.by = 'experiment', cols = col_experiment, order = o, pt.size = 2)
  print(h)
  plot.list = lapply(ma, function(n){
    FeaturePlot(srt, n, pt.size = 1, cols = c('grey',colors.module[n]))
  })
  print(CombinePlots(plot.list, ncol = 3))
}
srt.malignant = Reduce(merge2, srt.malignant.list)

srt.malignant$eci = factor(srt.malignant$eci, levels = c(
  "Adaptive-RagWT-F1", "Adaptive-RagWT-F2",
  "Adaptive-RagWT-M1", "Adaptive-RagWT-M2",
  "Adaptive-RagKO-F1", "Adaptive-RagKO-F2",
  "Adaptive-RagKO-M1", "Adaptive-RagKO-M2",
  "PancreasPeritoneum-Pancreas-1", "PancreasPeritoneum-Pancreas-2", 
  "PancreasPeritoneum-Peritoneum-1", "PancreasPeritoneum-Peritoneum-2",
  "LiverPeritoneum-Peritoneum-3", "LiverPeritoneum-Peritoneum-4",
  "LiverPeritoneum-Liver-5", "LiverPeritoneum-Liver-6"
))
srt.malignant$ec = factor(srt.malignant$ec, levels = c(
  'Adaptive-RagWT','Adaptive-RagKO',
  'PancreasPeritoneum-Pancreas','PancreasPeritoneum-Peritoneum',
  'LiverPeritoneum-Peritoneum','LiverPeritoneum-Liver'))
srt.malignant$experiment = factor(srt.malignant$experiment, levels = c(
  'Adaptive','PancreasPeritoneum','LiverPeritoneum'))

srt.malignant = RunPCA(srt.malignant)
srt.malignant = RunUMAP(srt.malignant, dims = 1:10)

for (m in ma){
  n = paste(m, 'cat', sep = '_')
  srt.malignant@meta.data[,n] = factor((srt.malignant@meta.data[,m] > 0.5)*as.numeric(srt.malignant$ec), labels = c(0, levels(srt.malignant$ec)))
}
col_cat_ec = c('lightgrey',col_ec)
names(col_cat_ec) = c('0',names(col_ec))

h = VlnPlot(srt.malignant, features = ma, group.by = 'eci', cols = col_eci, pt.size = 0, ncol = 1)
print(h)

o = sample(colnames(srt.malignant), size = ncol(srt.malignant), replace = FALSE)
h = DimPlot(srt.malignant, group.by = 'ec', cols = col_ec, order = o, pt.size = 2)
print(h)
h = DimPlot(srt.malignant, group.by = 'experiment', cols = col_experiment, order = o, pt.size = 2)
print(h)
print(h)
plot.list = lapply(ma, function(n){
  FeaturePlot(srt.malignant, n, pt.size = 1, cols = c('grey',colors.module[n]))
})
print(CombinePlots(plot.list, ncol = 3))
plot.list = lapply(paste(ma, 'cat', sep = '_'), function(n){
  h = DimPlot(srt.malignant, group.by = n, cols = col_cat_ec, order = o, pt.size = 0.1)
})
print(CombinePlots(plot.list, ncol = 2))

o = sample(1:ncol(srt.malignant), size = ncol(srt.malignant), replace = FALSE)
h = DimPlot(srt.malignant, group.by = 'ec', cols = col_ec, order = o, pt.size = 1, split.by = 'ec')
print(h)
h = DimPlot(srt.malignant, group.by = 'experiment', cols = col_experiment, order = o, pt.size = 1, split.by = 'ec')
print(h)
print(h)
plot.list = lapply(ma, function(n){
  FeaturePlot(srt.malignant, n, split.by = 'ec', pt.size = 1, cols = c('grey',colors.module[n]))
})
print(CombinePlots(plot.list, ncol = 1))
plot.list = lapply(paste(ma, 'cat', sep = '_'), function(n){
  h = DimPlot(srt.malignant, group.by = n, cols = col_cat_ec, order = (srt.malignant@meta.data[,n] == 0), pt.size = 0.1, split.by = 'ec')
})
print(CombinePlots(plot.list, ncol = 1))

meta = srt.malignant@meta.data
save(meta, file = 'meta.malignant.RData')
save(srt.malignant, file = 'srt.malignant.RData')

#### Plotting ####

mat = Reduce(cbind, lapply(ma, function(m){
  freq = (meta[,m] > 0.5)
  mat = data.frame('value' = as.matrix(by(freq, meta$eci, mean)))
  return(mat)
}))
colnames(mat) = ma
print(mat)*100

# Barplot

plot.list = lapply(ma, function(m){
  freq = (meta[,m] > 0.5)
  mat = data.frame('value' = as.matrix(by(freq, meta$eci, mean)))
  mat$eci = rownames(mat)
  mat$experiment = factor(sapply(strsplit(mat$eci, '-'), '[', 1))
  mat$condition = factor(sapply(strsplit(mat$eci, '-'), '[', 2))
  mat$individual = factor(sapply(strsplit(mat$eci, '-'), '[', 3))
  mat$ec = factor(paste(mat$experiment, mat$condition, sep = '-'), 
                  levels = c('Adaptive-RagWT','Adaptive-RagKO',
                             'PancreasPeritoneum-Pancreas','PancreasPeritoneum-Peritoneum',
                             'LiverPeritoneum-Peritoneum','LiverPeritoneum-Liver'))
  p = ggplot(mat, aes(x = ec, y = value, fill = ec, color = ec)) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = 'bar', width = 0.5) +
    geom_point(color = 'black') +
    scale_color_manual(values = col_ec) +
    scale_fill_manual(values = col_ec) +
    geom_signif(comparisons = list(levels(mat$ec)[1:2]), 
                test = 'wilcox.test', map_signif_level = TRUE, margin_top = -0.1) +
    geom_signif(comparisons = list(levels(mat$ec)[3:4],levels(mat$ec)[5:6]), 
                test = 't.test', map_signif_level = TRUE, margin_top = -0.1) +
    ylim(0, max(mat$value)+0.1) +
    ggtitle(m) +
    theme_classic() + 
    theme(
      legend.position = 'none',
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 0), 
      axis.title.x = element_text(size = 0), 
      axis.title.y = element_text(size = 0))
  p
})
print(CombinePlots(plot.list, ncol = 3))

# Diagonal plots

mat = data.frame(Reduce(rbind, lapply(ma, function(m){
  freq = (meta[,m] > 0.5)
  sam = data.frame('value' = as.matrix(by(freq, meta$eci, mean)))
  sam$eci = rownames(sam)
  sam$experiment = factor(sapply(strsplit(sam$eci, '-'), '[', 1))
  sam$condition = factor(sapply(strsplit(sam$eci, '-'), '[', 2))
  sam$individual = factor(sapply(strsplit(sam$eci, '-'), '[', 3))
  sam$ec = factor(paste(sam$experiment, sam$condition, sep = '-'), 
                  levels = c('Adaptive-RagWT','Adaptive-RagKO',
                             'PancreasPeritoneum-Pancreas','PancreasPeritoneum-Peritoneum',
                             'LiverPeritoneum-Peritoneum','LiverPeritoneum-Liver'))
  stats_mean = as.matrix(by(sam$value, sam$ec, mean))
  rownames(stats_mean) = paste('mean', rownames(stats_mean), sep = '_')
  stats_sd = as.matrix(by(sam$value, sam$ec, sd))
  rownames(stats_sd) = paste('sd', rownames(stats_sd), sep = '_')
  stats_se = stats_sd/ sqrt(length(unique(sam$eci))/length(unique(sam$ec)))
  rownames(stats_se) = gsub('sd','se',rownames(stats_se))
  stats = t(rbind(stats_mean, stats_sd, stats_se))
  rownames(stats) = m
  return(stats)
})))
mat$module = rownames(mat)
plot.list = list()
plot.list[[1]] = ggplot(data = mat, aes(x = mean_Adaptive.RagWT, y = mean_Adaptive.RagKO, color = module)) + 
  geom_point() +
  scale_fill_manual(values = colors.module) +
  scale_color_manual(values = colors.module) +
  geom_errorbarh(aes(xmax = mean_Adaptive.RagWT + se_Adaptive.RagWT, xmin = mean_Adaptive.RagWT - se_Adaptive.RagWT), height = 0.002) +
  geom_errorbar(aes(ymax = mean_Adaptive.RagKO + se_Adaptive.RagKO, ymin = mean_Adaptive.RagKO - se_Adaptive.RagKO), width = 0.002) + 
  geom_text(aes(label = module), nudge_y = 0.01) +
  geom_abline(slope = 1, intercept = 0, color = 'grey') + 
  ggtitle('') +
  xlim(0,0.5) +
  ylim(0,0.5) +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1), 
    axis.text.y = element_text(angle = 0, hjust = 1))
plot.list[[2]] = ggplot(data = mat, aes(x = mean_PancreasPeritoneum.Pancreas, y = mean_PancreasPeritoneum.Peritoneum, color = module)) + 
  geom_point() +
  scale_fill_manual(values = colors.module) +
  scale_color_manual(values = colors.module) +
  geom_errorbarh(aes(xmax = mean_PancreasPeritoneum.Pancreas + se_PancreasPeritoneum.Pancreas, xmin = mean_PancreasPeritoneum.Pancreas - se_PancreasPeritoneum.Pancreas), height = 0.002) +
  geom_errorbar(aes(ymax = mean_PancreasPeritoneum.Peritoneum + se_PancreasPeritoneum.Peritoneum, ymin = mean_PancreasPeritoneum.Peritoneum - se_PancreasPeritoneum.Peritoneum), width = 0.002) + 
  geom_text(aes(label = module), nudge_y = 0.01) +
  geom_abline(slope = 1, intercept = 0, color = 'grey') + 
  ggtitle('') +
  xlim(0,0.5) +
  ylim(0,0.5) +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1), 
    axis.text.y = element_text(angle = 0, hjust = 1))
plot.list[[3]] = ggplot(data = mat, aes(x = mean_LiverPeritoneum.Peritoneum, y = mean_LiverPeritoneum.Liver, color = module)) + 
  geom_point() +
  scale_fill_manual(values = colors.module) +
  scale_color_manual(values = colors.module) +
  geom_errorbarh(aes(xmax = mean_LiverPeritoneum.Peritoneum + se_LiverPeritoneum.Peritoneum, xmin = mean_LiverPeritoneum.Peritoneum - se_LiverPeritoneum.Peritoneum), height = 0.002) +
  geom_errorbar(aes(ymax = mean_LiverPeritoneum.Liver + se_LiverPeritoneum.Liver, ymin = mean_LiverPeritoneum.Liver - se_LiverPeritoneum.Liver), width = 0.002) + 
  geom_text(aes(label = module), nudge_y = 0.01) +
  geom_abline(slope = 1, intercept = 0, color = 'grey') + 
  ggtitle('') +
  xlim(0,0.5) +
  ylim(0,0.5) +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1), 
    axis.text.y = element_text(angle = 0, hjust = 1))
print(CombinePlots(plot.list, ncol = 2))

# # TME
# 
# state = data.frame(Reduce(cbind, lapply(ma, function(m){
#   freq = (meta[,m] > 0.5)
#   sam = data.frame('value' = as.matrix(by(freq, meta$eci, mean)))
#   sam$eci = rownames(sam)
#   sam$experiment = factor(sapply(strsplit(sam$eci, '-'), '[', 1))
#   sam$condition = factor(sapply(strsplit(sam$eci, '-'), '[', 2))
#   sam$individual = factor(sapply(strsplit(sam$eci, '-'), '[', 3))
#   sam$ec = factor(paste(sam$experiment, sam$condition, sep = '-'), 
#                   levels = c('Adaptive-RagWT','Adaptive-RagKO',
#                              'PancreasPeritoneum-Pancreas','PancreasPeritoneum-Peritoneum',
#                              'LiverPeritoneum-Peritoneum','LiverPeritoneum-Liver'))
#   stats_mean = as.matrix(by(sam$value, sam$eci, mean))
#   colnames(stats_mean) = m
#   return(stats_mean)
# })))
# 
# # Multiple regression
# exclude = c('Malignant')
# tme = tme[, setdiff(colnames(tme), exclude)]
# 
# test = sapply(ma, function(m){
#   fit = summary(lm(state[,m] ~ tme))$coefficients
#   pval = -log10(fit[,4])
#   pval[pval < pval_thresh] = 0
#   names(pval) = gsub('tme','',names(pval))
#   est = fit[,1]
#   names(est) = gsub('tme','',names(est))
#   val = pval*sign(est)
#   return(val)
# })
# h = Heatmap(test, 
#             show_row_names = TRUE, show_column_names = TRUE,
#             cluster_rows = FALSE, cluster_columns = FALSE,
#             breaks = seq(-3, 3, length = 11), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 11))
# print(h)
# 
# # Plot hits
# par(mfrow = c(4,4))
# for (ct in colnames(tme)){
#   for (m in ma){
#     if (!test[ct,m] == 0){
#       fit = summary(lm(state[,m] ~ tme[,ct]))$coefficients
#       pval = -log10(fit[2,4])
#       plot(state[,m] ~ tme[,ct], col = col_eci[rownames(state)], pch = 20, xlab = ct, ylab = m, main = paste(round(pval, digit = 2), round(test[ct,m], digit = 2)))
#       abline(a = fit[1,1], b = fit[2,1], lty = 2)
#     }
#   }
# }
# par(mfrow = c(1,1))
# 
# # Plot hits, log
# par(mfrow = c(4,4))
# for (ct in colnames(tme)){
#   for (m in ma){
#     if (!test[ct,m] == 0){
#       fit = summary(lm(state[,m] ~ tme[,ct]))$coefficients
#       pval = -log10(fit[2,4])
#       plot(10^-4+tme[,ct], 10^-4+state[,m], col = col_eci[rownames(state)], pch = 20, xlab = ct, ylab = m, main = paste(round(pval, digit = 2), round(test[ct,m], digit = 2)), log = 'xy')
#     }
#   }
# }
# par(mfrow = c(1,1))

# Individual genes

srt = loadRData('Adaptive.malignant.RData')
srt$combined = paste(srt$condition, srt$Interferon > 0.5, sep = '_')
Idents(srt) = 'combined'
avg = AverageExpression(srt, slot = 'data')$SCT

MHCI = c('B2m','H2-D1','H2-K1')

par(mfrow = c(5,5))

plot(avg[modules$Interferon, 'RagWT_TRUE'], avg[modules$Interferon, 'RagKO_TRUE'], 
     pch = 20, xlab = 'RagWT', ylab = 'RagKO', main = 'Interferon', xlim = c(0,18), ylim = c(0,18))
text(avg[MHCI, 'RagWT_TRUE'], avg[MHCI, 'RagKO_TRUE'], MHCI, cex = 0.7, pos = 3, offset = 0.2)
abline(a = 0, b = 1, col = 'grey')

plot(avg[modules$Interferon, 'RagWT_FALSE'], avg[modules$Interferon, 'RagKO_FALSE'], 
     pch = 20, xlab = 'RagWT', ylab = 'RagKO', main = 'Other', xlim = c(0,18), ylim = c(0,18))
text(avg[MHCI, 'RagWT_FALSE'], avg[MHCI, 'RagKO_FALSE'], MHCI, cex = 0.7, pos = 3, offset = 0.2)
abline(a = 0, b = 1, col = 'grey')

plot(avg[modules$Interferon, 'RagWT_FALSE'], avg[modules$Interferon, 'RagWT_TRUE'], 
     pch = 20, xlab = 'Other', ylab = 'Interferon', main = 'RagWT', xlim = c(0,18), ylim = c(0,18))
text(avg[MHCI, 'RagWT_FALSE'], avg[MHCI, 'RagWT_TRUE'], MHCI, cex = 0.7, pos = 3, offset = 0.2)
abline(a = 0, b = 1, col = 'grey')

plot(avg[modules$Interferon, 'RagKO_FALSE'], avg[modules$Interferon, 'RagKO_TRUE'], 
     pch = 20, xlab = 'Other', ylab = 'Interferon', main = 'RagKO', xlim = c(0,18), ylim = c(0,18))
text(avg[MHCI, 'RagKO_FALSE'], avg[MHCI, 'RagKO_TRUE'], MHCI, cex = 0.7, pos = 3, offset = 0.2)
abline(a = 0, b = 1, col = 'grey')

par(mfrow = c(5,1))

o = order(avg[modules$Interferon,'RagWT_TRUE'], decreasing = TRUE)
barplot(t(avg[modules$Interferon[o],]), beside = TRUE, col = col_ec[c(1,1,2,2)], las = 2)
barplot(t(avg[modules$Interferon[o],]), beside = TRUE, col = colors.module['Interferon'], density = c(0,20,0,20), las = 2, add = TRUE)

par(mfrow = c(1,1))

dev.off()

pvals = sapply(ma, function(m){
  sapply(levels(meta$experiment), function(e){
    sub = meta[meta$experiment == e,]
    sub$ec = factor(sub$ec)
    i = levels(sub$ec)[1]
    j = levels(sub$ec)[2]
    ks.test(sub[sub$ec == i,m], sub[sub$ec == j,m])$p.value
  })
})
pvals = -log10(pvals)
pvals[is.infinite(pvals)] = max(pvals[is.finite(pvals)])
barplot(pvals, beside = TRUE, las= 2, col = col_experiment)

pvals = sapply(ma, function(m){
  sapply(levels(meta$experiment), function(e){
    sub = meta[meta$experiment == e,]
    sub$ec = factor(sub$ec)
    i = levels(sub$ec)[1]
    j = levels(sub$ec)[2]
    li = unique(sub$eci[sub$ec == i])
    lj = unique(sub$eci[sub$ec == j])
    max(sapply(li, function(i){
      sapply(lj, function(j){
          ks.test(sub[sub$eci == li,m], sub[sub$eci == lj,m])$p.value
      })
    }))
  })
})
pvals = -log10(pvals)
pvals[is.infinite(pvals)] = max(pvals[is.finite(pvals)])
barplot(pvals, beside = TRUE, las= 2, col = col_experiment)



