#!/usr/bin/env Rscript
source('~/Documents/R/seurat_functions.R')

#### Arguments ####
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 3){
  cancer = args[1]
  species = args[2]
  author = ''
  sample = args[3]
} else if (length(args) == 4){
  cancer = args[1]
  species = args[2]
  author = args[3]
  sample = args[4]
}
numi.min = -Inf
ngene.min = -Inf
mt.max = Inf
range = 2:25
gmin = 5

# Combining
orig.ident = paste0(cancer, species, author, sample)
print(orig.ident)

#### Load ####
st = Load10X_Spatial(getwd())
st$cancer = cancer
st$author = author
st$species = species
st$technique = 'Visium'
st$sample = as.character(sample)
st$orig.ident = orig.ident
Idents(st) = 'orig.ident'

#### Filter ####
pdf('filter.pdf')
#st = FilterCells(st, do.plot = TRUE, assay = 'Spatial', 
#                 numi.min = numi.min, ngene.min = ngene.min, mt.max = mt.max)
st = FilterGenes(st)
st = SCTransform(st, assay = 'Spatial', return.only.var.genes = FALSE)
st = FindSpatiallyVariableFeatures(st, selection.method = 'markvariogram')
st = RunPCA(st, features = SpatiallyVariableFeatures(st))
st = RunUMAP(st, dims = 1:10)
dev.off()

#### CNV ####
if (file.exists('Clones/labels_Clones.csv')){
  cnv = read.csv('Clones/labels_Clones.csv')
  rownames(cnv) = colnames(st)
  st = AddMetaData(st, cnv[,2], col.name = 'clone')
} else {
  st$clone = 0
}
st$clone = as.factor(st$clone)
col_clone = brewer_pal(palette = 'Reds')(1+nlevels(st$clone))[-1]
names(col_clone) = levels(st$clone)

#### Cluster ####
st = FindNeighbors(st)
st = FindClusters(st)
st$cluster = Idents(st)
col_cluster = hue_pal()(nlevels(st$cluster))
names(col_cluster) = levels(st$cluster)

#### NMF ####
# Run NMF
data = as.matrix(GetAssayData(st, assay = 'SCT', slot = 'scale.data'))
data = data[SpatiallyVariableFeatures(st),]
data[data < 0] = 0
data = data[apply(data, 1, var) > 0, ]
print(dim(data))
res.list = mclapply(range, function(r){
  nmf(data, rank = r, nrun = 1, seed = 'ica', method = 'nsNMF')
}, mc.cores = ncores)
names(res.list) = range
save(res.list, file = 'res.list.RData')
# Select rank
modules.list = lapply(res.list, NMFToModules, gmin = gmin)
print(sapply(modules.list,length))
comp = as.numeric(names(modules.list)) - sapply(modules.list, length)
mi = min(comp)
r = names(which(comp == mi))
r = r[length(r)]
print(r)
res = res.list[[r]]
save(res, file = 'res.RData')
# Process output
tryCatch(expr = {
  modules = NMFToModules(res, gmin = gmin)
  scores = basis(res)
  colnames(scores) = names(modules)
  coefs = coefficients(res)
  rownames(coefs) = names(modules)
  # Order modules
  h = Heatmap(coefs, clustering_distance_columns = 'euclidean')
  o = row_order(h)
  scores = scores[, o]
  coefs = coefs[o, ]
  modules = modules[o]
  print(modules)
  st = AddMetaData(st, t(coefs), col.name = rownames(coefs))
  # Cluster NMF output
  h = Heatmap(coefs, clustering_distance_columns = 'euclidean')
  hcl = as.hclust(column_dend(h))
  sig = cutree(hcl, k = length(modules))
  nmf = c(by(t(coefs), INDICES = sig, FUN = function(x){names(modules)[which.max(colMeans(x))]}))[sig]
  st$nmf = factor(nmf, levels = names(modules))
  col_nmf = c(brewer.pal(12, 'Set3'), brewer.pal(8, 'Set1'))[1:nlevels(st$nmf)]
  names(col_nmf) = levels(st$nmf)
}, error = function(e){c()})
save(st, file = 'st.RData')

#### Saving ####

#### Plotting ####
pdf('dimplots.pdf', height = 15, width = 15)
for (reduction in c('pca','umap')){
  h = DimPlot(st, pt.size = 2, reduction = reduction, group.by = 'cluster', cols = col_cluster, label = TRUE)
  print(h)
  h = DimPlot(st, pt.size = 2, reduction = reduction, group.by = 'clone', cols = col_clone, label = FALSE)
  print(h)
  h = DimPlot(st, pt.size = 2, reduction = reduction, group.by = 'nmf', cols = col_nmf)
  print(h)
  h = FeaturePlot(st, reduction = reduction, features = names(modules))
  print(h)
}
h = SpatialDimPlot(st, pt.size = 2, group.by = 'cluster', cols = col_cluster)
print(h)
h = SpatialDimPlot(st, pt.size = 2, group.by = 'cluster', cols = col_cluster)
print(h)
h = SpatialDimPlot(st, pt.size = 2, group.by = 'nmf', cols = col_nmf)
print(h)
h = SpatialFeaturePlot(st, features = names(modules))
print(h)
dev.off()
pdf('heatmaps.pdf', height = 20, width = 10)
sink('heatmaps.txt')
# Clusters
markers = FindAllMarkers2(st, do.print = TRUE, do.plot = TRUE, enrichment.type = 'GO', group.by = 'cluster', cols = col_cluster, print.bar = FALSE)
# NMF
markers = FindAllMarkers2(st, do.print = TRUE, do.plot = TRUE, enrichment.type = 'GO', group.by = 'nmf', cols = col_nmf, print.bar = FALSE)
# Modules
top_ann = HeatmapAnnotation(df = data.frame('cluster' = st$cluster,
                                            'nmf' = st$nmf), 
                            col = list('cluster' = col_cluster,
                                       'nmf' = col_nmf), which = 'column')
h = Heatmap(coefs, name = 'coefs', top_ann = top_ann,
            show_row_names = TRUE, 
            cluster_rows = FALSE,
            cluster_columns = FALSE, column_order = order(st$nmf, st$cluster),
            breaks = c(0,max(coefs)/2), colors = c('white','red'))
print(h)
h = Heatmap(scores[unlist(modules),], name = 'scores', 
            show_row_names = TRUE, row_names_gp = gpar(cex = 0.5),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            split = factor(unlist(mapply(rep, names(modules), sapply(modules, length))), levels = names(modules)))
print(h)
print(modules)
enrichment.matrix = BuildEnrichmentMatrix(rownames(st), type = 'GO')
go_10 = sapply(modules, function(tbl){
  e = Enrichment(geneset = tbl, enrichment.matrix = enrichment.matrix)
  names(sort(e))[1:10]
})
colnames(go_10) = names(modules)
print(go_10)
sink()
dev.off()



