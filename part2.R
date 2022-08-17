#!/usr/bin/env Rscript
source('~/Documents/R/seurat_functions.R')
load('srt.RData')

##### Arguments #####
args <- commandArgs(trailingOnly = TRUE)
cluster.malignant = eval(parse(text = args[1]))
clone.malignant = eval(parse(text = args[2]))
k_obs_groups = 2
#if (unique(srt$author) == ''){
  range = 2:25
#} else {
#  range = 5:35
#}
gmin = 5
print(unique(srt$orig.ident))

##### Determining type ####
# Determining criteria
singler.malignant = c('Embryonic_stem_cells','iPS_cells','Tissue_stem_cells','Chondrocytes','Fibroblasts','Smooth_muscle_cells','Epithelial_cells','Hepatocytes','Keratinocytes','Astrocyte','Neurons','Neuroepithelial_cell')
cluster.normal = setdiff(levels(srt$cluster), cluster.malignant)
clone.normal = setdiff(levels(srt$clone), clone.malignant)
singler.normal = setdiff(levels(srt$singler), singler.malignant)
# Finding cells
cells.malignant = WhichCells(srt, expression = singler %in% singler.malignant
                             & cluster %in% cluster.malignant
                             & clone %in% clone.malignant)
if (length(cells.malignant) == ncol(srt)){
  donormal = FALSE
} else {
  donormal = TRUE
}
if (donormal){
  cells.normal = WhichCells(srt, expression = cluster %in% cluster.normal)
  cells.normal = setdiff(cells.normal, cells.malignant)
} else {
  cells.normal = c()
}
# Filtering
test = intersect(c('B_cell','DC','Macrophage','Monocyte','Neutrophils','NK_cell','T_cells'), levels(srt$singler))
print(length(cells.malignant))
for (i in test){
  s = paste('scores',i ,sep = '.')
  thresh = min(srt@meta.data[srt$singler == i,s])
  #thresh = quantile(srt@meta.data[srt$singler == i,s], seq(0,1,by=0.01))['5%']
  cells.malignant = cells.malignant[srt@meta.data[cells.malignant, s] <= thresh]
  print(i)
  print(length(cells.malignant))
}
Idents(srt) = 'singler'
avg = AverageExpression(srt, assay = 'SCT', slot = 'data')$SCT
data = GetData(srt, assay = 'SCT', slot = 'data')[rownames(avg),cells.malignant]
for (i in test){
  corr = apply(data, 2, function(x){
    cor(x, avg[,i], method = 'pearson')
  })
  cells.malignant = setdiff(cells.malignant, colnames(data)[corr > mean(corr) + 2*sd(corr)])
  print(i)
  print(length(cells.malignant))
}
if (unique(srt$author) == 'Puram'){
  cells.malignant = cells.malignant[srt@meta.data[cells.malignant, 'malignant'] == 1]
}
# Assigning
srt$type = 'undetermined'
srt@meta.data[cells.malignant, 'type'] = 'malignant'
srt@meta.data[cells.normal, 'type'] = 'normal'
srt$type = factor(srt$type, levels = c('malignant','normal','undetermined'))
print(table(srt$type))
print(table(srt$singler, srt$cluster, srt$type))
# Saving
save(srt, file = 'srt.RData')

##### Removing undetermined #####
srt = srt[, srt$type %in% c('malignant','normal')]
srt$type = factor(srt$type, levels = c('malignant','normal'))
col_type = c(brewer_pal(palette = 'Blues')(5)[4], brewer_pal(palette = 'Greens')(5)[4], 'grey')[1:nlevels(srt$type)]
names(col_type) = levels(srt$type)
srt = SCTransform(srt, return.only.var.genes = FALSE)
srt = RunPCA(srt, features = VariableFeatures(srt))
srt = RunUMAP(srt, dims = 1:10)
# Cluster
col_cluster = hue_pal()(nlevels(srt$cluster))
names(col_cluster) = levels(srt$cluster)
# Singler
col_singler = c(brewer.pal(8, 'Dark2'), brewer.pal(8, 'Accent'))[1:nlevels(srt$singler)]
names(col_singler) = levels(srt$singler)
# Pop
srt$pop = as.character(srt$type)
srt$pop[srt$type == 'normal'] = as.character(srt$singler[srt$type == 'normal'])
srt$pop = gsub('malignant','Malignant',srt$pop)
srt$pop = factor(srt$pop)
srt$pop = relevel(srt$pop, ref = 'Malignant')
col_pop = col_singler[intersect(levels(srt$pop),names(col_singler))]
col_pop = c(col_pop, 'Malignant' = 'grey20')

#### InferCNV ####
# Finding subclones
dir.create('Subclones', showWarnings = FALSE, recursive = TRUE, mode = "0777")
setwd('Subclones')
srt$subclone = 0
srt$subclone = factor(srt$subclone)
col_subclone = brewer_pal(palette = 'Purples')(1+nlevels(srt$subclone))[-1]
names(col_subclone) = levels(srt$subclone)
tryCatch(expr = {
annotations_file = data.frame(srt@meta.data$type, row.names = rownames(srt@meta.data), stringsAsFactors = FALSE)
raw_counts_matrix = as.matrix(GetAssayData(srt, assay = 'SCT', slot = 'counts'))
if (unique(srt$species) == 'Hs'){
  gene_order_file = '~/Documents/CNV/Hs/gene_order.tsv'
}
if (unique(srt$species) == 'Mm'){
  gene_order_file = '~/Documents/CNV/Mm/gene_order.tsv'
}
if (donormal){
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = raw_counts_matrix,
                                      annotations_file = annotations_file,
                                      gene_order_file = gene_order_file,
                                      delim = "\t",
                                      ref_group_names = 'normal',
                                      chr_exclude = NULL)
} else {
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = raw_counts_matrix,
                                      annotations_file = annotations_file,
                                      gene_order_file = gene_order_file,
                                      delim = "\t",
                                      ref_group_names = NULL,
                                      chr_exclude = NULL)
}

infercnv_obj = infercnv::run(infercnv_obj,
                             resume_mode = FALSE,
                             out_dir = getwd(),
                             cutoff = 0.1,
                             cluster_by_groups = FALSE,
                             analysis_mode = 'subclusters',
                             tumor_subcluster_partition_method = 'qnorm',
                             denoise = FALSE,
                             HMM = FALSE,
                             HMM_type = 'i6',
                             num_threads = ncores,
                             no_plot = TRUE,
                             plot_steps = FALSE)
save(infercnv_obj, file = 'infer_cnv.RData')
obs_hcl = hclust(dist(t(infercnv_obj@expr.data[, unlist(infercnv_obj@observation_grouped_cell_indices)]), 'euclidean'), 'ward.D2')
subclone = cutree(tree = as.hclust(obs_hcl), k = k_obs_groups)
srt$subclone = 0
srt$subclone[names(subclone)] = subclone[names(subclone)]
srt$subclone = factor(srt$subclone)
col_subclone = brewer_pal(palette = 'Purples')(1+nlevels(srt$subclone))[-1]
names(col_subclone) = levels(srt$subclone)
# Plotting heatmap with subclones
annotations_file = data.frame(srt$subclone, row.names = rownames(srt@meta.data))
raw_counts_matrix = as.matrix(GetAssayData(srt, assay = 'SCT', slot = 'counts'))
if (unique(srt$species) == 'Hs'){
  gene_order_file = '~/Documents/CNV/Hs/gene_order.tsv'
}
if (unique(srt$species) == 'Mm'){
  gene_order_file = '~/Documents/CNV/Mm/gene_order.tsv'
}
if (donormal){
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = raw_counts_matrix,
                                      annotations_file = annotations_file,
                                      gene_order_file = gene_order_file,
                                      delim = "\t",
                                      ref_group_names = '0',
                                      chr_exclude = NULL)
} else {
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = raw_counts_matrix,
                                      annotations_file = annotations_file,
                                      gene_order_file = gene_order_file,
                                      delim = "\t",
                                      ref_group_names = NULL,
                                      chr_exclude = NULL)
}
infercnv_obj = infercnv::run(infercnv_obj,
                             resume_mode = FALSE,
                             out_dir = getwd(),
                             cutoff = 0.1,
                             cluster_by_groups = TRUE,
                             k_obs_groups = k_obs_groups,
                             analysis_mode = 'subclusters',
                             tumor_subcluster_partition_method = 'qnorm',
                             denoise = TRUE,
                             sd_amplifier=2,  # sets midpoint for logistic
                             noise_logistic=TRUE, # turns gradient filtering on
                             HMM = FALSE,
                             HMM_type = 'i6',
                             num_threads = ncores,
                             no_plot = FALSE,
                             plot_steps = FALSE)
infercnv_obj_medianfiltered = infercnv::apply_median_filtering(infercnv_obj)
infercnv::plot_cnv(infercnv_obj_medianfiltered,
                   output_filename='infercnv.median_filtered',
                   x.range=c(0.9,1.1),
                   x.center=1,
                   title = "infercnv",
                   color_safe_pal = FALSE)
}, error = function(e){c()})
save(infercnv_obj, file = 'infer_cnv.RData')
setwd('../')

#### Plotting ####
pdf('dimplots.all.pdf', height = 15, width = 15)
sumsq = apply(infercnv_obj@expr.data, 2, function(prof){
  sum(prof^2)
})
srt$sumsq = sumsq
top = apply(infercnv_obj@expr.data[,names(sort(sumsq, decreasing = TRUE))[1:10]], 1, mean)
topcor = apply(infercnv_obj@expr.data, 2, function(prof){
  cor(prof, top)
})
srt$topcor = topcor
FeatureScatter(srt, pt.size = 2, 'sumsq', 'topcor', group.by = 'type', col = col_type)
FeatureScatter(srt, pt.size = 2, 'sumsq', 'topcor', group.by = 'subclone', col = col_subclone)
FeatureScatter(srt, pt.size = 2, 'sumsq', 'topcor', group.by = 'cluster', col = col_cluster)
for (reduction in c('pca','umap')){
  h = DimPlot(srt, pt.size = 2, reduction = reduction, group.by = 'type', cols = col_type)
  print(h)
  h = DimPlot(srt, pt.size = 2, reduction = reduction, group.by = 'cluster', cols = col_cluster, label = TRUE)
  print(h)
  h = DimPlot(srt, pt.size = 2, reduction = reduction, group.by = 'singler', cols = col_singler)
  print(h)
  h = DimPlot(srt, pt.size = 2, reduction = reduction, group.by = 'pop', cols = col_pop)
  print(h)
  h = DimPlot(srt, pt.size = 2, reduction = reduction, group.by = 'subclone', cols = col_subclone)
  print(h)
  h = FeaturePlot(srt, pt.size = 2, reduction = reduction, features = 'sumsq')
  print(h)
  h = FeaturePlot(srt, pt.size = 2, reduction = reduction, features = 'topcor')
  print(h)
}
dev.off()

#### Spliting ####
srt.all = srt
save(srt.all, file = 'srt.all.RData')
srt.malignant = srt[, srt$type == 'malignant']
save(srt.malignant, file = 'srt.malignant.RData')
if (donormal){
  srt.normal = srt[, srt$type == 'normal']
  save(srt.normal, file = 'srt.normal.RData')
}

#### Malignant ####
srt = srt.malignant
srt = SCTransform(srt, return.only.var.genes = FALSE)
srt = RunPCA(srt, features = VariableFeatures(srt))
srt = RunUMAP(srt, dims = 1:10)
## NMF
# Run NMF
data = as.matrix(GetAssayData(srt, assay = 'SCT', slot = 'scale.data'))
#if (unique(srt$author) == ''){
  data = data[VariableFeatures(srt),]
#}
data[data < 0] = 0
data = data[apply(data, 1, var) > 0, ]
print(dim(data))
res.list = mclapply(range, function(r){
  nmf(data, rank = r, nrun = 1, seed = 'ica', method = 'nsNMF')
}, mc.cores = ncores)
names(res.list) = range
# Select rank
modules.list = lapply(res.list, NMFToModules, gmin = gmin)
print(sapply(modules.list,length))
comp = as.numeric(names(modules.list)) - sapply(modules.list, length)
mi = min(comp)
r = names(which(comp == mi))
r = r[length(r)]
print(r)
res = res.list[[r]]
# Process output
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
srt = AddMetaData(srt, t(coefs), col.name = rownames(coefs))
# Cluster NMF output
h = Heatmap(coefs, clustering_distance_columns = 'euclidean')
hcl = as.hclust(column_dend(h))
sig = cutree(hcl, k = length(modules))
nmf = c(by(t(coefs), INDICES = sig, FUN = function(x){names(modules)[which.max(colMeans(x))]}))[sig]
srt$nmf = factor(nmf, levels = names(modules))
col_nmf = c(brewer.pal(12, 'Set3'), brewer.pal(8, 'Set1'))[1:nlevels(srt$nmf)]
names(col_nmf) = levels(srt$nmf)
## Saving
srt.malignant = srt
save(srt.malignant, file = 'srt.malignant.RData')
res.list.malignant = res.list
save(res.list.malignant, file = 'res.list.malignant.RData')
res.malignant = res
save(res.malignant, file = 'res.malignant.RData')
## Plotting
pdf('dimplots.malignant.pdf', height = 15, width = 15)
for (reduction in c('pca','umap')){
  h = DimPlot(srt, pt.size = 2, reduction = reduction, group.by = 'nmf', cols = col_nmf)
  print(h)
  h = DimPlot(srt, pt.size = 2, reduction = reduction, group.by = 'cluster', cols = col_cluster)
  print(h)
  h = FeaturePlot(srt, reduction = reduction, features = names(modules))
  print(h)
}
dev.off()
pdf('heatmaps.malignant.pdf', height = 20, width = 10)
sink('heatmaps.malignant.txt')
# NMF
# Modules
top_ann = HeatmapAnnotation(df = data.frame('cluster' = srt$cluster, 'singler' = srt$singler,
                                            'nmf' = srt$nmf, 'subclone' = srt$subclone), 
                            col = list('cluster' = col_cluster, 'singler' = col_singler, 
                                       'nmf' = col_nmf, 'subclone' = col_subclone), which = 'column')
h = Heatmap(coefs, name = 'coefs', top_ann = top_ann,
            show_row_names = TRUE, 
            cluster_rows = FALSE,
            cluster_columns = FALSE, column_order = order(srt$nmf, srt$cluster, srt$singler),
            breaks = c(0,max(coefs)/2), colors = c('white','red'))
print(h)
h = Heatmap(scores[unlist(modules),], name = 'scores', 
            show_row_names = TRUE, row_names_gp = gpar(cex = 0.5),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            split = factor(unlist(mapply(rep, names(modules), sapply(modules, length))), levels = names(modules)))
print(h)
print(modules)
enrichment.matrix = BuildEnrichmentMatrix(rownames(srt), type = 'GO')
go_10 = sapply(modules, function(tbl){
  e = Enrichment(geneset = tbl, enrichment.matrix = enrichment.matrix)
  names(sort(e))[1:10]
})
colnames(go_10) = names(modules)
print(go_10)
sink()
dev.off()

#### Normal ####
if (donormal){
  srt = srt.normal
  srt = SCTransform(srt, return.only.var.genes = FALSE)
  srt = RunPCA(srt, features = VariableFeatures(srt))
  srt = RunUMAP(srt, dims = 1:10)
  ## Saving
  srt.normal = srt
  save(srt.normal, file = 'srt.normal.RData')
}

