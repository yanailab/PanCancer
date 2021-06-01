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
numi.min = 500
ngene.min = 0
mt.max = 0.3
k_obs_groups = 2

# Combining
orig.ident = paste0(cancer, species, author, sample)
print(orig.ident)

#### Load ####
if (author == ''){
  files = list.files(recursive = TRUE, pattern = 'mtx')
  srt.list = lapply(1:length(files), function(lib){
    data = readMM(files[lib])
    colnames(data) = paste(orig.ident, lib, 1:dim(data)[2], sep = '_')
    rownames(data) = read.csv('genes.tsv', sep = '\t', stringsAsFactors = FALSE)$name
    srt = CreateSeuratObject(counts = data)
    meta.data = srt@meta.data
    meta.data$cancer = cancer
    meta.data$author = author
    meta.data$species = species
    meta.data$sample = as.character(sample)
    meta.data$library = lib
    meta.data$technique = 'inDrop'
    meta.data$orig.ident = orig.ident
    srt = AddMetaData(srt, meta.data)
    return(srt)
  })
  srt = Reduce(merge, srt.list)
  srt@meta.data[] = lapply(srt@meta.data, function(x){
    if (is.character(x)){ x = as.factor(x) }
    return(x)
  })
  Idents(srt) = 'orig.ident' 
} else {
  load('srt.RData')
}

#### Filter ####
pdf('filter.pdf')
if (!unique(srt$author) == 'Neftel' & !unique(srt$author) == 'Sharma'){
  srt = FilterCells(srt, do.plot = TRUE, numi.min = numi.min, ngene.min = ngene.min, mt.max = mt.max)
}
srt = FilterGenes(srt)
if (!unique(srt$author) == 'Neftel' & !unique(srt$author) == 'Sharma'){
  srt = SCTransform(srt, return.only.var.genes = FALSE)
} else {
  srt[['SCT']] = CreateAssayObject(counts = GetData(srt, assay = 'RNA', slot = 'counts'))
  srt = SetAssayData(srt, assay = 'SCT', slot = 'data', GetData(srt, assay = 'RNA', slot = 'counts'))
  DefaultAssay(srt) = 'SCT'
  srt = ScaleData(srt, do.center = TRUE, do.scale = FALSE, features = rownames(srt))
  srt = FindVariableFeatures(srt, selection.method = 'vst', nfeatures = 3000)
}
srt = RunPCA(srt, features = VariableFeatures(srt))
srt = RunUMAP(srt, dims = 1:10)
dev.off()
#### Cluster ####
srt = FindNeighbors(srt)
srt = FindClusters(srt)
srt$cluster = Idents(srt)
col_cluster = hue_pal()(nlevels(srt$cluster))
names(col_cluster) = levels(srt$cluster)

#### SingleR ####
# Loading expression data
exp_mat = GetAssayData(srt, assay = "SCT", slot = "data")
exp_mat = as.matrix(exp_mat)
# Loading singleR data
singler_se = celldex::HumanPrimaryCellAtlasData()
keep = c('Embryonic_stem_cells','iPS_cells','Tissue_stem_cells','Chondrocytes','Fibroblasts','Smooth_muscle_cells','Epithelial_cells','Hepatocytes','Keratinocytes','Astrocyte','Neurons','Neuroepithelial_cell',
  'B_cell','T_cells','NK_cell','DC','Neutrophils','Macrophage','Endothelial_cells')
singler_se = singler_se[, colData(singler_se)$label.main %in% keep]
# Keeping only common genes
common_genes = intersect(rownames(exp_mat), rownames(singler_se))
common_genes = sort(common_genes)
exp_common_mat = exp_mat[common_genes, ]
singler_se = singler_se[common_genes, ]
# Predicting, round 1
singler_pred = SingleR(
  test = exp_common_mat,
  ref = singler_se,
  labels = singler_se$label.main
)
keep = names(which(table(singler_pred$labels) >= 10))
# Predicting, round 2
singler_se = singler_se[, colData(singler_se)$label.main %in% keep]
singler_pred = SingleR(
  test = exp_common_mat,
  ref = singler_se,
  labels = singler_se$label.main
)
# Adding to srt
srt = AddMetaData(srt, as.data.frame(singler_pred))
srt$singler = factor(srt$labels)
col_singler = c(brewer.pal(8, 'Dark2'), brewer.pal(8, 'Accent'))[1:nlevels(srt$singler)]
names(col_singler) = levels(srt$singler)
print(table(srt$singler))

#### InferCNV ####
# Finding clones
dir.create('Clones', showWarnings = FALSE, recursive = TRUE, mode = "0777")
setwd('Clones')
annotations_file = data.frame(srt$orig.ident, row.names = rownames(srt@meta.data), stringsAsFactors = FALSE)
raw_counts_matrix = as.matrix(GetAssayData(srt, assay = 'SCT', slot = 'counts'))
if (unique(srt$species) == 'Hs'){
  gene_order_file = '~/Documents/CNV/Hs/gene_order.tsv'
}
if (unique(srt$species) == 'Mm'){
  gene_order_file = '~/Documents/CNV/Mm/gene_order.tsv'
}
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = raw_counts_matrix,
                                    annotations_file = annotations_file,
                                    gene_order_file = gene_order_file,
                                    delim = "\t",
                                    ref_group_names = NULL,
                                    chr_exclude = NULL)
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
obs_hcl = hclust(dist(t(infercnv_obj@expr.data), 'euclidean'), 'ward.D2')
clone = cutree(tree = as.hclust(obs_hcl), k = k_obs_groups)
srt$clone = 0
srt$clone[names(clone)] = clone[names(clone)]
srt$clone = factor(srt$clone)
col_clone = brewer_pal(palette = 'Reds')(1+nlevels(srt$clone))[-1]
names(col_clone) = levels(srt$clone)
# Plotting heatmap with clones
annotations_file = data.frame(srt$clone, row.names = rownames(srt@meta.data))
raw_counts_matrix = as.matrix(GetAssayData(srt, assay = 'SCT', slot = 'counts'))
if (unique(srt$species) == 'Hs'){
  gene_order_file = '~/Documents/CNV/Hs/gene_order.tsv'
}
if (unique(srt$species) == 'Mm'){
  gene_order_file = '~/Documents/CNV/Mm/gene_order.tsv'
}
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = raw_counts_matrix,
                                    annotations_file = annotations_file,
                                    gene_order_file = gene_order_file,
                                    delim = "\t",
                                    ref_group_names = NULL,
                                    chr_exclude = NULL)
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
save(infercnv_obj, file = 'infer_cnv.RData')
setwd('../')

#### Saving ####
save(srt, file = 'srt.RData')

#### Plotting ####
pdf('dimplots.pdf', height = 15, width = 15)
for (reduction in c('pca','umap')){
  h = DimPlot(srt, pt.size = 2, reduction = reduction, group.by = 'cluster', cols = col_cluster, label = TRUE)
  print(h)
  h = DimPlot(srt, pt.size = 2, reduction = reduction, group.by = 'singler', cols = col_singler)
  print(h)
  h = DimPlot(srt, pt.size = 2, reduction = reduction, group.by = 'clone', cols = col_clone)
  print(h)
}
dev.off()
pdf('heatmaps.pdf', height = 20, width = 10)
sink('heatmaps.txt')
# Clusters
markers = FindAllMarkers2(srt, do.print = TRUE, do.plot = TRUE, enrichment.type = 'GO', group.by = 'cluster', cols = col_cluster, print.bar = FALSE)
# SingleR
markers = FindAllMarkers2(srt, do.print = TRUE, do.plot = TRUE, enrichment.type = 'GO', group.by = 'singler', cols = col_singler, print.bar = FALSE)
sink()
dev.off()

