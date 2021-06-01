#!/usr/bin/env Rscript
source('~/Documents/R/seurat_functions.R')
setwd('~/Documents/Analysis/Tumors/')

pdf('Finding.pdf', height = 10, width = 12)
sink('Finding.txt')

load('parameters.RData')

## Data ##

files = files.primary

res.list = lapply(files, function(f){
  loadRData(paste0(f,'res.malignant.RData'))
})
modules.list = lapply(res.list, NMFToModules)
genes.all = sort(unique(unlist(modules.list)))

srt.list = lapply(files, function(f){
  loadRData(paste0(f,'srt.malignant.RData'))
})
genes.all = sort(unique(unlist(lapply(srt.list, rownames))))
names(genes.all) = genes.all

types = c('GO','KEGG','REACTOME','HALLMARK','MOTIF')
enrichment.matrix.list = lapply(types, function(type){
  BuildEnrichmentMatrix(genes = genes.all, type = type)
})
names(enrichment.matrix.list) = types

ref = read.delim("PanglaoDB_markers_27_Mar_2020.tsv.gz")
ref$germ.layer = gsub('mesoderm','Mesoderm',ref$germ.layer)
byct = lapply(split(ref$official.gene.symbol, ref$cell.type), unique)
bygl = lapply(split(ref$official.gene.symbol, ref$germ.layer), unique)
byo = lapply(split(ref$official.gene.symbol, ref$organ), unique)
enrichment.matrix.list[['Cell.type']] = BuildEnrichmentMatrix(genes = genes.all, db = byct)
enrichment.matrix.list[['Germ.layer']] = BuildEnrichmentMatrix(genes = genes.all, db = bygl)
enrichment.matrix.list[['Organ']] = BuildEnrichmentMatrix(genes = genes.all, db = byo)
types = names(enrichment.matrix.list)

gbm = read.csv('gbm.csv', header = TRUE, stringsAsFactors = FALSE)
gbm = as.list(gbm)
gbm = lapply(gbm, function(x){x[!x=='']})
names(gbm) = paste0('Neftel_',names(gbm))
enrichment.matrix.list[['gbm']] = BuildEnrichmentMatrix(genes = genes.all, db = gbm)
types = names(enrichment.matrix.list)

hn = read.csv('hn.csv', header = TRUE, stringsAsFactors = FALSE)
hn = as.list(hn)
hn = lapply(hn, function(x){x[!x=='']})
names(hn) = paste0('Puram_',names(hn))
enrichment.matrix.list[['hn']] = BuildEnrichmentMatrix(genes = genes.all, db = hn)
types = names(enrichment.matrix.list)

sk = read.csv('sk.csv', header = TRUE, stringsAsFactors = FALSE)
sk = as.list(sk)
sk = lapply(sk, function(x){x[!x=='']})
names(sk) = paste0('Ji_',names(sk))
enrichment.matrix.list[['sk']] = BuildEnrichmentMatrix(genes = genes.all, db = sk)
types = names(enrichment.matrix.list)

cancer = factor(apply(sapply(strsplit(names(res.list), ''), '[', 1:4), 2, paste0, collapse = ''))
system = factor(corr[as.character(cancer),'system'])

## Modules from graph ##

modules.list = lapply(res.list, NMFToModules)
all = unlist(modules.list, recursive = FALSE, use.names = FALSE)
names(all) = unlist(sapply(modules.list, names))
ta = table(unlist(all))
genes.use = names(ta)[ta > 1]

# Filter non-overlapping modules
for (i in 1:5){
  all = unlist(modules.list, recursive = FALSE, use.names = TRUE)
  all = lapply(all, intersect, genes.all)
  sim = sapply(all, function(x){
    sapply(all, function(y){
      length(intersect(x,y))/length(union(x,y))
    })
  })
  keep = rownames(sim)[apply(sim, 1, function(x){
    sum(x > 0.05) >= 3
  })]
  all = all[keep]
  modules.list = lapply(names(modules.list), function(x){
    li = modules.list[[x]]
    li[names(li)[paste(x,names(li),sep='.') %in% keep]]
  })
  names(modules.list) = names(res.list)
  ta = table(unlist(all))
  genes.use = names(ta)[ta > 1] 
  print(length(all))
}

# Adjacency matrix, list by cancer
adj = matrix(0, nrow = length(genes.use), ncol = length(genes.use))
adj.list = list()
for (can in levels(cancer)){
  sub = matrix(0, nrow = length(genes.use), ncol = length(genes.use))
  rownames(sub) = genes.use
  colnames(sub) = genes.use
  for (s in names(modules.list)[cancer == can]){
    for (mod in modules.list[[s]]){
      mod = intersect(mod, genes.use)
      for (x in mod){
        for (y in mod){
          sub[x,y] = sub[x,y] + 1
        }
      }
    }
  }
  diag(sub) = 0
  adj.list[[can]] = sub
  #adj = adj + (sub > 0)
  adj = adj + sub
}
adj_keep = adj
adj = adj_keep

# Remove low connections
adj = adj_keep
adj[] = (adj >= v_min)
#adj[adj <= 1] = 0
for (i in 1:5){
  keep = names(which(rowSums(adj) >= s_min))
  adj = adj[keep,keep]
  print(dim(adj))
}

# Cluster
g = graph_from_adjacency_matrix(adj, diag = FALSE, mode = 'undirected', weighted = TRUE)
modules = communities(cluster_infomap(g, nb.trials = 100))
names(modules) = paste0('m_', sapply(modules, '[', 1))
modules = modules[c('m_ACYP1','m_ABHD3','m_ACTN1',
                    'm_ADM','m_ATOX1','m_FKBP5',
                    'm_AQP1','m_ACTB',
                    'm_AGER','m_ALDH3A1','m_AKR1B10','m_AGR2','m_AC007906.2',
                    'm_AGT','m_ASCL1','m_ARL4D')]
names(modules) = c('Cycle','Stress','Interferon',
                   'Hypoxia','Oxphos','Metal',
                   'Mesenchymal','pEMT',
                   'Alveolar','Basal','Squamous','Glandular','Ciliated',
                   'AC','OPC','NPC')
modules
save(modules, file = 'modules.RData')

colors.module = c(hue_pal()(length(modules)-3),
                  'thistle3',
                  'thistle2',
                  'thistle1')
names(colors.module) = names(modules)
save(colors.module, file = 'colors.module.RData')

types.module = c(rep('Process',6), rep('Type',length(modules)-6))
names(types.module) = names(modules)
save(types.module, file = 'types.module.RData')

fi = data.frame(sapply(modules,'[', 1:max(sapply(modules, length))), stringsAsFactors = FALSE)
fi[is.na(fi)] = ''
write.table(fi, file = 'modules.csv', sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE)

## Enrichment ##

# LUAD vs LUSC

data = read.csv('SCCvsAC.csv', row.names=1)
data$PValue.max = apply(data[,c('PValue.x', 'PValue.y')], 1, max)
data$logFC.mean = apply(data[,c('logFC.x', 'logFC.y')], 1, mean)
data$type = 'Other'
for (m in names(modules)){
  data[intersect(rownames(data),modules[[m]]), 'type'] = m
}
plot(data$logFC.mean, -log10(data$PValue.max), ylim = c(0,50), pch = 20, cex = 0.5, col = colors.module[data$type], xlab = 'LogFC', ylab = 'p-value (-log10)')
abline(v = 0, lty = 2)
text(data$logFC.mean[data$type %in% names(modules)], -log10(data$PValue.max)[data$type %in% names(modules)], labels = rownames(data)[data$type %in% names(modules)], cex = 1, col = colors.module[data$type[data$type %in% names(modules)]])
for (m in names(modules)){
  plot(data$logFC.mean, -log10(data$PValue.max), ylim = c(0,50), pch = 20, cex = 0.5, col = colors.module[data$type], xlab = 'LogFC', ylab = 'p-value (-log10)', main = m)
  abline(v = 0, lty = 2)
  tryCatch(expr = {
    text(data$logFC.mean[data$type == m], -log10(data$PValue.max)[data$type == m], labels = rownames(data)[data$type == m], cex = 1, col = colors.module[data$type[data$type == m]])
  }, error = function(e){c()})
}
pvals = sapply(names(modules), function(m){
  a = mHG.test((data$type == m)[order(data$logFC.mean, decreasing = FALSE)])$p.value
  b = mHG.test((data$type == m)[order(data$logFC.mean, decreasing = TRUE)])$p.value
  return(-log10(c(a,b)))
})
rownames(pvals) = c('Squamous','Adenocarcinoma')
h = Heatmap(pvals, name = 'mHG p-pvalue', 
            row_names_gp = gpar(cex = 0.3), column_names_gp = gpar(cex = 0.3, col = colors.module),
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = TRUE, show_column_names = TRUE,
            breaks = seq(0, 8, length = 9), colors = brewer_pal(palette = 'Greens')(n = 9))
print(h)

# Primary vs Met

data.list = lapply(c('5K','10K','30K','100K'), function(x){
  data = read.csv(paste0("~/Documents/Analysis/Tumors/PrimaryMet",x,".csv"))
  data = data[!duplicated(data[,1]),]
  rownames(data) = data[,1]
  colnames(data) = paste0(colnames(data),'.',x)
  return(data)
})
keep = sort(Reduce(intersect, lapply(data.list, rownames)))
data = Reduce(cbind, lapply(data.list, function(x){x[keep,]}))
data$pvalues.max = apply(data[,grep('pvalues',colnames(data))], 1, max)
data$log2fc.mean = apply(data[,grep('log2fc',colnames(data))], 1, mean)
data = data[data$pvalues.max < 0.01,]
data$type = 'Other'
for (m in names(modules)){
  data[intersect(rownames(data),modules[[m]]), 'type'] = m
}
plot(data$log2fc.mean, -log10(data$pvalues.max), xlim = c(-2,2), ylim = c(0,50), pch = 20, cex = 0.5, col = colors.module[data$type], xlab = 'LogFC', ylab = 'p-value (-log10)')
abline(v = 0, lty = 2)
text(data$log2fc.mean[data$type %in% names(modules)], -log10(data$pvalues.max)[data$type %in% names(modules)], labels = rownames(data)[data$type %in% names(modules)], cex = 1, col = colors.module[data$type[data$type %in% names(modules)]])
for (m in names(modules)){
  plot(data$log2fc.mean, -log10(data$pvalues.max), xlim = c(-2,2), ylim = c(0,50), pch = 20, cex = 0.5, col = colors.module[data$type], xlab = 'LogFC', ylab = 'p-value (-log10)', main = m)
  abline(v = 0, lty = 2)
  tryCatch(expr = {
    text(data$log2fc.mean[data$type == m], -log10(data$pvalues.max)[data$type == m], labels = rownames(data)[data$type == m], cex = 1, col = colors.module[data$type[data$type == m]])
  }, error = function(e){c()})
}
pvals = sapply(names(modules), function(m){
  a = mHG.test((data$type == m)[order(data$log2fc.mean, decreasing = FALSE)])$p.value
  b = mHG.test((data$type == m)[order(data$log2fc.mean, decreasing = TRUE)])$p.value
  return(-log10(c(a,b)))
})
rownames(pvals) = c('Primary','Metastasis')
h = Heatmap(pvals, name = 'mHG p-pvalue', 
            row_names_gp = gpar(cex = 0.3), column_names_gp = gpar(cex = 0.3, col = colors.module),
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = TRUE, show_column_names = TRUE,
            breaks = seq(0, 8, length = 9), colors = brewer_pal(palette = 'Greens')(n = 9))
print(h)

# Enrichment

enrichment = lapply(types, function(type){
  enrichment.matrix = enrichment.matrix.list[[type]]
  e = sapply(modules, function(mod){
    Enrichment(mod, enrichment.matrix = enrichment.matrix)
  })
  e = -log10(e)
  e[is.infinite(e)] = max(e[is.finite(e)])
  e = e[!apply(is.na(e), 1, any), ]
  rownames(e) = gsub(paste0(type,'_'), '', rownames(e))
  rownames(e) = gsub('_', ' ', rownames(e))
  rownames(e) = gsub('^ ', '', rownames(e))
  return(e)
})
names(enrichment) = types

for (e in enrichment){
  keep = unique(rownames(e)[apply(e, 2, function(x){order(x, decreasing = TRUE)[1:5]})])
  e = e[keep, ]
  h = Heatmap(e, name = 'Enrichment p-pvalue', 
              row_names_gp = gpar(cex = 0.6), column_names_gp = gpar(cex = 0.6, col = colors.module),
              cluster_rows = TRUE, cluster_columns = FALSE, row_dend_reorder = 1:length(keep),              
              show_row_names = TRUE, show_column_names = TRUE,
              breaks = seq(0, 8, length = 9), colors = brewer_pal(palette = 'Greens')(n = 9))
  print(h)
}

for (e in enrichment){
  keep = unique(rownames(e)[apply(e, 2, function(x){order(x, decreasing = TRUE)[1:5]})])
  e = e[keep, ]
  h = Heatmap(e, name = 'Enrichment p-pvalue', 
              row_names_gp = gpar(cex = 0.6), column_names_gp = gpar(cex = 0.6, col = colors.module),
              cluster_rows = FALSE, cluster_columns = FALSE,         
              show_row_names = TRUE, show_column_names = TRUE,
              breaks = seq(0, 8, length = 9), colors = brewer_pal(palette = 'Greens')(n = 9))
  print(h)
}

for (e in enrichment){
  keep = unique(rownames(e)[apply(e, 2, function(x){order(x, decreasing = TRUE)[1:2]})])
  e = e[keep, ]
  h = Heatmap(e, name = 'Enrichment p-pvalue', 
              row_names_gp = gpar(cex = 0.6), column_names_gp = gpar(cex = 0.6, col = colors.module),
              cluster_rows = FALSE, cluster_columns = FALSE,               
              show_row_names = TRUE, show_column_names = TRUE,
              breaks = seq(0, 8, length = 9), colors = brewer_pal(palette = 'Greens')(n = 9))
  print(h)
}

for (e in enrichment){
  keep = unique(rownames(e)[apply(e, 2, function(x){order(x, decreasing = TRUE)[1]})])
  e = e[keep, ]
  h = Heatmap(e, name = 'Enrichment p-pvalue', 
              row_names_gp = gpar(cex = 0.6), column_names_gp = gpar(cex = 0.6, col = colors.module),
              cluster_rows = FALSE, cluster_columns = FALSE,               
              show_row_names = TRUE, show_column_names = TRUE,
              breaks = seq(0, 8, length = 9), colors = brewer_pal(palette = 'Greens')(n = 9))
  print(h)
}

## Graph plot ##

adj = adj_keep
adj = adj[unlist(modules),unlist(modules)]
adj[adj <= 1] = 0
g = graph_from_adjacency_matrix(adj, diag = FALSE, mode = 'undirected', weighted = TRUE)

for (coords in list(layout.kamada.kawai(g), layout.fruchterman.reingold(g))){
  E(g)$width = E(g)$weight
  col = colors.module[unlist(mapply(rep, 1:length(modules), sapply(modules, length)))]
  names(col) = unlist(modules)
  col = col[!duplicated(names(col))]
  V(g)$color = col[names(V(g))]
  E(g)$color = 'lightgrey'
  plot(g, 
       layout = coords,
       vertex.label.color = V(g)$color, vertex.label.cex = 0.5,
       vertex.shape = 'none', vertex.size = 0,
       frame=FALSE)
  lab = sapply(names(V(g)), function(x){
    if (x %in% c('TOP2A','PCNA',
                 'JUN', 'FOS',
                 'IFIT1', 'STAT1','HLA-A','HLA-DRA',
                 'HILDPA','VEGFA',
                 'ATP5H','LAMTOR2',
                 'MT1G','MT1E',
                 'COL1A1','FN1',
                 'LAMC2','VIM',
                 'KRT5','KRT15',
                 'KLK10','LY6D',
                 'CLU','MUC5B','TFF1','CEACAM6',
                 'AGER','CAV1',
                 'FOXJ1','PIFO',
                 'ALDOC','APOE',
                 'OLIG1','OLIG2','SOX8',
                 'DLX1','DLX5','SOX11'
                 )){
      return(x)
    } else {
      return('')
    }
  })
  plot(g, 
       layout = coords,
       vertex.label.color = V(g)$color, vertex.label.cex = 1, vertex.label = lab,
       vertex.shape = 'circle', vertex.size = 1, vertex.frame.color = NA,
       frame=FALSE)
}

## Overlap plots ##

# Similarity
modules.list = lapply(res.list, NMFToModules)
all = unlist(modules.list, recursive = FALSE, use.names = FALSE)
names(all) = unlist(sapply(modules.list, names))
ta = table(unlist(all))
genes.use = names(ta)[ta > 1]
for (i in 1:5){
  all = unlist(modules.list, recursive = FALSE, use.names = TRUE)
  all = lapply(all, intersect, genes.all)
  sim = sapply(all, function(x){
    sapply(all, function(y){
      length(intersect(x,y))/length(union(x,y))
    })
  })
  keep = rownames(sim)[apply(sim, 1, function(x){
    sum(x > 0.05) >= 3
  })]
  all = all[keep]
  modules.list = lapply(names(modules.list), function(x){
    li = modules.list[[x]]
    li[names(li)[paste(x,names(li),sep='.') %in% keep]]
  })
  names(modules.list) = names(res.list)
  ta = table(unlist(all))
  genes.use = names(ta)[ta > 1] 
  print(length(all))
}
sim = sapply(all, function(x){
  sapply(modules, function(y){
    pval = phyper(length(intersect(x, y)), length(x), 2*10^4 - length(x), length(y), lower.tail = FALSE)
    return(-log10(pval))
  })
})
sim[is.infinite(sim)] = max(sim[is.finite(sim)])
sim = sim[,apply(sim, 2, function(x){any(x > 3)})]
df = data.frame('sample' = sapply(colnames(sim), function(x){
  y = sapply(strsplit(x, '.', fixed = TRUE), '[', 1)
}))
df$cancer = sapply(as.character(df$sample), function(x){
  y = sapply(strsplit(x, ''), '[', 1:4)
  paste0(y, collapse = '')
})
df$system = factor(corr[as.character(df$cancer),'system'])
df$top = apply(sim, 2, which.max)
df$top = factor(names(modules)[df$top], levels = names(modules))
ta = table(df[c('top','system')])
pval = sapply(colnames(ta), function(y){
  sapply(rownames(ta), function(x){
    phyper(ta[x,y], sum(ta[,y]), sum(ta) - sum(ta[,y]), sum(ta[x,]), lower.tail = FALSE)
  })
})
pval = -log10(pval)
pval[is.infinite(pval)] = max(pval[is.finite(pval)])
pval = (pval > max(3,max(apply(pval, 1, function(x){sort(x, decreasing = TRUE)[2]}))))
category = factor(apply(pval, 1, function(x){
  if (any(x)){
    return(names(which(x)))
  } else {
    return('Universal')
  }
}))
df$category = category[df$top]
top_ann = HeatmapAnnotation(df = df[, c('sample','category','top')], 
                            col = list('sample' = colors, 'category' = colors.system, 'top' = colors.module), 
                            which = 'column')
side_ann = HeatmapAnnotation(df = df[, c('sample','category','top')], 
                             col = list('sample' = colors, 'category' = colors.system, 'top' = colors.module), 
                             which = 'row')
h = Heatmap(name = 'Module overlap p-value', sim, 
            top_annotation = top_ann,
            cluster_columns = FALSE, cluster_rows = FALSE,
            column_order = order(df$top),
            show_row_names = TRUE, row_names_gp = gpar(col = colors.module),
            breaks = seq(0, 6, length = 7), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 7))
print(h)
decorate_heatmap_body('Module overlap p-value', {
  l = c(0,cumsum(table(df$top))[1:length(table(df$top))])
  for (k in 1:length(l)){
    i = unit(l[k]/ncol(sim), 'npc')
    j = unit(k/nrow(sim), 'npc')
    grid.lines(i, c(0,1), gp = gpar(lwd = 1, col = 'black'))
    grid.lines(c(0,1),j, gp = gpar(lwd = 1, col = 'black'))
  }
})
# Overlap
ovlp = sapply(all, function(x){
  sapply(all, function(y){
    pval = phyper(length(intersect(x, y)), length(x), 2*10^4 - length(x), length(y), lower.tail = FALSE)
    return(-log10(pval))
  })
})
ovlp[is.infinite(ovlp)] = max(ovlp[is.finite(ovlp)])
ovlp = ovlp[colnames(sim),colnames(sim)]
h = Heatmap(name = 'Module overlap p-value', ovlp, 
            top_annotation = top_ann,
            is.symmetric = TRUE, row_order = order(df$top),
            #show_row_names = TRUE, row_names_gp = gpar(cex = 0.4),
            breaks = seq(0, 6, length = 7), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 7)) + side_ann
print(h)
# decorate_heatmap_body('Module overlap p-value', {
#   l = c(0,cumsum(table(df$top))[1:length(table(df$top))])
#   for (k in 1:length(l)){
#     i = unit(l[k]/ncol(ovlp), 'npc')
#     j = unit(1-l[k]/nrow(ovlp), 'npc')
#     grid.lines(i, c(0,1), gp = gpar(lwd = 1, col = 'black'))
#     grid.lines(c(0,1), j, gp = gpar(lwd = 1, col = 'black'))
#   }
# })
decorate_heatmap_body('Module overlap p-value', {
  l = c(0,cumsum(table(df$top))[1:length(table(df$top))])
  for (k in 1:length(l)){
    i = unit(l[k:(k+1)]/ncol(ovlp), 'npc')
    j = unit(1-l[k:(k+1)]/nrow(ovlp), 'npc')
    grid.lines(i, j[1], gp = gpar(lwd = 1, col = 'black'))
    grid.lines(i, j[2], gp = gpar(lwd = 1, col = 'black'))
    grid.lines(i[1], j, gp = gpar(lwd = 1, col = 'black'))
    grid.lines(i[2], j, gp = gpar(lwd = 1, col = 'black'))
  }
})

# Schematic
for (s in names(srt.list)){
  mod = modules.list[[s]]
  srt = srt.list[[s]]
  mod = lapply(mod, sort)
  labels = unlist(lapply(names(mod), function(m){
    a = intersect(mod[[m]], modules[[df[paste(s, m, sep = '.'),'top']]])
    if (length(a) > 5){
      a = sample(a, size = 5, replace = FALSE)
    }
    return(a)
  }))
  data = GetData(srt, slot = 'scale.data')[unlist(mod), ]
  h = Heatmap(name = s, data,
              cluster_rows = FALSE,
              cluster_columns = TRUE, show_column_dend = FALSE,
              breaks = c(-3,0,3), colors = c('blue','white','red'))
  col = colors.module[unlist(mapply(rep, 1:length(modules), sapply(modules, length)))]
  names(col) = unlist(modules)
  col = col[!duplicated(names(col))]
  r = rowAnnotation(link = anno_mark(at = which(rownames(data) %in% labels), labels = labels,
                                     labels_gp = gpar(cex = 1.5, col = col[labels])))
  print(h+r)
  decorate_heatmap_body(s, {
    i = c(0,cumsum(sapply(rev(mod), length)))
    x = unit(i/max(i), 'npc')
    for (k in x){
      grid.lines(c(0,1), c(k,k), gp = gpar(lwd = 1, col = 'black'))
    }
    grid.lines(c(0,0), c(0,1), gp = gpar(lwd = 1, col = 'black'))
    grid.lines(c(1,1), c(0,1), gp = gpar(lwd = 1, col = 'black'))
  })
}

# Scenic ##

regulons.list = lapply(files, function(f){
  loadRData(paste0(f,'regulons.malignant.RData'))
})
regulons = unlist(regulons.list, recursive = FALSE, use.names = TRUE)
# Similarity
sim = sapply(regulons, function(x){
  sapply(modules[setdiff(names(modules), c('Stress','Cycle'))], function(y){
    pval = phyper(length(intersect(x, y)), length(x), 2*10^4 - length(x), length(y), lower.tail = FALSE)
    return(-log10(pval))
  })
})
sim[is.infinite(sim)] = max(sim[is.finite(sim)])
sim = sim[,apply(sim, 2, function(x){any(x > 3)})]
df = data.frame('sample' = sapply(colnames(sim), function(x){
  y = sapply(strsplit(x, '.', fixed = TRUE), '[', 1)
}))
df$tf = sapply(colnames(sim), function(x){
  y = sapply(strsplit(x, '.', fixed = TRUE), '[', 2)
})
df$cancer = sapply(as.character(df$sample), function(x){
  y = sapply(strsplit(x, ''), '[', 1:4)
  paste0(y, collapse = '')
})
df$system = factor(corr[as.character(df$cancer),'system'])
df$top = apply(sim, 2, which.max)
df$top = factor(rownames(sim)[df$top], levels = rownames(sim))
ta = table(df[c('top','system')])
top_ann = HeatmapAnnotation(df = df[, c('sample','top')],
                            col = list('sample' = colors, 'top' = colors.module),
                            which = 'column')
h = Heatmap(name = 'Scenic overlap p-value', sim,
            top_annotation = top_ann,
            cluster_columns = FALSE, cluster_rows = FALSE,
            column_order = order(df$top),
            show_row_names = TRUE, row_names_gp = gpar(col = colors.module),
            breaks = seq(0, 6, length = 7), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 7))
print(h)
decorate_heatmap_body('Scenic overlap p-value', {
  l = c(0,cumsum(table(df$top))[1:length(table(df$top))])
  for (k in 1:length(l)){
    i = unit(l[k]/ncol(sim), 'npc')
    j = unit(k/nrow(sim), 'npc')
    grid.lines(i, c(0,1), gp = gpar(lwd = 1, col = 'black'))
    grid.lines(c(0,1),j, gp = gpar(lwd = 1, col = 'black'))
  }
})
print(sapply(split(df$tf, df$top), function(tf){sort(table(tf), decreasing = TRUE)}))

regulons.list = lapply(files, function(f){
  loadRData(paste0(f,'regulons.malignant.RData'))
})
regulons = unlist(regulons.list, recursive = FALSE, use.names = TRUE)
# Similarity
sim = sapply(regulons, function(x){
  sapply(modules, function(y){
    pval = phyper(length(intersect(x, y)), length(x), 2*10^4 - length(x), length(y), lower.tail = FALSE)
    return(-log10(pval))
  })
})
sim[is.infinite(sim)] = max(sim[is.finite(sim)])
sim = sim[,apply(sim, 2, function(x){any(x > 3)})]
df = data.frame('sample' = sapply(colnames(sim), function(x){
  y = sapply(strsplit(x, '.', fixed = TRUE), '[', 1)
}))
df$tf = sapply(colnames(sim), function(x){
  y = sapply(strsplit(x, '.', fixed = TRUE), '[', 2)
})
df$cancer = sapply(as.character(df$sample), function(x){
  y = sapply(strsplit(x, ''), '[', 1:4)
  paste0(y, collapse = '')
})
df$system = factor(corr[as.character(df$cancer),'system'])
df$top = apply(sim, 2, which.max)
df$top = factor(rownames(sim)[df$top], levels = rownames(sim))
ta = table(df[c('top','system')])
top_ann = HeatmapAnnotation(df = df[, c('sample','top')],
                            col = list('sample' = colors, 'top' = colors.module),
                            which = 'column')
h = Heatmap(name = 'Scenic overlap p-value', sim,
            top_annotation = top_ann,
            cluster_columns = FALSE, cluster_rows = FALSE,
            column_order = order(df$top),
            show_row_names = TRUE, row_names_gp = gpar(col = colors.module),
            breaks = seq(0, 6, length = 7), colors = brewer_pal(palette = 'RdBu', direction = -1)(n = 7))
print(h)
decorate_heatmap_body('Scenic overlap p-value', {
  l = c(0,cumsum(table(df$top))[1:length(table(df$top))])
  for (k in 1:length(l)){
    i = unit(l[k]/ncol(sim), 'npc')
    j = unit(k/nrow(sim), 'npc')
    grid.lines(i, c(0,1), gp = gpar(lwd = 1, col = 'black'))
    grid.lines(c(0,1),j, gp = gpar(lwd = 1, col = 'black'))
  }
})
print(sapply(split(df$tf, df$top), function(tf){sort(table(tf), decreasing = TRUE)}))

sink()
dev.off()
