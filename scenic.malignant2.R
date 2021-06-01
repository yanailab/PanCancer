#!/usr/bin/env Rscript
source('~/Documents/R/seurat_functions.R')
load('srt.malignant.RData')
srt = srt.malignant

# Running RScenic
setwd('Scenic.Malignant')
expr = as.matrix(GetAssayData(srt, assay = 'SCT', slot = 'data'))
expr.scale = as.matrix(GetAssayData(srt, assay = 'SCT', slot = 'scale.data'))
minCountsPerGene = max(10, 0.05*ncol(expr))
minSamples = max(10, ncol(expr)*0.05)
scenicOptions <- readRDS("scenicOptions.Rds")
GRNBoost_linkList = importArboreto("grnboost_output.csv")
colnames(GRNBoost_linkList) = c('TF','Target','weight')
saveRDS(GRNBoost_linkList, file = getIntName(scenicOptions, "genie3ll"))
genes.sub = geneFiltering(expr, scenicOptions=scenicOptions,
                          minCountsPerGene=minCountsPerGene,
                          minSamples=minSamples)
expr.sub = expr[genes.sub, ]
expr.scale.sub = expr.scale[genes.sub, ]
corrMat = cor(t(expr.scale.sub), method = "spearman")
saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
regulons = read.table('output/Step2_regulonTargetsInfo.tsv', header = TRUE, stringsAsFactors = FALSE)
regulons = split(regulons$gene, regulons$TF)
setwd('../')

# Saving
save(regulons, file = 'regulons.malignant.RData')
