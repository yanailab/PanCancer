#!/usr/bin/env Rscript
source('~/Documents/R/seurat_functions.R')
load('srt.malignant.RData')
srt = srt.malignant

# Exporting for Scenic
dir.create('Scenic.Malignant', showWarnings = FALSE, recursive = TRUE, mode = "0777")
setwd('Scenic.Malignant')
expr = as.matrix(GetAssayData(srt, assay = 'SCT', slot = 'data'))
expr.scale = as.matrix(GetAssayData(srt, assay = 'SCT', slot = 'scale.data'))
minCountsPerGene = max(10, 0.05*ncol(expr))
minSamples = max(10, ncol(expr)*0.05)
if (unique(srt$species) == 'Hs'){
  org = "hgnc"
}
if (unique(srt$species) == 'Mm'){
  org = "mgi"
}
dbDir= "~/Documents/Scenic"
data(defaultDbNames)
dbs = defaultDbNames[[org]]
scenicOptions = initializeScenic(org=org, dbDir=dbDir, nCores=10) 
scenicOptions@settings$verbose = TRUE
scenicOptions@settings$nCores = 10
scenicOptions@settings$seed = 123
scenicOptions@settings$defaultTsne$aucType = "AUC"
scenicOptions@settings$defaultTsne$dims = 10
scenicOptions@settings$defaultTsne$perpl = round(sqrt(dim(expr)[2]))
saveRDS(scenicOptions, file="scenicOptions.Rds")
genes.sub = geneFiltering(expr, scenicOptions=scenicOptions,
                          minCountsPerGene=minCountsPerGene,
                          minSamples=minSamples)
expr.sub = expr[genes.sub, ]
expr.scale.sub = expr.scale[genes.sub, ]
exportsForArboreto(expr.scale.sub, scenicOptions, dir = getwd())
setwd('../')

# Saving
srt.malignant = srt
save(srt.malignant, file = 'srt.malignant.RData')