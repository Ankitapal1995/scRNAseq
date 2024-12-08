library(Seurat)
library(dplyr)
library(Matrix)
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size
sparse.size <- object.size(x = pbmc.data) 
sparse.size

dense.size/sparse.size
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200, project = "10X_PBMC")

# Access raw count data
counts <- GetAssayData(object = pbmc, layer = "counts", assay = "RNA")

nUMI <- Matrix::colSums(GetAssayData(object = pbmc, layer = "counts", assay = "RNA"))
nGene <- Matrix::colSums(GetAssayData(object = pbmc, layer = "counts", assay = "RNA") > 0)
pbmc <- AddMetaData(object = pbmc, metadata = nGene, col.name = "nGene")
pbmc <- AddMetaData(object = pbmc, metadata = nUMI, col.name = "nUMI")
# Identify mitochondrial genes
mito.genes <- grep("^MT-", rownames(counts), value = TRUE)

# Calculate percentage of mitochondrial gene expression
percent.mito <- Matrix::colSums(counts[mito.genes, ]) / Matrix::colSums(counts) * 100

# Add as metadata to the Seurat object
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

VlnPlot(object = pbmc, features = c("nGene", "nUMI", "percent.mito"), ncol = 3)

# content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = pbmc, feature1 = "nUMI", feature2 = "percent.mito")
FeatureScatter(object = pbmc, feature1 = "nUMI", feature2 = "nGene")
pbmc <- subset(pbmc, subset = nGene > 200 & nGene < 2500 & percent.mito < 0.05)

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                             selection.method = "vst", x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))


DimPlot(pbmc, reduction = "pca", dims = c(1, 2))

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# PCElbowPlot() can be used to reduce computation time
pbmc <- JackStraw(object = pbmc, num.replicate = 100, do.print = FALSE)

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = pbmc)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = pbmc)
save(pbmc, file = "./pbmc.Robj")
# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = c(0, 3), 
                                min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, thresh.use = 0.25, 
                                test.use = "roc", only.pos = TRUE)

VlnPlot(object = pbmc, features.plot = c("MS4A1", "CD79A"))

VlnPlot(object = pbmc, features.plot = c("NKG7", "PF4"), use.raw = TRUE, y.log = TRUE)

FeaturePlot(object = pbmc, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", 
                                             "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", 
                     "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5)

# First lets stash our identities for later
pbmc <- StashIdent(object = pbmc, save.name = "ClusterNames_0.6")

# Note that if you set save.snn=T above, you don't need to recalculate the
# SNN, and can simply put: pbmc <- FindClusters(pbmc,resolution = 0.8)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.8, print.output = FALSE)

# points based on different criteria
plot1 <- TSNEPlot(object = pbmc, do.return = TRUE, no.legend = TRUE, do.label = TRUE)
plot2 <- TSNEPlot(object = pbmc, do.return = TRUE, group.by = "ClusterNames_0.6", 
                  no.legend = TRUE, do.label = TRUE)
plot_grid(plot1, plot2)

# Find discriminating markers
tcell.markers <- FindMarkers(object = pbmc, ident.1 = 0, ident.2 = 1)

# palette from low to high expression
FeaturePlot(object = pbmc, features.plot = c("S100A4", "CCR7"), cols.use = c("green", 
                                                                             "blue"))
pbmc <- SetAllIdent(object = pbmc, id = "ClusterNames_0.6")
save(pbmc, file = "./pbmc3k_final.Rda")
pbmc <- readRDS("pbmc_seurat_object.rds")