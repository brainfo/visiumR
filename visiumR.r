set.seed(seed = 404)
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(SeuratDisk)

rm_cache <- function() {
    rm(de_markers)
    rm(top2)
    rm(top5)
    rm(top5.heatmap)
    rm(feature.plot)
    rm(p1)
    rm(p2)
    rm(p3)
    rm(p4)
    rm(unique.top2)
    rm(top.features)
}


run_single <- function(sample_id){
data_dir <- paste("testdata", sample_id,sep = '/')
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
manual <- Load10X_Spatial(data_dir, filename = paste(sample_id, "filtered_feature_bc_matrix.h5", sep='_'), assay = "Spatial", slice = sample_id, filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

manual <- RenameGenesSeurat(obj = manual, newnames = new.names)

## Save a table of all.genes
all.genes <- row.names(manual)
write.table(all.genes, "all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
manual <- PercentageFeatureSet(manual, pattern = "^MT-", col.name = "percent.mt")
manual <- PercentageFeatureSet(manual, pattern = "^RP[SL]", col.name = "percent.ribo")
manual <- PercentageFeatureSet(manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
manual<- NormalizeData(manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
manual <- CellCycleScoring(manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
manual$CC.Difference <- manual$S.Score - manual$G2M.Score

acute.viral <- subset(manual, features = rownames(manual)[rownames(manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

# ncol(acute.viral)   # 1505
# acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(manual)[rownames(manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

manual$orig.ident <- sample_id

## plot UMI counts/spot 
plot1 <- VlnPlot(manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave(paste0( sample_id, "manual_unfiltered_nCount_Spatial.pdf"), plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave(paste0( sample_id, "manual_unfiltered_nFeature_Spatial.pdf"), plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave(paste0( sample_id, "manual_unfiltered_QC_vlnplot.pdf"), plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(manual)
  #  1505
manual2 <- subset(manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(manual2)
  #   1480
qc.vlnplot <- VlnPlot(manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave(paste0( sample_id, "manual_filtered_QC_vlnplot.pdf"), plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
manual <- manual2
rm(manual2)
qc.vlnplot <- VlnPlot(manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave(paste0( sample_id, "manual_filtered_QC_vlnplot.pdf"), plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave(paste0( sample_id, "manual_filtered__Spatial.pdf"), plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S15: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
manual <- SCTransform(manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
manual <- RunPCA(manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(manual)
ggsave(paste0( sample_id, "manual_elbow_plot.pdf"), plot = elbow.plot, device = "pdf")
rm(elbow.plot)

manual <- FindNeighbors(manual, reduction = "pca", dims = 1:30)
manual <- FindClusters(manual, verbose = T)
manual <- RunUMAP(manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave(paste0( sample_id, "manual_UMAP.pdf"), plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(manual, cells.highlight = CellsByIdentities(object = manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave(paste0( sample_id, "manual_UMAP_Spatial-split.by.cluster.pdf"), plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
manual <- FindSpatiallyVariableFeatures(manual, assay = "SCT", features = VariableFeatures(manual)[1:1000], selection.method = "markvariogram")
manual.var <- manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave(paste0( sample_id, "manual_spatial.var.features.pdf"), plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
manual <- NormalizeData(manual, normalization.method = "LogNormalize", scale.factor = 10000)
manual <- ScaleData(manual, features = manual.var, assay="Spatial")

Idents(manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(manual, features = manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S15_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(manual, features = manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# View(top5)

top5.heatmap <- DoHeatmap(manual, features = top5$gene, raster = FALSE)
ggsave(paste0( sample_id, "manual_top5_markers_logfc_heatmaptop3k.pdf"), plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave(paste0( sample_id, "manual_top5_markers_counts_heatmaptop3k.pdf"), plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(manual, features = unique.top2)
png(file=paste0( sample_id, "manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(manual, features = unique.top2)
png(file=paste0( sample_id, "manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(manual, features = unique.top2)
png(file=paste0( sample_id, "manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm_cache()
}