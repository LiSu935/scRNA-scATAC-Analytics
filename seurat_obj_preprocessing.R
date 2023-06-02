# this is for stowers HSC project.
library(ggplot2)

install.packages('BiocManager')
BiocManager::install('limma')


hsc = readRDS("multiome.motif.rds")
colnames(hsc[[]])
#  [1] "orig.ident"        "nCount_RNA"        "nFeature_RNA"
# [4] "nCount_ATAC"       "nFeature_ATAC"     "rna_clusters"
# [7] "nCount_SCT"        "nFeature_SCT"      "SCT_snn_res.0.5"
#[10] "seurat_clusters"   "ATAC_snn_res.0.5"  "nCount_ACTIVITY"
#[13] "nFeature_ACTIVITY"

#hsc@reductions$
#hsc@reductions$pca           hsc@reductions$harmony_atac
#hsc@reductions$harmony_rna   hsc@reductions$umap_atac
#hsc@reductions$umap_rna      hsc@reductions$wnn.umap
#hsc@reductions$lsi


hsc[["percent.mt"]] <- PercentageFeatureSet(hsc, pattern = "^MT-")
Idents(hsc) = hsc$"orig.ident"
VlnPlot(hsc, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()
  
hsc <- FindMultiModalNeighbors(hsc, reduction.list = list("harmony_rna", "harmony_atac"), dims.list = list(1:50, 1:49))
hsc <- RunUMAP(hsc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
hsc <- FindClusters(hsc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 = DimPlot(hsc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p2 = DimPlot(hsc, reduction = "wnn.umap", label = TRUE, group.by = "orig.ident", label.size = 2.5, repel = TRUE) + ggtitle("orig.ident")
p1 | p2

# unscaled :
hsc@assays$SCT@data
# scaled:
hsc@assays$SCT@scale.data

# find the cluster name. print it out. Then
print(colnames(hsc[[]]))

# output the followings as the seurat to scenic:
# unscaled :
hsc@assays$SCT@data
# scaled:
hsc@assays$SCT@scale.data 

DefaultAssay = "RNA"
markers <- FindAllMarkers(hsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    
markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(hsc, features = top10$gene) + NoLegend()  
ggsave("hsc_wnn_0.5_sct_heatmap.png"), units = "in", width = 15, height = 15, dpi = 100)





# just for checking:

day0_10x = Read10X_h5("/storage/htc/joshilab/Su_Li/StowersHSC/data/round2/Day0.rna_atac.filtered_feature_barcode_matrix.h5")

# extract RNA and ATAC data
rna_counts <- day0_10x$`Gene Expression`
atac_counts <- day0_10x$Peaks

# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "../data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
pbmc[["ATAC"]] <- chrom_assay
