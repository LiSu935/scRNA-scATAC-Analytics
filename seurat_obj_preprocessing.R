# this is for stowers HSC project. """conda activate r-mofa"""
library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(Signac)
library(cowplot)
library(dplyr)

hsc = readRDS("/storage/htc/joshilab/Su_Li/StowersHSC/from_ningzhang/multiome.motif.rds")
prefix = "hsc"
output_dir = '/storage/htc/joshilab/Su_Li/StowersHSC/scenic_application/seurat_scenicInput/'

print(colnames(hsc[[]]))
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

# find the cluster name. print it out. Then
print(colnames(hsc[[]]))

# output the followings as the seurat to scenic:
# unscaled :
DefaultAssay(object = hsc) <- "SCT"
c1_1 = hsc

c1_1 = DietSeurat(c1_1, assays = "SCT",  scale.data = FALSE, dimreducs = c("harmony_rna", "pca"), graphs = TRUE)
# Already saved once from interactive job, skip now
#SaveH5Seurat(c1_1, filename = paste0(output_dir, prefix, "_slim.h5Seurat"),overwrite = TRUE)
#Convert(paste0(output_dir, prefix, "_slim.h5Seurat"), dest = "h5ad", overwrite = TRUE)


# scaled:
#hsc@assays$SCT@scale.data 

DefaultAssay = "RNA"
markers <- FindAllMarkers(hsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    
markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

v=colnames(hsc[[]])
suffix = v[length(v)-1]

write.table(top10, file=paste(output_dir, prefix, "_", suffix, "_goodmarkers.txt",sep=""), row.names = TRUE, col.names = NA, sep="\t")

DoHeatmap(hsc, features = top10$gene) + NoLegend()  

ggsave(paste0(output_dir, prefix, "_", suffix,"_heatmap.png"), units = "in", width = 15, height = 15, dpi = 100)

#saveRDS(hsc, paste0(output_dir,prefix,"_ob.rds"))

output_dir = '/storage/htc/joshilab/Su_Li/StowersHSC/scenic_application/seurat_scenicInput/'
prefix = "hsc"
#hsc = readRDS(paste0(output_dir,prefix,"_ob.rds"))
print(colnames(hsc[[]]))

hsc <- FindClusters(hsc, graph.name = "wsnn", algorithm = 3, verbose = FALSE, res = 0.6)
hsc <- FindClusters(hsc, graph.name = "wsnn", algorithm = 3, verbose = FALSE, res = 0.4)

v=colnames(hsc[[]])
suffix = v[length(v)]

p1 = DimPlot(hsc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p2 = DimPlot(hsc, reduction = "wnn.umap", label = TRUE, group.by = "orig.ident", label.size = 2.5, repel = TRUE) + ggtitle("orig.ident")
p1 | p2
ggsave(paste0(output_dir,prefix,"wnn_umap_ident.png"), units = "in", width = 10, height = 5, dpi = 100)

DefaultAssay = "RNA"
markers <- FindAllMarkers(hsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    
markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

write.table(top10, file=paste(output_dir, prefix, "_", suffix, "_goodmarkers.txt",sep=""), row.names = TRUE, col.names = NA, sep="\t")

DoHeatmap(hsc, features = top10$gene) + NoLegend()  
ggsave(paste0(output_dir, prefix, "_", suffix,"_heatmap.png"), units = "in", width = 15, height = 15, dpi = 100)

saveRDS(hsc, paste0(output_dir,prefix,"_ob.rds"))

#DefaultAssay(object = hsc) <- "SCT"
#c1_1 = hsc

#c1_1 = DietSeurat(c1_1, assays = "SCT",  scale.data = TRUE, dimreducs = c("wnn.umap","harmony_rna"), graphs = TRUE)

#SaveH5Seurat(c1_1, filename = paste0(output_dir, prefix, "_integrated.h5Seurat"),overwrite = TRUE)
#Convert(paste0(output_dir, prefix, "_integrated.h5Seurat"), dest = "h5ad", overwrite = TRUE)


