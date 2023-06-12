# 20230611
# subclustering of cluster 4&11

library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(Signac)
library(cowplot)
library(dplyr)

output_dir = '/storage/htc/joshilab/Su_Li/StowersHSC/scenic_application/seurat_scenicInput/'
prefix = "hsc"
hsc = readRDS(paste0(output_dir,prefix,"_ob.rds"))
print(colnames(hsc[[]]))

hsc <- FindSubCluster(hsc, cluster = 4, graph.name = "wsnn", algorithm = 3, subcluster.name = "sub.cluster4")
print(colnames(hsc[[]]))

hsc <- FindSubCluster(hsc, cluster = 11, graph.name = "wsnn", algorithm = 3, subcluster.name = "sub.cluster11")
print(colnames(hsc[[]]))

# Generate new column based on conditions
data = hsc[[]]
hsc[["sub.cluster"]] <- ifelse(data$sub.cluster4 == data$sub.cluster11, data$sub.cluster4,
                          ifelse(grepl("^4_", data$sub.cluster4), data$sub.cluster4,
                                 ifelse(grepl("^11_", data$sub.cluster11), data$sub.cluster11, NA)))

hsc <- SetIdent(hsc, value = hsc@meta.data$sub.cluster)
markers <- FindAllMarkers(hsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

write.table(top10, file=paste(output_dir, prefix, "_", "sub.cluster", "_goodmarkers.txt",sep=""), row.names = TRUE, col.names = NA, sep="\t")
p1 = DimPlot(hsc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p2 = DimPlot(hsc, reduction = "wnn.umap", label = TRUE, group.by = "orig.ident", label.size = 2.5, repel = TRUE) + ggtitle("orig.ident")
p1 | p2
ggsave(paste0(output_dir,prefix,"_", "sub.cluster", "_wnn_umap_ident.png"), units = "in", width = 10, height = 5, dpi = 100)

saveRDS(hsc, paste0(output_dir,prefix,"_ob_sub_4_11.rds"))

DefaultAssay(object = hsc) <- "SCT"
c1_1 = hsc

c1_1 = DietSeurat(c1_1, assays = "SCT",  scale.data = TRUE, dimreducs = c("wnn.umap","tsne"), graphs = TRUE)

SaveH5Seurat(c1_1, filename = paste0(output_dir, prefix, "_integrated_4_11.h5Seurat"),overwrite = TRUE)
Convert(paste0(output_dir, prefix, "_integrated_4_11.h5Seurat"), dest = "h5ad", overwrite = TRUE)


