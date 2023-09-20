# libraries 
library(Seurat)
library(tidyverse)
library(patchwork)

# Perform analysis of various cluster resolutions
# pull data from the previous script (multimodal analysis, clustered, with TCR data)
bri.integrated <- readRDS("2-multimodal-sjs-clustered-tcr.rds")

## Using resolution 0.6, we combined cluster 0, cluster 3, and cluster 2 and then sub-clustered the resulting combination for annotation####

bri.integrated_cluster0 <- subset(bri.integrated, subset = (wsnn_res.0.6 %in% c(0,2,3)))
DefaultAssay(bri.integrated_cluster0) <- "integratedRNA"

## scale and pca ####
DefaultAssay(bri.integrated_cluster0) <- "integratedRNA"
bri.integrated_cluster0 <- ScaleData(bri.integrated_cluster0) %>% RunPCA(features=rownames(bri.integrated_cluster0))
DefaultAssay(bri.integrated_cluster0) <- 'integratedADT'
VariableFeatures(bri.integrated_cluster0) <- rownames(bri.integrated_cluster0[["ADT"]])
bri.integrated_cluster0 <- ScaleData(bri.integrated_cluster0) %>% RunPCA(reduction.name = 'apca')

# multimodal analysis
bri.integrated_cluster0 <- FindMultiModalNeighbors(bri.integrated_cluster0, reduction.list = list("pca", "apca"), 
                                                   dims.list = list(1:30, 1:12), modality.weight.name = "RNA.weight")

# UMAP for WNN
bri.integrated_cluster0 <- RunUMAP(bri.integrated_cluster0, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p1 <- DimPlot(bri.integrated_cluster0, reduction = "wnn.umap", group.by = "hash.ID")
p1
p2 <- DimPlot(bri.integrated_cluster0, reduction = 'rnaintegrated.umap', group.by = "hash.ID")
p3 <- DimPlot(bri.integrated_cluster0, reduction = 'adtintegrated.umap', group.by = "hash.ID")
p2 + p3

# Clustering for resolutions c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
bri.integrated_cluster0 <- FindClusters(bri.integrated_cluster0, graph.name = "wsnn", algorithm = 3, 
                                        resolution = c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0), verbose = FALSE)

# Visualization for resolutions c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
Idents(bri.integrated_cluster0) <- "wsnn_res.0.4"
p1 <- DimPlot(bri.integrated_cluster0, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.4")
Idents(bri.integrated_cluster0) <- "wsnn_res.0.6"
p2 <- DimPlot(bri.integrated_cluster0, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.6")
Idents(bri.integrated_cluster0) <- "wsnn_res.0.8"
p3 <- DimPlot(bri.integrated_cluster0, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.8")
Idents(bri.integrated_cluster0) <- "wsnn_res.1"
p4 <- DimPlot(bri.integrated_cluster0, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.0")
Idents(bri.integrated_cluster0) <- "wsnn_res.1.2"
p5 <- DimPlot(bri.integrated_cluster0, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.2")
Idents(bri.integrated_cluster0) <- "wsnn_res.1.4"
p6 <- DimPlot(bri.integrated_cluster0, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.4")
Idents(bri.integrated_cluster0) <- "wsnn_res.1.6"
p7 <- DimPlot(bri.integrated_cluster0, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.6")
Idents(bri.integrated_cluster0) <- "wsnn_res.1.8"
p8 <- DimPlot(bri.integrated_cluster0, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.8")
p1 + p2 + p3 + p4
p5 + p6 + p7 + p8


####
## Find ALL Markers #######
####

res <- c("0.4","0.6","0.8","1","1.2","1.4","1.6","1.8","2")
idents <- paste0("wsnn_res.",res)
rnamarkers <- adtmarkers <- top_rnamarkers <- top_adtmarkers <- list()

for(i in 1:length(idents)){
  
  # set ident
  Idents(bri.integrated_cluster0) <- idents[i]
  
  # RNA markers
  DefaultAssay(bri.integrated_cluster0) <- "integratedRNA"
  rnamarkers[[idents[i]]] <- FindAllMarkers(bri.integrated_cluster0, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  top_rnamarkers[[idents[i]]] <- rnamarkers[[idents[i]]] %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 12)
  
  # ADT markers
  DefaultAssay(bri.integrated_cluster0) <- "integratedADT"
  adtmarkers[[idents[i]]] <- FindAllMarkers(bri.integrated_cluster0, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  top_adtmarkers[[idents[i]]] <- adtmarkers[[idents[i]]] %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 12)
}

## We subclustered cluster 1 of res 0.6 for annotation####

## split and filter ####
bri.integrated_cluster1 <- subset(bri.integrated, subset = (wsnn_res.0.6 == "1"))
DefaultAssay(bri.integrated_cluster1) <- "integratedRNA"

## scale and pca ####
DefaultAssay(bri.integrated_cluster1) <- "integratedRNA"
bri.integrated_cluster1 <- ScaleData(bri.integrated_cluster1) %>% RunPCA(features=rownames(bri.integrated_cluster1))
DefaultAssay(bri.integrated_cluster1) <- 'integratedADT'
VariableFeatures(bri.integrated_cluster1) <- rownames(bri.integrated_cluster1[["ADT"]])
bri.integrated_cluster1 <- ScaleData(bri.integrated_cluster1) %>% RunPCA(reduction.name = 'apca')

# multimodal analysis
bri.integrated_cluster1 <- FindMultiModalNeighbors(bri.integrated_cluster1, reduction.list = list("pca", "apca"), 
                                                   dims.list = list(1:30, 1:12), modality.weight.name = "RNA.weight")

# UMAP for WNN
bri.integrated_cluster1 <- RunUMAP(bri.integrated_cluster1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p1 <- DimPlot(bri.integrated_cluster1, reduction = "wnn.umap", group.by = "hash.ID")
p1

# Clustering for resolutions c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
bri.integrated_cluster1 <- FindClusters(bri.integrated_cluster1, graph.name = "wsnn", algorithm = 3, 
                                        resolution = c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0), verbose = FALSE)

# # Visualization for resolutions c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
Idents(bri.integrated_cluster1) <- "wsnn_res.0.4"
p1 <- DimPlot(bri.integrated_cluster1, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.4")
Idents(bri.integrated_cluster1) <- "wsnn_res.0.6"
p2 <- DimPlot(bri.integrated_cluster1, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.6")
Idents(bri.integrated_cluster1) <- "wsnn_res.0.8"
p3 <- DimPlot(bri.integrated_cluster1, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.8")
Idents(bri.integrated_cluster1) <- "wsnn_res.1"
p4 <- DimPlot(bri.integrated_cluster1, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.0")
Idents(bri.integrated_cluster1) <- "wsnn_res.1.2"
p5 <- DimPlot(bri.integrated_cluster1, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.2")
Idents(bri.integrated_cluster1) <- "wsnn_res.1.4"
p6 <- DimPlot(bri.integrated_cluster1, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.4")
Idents(bri.integrated_cluster1) <- "wsnn_res.1.6"
p7 <- DimPlot(bri.integrated_cluster1, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.6")
Idents(bri.integrated_cluster1) <- "wsnn_res.1.8"
p8 <- DimPlot(bri.integrated_cluster1, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.8")
p1 + p2 + p3 + p4
p5 + p6 + p7 + p8


DefaultAssay(bri.integrated_cluster1) <- "integratedRNA"
features_RNA_CD8 <- c("IFNG", "NKG7", "GNLY", "GZMA","GZMB", "GZMH", "GZMK")
FeaturePlot(bri.integrated_cluster1, reduction="wnn.umap", features = features_RNA_CD8)

####
## Find ALL Markers #######
####

res <- c("0.4","0.6","0.8","1","1.2","1.4","1.6","1.8","2")
idents <- paste0("wsnn_res.",res)
rnamarkers <- adtmarkers <- top_rnamarkers <- top_adtmarkers <- list()

for(i in 1:length(idents)){
  
  # set ident
  Idents(bri.integrated_cluster1) <- idents[i]
  
  # RNA markers
  DefaultAssay(bri.integrated_cluster1) <- "integratedRNA"
  rnamarkers[[idents[i]]] <- FindAllMarkers(bri.integrated_cluster1, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  top_rnamarkers[[idents[i]]] <- rnamarkers[[idents[i]]] %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 12)
  
  # ADT markers
  DefaultAssay(bri.integrated_cluster1) <- "integratedADT"
  adtmarkers[[idents[i]]] <- FindAllMarkers(bri.integrated_cluster1, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  top_adtmarkers[[idents[i]]] <- adtmarkers[[idents[i]]] %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 12)
}


## Cluster 1 (of res0.6 from original cluster) and use res 0.4
Idents(bri.integrated_cluster1) <- "wsnn_res.0.4"
bri.integrated_cluster1 <- FindSubCluster(bri.integrated_cluster1, 
                                          cluster = 1, resolution = 0.5, graph.name = "wsnn")
Idents(bri.integrated_cluster1) <- "sub.cluster"
DimPlot(bri.integrated_cluster1, reduction = "wnn.umap", label = TRUE)

## We subclustered cluster 9 of res 0.6 for annotation ####

## split and filter ####
bri.integrated_cluster9 <- subset(bri.integrated, subset = (wsnn_res.0.6 == "9"))
DefaultAssay(bri.integrated_cluster9) <- "integratedRNA"

## scale and pca ####
DefaultAssay(bri.integrated_cluster9) <- "integratedRNA"
bri.integrated_cluster9 <- ScaleData(bri.integrated_cluster9) %>% RunPCA(features=rownames(bri.integrated_cluster9))
DefaultAssay(bri.integrated_cluster9) <- 'integratedADT'
VariableFeatures(bri.integrated_cluster9) <- rownames(bri.integrated_cluster9[["ADT"]])
bri.integrated_cluster9 <- ScaleData(bri.integrated_cluster9) %>% RunPCA(reduction.name = 'apca')

# multimodal analysis
bri.integrated_cluster9 <- FindMultiModalNeighbors(bri.integrated_cluster9, reduction.list = list("pca", "apca"), 
                                                   dims.list = list(1:30, 1:12), modality.weight.name = "RNA.weight")

# UMAP for WNN
bri.integrated_cluster9 <- RunUMAP(bri.integrated_cluster9, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p1 <- DimPlot(bri.integrated_cluster9, reduction = "wnn.umap", group.by = "hash.ID")
p1

# Clustering for resolutions c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
bri.integrated_cluster9 <- FindClusters(bri.integrated_cluster9, graph.name = "wsnn", algorithm = 3, 
                                        resolution = c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0), verbose = FALSE)

# # Visualization for resolutions c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
Idents(bri.integrated_cluster9) <- "wsnn_res.0.4"
p1 <- DimPlot(bri.integrated_cluster9, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.4")
Idents(bri.integrated_cluster9) <- "wsnn_res.0.6"
p2 <- DimPlot(bri.integrated_cluster9, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.6")
Idents(bri.integrated_cluster9) <- "wsnn_res.0.8"
p3 <- DimPlot(bri.integrated_cluster9, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.8")
Idents(bri.integrated_cluster9) <- "wsnn_res.1"
p4 <- DimPlot(bri.integrated_cluster9, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.0")
Idents(bri.integrated_cluster9) <- "wsnn_res.1.2"
p5 <- DimPlot(bri.integrated_cluster9, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.2")
Idents(bri.integrated_cluster9) <- "wsnn_res.1.4"
p6 <- DimPlot(bri.integrated_cluster9, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.4")
Idents(bri.integrated_cluster9) <- "wsnn_res.1.6"
p7 <- DimPlot(bri.integrated_cluster9, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.6")
Idents(bri.integrated_cluster9) <- "wsnn_res.1.8"
p8 <- DimPlot(bri.integrated_cluster9, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.8")
p1 + p2 + p3 + p4
p5 + p6 + p7 + p8

####
## Find ALL Markers #######
####

res <- c("0.4","0.6","0.8","1","1.2","1.4","1.6","1.8","2")
idents <- paste0("wsnn_res.",res)
rnamarkers <- adtmarkers <- top_rnamarkers <- top_adtmarkers <- list()

for(i in 1:length(idents)){
  
  # set ident
  Idents(bri.integrated_cluster9) <- idents[i]
  
  # RNA markers
  DefaultAssay(bri.integrated_cluster9) <- "integratedRNA"
  rnamarkers[[idents[i]]] <- FindAllMarkers(bri.integrated_cluster9, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  top_rnamarkers[[idents[i]]] <- rnamarkers[[idents[i]]] %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 12)
  
  # ADT markers
  DefaultAssay(bri.integrated_cluster9) <- "integratedADT"
  adtmarkers[[idents[i]]] <- FindAllMarkers(bri.integrated_cluster9, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  top_adtmarkers[[idents[i]]] <- adtmarkers[[idents[i]]] %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 12)
}

## We subclustered cluster 8 of res 0.6 for annotation ####

## split and filter ####
bri.integrated_cluster8 <- subset(bri.integrated, subset = (wsnn_res.0.6 == "8"))

## scale and pca ####
DefaultAssay(bri.integrated_cluster8) <- "integratedRNA"
bri.integrated_cluster8 <- ScaleData(bri.integrated_cluster8) %>% RunPCA(features=rownames(bri.integrated_cluster8))
DefaultAssay(bri.integrated_cluster8) <- 'integratedADT'
VariableFeatures(bri.integrated_cluster8) <- rownames(bri.integrated_cluster8[["ADT"]])
bri.integrated_cluster8 <- ScaleData(bri.integrated_cluster8) %>% RunPCA(reduction.name = 'apca')

# multimodal analysis
bri.integrated_cluster8 <- FindMultiModalNeighbors(bri.integrated_cluster8, reduction.list = list("pca", "apca"), 
                                                   dims.list = list(1:30, 1:12), modality.weight.name = "RNA.weight")

# UMAP for WNN
bri.integrated_cluster8 <- RunUMAP(bri.integrated_cluster8, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p1 <- DimPlot(bri.integrated_cluster8, reduction = "wnn.umap", group.by = "hash.ID")
p1

# Clustering for resolutions c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
bri.integrated_cluster8 <- FindClusters(bri.integrated_cluster8, graph.name = "wsnn", algorithm = 3, 
                                        resolution = c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0), verbose = FALSE)

# # Visualization for resolutions c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
Idents(bri.integrated_cluster8) <- "wsnn_res.0.4"
p1 <- DimPlot(bri.integrated_cluster8, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.4")
Idents(bri.integrated_cluster8) <- "wsnn_res.0.6"
p2 <- DimPlot(bri.integrated_cluster8, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.6")
Idents(bri.integrated_cluster8) <- "wsnn_res.0.8"
p3 <- DimPlot(bri.integrated_cluster8, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.8")
Idents(bri.integrated_cluster8) <- "wsnn_res.1"
p4 <- DimPlot(bri.integrated_cluster8, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.0")
Idents(bri.integrated_cluster8) <- "wsnn_res.1.2"
p5 <- DimPlot(bri.integrated_cluster8, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.2")
Idents(bri.integrated_cluster8) <- "wsnn_res.1.4"
p6 <- DimPlot(bri.integrated_cluster8, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.4")
Idents(bri.integrated_cluster8) <- "wsnn_res.1.6"
p7 <- DimPlot(bri.integrated_cluster8, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.6")
Idents(bri.integrated_cluster8) <- "wsnn_res.1.8"
p8 <- DimPlot(bri.integrated_cluster8, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.8")
p1 + p2 + p3 + p4
p5 + p6 + p7 + p8

####
## Find ALL Markers #######
####

res <- c("0.4","0.6","0.8","1","1.2","1.4","1.6","1.8","2")
idents <- paste0("wsnn_res.",res)
rnamarkers <- adtmarkers <- top_rnamarkers <- top_adtmarkers <- list()

for(i in 1:length(idents)){
  
  # set ident
  Idents(bri.integrated_cluster8) <- idents[i]
  
  # RNA markers
  DefaultAssay(bri.integrated_cluster8) <- "integratedRNA"
  rnamarkers[[idents[i]]] <- FindAllMarkers(bri.integrated_cluster8, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  top_rnamarkers[[idents[i]]] <- rnamarkers[[idents[i]]] %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 12)
  
  # ADT markers
  DefaultAssay(bri.integrated_cluster8) <- "integratedADT"
  adtmarkers[[idents[i]]] <- FindAllMarkers(bri.integrated_cluster8, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  top_adtmarkers[[idents[i]]] <- adtmarkers[[idents[i]]] %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 12)
}

## Save all subclusterings ####
save(bri.integrated_cluster0, bri.integrated_cluster1, bri.integrated_cluster9, bri.integrated_cluster8, file = "3-multimodal-sjs-clustered-tcr-subclusters.Rdata")


