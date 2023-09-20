# For supplemental figure, testing independent clustering of pbmcs and independent clustering of skin vs combined


# libraries
library(tidyverse)
library(Seurat)
library(patchwork)

library(ComplexHeatmap)
library(circlize)
library(speckle)
library(scales)

#### 
## Data Import ##############
####

# # read cellranger 10x filtered feature matrices 
bri793 <- Read10X(data.dir="BRI-793/cite/filtered_feature_bc_matrix/")
bri793hash <- Read10X(data.dir="BRI-793/hash/filtered_feature_bc_matrix/")
bri817 <- Read10X(data.dir="BRI-817/cite/filtered_feature_bc_matrix/")
bri817hash <- Read10X(data.dir="BRI-817/hash/filtered_feature_bc_matrix/")
bri820 <- Read10X(data.dir="BRI-820/cite/filtered_feature_bc_matrix/")
bri820hash <- Read10X(data.dir="BRI-820/hash/filtered_feature_bc_matrix/")

# appoint ADT, HTO and GEX to Seurat objects
bri793.joint.bcs <- intersect(colnames(bri793[["Antibody Capture"]]),colnames(bri793hash[["Antibody Capture"]]))
bri793.obj <- CreateSeuratObject(counts=bri793[["Gene Expression"]][, bri793.joint.bcs])
bri793.obj[["HTO"]] <- CreateAssayObject(bri793hash[["Antibody Capture"]][, bri793.joint.bcs])
bri793.obj[["ADT"]] <- CreateAssayObject(bri793[["Antibody Capture"]][, bri793.joint.bcs])

bri817.joint.bcs <- intersect(colnames(bri817[["Antibody Capture"]]),colnames(bri817hash[["Antibody Capture"]]))
bri817.obj <- CreateSeuratObject(counts=bri817[["Gene Expression"]][, bri817.joint.bcs])
bri817.obj[["HTO"]] <- CreateAssayObject(bri817hash[["Antibody Capture"]][, bri817.joint.bcs])
bri817.obj[["ADT"]] <- CreateAssayObject(bri817[["Antibody Capture"]][, bri817.joint.bcs])

bri820.joint.bcs <- intersect(colnames(bri820[["Antibody Capture"]]),colnames(bri820hash[["Antibody Capture"]]))
bri820.obj <- CreateSeuratObject(counts=bri820[["Gene Expression"]][, bri820.joint.bcs])
bri820.obj[["HTO"]] <- CreateAssayObject(bri820hash[["Antibody Capture"]][, bri820.joint.bcs])
bri820.obj[["ADT"]] <- CreateAssayObject(bri820[["Antibody Capture"]][, bri820.joint.bcs])

rownames(bri793.obj[["ADT"]])
rownames(bri817.obj[["ADT"]])
rownames(bri820.obj[["ADT"]])

#### 
## HTO Analysis ###########
####

# Normalize HTO with central log transformation, use 99th quantile for cutoff
bri793.obj <- NormalizeData(bri793.obj, assay = "HTO", normalization.method = "CLR")
bri793.obj <- HTODemux(bri793.obj, assay = "HTO", positive.quantile = 0.99)
bri817.obj <- NormalizeData(bri817.obj, assay = "HTO", normalization.method = "CLR")
bri817.obj <- HTODemux(bri817.obj, assay = "HTO", positive.quantile = 0.99)
bri820.obj <- NormalizeData(bri820.obj, assay = "HTO", normalization.method = "CLR")
bri820.obj <- HTODemux(bri820.obj, assay = "HTO", positive.quantile = 0.99)

# show how many singlets & doublets by HTO counts (still possibly for two of the same HTO to be counted together as 'singlet')
p1 <- table(bri793.obj@meta.data$HTO_classification.global) %>%
  as.data.frame()  %>%
  ggplot(., aes(x=Var1, y = Freq, label = as.character(Freq))) + 
  geom_bar(aes(fill=Var1), position = "dodge", stat = "identity", show.legend = FALSE) +
  geom_label(label.size = 2) + 
  labs(title = "BRI793: # of Cells per HTO Class", x = "HTO Class", y = "# of Cells") 

p2 <- table(bri817.obj@meta.data$HTO_classification.global) %>%
  as.data.frame()  %>%
  ggplot(., aes(x=Var1, y = Freq, label = as.character(Freq))) + 
  geom_bar(aes(fill=Var1), position = "dodge", stat = "identity", show.legend = FALSE) +
  geom_label(label.size = 2) + 
  labs(title = "BRI817: # of Cells per HTO Class", x = "HTO Class", y = "# of Cells") 

p3 <- table(bri820.obj@meta.data$HTO_classification.global) %>%
  as.data.frame()  %>%
  ggplot(., aes(x=Var1, y = Freq, label = as.character(Freq))) + 
  geom_bar(aes(fill=Var1), position = "dodge", stat = "identity", show.legend = FALSE) +
  geom_label(label.size = 2) + 
  labs(title = "BRI820: # of Cells per HTO Class", x = "HTO Class", y = "# of Cells") 

p1 + p2 + p3

# Visualize
FeatureScatter(bri820.obj, feature1 = "SJS001SKIN", feature2 = "SJS001PBMC") + 
  labs(title = "HTO Clusters for Singlet Detection for BRI820", x = "Norm. Exp. of HTO SJS001SKIN", y = "Norm. Exp. of HTO SJS001PBMC") + 
  NoLegend()

FeatureScatter(bri820.obj, feature1 = "SJS001SKIN", feature2 = "SJS001PBMC", pt.size = 2) + 
  labs(title = "", x = "Hashtag 1 (HTO)", y = "Hashtag 2 (HTO)") + 
  NoLegend() + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                     axis.text.x = element_blank(), axis.text.y = element_blank()) + 
  ylim(-0.1, 4.3)

# filter out doublets and negatives
bri793.obj <- subset(bri793.obj, subset=HTO_classification.global == 'Singlet')
bri817.obj <- subset(bri817.obj, subset=HTO_classification.global == 'Singlet')
bri820.obj <- subset(bri820.obj, subset=HTO_classification.global == 'Singlet')

# define lists of Seurat Objects based on the the hash.id, 
bri793.list <- SplitObject(bri793.obj, split.by="hash.ID")
bri793.list <- bri793.list[-3] # remove Control030Skin which has very low cell #'s based on VlnPlots from downstream steps
bri817.list <- SplitObject(bri817.obj, split.by="hash.ID")
bri820.list <- SplitObject(bri820.obj, split.by="hash.ID")
bri.list <- do.call(c, list(bri793.list,bri817.list,bri820.list))

####
# Filtering RNA and ADT, Normalize, Find variable features
# low RNA 1000, high RNA 10000 
# remove ribosomal proteins (RPS and RPL) and remove CD4 and CD8 doublets
# Normalize, filter on UMIs based on previous Ridgeplot
# then Integrate on the RNA expression
names(bri.list) #above has removed Ctrl030Skin which only had few cells.
adt_features <- paste0("adt_",rownames(bri.list[[1]][['ADT']]))

for (i in 1:length(bri.list)) {
  DefaultAssay(bri.list[[i]]) <- "RNA"
  bri.list[[i]][["percent.mt"]] <- PercentageFeatureSet(bri.list[[i]], pattern = "^MT-")
  bri.list[[i]] <- subset(bri.list[[i]], subset = nCount_RNA > 1000 & nCount_RNA < 20000 & percent.mt < 20)
  bri.list[[i]] <- subset(bri.list[[i]], features = c(rownames(bri.list[[i]][['RNA']])[!grepl("^RP[SL]|^MRP[SL]",rownames(bri.list[[i]][['RNA']]))],adt_features))
  bri.list[[i]] <- NormalizeData(bri.list[[i]], verbose = FALSE)
  bri.list[[i]] <- FindVariableFeatures(bri.list[[i]], selection.method = "vst",
                                        nfeatures = 2000, verbose = FALSE)
}

####
## Integration for RNA and ADT: #######
####

bri.list.pbmc <- bri.list[grepl("PBMC",names(bri.list))]
bri.list.skin <- bri.list[!grepl("PBMC",names(bri.list))]

### PBMC: #######

# Integration for RNA 
bri.anchors <- FindIntegrationAnchors(object.list=bri.list.pbmc, dims=1:30)
bri.integrated.pbmc <- IntegrateData(anchorset = bri.anchors, new.assay.name = "integratedRNA", dims=1:30, k.weight=30)

# split objects (which are integratedRNA) to RE-integrate on ADT
# remove cd4 and cd8 doublets
bri.integrated.pbmc.list <- SplitObject(bri.integrated.pbmc, split.by = "hash.ID")

for (i in 1:length(bri.integrated.pbmc.list)) {
  DefaultAssay(bri.integrated.pbmc.list[[i]]) <- "ADT"
  bri.integrated.pbmc.list[[i]] <- NormalizeData(bri.integrated.pbmc.list[[i]], normalization.method = "CLR")
  bri.integrated.pbmc.list[[i]] <- subset(bri.integrated.pbmc.list[[i]], subset = adt_CD4.1 < 0.75 | adt_CD8 < 1)
}

# Integration for ADT
bri.anchors <- FindIntegrationAnchors(object.list=bri.integrated.pbmc.list, dims=1:12)
bri.integrated.pbmc <- IntegrateData(anchorset = bri.anchors, new.assay.name = "integratedADT", dims=1:12, k.weight = 50)


### Skin: #######

# Integration for RNA 
bri.anchors <- FindIntegrationAnchors(object.list=bri.list.skin, dims=1:30)
bri.integrated.skin <- IntegrateData(anchorset = bri.anchors, new.assay.name = "integratedRNA", dims=1:30, k.weight=30) #k.weight 100 by default, lower it if any sample has too few cells

# split objects (which are integratedRNA) to RE-integrate on ADT
bri.integrated.skin.list <- SplitObject(bri.integrated.skin, split.by = "hash.ID")

# remove cd4 and cd8 doublets
for (i in 1:length(bri.integrated.skin.list)) {
  DefaultAssay(bri.integrated.skin.list[[i]]) <- "ADT"
  bri.integrated.skin.list[[i]] <- NormalizeData(bri.integrated.skin.list[[i]], normalization.method = "CLR")
  bri.integrated.skin.list[[i]] <- subset(bri.integrated.skin.list[[i]], subset = adt_CD4.1 < 0.75 | adt_CD8 < 1)
}

# Integration for ADT
bri.anchors <- FindIntegrationAnchors(object.list=bri.integrated.skin.list, dims=1:12)
bri.integrated.skin <- IntegrateData(anchorset = bri.anchors, new.assay.name = "integratedADT", dims=1:12, k.weight = 50)

####
## MultiModal Analysis for RNA and ADT: #######
####

### PBMC: #######
DefaultAssay(bri.integrated.pbmc) <- "integratedRNA"
bri.integrated.pbmc <- ScaleData(bri.integrated.pbmc) %>% RunPCA(features=rownames(bri.integrated.pbmc))

# Examine and visualize PCA results
print(bri.integrated.pbmc@reductions$pca)
VizDimLoadings(bri.integrated.pbmc, dims = 1:20, projected = FALSE, nfeatures = 15)
ElbowPlot(bri.integrated.pbmc, ndims = 30)

DefaultAssay(bri.integrated.pbmc) <- 'integratedADT'
VariableFeatures(bri.integrated.pbmc) <- rownames(bri.integrated.pbmc[["ADT"]])
bri.integrated.pbmc <- ScaleData(bri.integrated.pbmc) %>% RunPCA(reduction.name = 'apca')

# wnn analysis 
bri.integrated.pbmc <- FindMultiModalNeighbors(bri.integrated.pbmc, reduction.list = list("pca", "apca"),
                                               dims.list = list(1:30, 1:12), modality.weight.name = "RNA.weight")

# UMAP for WNN
bri.integrated.pbmc <- RunUMAP(bri.integrated.pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p.pbmc <- DimPlot(bri.integrated.pbmc, reduction = "wnn.umap", group.by = "hash.ID")
p.pbmc

### Skin: #######

# DefaultAssay(bri.integrated) <- 'RNA'
DefaultAssay(bri.integrated.skin) <- "integratedRNA"
bri.integrated.skin <- NormalizeData(bri.integrated.skin) %>% FindVariableFeatures()
bri.integrated.skin <- ScaleData(bri.integrated.skin) %>% RunPCA(features=rownames(bri.integrated.skin))

# Examine and visualize PCA results
print(bri.integrated.skin@reductions$pca)
VizDimLoadings(bri.integrated.skin, dims = 1:20, projected = FALSE, nfeatures = 15) #default is 30, but too much text overlap
ElbowPlot(bri.integrated.skin, ndims = 30)

DefaultAssay(bri.integrated.skin) <- 'integratedADT'
VariableFeatures(bri.integrated.skin) <- rownames(bri.integrated.skin[["ADT"]])
bri.integrated.skin <- ScaleData(bri.integrated.skin) %>% RunPCA(reduction.name = 'apca')

# wnn analysis 
bri.integrated.skin <- FindMultiModalNeighbors(bri.integrated.skin, reduction.list = list("pca", "apca"),
                                               dims.list = list(1:30, 1:12), modality.weight.name = "RNA.weight")

# UMAP for WNN
bri.integrated.skin <- RunUMAP(bri.integrated.skin, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p.skin <- DimPlot(bri.integrated.skin, reduction = "wnn.umap", group.by = "hash.ID")
p.skin


####
## Match barcodes of combined integrated data vs the new independent skin and blood clusters #######
####

oldintegrateddata <- readRDS("final-multimodal-sjs-clustered-tcr-annotated.rds")

# add cytotoxicity
cytotoxic_markers <- c("NKG7","GZMA", "GZMB", "GNLY", "PRF1")
rna_counts <- bri.integrated.skin@assays$integratedRNA@data
agg_cytotoxic <- colMeans(rna_counts[cytotoxic_markers,])
bri.integrated.skin$cytotoxic_expression <- agg_cytotoxic
rna_counts <- bri.integrated.pbmc@assays$integratedRNA@data
agg_cytotoxic <- colMeans(rna_counts[cytotoxic_markers,])
bri.integrated.pbmc$cytotoxic_expression <- agg_cytotoxic

# correct barcodes and hash names
temp <- oldintegrateddata$hash.ID
temp <- gsub("-","",temp)
oldintegrateddata$hash.ID <- temp

temp <- bri.integrated.skin$hash.ID
temp <- gsub("-","",temp)
bri.integrated.skin$hash.ID <- temp

temp <- bri.integrated.pbmc$hash.ID
temp <- gsub("-","",temp)
bri.integrated.pbmc$hash.ID <- temp

temp <- sapply(colnames(bri.integrated.skin), function(x) strsplit(x, split = "-")[[1]][1], USE.NAMES = FALSE)
temp <- paste(temp, bri.integrated.skin$hash.ID, sep = "-")
bri.integrated.skin <- RenameCells(bri.integrated.skin, old.names = colnames(bri.integrated.skin), new.names = temp)

temp <- sapply(colnames(bri.integrated.pbmc), function(x) strsplit(x, split = "-")[[1]][1], USE.NAMES = FALSE)
temp <- paste(temp, bri.integrated.pbmc$hash.ID, sep = "-")
bri.integrated.pbmc <- RenameCells(bri.integrated.pbmc, old.names = colnames(bri.integrated.pbmc), new.names = temp)

temp <- sapply(colnames(oldintegrateddata), function(x) strsplit(x, split = "-")[[1]][1], USE.NAMES = FALSE)
temp <- paste(temp, oldintegrateddata$hash.ID, sep = "-")
oldintegrateddata <- RenameCells(oldintegrateddata, old.names = colnames(oldintegrateddata), new.names = temp)

bri.integrated.skin$CellType <- oldintegrateddata$CellType[match(colnames(bri.integrated.skin), colnames(oldintegrateddata))]
bri.integrated.pbmc$CellType <- oldintegrateddata$CellType[match(colnames(bri.integrated.pbmc), colnames(oldintegrateddata))]

bri.integrated.skin.nonna <- subset(bri.integrated.skin, subset = CellType != "NA")
bri.integrated.skin.nonna <- subset(bri.integrated.skin.nonna, subset = CellType != "Dead")
bri.integrated.skin.nonna <- subset(bri.integrated.skin.nonna, subset = CellType != "Proliferating")
Idents(bri.integrated.skin.nonna) <- "CellType"
colors <- hue_pal()(length(unique(bri.integrated.skin.nonna$CellType)))
color_ind1 <- 1:12
color_ind2 <- c(13:23,1)
color_ind <- c(color_ind1, color_ind2)
color_ind <- color_ind[order(c(seq_along(color_ind1)*2-1,seq_along(color_ind2)*2))]
show_col(colors[color_ind])
DimPlot(bri.integrated.skin.nonna, reduction = 'wnn.umap', label = TRUE, label.size = 3, pt.size = 2, cols = colors)
FeaturePlot(bri.integrated.skin.nonna, reduction = 'wnn.umap', features = "CD45RA")

bri.integrated.pbmc.nonna <- subset(bri.integrated.pbmc, subset = CellType != "NA")
bri.integrated.pbmc.nonna <- subset(bri.integrated.pbmc.nonna, subset = CellType != "Dead")
bri.integrated.pbmc.nonna <- subset(bri.integrated.pbmc.nonna, subset = CellType != "Proliferating")
Idents(bri.integrated.pbmc.nonna) <- "CellType"
colors <- hue_pal()(length(unique(bri.integrated.skin.nonna$CellType)))
color_ind <- c(rbind(1:11,12:22))
show_col(colors[color_ind])
DimPlot(bri.integrated.pbmc.nonna, reduction = 'wnn.umap', label = TRUE, label.size = 3, pt.size = 2, cols = colors[color_ind]) + NoLegend()


####
## Clustering and compare #######
####

Idents(bri.integrated.pbmc) <- "CellType"
DimPlot(bri.integrated.pbmc, reduction = 'wnn.umap', label = TRUE, label.size = 5, pt.size = 2, cols = colors) +  guides(color = guide_legend(ncol=1))
Idents(bri.integrated.skin) <- "CellType"
DimPlot(bri.integrated.skin, reduction = 'wnn.umap', label = TRUE, label.size = 5, pt.size = 2, cols = colors) + guides(color = guide_legend(ncol=1))

### skin ####
bri.integrated.skin <- FindClusters(bri.integrated.skin, graph.name = "wsnn", algorithm = 3,
                                    resolution = c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0), verbose = FALSE)

DimPlot(bri.integrated.skin, group.by = "hash.ID")
DimPlot(bri.integrated.skin, group.by = "seurat_clusters")

# Visualization for resolutions c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
Idents(bri.integrated.skin) <- "wsnn_res.0.4"
p1 <- DimPlot(bri.integrated.skin, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.4")
Idents(bri.integrated.skin) <- "wsnn_res.0.6"
p2 <- DimPlot(bri.integrated.skin, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.6")
Idents(bri.integrated.skin) <- "wsnn_res.0.8"
p3 <- DimPlot(bri.integrated.skin, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.8")
Idents(bri.integrated.skin) <- "wsnn_res.1"
p4 <- DimPlot(bri.integrated.skin, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.0")
Idents(bri.integrated.skin) <- "wsnn_res.1.2"
p5 <- DimPlot(bri.integrated.skin, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.2")
Idents(bri.integrated.skin) <- "wsnn_res.1.4"
p6 <- DimPlot(bri.integrated.skin, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.4")
Idents(bri.integrated.skin) <- "wsnn_res.1.6"
p7 <- DimPlot(bri.integrated.skin, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.6")
Idents(bri.integrated.skin) <- "wsnn_res.1.8"
p8 <- DimPlot(bri.integrated.skin, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.8")
p1 + p2 + p3 + p4
p5 + p6 + p7 + p8

bri.integrated.pbmc <- FindClusters(bri.integrated.pbmc, graph.name = "wsnn", algorithm = 3,
                                    resolution = c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0), verbose = FALSE)

DimPlot(bri.integrated.pbmc, group.by = "hash.ID")
DimPlot(bri.integrated.pbmc, group.by = "seurat_clusters")

# Visualization for resolutions c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
Idents(bri.integrated.pbmc) <- "wsnn_res.0.4"
p1 <- DimPlot(bri.integrated.pbmc, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.4")
Idents(bri.integrated.pbmc) <- "wsnn_res.0.6"
p2 <- DimPlot(bri.integrated.pbmc, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.6")
Idents(bri.integrated.pbmc) <- "wsnn_res.0.8"
p3 <- DimPlot(bri.integrated.pbmc, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.8")
Idents(bri.integrated.pbmc) <- "wsnn_res.1"
p4 <- DimPlot(bri.integrated.pbmc, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.0")
Idents(bri.integrated.pbmc) <- "wsnn_res.1.2"
p5 <- DimPlot(bri.integrated.pbmc, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.2")
Idents(bri.integrated.pbmc) <- "wsnn_res.1.4"
p6 <- DimPlot(bri.integrated.pbmc, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.4")
Idents(bri.integrated.pbmc) <- "wsnn_res.1.6"
p7 <- DimPlot(bri.integrated.pbmc, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.6")
Idents(bri.integrated.pbmc) <- "wsnn_res.1.8"
p8 <- DimPlot(bri.integrated.pbmc, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.8")
p1 + p2 + p3 + p4
p5 + p6 + p7 + p8

####
## Relabeling of annotations with markers #######
####

bri.integrated.pbmc <- bri.integrated.pbmc[,bri.integrated.pbmc$CellType != "NA"]
bri.integrated.pbmc <- bri.integrated.pbmc[,bri.integrated.pbmc$CellType != "Dead"]

bri.integrated.skin <- bri.integrated.skin[,bri.integrated.skin$CellType != "NA"]
bri.integrated.skin <- bri.integrated.skin[,bri.integrated.skin$CellType != "Dead"]

# Rename PBMC
Idents(bri.integrated.pbmc) <- "CellType"
bri.integrated.pbmc <- RenameIdents(object = bri.integrated.pbmc,
                                    "CD8+ CD56+ T cells" = "CD8+ CD56+ T cells",
                                    "CD8+ CD103- CD69+ TRM (nf)" = "CD8+ CD103- CD69+ (nf)",
                                    "CD4+ CD103+ CD69+ TRM" = "CD4+ CD103+",
                                    "Proliferating" = "Proliferating",
                                    "CD4+ TCM" = "CD4+ TCM",
                                    "CD4+ Treg 1" = "CD4+ Treg 1",
                                    "CD4+ Naive" = "CD4+ Naive",
                                    "CD4+ T Effector" = "CD4+ T Effector",
                                    "CD8+ Naive" = "CD8+ Naive",
                                    "CD4+ CD45RA- CD62L-low" = "CD4+ CD45RA- CD62L-low",
                                    "CD4+ Treg 2" = "CD4+ Treg 2",
                                    "CD8+ CD103+ TEMRA" = "CD8+ CD103+ TEMRA",
                                    "CD8+ CD103- TEMRA" = "CD8+ CD103- TEMRA",
                                    "CD4+ CD45RA+ CD62L-low" = "CD4+ CD45RA+ CD62L-low",
                                    "CD8+ CD103- CD69+ TRM" = "CD8+ CD103- CD69+",
                                    "CD8+ TMM" = "CD8+ TMM",
                                    "CD8+ TCM" = "CD8+ TCM",
                                    "CD8+ TEM" = "CD8+ TEM",
                                    "CD4+ CD103- CD69+ TRM" = "CD4+ CD103- CD69+",
                                    "CD8+ T Effector" = "CD8+ T Effector",
                                    "CD8+ CD103+ CD69+ TRM" = "CD8+ CD103+ CD69+",
                                    "GD T cells" = "GD T cells")
bri.integrated.pbmc$CellTypePBMCmarkers <- Idents(bri.integrated.pbmc)

# Rename Skin
Idents(bri.integrated.skin) <- "CellType"
bri.integrated.skin <- RenameIdents(object = bri.integrated.skin,
                                    "CD8+ CD56+ T cells" = "CD8+ CD56+ T cells",
                                    "CD8+ CD103- CD69+ TRM (nf)" = "CD8+ CD103- CD69+ TRM (nf)",
                                    "CD4+ CD103+ CD69+ TRM" = "CD4+ CD103+ CD69+ TRM",
                                    "Proliferating" = "Proliferating",
                                    "CD4+ TCM" = "CD4 CD62L- IL7RA+ RO+",
                                    "CD4+ Treg 1" = "CD4+ Treg 1",
                                    "CD4+ Naive" = "CD4+ CD62L+ RA-",
                                    "CD4+ T Effector" = "CD4+ T Effector",
                                    "CD8+ Naive" = "CD8+ IL7RA+ RA+",
                                    "CD4+ CD45RA- CD62L-low" = "CD4+ IL7RA+ CD62L-low",
                                    "CD4+ Treg 2" = "CD4+ Treg 2",
                                    "CD8+ CD103+ TEMRA" = "CD8+ CD103+ TEMRA",
                                    "CD8+ CD103- TEMRA" = "CD8+ CD103- TEMRA",
                                    "CD4+ CD45RA+ CD62L-low" = "CD4+ IL7RA+ CD62L-low",
                                    "CD8+ CD103- CD69+ TRM" = "CD8+ CD103- CD69+",
                                    "CD8+ TMM" = "CD8+ CD62L- IL7RA+ RA-",
                                    "CD8+ TCM" = "CD8+ CD62L+ IL7RA+ RO+",
                                    "CD8+ TEM" = "CD8+ TEM",
                                    "CD4+ CD103- CD69+ TRM" = "CD4+ CD103- CD69+ TRM",
                                    "CD8+ T Effector" = "CD8+ T Effector",
                                    "CD8+ CD103+ CD69+ TRM" = "CD8+ CD103+ CD69+ TRM",
                                    "GD T cells" = "GD T cells")
bri.integrated.skin$CellTypeSkinmarkers <- Idents(bri.integrated.skin)


####
## Match Clusters with annotations #######
####

### Skin ####
# match columns and rows
col_fun = colorRamp2(c(0, 4), c("white", "red"))
compare_table <- unclass(table(bri.integrated.skin$CellTypeSkin,bri.integrated.skin$wsnn_res.1.6))
compare_table_rev <- max(compare_table) - compare_table
# compare_table_ind <- apply(compare_table, 1, function(x) colnames(compare_table)[which.max(x)])
# compare_table <- compare_table[,c(unique(compare_table_ind),colnames(compare_table)[!colnames(compare_table) %in% unique(compare_table_ind)])]
compare_table_alg <- HungarianSolver(compare_table_rev)
compare_table_reorder<- compare_table[compare_table_alg$pairs[,1], compare_table_alg$pairs[,2]]
compare_table_empty <- compare_table[which(compare_table_alg$pairs[,2] == 0),]
compare_table_empty_rev <- max(compare_table_empty) - compare_table_empty
compare_table_empty_alg <- HungarianSolver(compare_table_empty_rev)
compare_table_alg$pairs[which(compare_table_alg$pairs[,2] == 0),2] <- compare_table_empty_alg$pairs[,2]
compare_table_reorder<- compare_table[compare_table_alg$pairs[,1], unique(compare_table_alg$pairs[,2])]
colnames_table <- colnames(compare_table_reorder)
compare_table_reorder <- t(apply(compare_table_reorder,1,scale))
colnames(compare_table_reorder) <- colnames_table
Heatmap(unclass(compare_table_reorder), 
        col = col_fun,
        cluster_columns = FALSE, cluster_rows = FALSE,
        show_row_dend  = FALSE, column_title_gp = grid::gpar(fontsize = 16),
        show_column_dend  = FALSE, column_title = "Skin Clusters vs Cell Types", 
        heatmap_legend_param = list(title = "# of Cells \n (Scaled)"))

### PBMC ####

# remove high mt 
bri.integrated.pbmc_nonmt <- bri.integrated.pbmc[,bri.integrated.pbmc$wsnn_res.1.6 != "11"]

# match columns and rows
col_fun = colorRamp2(c(0, 4), c("white", "red"))
compare_table <- unclass(table(bri.integrated.pbmc_nonmt$CellTypePBMC,bri.integrated.pbmc_nonmt$wsnn_res.1.4))
compare_table_rev <- max(compare_table) - compare_table
# compare_table_ind <- apply(compare_table, 1, function(x) colnames(compare_table)[which.max(x)])
# compare_table <- compare_table[,c(unique(compare_table_ind),colnames(compare_table)[!colnames(compare_table) %in% unique(compare_table_ind)])]
compare_table_alg <- HungarianSolver(compare_table_rev)
compare_table_reorder<- compare_table[compare_table_alg$pairs[,1], compare_table_alg$pairs[,2]]
compare_table_empty <- compare_table[which(compare_table_alg$pairs[,2] == 0),]
compare_table_empty_rev <- max(compare_table_empty) - compare_table_empty
compare_table_empty_alg <- HungarianSolver(compare_table_empty_rev)
compare_table_alg$pairs[which(compare_table_alg$pairs[,2] == 0),2] <- compare_table_empty_alg$pairs[,2]
compare_table_reorder<- compare_table[compare_table_alg$pairs[,1], unique(compare_table_alg$pairs[,2])]
colnames_table <- colnames(compare_table_reorder)
compare_table_reorder <- t(apply(compare_table_reorder,1,scale))
colnames(compare_table_reorder) <- colnames_table
Heatmap(unclass(compare_table_reorder), 
        col = col_fun,
        cluster_columns = FALSE, cluster_rows = FALSE,
        show_row_dend  = FALSE, column_title_gp = grid::gpar(fontsize = 16),
        show_column_dend  = FALSE, column_title = "PBMC Clusters vs Integrated Clusters", 
        heatmap_legend_param = list(title = "# of Cells \n (Scaled)"))


# save RDS files for independently clustered skin and blood
saveRDS(bri.integrated.skin, "8-multimodal-sjs-skin-only.rds")
saveRDS(bri.integrated.pbmc_nonmt, "8-multimodal-sjs-pbmc-only.rds")


















####
## Compare Markers acros skin and blood #######
####

# compare blood markers
skin.bloodmarkers <- VlnPlot(bri.integrated.skin, features = c("adt_CD62L", "adt_IL-7Ra", "adt_CCR7.1", "adt_CD45RA", "adt_CD45RO", "cytotoxic_expression"), ncol = 6, pt.size = 0.5, combine = FALSE) 
pbmc.bloodmarkers <- VlnPlot(bri.integrated.pbmc, features = c("adt_CD62L", "adt_IL-7Ra", "adt_CCR7.1", "adt_CD45RA", "adt_CD45RO","cytotoxic_expression"), ncol = 6, pt.size = 0.5, combine = FALSE) 
bloodmarkers <- c(skin.bloodmarkers, pbmc.bloodmarkers)
for(i in 1:length(bloodmarkers)){
  bloodmarkers[[i]] <- bloodmarkers[[i]] + NoLegend() + theme(axis.text.x = element_text(size = 5)) 
}
CombinePlots(bloodmarkers, ncol = 6)

all.bloodmarkers <- VlnPlot(oldintegrateddata, features = c("adt_CD62L", "adt_IL-7Ra", "CCR7.1", "adt_CD45RA", "adt_CD45RO"), ncol = 5, pt.size = 0.5, group.by = "CellType", combine = FALSE)
for(i in 1:length(all.bloodmarkers)){
  all.bloodmarkers[[i]] <- all.bloodmarkers[[i]] + NoLegend() + theme(axis.text.x = element_text(size = 5)) 
}
CombinePlots(all.bloodmarkers, ncol = 5)

# compare skin markers
skin.skinmarkers <- VlnPlot(bri.integrated.skin, features = c("adt_CD103", "adt_IL-7Ra", "adt_CD69.1", "adt_CD45RA", "adt_CD45RO", "cytotoxic_expression"), ncol = 6, pt.size = 0.5, combine = FALSE) 
pbmc.skinmarkers <- VlnPlot(bri.integrated.pbmc, features = c("adt_CD103", "adt_IL-7Ra", "adt_CD69.1", "adt_CD45RA", "adt_CD45RO", "cytotoxic_expression"), ncol = 6, pt.size = 0.5, combine = FALSE) 
bloodmarkers <- c(skin.skinmarkers, pbmc.skinmarkers)
for(i in 1:length(bloodmarkers)){
  bloodmarkers[[i]] <- bloodmarkers[[i]] + NoLegend() + theme(axis.text.x = element_text(size = 5)) 
}
CombinePlots(bloodmarkers, ncol = 6)

# CD56 and gd markers
skin.skinmarkers <- VlnPlot(bri.integrated.skin, features = c("adt_CD45RA", "adt_CD45RO", "adt_CD4.1", "adt_CD8", "adt_CD56","rna_GNLY","rna_TRDV2"), ncol = 7, pt.size = 0.5, combine = FALSE) 
pbmc.skinmarkers <- VlnPlot(bri.integrated.pbmc, features = c("adt_CD45RA", "adt_CD45RO","adt_CD4.1", "adt_CD8", "adt_CD56", "GNLY","rna_TRDV2"), ncol = 7, pt.size = 0.5, combine = FALSE) 
bloodmarkers <- c(skin.skinmarkers, pbmc.skinmarkers)
for(i in 1:length(bloodmarkers)){
  bloodmarkers[[i]] <- bloodmarkers[[i]] + NoLegend() + theme(axis.text.x = element_text(size = 5)) 
}
CombinePlots(bloodmarkers, ncol = 7)

# cytotoxic
cytotoxic_markers <- c("rna_NKG7","rna_GZMA", "rna_GZMB", "rna_GNLY", "rna_PRF1")
DefaultAssay(bri.integrated.skin) <- "RNA"
DefaultAssay(bri.integrated.pbmc) <- "RNA"
skin.skinmarkers <- VlnPlot(bri.integrated.skin, features = c(cytotoxic_markers,"cytotoxic_expression"), ncol = 6, pt.size = 0.5, combine = FALSE) 
pbmc.skinmarkers <- VlnPlot(bri.integrated.pbmc, features = c(cytotoxic_markers,"cytotoxic_expression"), ncol = 6, pt.size = 0.5, combine = FALSE) 
bloodmarkers <- c(skin.skinmarkers, pbmc.skinmarkers)
for(i in 1:length(bloodmarkers)){
  bloodmarkers[[i]] <- bloodmarkers[[i]] + NoLegend() + theme(axis.text.x = element_text(size = 5)) 
}
CombinePlots(bloodmarkers, ncol = 6)

### Heatmap of Expression #### 

#### PBMC #### 

##### full template ####
features_all <- c("CD8","CD8A","CD4.1","CD4", "CD45RA","CD45RO",
                  "CD69.1", "CD69", "CD103", "ITGAE", "KLF2", "S1PR1", "CD62L", "SELL", "CCR7.1", "CCR7", "IL-7Ra","IL7R", 
                  "CD56",'NCAM1', "CD335", "NCR1", "IFNG","GNLY", 'GZMA', 'GZMB', 'PRF1', 'NKG7', "CX3CR1",
                  "TRDV2", "TRGV10",
                  "FAS.1","FAS","LAMP-1","LAMP1", 'IL10', 'FOXP3', 'CTLA4', 'IL2RA', "GATA3","IL17A","KLRG1","B3GAT1",
                  "TUBA1B","STMN1")

features_RNA <- c('CD8A','CD4','NCAM1','ITGAE','SELL','CD69',"S1PR1", "TRDV2", "TRGV10", "KLF2",
                  "IFNG","GNLY", 'GZMA', 'GZMB', 'PRF1', 'NKG7', "FAS",
                  'IL10', 'FOXP3', 'CTLA4', 'IL2RA',"GATA3","IL17A","KLRG1","B3GAT1","CX3CR1", "CCR7","IL7R", "NCR1", "LAMP1","TUBA1B","STMN1")

features_ADT <- c("CD8","CD4.1","CD45RA","CD45RO","IL-7Ra","CD103",'CD62L', "CD69.1","CD56","FAS.1","CCR7.1", "CD335","LAMP-1")

label_order <- c("CD8+ Naive","CD8+ TCM","CD8+ TMM","CD8+ CD103- CD69+ (nf)","CD8+ CD103- TEMRA","CD8+ CD103+ TEMRA",
                 "CD8+ T Effector", "CD8+ TEM", "CD8+ CD103+ CD69+", "CD8+ CD103- CD69+", "CD8+ CD56+ T cells", "GD T cells", 
                 "CD4+ Naive", "CD4+ CD45RA+ CD62L-low", "CD4+ CD45RA- CD62L-low", "CD4+ TCM", 
                 "CD4+ CD103- CD69+", "CD4+ CD103+", "CD4+ T Effector", "CD4+ Treg 1", "CD4+ Treg 2", "Proliferating")

mean_ADT <- AverageExpression(bri.integrated.pbmc, assays = "ADT", features = features_ADT)
mean_ADT_scaled <- apply(mean_ADT$ADT, 1, scale)
rownames(mean_ADT_scaled) <- colnames(mean_ADT$ADT)
mean_ADT_scaled <- as.data.frame(mean_ADT_scaled)

mean_RNA <- AverageExpression(bri.integrated.pbmc, assays = "RNA", features = features_RNA)
mean_RNA_scaled <- apply(mean_RNA$RNA, 1, scale)
rownames(mean_RNA_scaled) <- colnames(mean_RNA$RNA)
mean_RNA_scaled <- as.data.frame(mean_RNA_scaled)

mean_all <- cbind(as.data.frame(mean_ADT_scaled), as.data.frame(mean_RNA_scaled))
mean_all <- mean_all[,features_all]
# mean_all <- mean_all[order(as.numeric(rownames(mean_all))),]
Heatmap(mean_all[label_order,], show_row_dend = FALSE, show_column_dend = FALSE,
        row_names_side = "left", column_names_side = "top",
        column_names_rot = 45, cluster_columns = FALSE, cluster_rows = FALSE,
        column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 10),
        column_split = rep(c("A","B","C","D","E","F"),c(6,12,11,2,12,2)), #change numbers to match number of features,
        heatmap_legend_param = list(title = ""), column_title = NULL)

##### for figure ####
Idents(bri.integrated.pbmc) <- "CellTypePBMC"
features_all <- c("CD8","CD8A","CD4.1","CD4", "CD45RA","CD45RO",
                  "CD69.1", "CD69", "CD103", "ITGAE", "KLF2", "S1PR1", "CD62L", "SELL", "CCR7.1", "CCR7", "IL-7Ra","IL7R", 
                  "CD56",'NCAM1',"IFNG","GNLY", 'GZMA', 'GZMB', 'PRF1', 'NKG7', "CX3CR1",
                  "TRDV2", "TRGV10",
                  "FAS.1","FAS", 'IL10', 'FOXP3', 'CTLA4', 'IL2RA',
                  "TUBA1B","STMN1")

features_RNA <- c('CD8A','CD4','NCAM1','ITGAE','SELL','CD69',"S1PR1", "TRDV2", "TRGV10", "KLF2",
                  "IFNG","GNLY", 'GZMA', 'GZMB', 'PRF1', 'NKG7', "FAS",
                  'IL10', 'FOXP3', 'CTLA4', 'IL2RA',"GATA3","IL17A","KLRG1","B3GAT1","CX3CR1", "CCR7","IL7R", "NCR1", "LAMP1","TUBA1B","STMN1")

features_ADT <- c("CD8","CD4.1","CD45RA","CD45RO","IL-7Ra","CD103",'CD62L', "CD69.1","CD56","FAS.1","CCR7.1", "CD335","LAMP-1")

label_order <- c("CD8+ Naive","CD8+ TCM","CD8+ TMM","CD8+ CD103- CD69+ TRM (nf)","CD8+ CD103- TEMRA","CD8+ CD103+ TEMRA",
                 "CD8+ T Effector", "CD8+ TEM", "CD8+ CD103+ TRM", "CD8+ CD103- CD69+ TRM", "CD8+ CD56+ T cells", "GD T cells", 
                 "CD4+ Naive", "CD4+ CD45RA+ CD62L-low", "CD4+ CD45RA- CD62L-low", "CD4+ TCM", 
                 "CD4+ CD103- CD69+ TRM", "CD4+ CD103+ TRM", "CD4+ T Effector", "CD4+ Treg 1", "CD4+ Treg 2", "Proliferating")

mean_ADT <- AverageExpression(bri.integrated.pbmc, assays = "ADT", features = features_ADT)
mean_ADT_scaled <- apply(mean_ADT$ADT, 1, scale)
rownames(mean_ADT_scaled) <- colnames(mean_ADT$ADT)
mean_ADT_scaled <- as.data.frame(mean_ADT_scaled)

mean_RNA <- AverageExpression(bri.integrated.pbmc, assays = "RNA", features = features_RNA)
mean_RNA_scaled <- apply(mean_RNA$RNA, 1, scale)
rownames(mean_RNA_scaled) <- colnames(mean_RNA$RNA)
mean_RNA_scaled <- as.data.frame(mean_RNA_scaled)

mean_all <- cbind(as.data.frame(mean_ADT_scaled), as.data.frame(mean_RNA_scaled))
mean_all <- mean_all[,features_all]
colnames(mean_all) <- features_all <- c("CD8 (ADT)","CD8A (RNA)","CD4.1 (ADT)","CD4 (RNA)", "CD45RA (ADT)","CD45RO (ADT)",
                                        "CD69.1 (ADT)", "CD69 (RNA)", "CD103 (ADT)", "ITGAE (RNA)", "KLF2", "S1PR1", "CD62L (ADT)", "SELL (RNA)", "CCR7.1 (ADT)", "CCR7 (RNA)", "IL-7Ra (ADT)","IL7R (RNA)", 
                                        "CD56 (ADT)",'NCAM1 (RNA)',"IFNG","GNLY", 'GZMA', 'GZMB', 'PRF1', 'NKG7', "CX3CR1",
                                        "TRDV2", "TRGV10",
                                        "FAS.1 (ADT)","FAS (RNA)", 'IL10', 'FOXP3', 'CTLA4', 'IL2RA',
                                        "TUBA1B","STMN1")
# mean_all <- mean_all[order(as.numeric(rownames(mean_all))),]
PBMC.heatmap <- Heatmap(mean_all[label_order,], show_row_dend = FALSE, show_column_dend = FALSE,
        row_names_side = "left", column_names_side = "top",
        column_names_rot = 45, cluster_columns = FALSE, cluster_rows = FALSE,
        column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 10),
        column_split = rep(c("A","B","C","D","E","F"),c(6,12,9,2,6,2)), #change numbers to match number of features
          column_title = "Average Expression of Cell Type Markers (Integrated Analysis) on PBMC Cells",
        heatmap_legend_param = list(title = ""), show_heatmap_legend = FALSE)
PBMC.heatmap <- draw(PBMC.heatmap, padding = unit(c(2, 2, 2, 6), "mm")) # add space for titles

#### Skin #### 

##### for figure ####
Idents(bri.integrated.skin) <- "CellTypeSkin"
features_all <- c("CD8","CD8A","CD4.1","CD4", "CD45RA","CD45RO",
                  "CD69.1", "CD69", "CD103", "ITGAE", "KLF2", "S1PR1", "CD62L", "SELL", "CCR7.1", "CCR7", "IL-7Ra","IL7R", 
                  "CD56",'NCAM1',"IFNG","GNLY", 'GZMA', 'GZMB', 'PRF1', 'NKG7', "CX3CR1",
                  "TRDV2", "TRGV10",
                  "FAS.1","FAS", 'IL10', 'FOXP3', 'CTLA4', 'IL2RA',
                  "TUBA1B","STMN1")

features_RNA <- c('CD8A','CD4','NCAM1','ITGAE','SELL','CD69',"S1PR1", "TRDV2", "TRGV10", "KLF2",
                  "IFNG","GNLY", 'GZMA', 'GZMB', 'PRF1', 'NKG7', "FAS",
                  'IL10', 'FOXP3', 'CTLA4', 'IL2RA',"GATA3","IL17A","KLRG1","B3GAT1","CX3CR1", "CCR7","IL7R", "NCR1", "LAMP1","TUBA1B","STMN1")

features_ADT <- c("CD8","CD4.1","CD45RA","CD45RO","IL-7Ra","CD103",'CD62L', "CD69.1","CD56","FAS.1","CCR7.1", "CD335","LAMP-1")

label_order <- c("CD8+ Naive","CD8+ TCM","CD8+ TMM","CD8+ CD103- CD69+ TRM (nf)","CD8+ CD103- TEMRA","CD8+ CD103+ TEMRA",
                 "CD8+ T Effector", "CD8+ TEM", "CD8+ CD103+ CD69+ TRM", "CD8+ CD103- CD69+ TRM", "CD8+ CD56+ T cells", "GD T cells", 
                 "CD4+ Naive", "CD4+ IL7RA+ CD62L-low", "CD4+ TCM", 
                 "CD4+ CD103- CD69+ TRM", "CD4+ CD103+ CD69+ TRM", "CD4+ T Effector", "CD4+ Treg 1", "CD4+ Treg 2", "Proliferating")

# label_order <- c("CD8+ IL7RA+ RA+","CD8+ CD62L+ IL7RA+ RO+","CD8+ CD62L- IL7RA+ RA-","CD8+ CD103- CD69+ TRM (nf)","CD8+ CD103- TEMRA","CD8+ CD103+ TEMRA",
#                  "CD8+ T Effector", "CD8+ TEM", "CD8+ CD103+ CD69+ TRM", "CD8+ CD103- CD69+", "CD8+ CD56+ T cells", "GD T cells", 
#                  "CD4+ CD62L+ RA-", "CD4+ IL7RA+ CD62L-low","CD4 CD62L- IL7RA+ RO+",
#                  "CD4+ CD103- CD69+ TRM", "CD4+ CD103+ CD69+ TRM", "CD4+ T Effector", "CD4+ Treg 1", "CD4+ Treg 2", "Proliferating")

mean_ADT <- AverageExpression(bri.integrated.skin, assays = "ADT", features = features_ADT)
mean_ADT_scaled <- apply(mean_ADT$ADT, 1, scale)
rownames(mean_ADT_scaled) <- colnames(mean_ADT$ADT)
mean_ADT_scaled <- as.data.frame(mean_ADT_scaled)

mean_RNA <- AverageExpression(bri.integrated.skin, assays = "RNA", features = features_RNA)
mean_RNA_scaled <- apply(mean_RNA$RNA, 1, scale)
rownames(mean_RNA_scaled) <- colnames(mean_RNA$RNA)
mean_RNA_scaled <- as.data.frame(mean_RNA_scaled)

mean_all <- cbind(as.data.frame(mean_ADT_scaled), as.data.frame(mean_RNA_scaled))
mean_all <- mean_all[,features_all]
colnames(mean_all) <- features_all <- c("CD8 (ADT)","CD8A (RNA)","CD4.1 (ADT)","CD4 (RNA)", "CD45RA (ADT)","CD45RO (ADT)",
                                        "CD69.1 (ADT)", "CD69 (RNA)", "CD103 (ADT)", "ITGAE (RNA)", "KLF2", "S1PR1", "CD62L (ADT)", "SELL (RNA)", "CCR7.1 (ADT)", "CCR7 (RNA)", "IL-7Ra (ADT)","IL7R (RNA)", 
                                        "CD56 (ADT)",'NCAM1 (RNA)',"IFNG","GNLY", 'GZMA', 'GZMB', 'PRF1', 'NKG7', "CX3CR1",
                                        "TRDV2", "TRGV10",
                                        "FAS.1 (ADT)","FAS (RNA)", 'IL10', 'FOXP3', 'CTLA4', 'IL2RA',
                                        "TUBA1B","STMN1")
# mean_all <- mean_all[order(as.numeric(rownames(mean_all))),]
Skin.Heatmap <- Heatmap(mean_all[label_order,], show_row_dend = FALSE, show_column_dend = FALSE,
        row_names_side = "left", column_names_side = "top",
        column_names_rot = 45, cluster_columns = FALSE, cluster_rows = FALSE,
        column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 10),
        column_split = rep(c("A","B","C","D","E","F"),c(6,12,9,2,6,2)), #change numbers to match number of features
        column_title = "Average Expression of Cell Type Markers (Integrated Analysis) on Skin Cells",
        heatmap_legend_param = list(title = ""), show_heatmap_legend = FALSE)
Skin.Heatmap <- draw(Skin.Heatmap, padding = unit(c(2, 2, 2, 6), "mm")) # add space for titles

#### Integrated #### 

##### for figure ####

#test to show rna cd8, cd4, ncam1, itgae, to go along with adt cd8, cd4, cd56, cd103 - good correlation avg expr level
### #test to show rna cd8, cd4, ncam1 to go along with adt cd8, cd4, cd56 - good correlation avg expr level
Idents(oldintegrateddata) <- "CellTypePranali"
features_all <- c("CD8","CD8A","CD4.1","CD4", "CD45RA","CD45RO",
                  "CD69.1", "CD69", "CD103", "ITGAE", "KLF2", "S1PR1", "CD62L", "SELL", "CCR7.1", "CCR7", "IL-7Ra","IL7R", 
                  "CD56",'NCAM1',"IFNG","GNLY", 'GZMA', 'GZMB', 'PRF1', 'NKG7', "CX3CR1",
                  "TRDV2", "TRGV10",
                  "FAS.1","FAS", 'IL10', 'FOXP3', 'CTLA4', 'IL2RA',
                  "TUBA1B","STMN1")

features_RNA <- c('CD8A','CD4','NCAM1','ITGAE','SELL','CD69',"S1PR1", "TRDV2", "TRGV10", "KLF2",
                  "IFNG","GNLY", 'GZMA', 'GZMB', 'PRF1', 'NKG7', "FAS",
                  'IL10', 'FOXP3', 'CTLA4', 'IL2RA',"GATA3","IL17A","KLRG1","B3GAT1","CX3CR1", "CCR7","IL7R", "NCR1", "LAMP1","TUBA1B","STMN1")

features_ADT <- c("CD8","CD4.1","CD45RA","CD45RO","IL-7Ra","CD103",'CD62L', "CD69.1","CD56","FAS.1","CCR7.1", "CD335","LAMP-1")

mean_ADT <- AverageExpression(oldintegrateddata, assays = "ADT", features = features_ADT)
mean_ADT_scaled <- apply(mean_ADT$ADT, 1, scale)
rownames(mean_ADT_scaled) <- colnames(mean_ADT$ADT)
mean_ADT_scaled <- as.data.frame(mean_ADT_scaled)

mean_RNA <- AverageExpression(oldintegrateddata, assays = "RNA", features = features_RNA)
mean_RNA_scaled <- apply(mean_RNA$RNA, 1, scale)
rownames(mean_RNA_scaled) <- colnames(mean_RNA$RNA)
mean_RNA_scaled <- as.data.frame(mean_RNA_scaled)

label_order <- c("CD8+ Naive","CD8+ TCM","CD8+ TMM","CD8+ CD103- CD69+ TRM (nf)","CD8+ CD103- TEMRA","CD8+ CD103+ TEMRA",
                 "CD8+ T Effector", "CD8+ TEM", "CD8+ CD103+ CD69+ TRM", "CD8+ CD103- CD69+ TRM", "CD8+ CD56+ T cells", "GD T cells", 
                 "CD4+ Naive", "CD4+ CD45RA+ CD62L-low", "CD4+ CD45RA- CD62L-low", "CD4+ TCM", 
                 "CD4+ CD103- CD69+ TRM", "CD4+ CD103+ CD69+ TRM", "CD4+ T Effector", 
                 "CD4+ Treg 1", "CD4+ Treg 2", "Proliferating")

mean_all <- cbind(as.data.frame(mean_ADT_scaled), as.data.frame(mean_RNA_scaled))
mean_all <- mean_all[,features_all]
colnames(mean_all) <- features_all <- c("CD8 (ADT)","CD8A (RNA)","CD4.1 (ADT)","CD4 (RNA)", "CD45RA (ADT)","CD45RO (ADT)",
                                        "CD69.1 (ADT)", "CD69 (RNA)", "CD103 (ADT)", "ITGAE (RNA)", "KLF2", "S1PR1", "CD62L (ADT)", "SELL (RNA)", "CCR7.1 (ADT)", "CCR7 (RNA)", "IL-7Ra (ADT)","IL7R (RNA)", 
                                        "CD56 (ADT)",'NCAM1 (RNA)',"IFNG","GNLY", 'GZMA', 'GZMB', 'PRF1', 'NKG7', "CX3CR1",
                                        "TRDV2", "TRGV10",
                                        "FAS.1 (ADT)","FAS (RNA)", 'IL10', 'FOXP3', 'CTLA4', 'IL2RA',
                                        "TUBA1B","STMN1")
Integrated.Heatmap <- Heatmap(mean_all[label_order,], show_row_dend = FALSE, show_column_dend = FALSE,
        row_names_side = "left", column_names_side = "top",
        column_names_rot = 45, cluster_columns = FALSE, cluster_rows = FALSE,
        column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 10),
        column_split = rep(c("A","B","C","D","E","F"),c(6,12,9,2,6,2)), #change numbers to match number of features
        column_title = "Average Expression of Cell Type Markers (Integrated Analysis) on PBMC and Skin (All) Cells",
        heatmap_legend_param = list(title = ""), show_heatmap_legend = FALSE)
Integrated.Heatmap <- draw(Integrated.Heatmap, padding = unit(c(2, 2, 2, 6), "mm")) # add space for titles

#### Combine Heatmaps #### 

PBMC.heatmap + Skin.Heatmap + Integrated.Heatmap

### Heatmap of counts #### 

#### PBMC #### 

cell_vs_number <- table(bri.integrated.pbmc$CellTypePranaliPBMC, bri.integrated.pbmc$hash.ID)
colnames(cell_vs_number) <- gsub("-","",colnames(cell_vs_number))
colnames(cell_vs_number) <- gsub("Control","Ctrl",colnames(cell_vs_number))
colnames(cell_vs_number) <- gsub("CtrlBCH","Ctrl",colnames(cell_vs_number))
colnames(cell_vs_number) <- gsub("BCH","Ctrl",colnames(cell_vs_number))
colnames(cell_vs_number) <- gsub("Skin","SKIN",colnames(cell_vs_number))

col_fun = colorRamp2(c(0, 300), c("white","red"))
# label_order <- c("CD8+ Naive","CD8+ TCM","CD8+ TMM","CD8+ CD103- CD69+ (nf)","CD8+ CD103- TEMRA","CD8+ CD103+ TEMRA",
#                  "CD8+ T Effector", "CD8+ TEM", "CD8+ CD103+ CD69+", "CD8+ CD103- CD69+", "CD8+ CD56+ T cells", "GD T cells", 
#                  "CD4+ Naive", "CD4+ CD45RA+ CD62L-low", "CD4+ CD45RA- CD62L-low", "CD4+ TCM", 
#                  "CD4+ CD103- CD69+", "CD4+ CD103+", "CD4+ T Effector", "CD4+ Treg 1", "CD4+ Treg 2", "Proliferating")
label_order <- c("CD8+ Naive","CD8+ TCM","CD8+ TMM","CD8+ CD103- CD69+ TRM (nf)","CD8+ CD103- TEMRA","CD8+ CD103+ TEMRA",
                 "CD8+ T Effector", "CD8+ TEM", "CD8+ CD103+ TRM", "CD8+ CD103- CD69+ TRM", "CD8+ CD56+ T cells", "GD T cells", 
                 "CD4+ Naive", "CD4+ CD45RA+ CD62L-low", "CD4+ CD45RA- CD62L-low", "CD4+ TCM", 
                 "CD4+ CD103- CD69+ TRM", "CD4+ CD103+ TRM", "CD4+ T Effector", "CD4+ Treg 1", "CD4+ Treg 2", "Proliferating")
colnames(cell_vs_number)[grepl("MDE006PBMC", colnames(cell_vs_number))] <- "MDE003PBMC"
colnames(cell_vs_number)[grepl("Ctrl001PBMC", colnames(cell_vs_number))] <- "Ctrl002PBMC"
colnames(cell_vs_number)[grepl("Ctrl026PBMC", colnames(cell_vs_number))] <- "Ctrl001PBMC"
colnames(cell_vs_number)[grepl("Ctrl048PBMC", colnames(cell_vs_number))] <- "Ctrl003PBMC"
cell_vs_number <- cell_vs_number[label_order,]
count_heatmap <- Heatmap(cell_vs_number, show_row_dend = FALSE, show_column_dend = FALSE, col = col_fun,
                         row_names_side = "left", column_names_side = "top", 
                         column_names_rot = 45, cluster_columns = FALSE, cluster_rows = FALSE,
                         heatmap_legend_param = list(title = "# of Cells"), column_split = rep(c("Control", "MDE","SJS"), c(3,3,3)),
                         column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8),
                         column_title_gp = grid::gpar(fontsize = 14),
                         column_title = "Cell Type Abundances in PBMC Samples",
                         cell_fun = function(j, i, x, y, width, height, fill) {
                           grid.text(paste0(sprintf("%.0f", cell_vs_number[i, j])), x, y, gp = gpar(fontsize = 10))
                         })
count_heatmap

##### Test for abundances #### 

# Assign tissue type, either broadly PBMC or Skin
hashIDs <- bri.integrated.pbmc@meta.data$hash.ID
hashIDs[grepl('PBMC',hashIDs)] <- "PBMC"
bri.integrated.pbmc@meta.data$Tissue <- hashIDs

# Assigng tissue type, specifically which type of PBMC or Skin
# care with the names - is case sensitive
hashIDs <- bri.integrated.pbmc@meta.data$hash.ID
hashIDs[grepl("PBMC",hashIDs) & grepl("SJS",hashIDs)] <- "SJS_PBMC"
hashIDs[grepl("PBMC",hashIDs) & grepl("MDE",hashIDs)] <- "MDE_PBMC"
hashIDs[grepl("PBMC",hashIDs) & grepl("BCH",hashIDs)] <- "BCH_PBMC"
hashIDs[grepl("Skin",hashIDs) & grepl("Control",hashIDs)] <- "Control_Skin"
bri.integrated.pbmc@meta.data$phenotypeTissue <- hashIDs

# phenotypetissue
bri.integrated.pbmc$phenotype <- sapply(bri.integrated.pbmc$phenotypeTissue, function(x) return(strsplit(x, split = "_")[[1]][1]))
bri.integrated.pbmc$phenotype <- gsub("BCH","Control", bri.integrated.pbmc$phenotype)

# test 1
bri.integrated.pbmc_test <- bri.integrated.pbmc[,grepl("SJS|MDE", bri.integrated.pbmc$phenotype)]
prop.list <- getTransformedProps(clusters = bri.integrated.pbmc_test$CellTypePranali, sample = bri.integrated.pbmc_test$hash.ID)
grp <- data.frame(cbind(bri.integrated.pbmc_test$phenotype, bri.integrated.pbmc_test$hash.ID)) %>% distinct()
rownames(grp) <- grp$X2
grp <- grp[colnames(prop.list$Proportions),"X1"]
design <- model.matrix(~0+grp)
contrasts <- c(-1,1)
propeller.ttest(prop.list = prop.list, design = design, contrasts = contrasts, trend = FALSE, sort = TRUE, robust=TRUE)

# test 2
bri.integrated.pbmc_test <- bri.integrated.pbmc[,grepl("SJS|MDE", bri.integrated.pbmc$phenotype)]
bri.integrated.pbmc_test <- bri.integrated.pbmc_test[,!bri.integrated.pbmc_test$hash.ID %in% "MDE002PBMC"]
prop.list <- getTransformedProps(clusters = droplevels(bri.integrated.pbmc_test$CellTypePranaliPBMC), sample = bri.integrated.pbmc_test$hash.ID)
grp <- data.frame(cbind(bri.integrated.pbmc_test$phenotype, bri.integrated.pbmc_test$hash.ID)) %>% distinct()
rownames(grp) <- grp$X2
grp <- grp[colnames(prop.list$Proportions),"X1"]
design <- model.matrix(~0+grp)
contrasts <- c(-1,1)
propeller.ttest(prop.list = prop.list, design = design, contrasts = contrasts, trend = FALSE, sort = TRUE, robust=TRUE)

#### Skin #### 

cell_vs_number <- table(bri.integrated.skin$CellTypePranaliSkin, bri.integrated.skin$hash.ID)
colnames(cell_vs_number) <- gsub("-","",colnames(cell_vs_number))
colnames(cell_vs_number) <- gsub("Control","Ctrl",colnames(cell_vs_number))
colnames(cell_vs_number) <- gsub("CtrlBCH","Ctrl",colnames(cell_vs_number))
colnames(cell_vs_number) <- gsub("BCH","Ctrl",colnames(cell_vs_number))
colnames(cell_vs_number) <- gsub("Skin","SKIN",colnames(cell_vs_number))

col_fun = colorRamp2(c(0, 200), c("white","red"))
# label_order <- c("CD8+ IL7RA+ RA+","CD8+ CD62L+ IL7RA+ RO+","CD8+ CD62L- IL7RA+ RA-","CD8+ CD103- CD69+ TRM (nf)","CD8+ CD103- TEMRA","CD8+ CD103+ TEMRA",
#                  "CD8+ T Effector", "CD8+ TEM", "CD8+ CD103+ CD69+ TRM", "CD8+ CD103- CD69+", "CD8+ CD56+ T cells", "GD T cells", 
#                  "CD4+ CD62L+ RA-", "CD4+ IL7RA+ CD62L-low","CD4 CD62L- IL7RA+ RO+",
#                  "CD4+ CD103- CD69+ TRM", "CD4+ CD103+ CD69+ TRM", "CD4+ T Effector", "CD4+ Treg 1", "CD4+ Treg 2", "Proliferating")
label_order <- c("CD8+ Naive","CD8+ TCM","CD8+ TMM","CD8+ CD103- CD69+ TRM (nf)","CD8+ CD103- TEMRA","CD8+ CD103+ TEMRA",
                 "CD8+ T Effector", "CD8+ TEM", "CD8+ CD103+ CD69+ TRM", "CD8+ CD103- CD69+ TRM", "CD8+ CD56+ T cells", "GD T cells", 
                 "CD4+ Naive", "CD4+ IL7RA+ CD62L-low", "CD4+ TCM", 
                 "CD4+ CD103- CD69+ TRM", "CD4+ CD103+ CD69+ TRM", "CD4+ T Effector", "CD4+ Treg 1", "CD4+ Treg 2", "Proliferating")
colnames(cell_vs_number)[grepl("MDE006SKIN", colnames(cell_vs_number))] <- "MDE003SKINC"
colnames(cell_vs_number)[grepl("Ctrl031SKIN", colnames(cell_vs_number))] <- "Ctrl001SKIN"
colnames(cell_vs_number)[grepl("Ctrl039SKIN", colnames(cell_vs_number))] <- "Ctrl002SKIN"
cell_vs_number <- cell_vs_number[label_order,]
count_heatmap <- Heatmap(cell_vs_number, show_row_dend = FALSE, show_column_dend = FALSE, col = col_fun,
                         row_names_side = "left", column_names_side = "top", 
                         column_names_rot = 45, cluster_columns = FALSE, cluster_rows = FALSE,
                         heatmap_legend_param = list(title = "# of Cells"), column_split = rep(c("Control", "MDE","SJS"), c(2,3,3)),
                         column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8),
                         column_title_gp = grid::gpar(fontsize = 14),
                         column_title = "Cell Type Abundances in Skin Samples",
                         cell_fun = function(j, i, x, y, width, height, fill) {
                           grid.text(paste0(sprintf("%.0f", cell_vs_number[i, j])), x, y, gp = gpar(fontsize = 10))
                         })
count_heatmap

##### Test for abundances #### 

# Assign tissue type, either broadly PBMC or Skin
hashIDs <- bri.integrated.skin@meta.data$hash.ID
hashIDs[grepl('Skin|SKIN',hashIDs)] <- "Skin"
bri.integrated.skin@meta.data$Tissue <- hashIDs

# Assigng tissue type, specifically which type of PBMC or Skin
# care with the names - is case sensitive
hashIDs <- bri.integrated.skin@meta.data$hash.ID
hashIDs[grepl("SKIN",hashIDs) & grepl("SJS",hashIDs)] <- "SJS_Skin"
hashIDs[grepl("SKIN",hashIDs) & grepl("MDE",hashIDs)] <- "MDE_Skin"
hashIDs[grepl("Skin",hashIDs) & grepl("MDE",hashIDs)] <- "MDE_Skin"
hashIDs[grepl("Skin",hashIDs) & grepl("Control",hashIDs)] <- "Control_Skin"
bri.integrated.skin@meta.data$phenotypeTissue <- hashIDs

# phenotypetissue
bri.integrated.skin$phenotype <- sapply(bri.integrated.skin$phenotypeTissue, function(x) return(strsplit(x, split = "_")[[1]][1]))
bri.integrated.skin$phenotype <- gsub("BCH","Control", bri.integrated.skin$phenotype)

# test 1
bri.integrated.skin_test <- bri.integrated.skin[,grepl("SJS|MDE", bri.integrated.skin$phenotype)]
prop.list <- getTransformedProps(clusters = bri.integrated.skin_test$CellTypePranali, sample = bri.integrated.skin_test$hash.ID)
grp <- data.frame(cbind(bri.integrated.skin_test$phenotype, bri.integrated.skin_test$hash.ID)) %>% distinct()
rownames(grp) <- grp$X2
grp <- grp[colnames(prop.list$Proportions),"X1"]
design <- model.matrix(~0+grp)
contrasts <- c(-1,1)
propeller.ttest(prop.list = prop.list, design = design, contrasts = contrasts, trend = FALSE, sort = TRUE, robust=TRUE)

# test 2
bri.integrated.pbmc_test <- bri.integrated.pbmc[,grepl("SJS|MDE", bri.integrated.pbmc$phenotype)]
bri.integrated.pbmc_test <- bri.integrated.pbmc_test[,!bri.integrated.pbmc_test$hash.ID %in% "MDE002SKIN"]
prop.list <- getTransformedProps(clusters = droplevels(bri.integrated.pbmc_test$CellTypePranali), sample = bri.integrated.pbmc_test$hash.ID)
grp <- data.frame(cbind(bri.integrated.pbmc_test$phenotype, bri.integrated.pbmc_test$hash.ID)) %>% distinct()
rownames(grp) <- grp$X2
grp <- grp[colnames(prop.list$Proportions),"X1"]
design <- model.matrix(~0+grp)
contrasts <- c(-1,1)
propeller.ttest(prop.list = prop.list, design = design, contrasts = contrasts, trend = FALSE, sort = TRUE, robust=TRUE)




bri.integrated.treg1 <- subset(bri.integrated, subset=CellTypePranali=='CD4+ Treg 1')
bri.integrated.treg1a <- subset(bri.integrated.treg1, subset=is.na(clonotype_id))

bri.integrated.treg1a <- bri.integrated.treg1@meta.data[is.na(bri.integrated.treg1$clonotype_id),]


## Modality weights

# VlnPlot(bri.integrated, features = "RNA.weight", group.by = 'CelltypePranali', sort = TRUE, pt.size = 0.1) +
#   NoLegend()

VlnPlot(bri.integrated, features = "integratedRNA.weight", group.by = 'CellTypePranali', sort = TRUE, pt.size = 0.1) +
  NoLegend()
VlnPlot(bri.integrated, features = "integratedADT.weight", group.by = 'CellTypePranali', sort = TRUE, pt.size = 0.1) +
  NoLegend()

#make heatmap
# display Rows and Columns
# make a dataframe
# row1:RNAnumbers, row2: ADTnumbers, columns are the CellTypePranalis
#dplyr will make col1 are all the cells, col2: celltypepranali, col3: RNA, col4: ADT
#then aggregate all the RNAs per celltypepranali, aggregate all the ADTs per celltypepranali
weightsdf <- bri.integrated@meta.data %>% select(CellTypePranali, integratedRNA.weight, integratedADT.weight) %>%
  group_by(CellTypePranali) %>% summarise(avgRNA.weight = mean(integratedRNA.weight), avgADT.weight = mean(integratedADT.weight)) %>% 
  filter(CellTypePranali != 'Dead')

weightsdf.long <- weightsdf %>% pivot_longer(cols=c('avgRNA.weight','avgADT.weight'), names_to='modality', values_to='weight') %>%
  spread(CellTypePranali,weight)
df.rowname <-c('ADTweight','RNAweight')
weightsdf.long <- weightsdf.long %>% select(-modality)
# rownames(weightsdf.long) <- c('ADTweight','RNAweight')

col_fun <- colorRamp2(c(0,1),c('white','blue'))
Heatmap(weightsdf.long, name = "weight",
        cluster_rows = FALSE, cluster_columns = FALSE, 
        right_annotation = rowAnnotation(text=anno_text(df.rowname)),
        col=col_fun,
        cell_fun = function(j,i,x,y,width,height,fill) {
          grid.text(sprintf("%.2f", weightsdf.long[i, j]), x, y, gp = gpar(fontsize = 10))
        })

#put the avgadtweight into the Weights col
# celltypepranali   modality  weight
# ...               RNA; ADT  score


