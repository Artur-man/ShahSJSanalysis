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

# change hash names
temp <- oldintegrateddata$hash.ID
temp <- gsub("-","",temp)
oldintegrateddata$hash.ID <- temp

temp <- bri.integrated.skin$hash.ID
temp <- gsub("-","",temp)
bri.integrated.skin$hash.ID <- temp

temp <- bri.integrated.pbmc$hash.ID
temp <- gsub("-","",temp)
bri.integrated.pbmc$hash.ID <- temp

#change cell names by pasting hash ID
temp <- sapply(colnames(bri.integrated.skin), function(x) strsplit(x, split = "-")[[1]][1], USE.NAMES = FALSE)
temp <- paste(temp, bri.integrated.skin$hash.ID, sep = "-")
bri.integrated.skin <- RenameCells(bri.integrated.skin, old.names = colnames(bri.integrated.skin), new.names = temp)

temp <- sapply(colnames(bri.integrated.pbmc), function(x) strsplit(x, split = "-")[[1]][1], USE.NAMES = FALSE)
temp <- paste(temp, bri.integrated.pbmc$hash.ID, sep = "-")
bri.integrated.pbmc <- RenameCells(bri.integrated.pbmc, old.names = colnames(bri.integrated.pbmc), new.names = temp)

temp <- sapply(colnames(oldintegrateddata), function(x) strsplit(x, split = "-")[[1]][1], USE.NAMES = FALSE)
temp <- paste(temp, oldintegrateddata$hash.ID, sep = "-")
oldintegrateddata <- RenameCells(oldintegrateddata, old.names = colnames(oldintegrateddata), new.names = temp)

# match cell names in tissue specific objects with the main integrated object, and transfer cell name names
bri.integrated.skin$CellType <- oldintegrateddata$CellType[match(colnames(bri.integrated.skin), colnames(oldintegrateddata))]
bri.integrated.pbmc$CellType <- oldintegrateddata$CellType[match(colnames(bri.integrated.pbmc), colnames(oldintegrateddata))]


####
## Clustering and compare #######

bri.integrated.skin <- FindClusters(bri.integrated.skin, graph.name = "wsnn", algorithm = 3,
                                    resolution = c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0), verbose = FALSE)

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


# # remove NA and dead cells
bri.integrated.skin <- subset(bri.integrated.skin, subset = CellTypePranali != "NA")
bri.integrated.skin <- subset(bri.integrated.skin, subset = CellTypePranali != "Dead")
bri.integrated.pbmc <- subset(bri.integrated.pbmc, subset = CellTypePranali != "NA")
bri.integrated.pbmc <- subset(bri.integrated.pbmc, subset = CellTypePranali != "Dead")
bri.integrated.pbmc <- bri.integrated.pbmc[,bri.integrated.pbmc$wsnn_res.1.6 != "11"]

# save RDS files for independently clustered skin and blood
saveRDS(bri.integrated.skin, "8-multimodal-sjs-skin-only.rds")
saveRDS(bri.integrated.pbmc_nonmt, "8-multimodal-sjs-pbmc-only.rds")





