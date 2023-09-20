library(tidyverse)
library(Seurat)
library(patchwork)

#### 
## Data Import ##############
####

# # read cellranger 10x filtered feature matrices
# # the filtered barcode matrices of 3 sequencing flowcells, which combine gene expression and antibody derived tags / hashtags (6 total barcode matrices)
bri793 <- Read10X(data.dir="BRI-793-/cite/filtered_feature_bc_matrix/")
bri793hash <- Read10X(data.dir="BRI-793-/hash/filtered_feature_bc_matrix/")
bri817 <- Read10X(data.dir="BRI-817-/cite/filtered_feature_bc_matrix/")
bri817hash <- Read10X(data.dir="BRI-817-/hash/filtered_feature_bc_matrix/")
bri820 <- Read10X(data.dir="BRI-820-/cite/filtered_feature_bc_matrix/")
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

#### 
## HTO Analysis ###########
####

# Normalize HTO with central log ratio transformation, use 99th quantile for cutoff
bri793.obj <- NormalizeData(bri793.obj, assay = "HTO", normalization.method = "CLR")
bri793.obj <- HTODemux(bri793.obj, assay = "HTO", positive.quantile = 0.99)
bri817.obj <- NormalizeData(bri817.obj, assay = "HTO", normalization.method = "CLR")
bri817.obj <- HTODemux(bri817.obj, assay = "HTO", positive.quantile = 0.99)
bri820.obj <- NormalizeData(bri820.obj, assay = "HTO", normalization.method = "CLR")
bri820.obj <- HTODemux(bri820.obj, assay = "HTO", positive.quantile = 0.99)

# show how many negatives, singlets and doublets by HTO counts
p1 <- table(bri793.obj@meta.data$HTO_classification.global) %>%
  as.data.frame()  %>%
  ggplot(., aes(x=Var1, y = Freq)) + 
  geom_bar(aes(fill=Var1), position = "dodge", stat = "identity", show.legend = FALSE) +
  labs(title = "BRI793 Hashtag Oligo Counts")

p2 <- table(bri817.obj@meta.data$HTO_classification.global) %>%
  as.data.frame()  %>%
  ggplot(., aes(x=Var1, y = Freq)) + 
  geom_bar(aes(fill=Var1), position = "dodge", stat = "identity", show.legend = FALSE) +
  labs(title = "BRI817 Hashtag Oligo Counts")

p3 <- table(bri820.obj@meta.data$HTO_classification.global) %>%
  as.data.frame()  %>%
  ggplot(., aes(x=Var1, y = Freq)) + 
  geom_bar(aes(fill=Var1), position = "dodge", stat = "identity", show.legend = FALSE) +
  labs(title = "BRI820 Hashtag Oligo Counts")

p1 + p2 + p3

# filter out doublets and negatives
bri793.obj <- subset(bri793.obj, subset=HTO_classification.global == 'Singlet')
bri817.obj <- subset(bri817.obj, subset=HTO_classification.global == 'Singlet')
bri820.obj <- subset(bri820.obj, subset=HTO_classification.global == 'Singlet')

# define lists of Seurat Objects based on the the hash.id, 
bri793.list <- SplitObject(bri793.obj, split.by="hash.ID")
bri793.list <- bri793.list[-3] # remove Control030Skin which has very low cells
bri817.list <- SplitObject(bri817.obj, split.by="hash.ID")
bri820.list <- SplitObject(bri820.obj, split.by="hash.ID")
bri.list <- do.call(c, list(bri793.list,bri817.list,bri820.list))

####
## Filter and Process RNA #######
# Filtering RNA and ADT counts, Normalize, Find variable features
# low RNA 1000, high RNA 20000, 20% percent.mt
# Remove ribosomal RNA transcripts (RPS and RPL) 
####

bri.list <- do.call(c, list(bri793.list,bri817.list,bri820.list))
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
## Integration for RNA #######
####

bri.anchors <- FindIntegrationAnchors(object.list=bri.list, dims=1:30)
bri.integrated <- IntegrateData(anchorset = bri.anchors, new.assay.name = "integratedRNA", dims=1:30)

####
## Filter and Process RNA #######
# Filtering ADT and Normalize 
# high ADT 10000 
# Remove CD4 and CD8 doublets
####

DefaultAssay(bri.integrated) <- "ADT"
Idents(bri.integrated) <- "hash.ID"
VlnPlot(bri.integrated, features = c("nCount_ADT"))
bri.integrated.filterADT <- subset(bri.integrated, nCount_ADT < 10000)
VlnPlot(bri.integrated.filterADT , features = c("nCount_ADT"))

bri.integrated.list <- SplitObject(bri.integrated.filterADT, split.by = "hash.ID")

for (i in 1:length(bri.integrated.list)) {
  DefaultAssay(bri.integrated.list[[i]]) <- "ADT"
  bri.integrated.list[[i]] <- NormalizeData(bri.integrated.list[[i]], normalization.method = "CLR")
  bri.integrated.list[[i]] <- subset(bri.integrated.list[[i]], subset = adt_CD4.1 < 0.75 | adt_CD8 < 1)
}

####
## Integration for ADT #######
####

bri.anchors <- FindIntegrationAnchors(object.list=bri.integrated.list, dims=1:12)
bri.integrated <- IntegrateData(anchorset = bri.anchors, new.assay.name = "integratedADT", dims=1:12)


####
## MultiModal Analysis for RNA and ADT: #######
####

DefaultAssay(bri.integrated) <- "integratedRNA"
bri.integrated <- ScaleData(bri.integrated) %>% RunPCA(features=rownames(bri.integrated))

DefaultAssay(bri.integrated) <- 'integratedADT'
VariableFeatures(bri.integrated) <- rownames(bri.integrated[["ADT"]])
bri.integrated <- ScaleData(bri.integrated) %>% RunPCA(reduction.name = 'apca')

# wnn analysis 
bri.integrated <- FindMultiModalNeighbors(bri.integrated, reduction.list = list("pca", "apca"),
                                          dims.list = list(1:30, 1:12), modality.weight.name = "RNA.weight")


# UMAP for WNN
bri.integrated <- RunUMAP(bri.integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p1 <- DimPlot(bri.integrated, reduction = "wnn.umap", group.by = "hash.ID")
p1

#UMAP for RNA and ADT
bri.integrated <- RunUMAP(bri.integrated, reduction = 'pca', dims = 1:30, assay = 'integratedRNA', 
                          reduction.name = 'rnaintegrated.umap', reduction.key = 'rnaUMAP_')

# setting apca to dims 1:13 makes for out of bounds
bri.integrated <- RunUMAP(bri.integrated, reduction = 'apca', dims = 1:12, assay = 'integratedADT', 
                          reduction.name = 'adtintegrated.umap', reduction.key = 'adtintegratedUMAP_')

p2 <- DimPlot(bri.integrated, reduction = 'rnaintegrated.umap', group.by = "hash.ID")
p3 <- DimPlot(bri.integrated, reduction = 'adtintegrated.umap', group.by = "hash.ID")
p2 + p3

####
## Additional Metadata ####
####

# Assign tissue type, either broadly PBMC or Skin
hashIDs <- bri.integrated@meta.data$hash.ID
hashIDs[grepl('PBMC',hashIDs)] <- "PBMC"
hashIDs[grepl("Skin|SKIN",hashIDs)] <- "Skin"
bri.integrated@meta.data$Tissue <- hashIDs

# Assigng tissue and sample type, specifically which type of PBMC or Skin
hashIDs <- bri.integrated@meta.data$hash.ID
hashIDs[grepl("PBMC",hashIDs) & grepl("SJS",hashIDs)] <- "SJS_PBMC"
hashIDs[grepl("PBMC",hashIDs) & grepl("MDE",hashIDs)] <- "MDE_PBMC"
hashIDs[grepl("PBMC",hashIDs) & grepl("BCH",hashIDs)] <- "BCH_PBMC"
hashIDs[grepl("SKIN",hashIDs) & grepl("SJS",hashIDs)] <- "SJS_Skin"
hashIDs[grepl("SKIN",hashIDs) & grepl("MDE",hashIDs)] <- "MDE_Skin"
hashIDs[grepl("Skin",hashIDs) & grepl("MDE",hashIDs)] <- "MDE_Skin"
hashIDs[grepl("Skin",hashIDs) & grepl("Control",hashIDs)] <- "Control_Skin"
bri.integrated@meta.data$phenotypeTissue <- hashIDs


####
## Clustering #######
####

# Clustering for resolutions c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
bri.integrated <- FindClusters(bri.integrated, graph.name = "wsnn", algorithm = 3, 
                               resolution = c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0), verbose = FALSE)


## Visualize clusters #### 
Idents(bri.integrated) <- "wsnn_res.0.4"
p1 <- DimPlot(bri.integrated, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.4")
Idents(bri.integrated) <- "wsnn_res.0.6"
p2 <- DimPlot(bri.integrated, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.6")
Idents(bri.integrated) <- "wsnn_res.0.8"
p3 <- DimPlot(bri.integrated, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 0.8")
Idents(bri.integrated) <- "wsnn_res.1"
p4 <- DimPlot(bri.integrated, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.0")
Idents(bri.integrated) <- "wsnn_res.1.2"
p5 <- DimPlot(bri.integrated, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.2")
Idents(bri.integrated) <- "wsnn_res.1.4"
p6 <- DimPlot(bri.integrated, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.4")
Idents(bri.integrated) <- "wsnn_res.1.6"
p7 <- DimPlot(bri.integrated, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.6")
Idents(bri.integrated) <- "wsnn_res.1.8"
p8 <- DimPlot(bri.integrated, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + labs(title = "resolution 1.8")
p1 + p2 + p3 + p4
p5 + p6 + p7 + p8

## save bri.integrated to a rds file ####
saveRDS(bri.integrated, "1-multimodal-sjs-clustered.rds")





