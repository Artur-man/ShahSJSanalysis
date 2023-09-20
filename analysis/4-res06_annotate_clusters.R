# libraries 
library(Seurat)
library(tidyverse)
library(patchwork)

####
## import Data #####
#### 
bri.integrated <- readRDS("2-multimodal-sjs-clustered-tcr.rds")
load("3-multimodal-sjs-clustered-tcr-subclusters.Rdata")

# import subclusters
bri.integrated$seurat_clusters <- as.character(bri.integrated$wsnn_res.0.6)
bri.integrated$seurat_clusters[bri.integrated$seurat_clusters==1] <- paste(1,bri.integrated_cluster1$sub.cluster, sep = "_")
bri.integrated$seurat_clusters[bri.integrated$seurat_clusters==9] <- paste(9,bri.integrated_cluster9$wsnn_res.0.8, sep = "_")
bri.integrated$seurat_clusters[bri.integrated$seurat_clusters %in% c(0,3,2)] <- paste("032",bri.integrated_cluster0$wsnn_res.0.6, sep = "_")
bri.integrated$seurat_clusters[bri.integrated$seurat_clusters==8] <- paste(8,bri.integrated_cluster8$wsnn_res.0.4, sep = "_")

# Delete a small artifact cluster with only 2 cells
bri.integrated <- subset(bri.integrated, subset = seurat_clusters != "14")

# Cell Types
Idents(bri.integrated) <- "seurat_clusters"
bri.integrated <- RenameIdents(object = bri.integrated,
                                  "032_0" = "CD4+ TCM",
                                  "032_1" = "CD4+ CD45RA+ CD62L-low",
                                  "032_2" = "CD4 Naive",
                                  "032_3" = "CD4+ CD45RA- CD62L-low",
                                  "032_4" = "CD8+ TCM",
                                  "032_5" = "CD4+ CD103- TRM",
                                  "1_0" = "CD8+ T Effector",
                                  "1_1_0" = "CD8+ CD103- TRM", 
                                  "1_1_1" = "CD8+ TEM",
                                  "1_1_2" = "CD8+ TMM",
                                  "1_2" = "CD8+ TEM", 
                                  "1_3" = "CD8+ TMM",
                                  "4" = "CD4+ Treg 2",   
                                  "5" = "CD8+ CD103- TEMRA",         
                                  "6" = "CD4+ T Effector",
                                  "7" = "CD8+ Naive",            
                                  "8_1" = "GD T cells",
                                  "8_0" = "CD8+ CD56+ T cells",
                                  "8_2" = "CD8+ CD56+ T cells",
                                  "8_3" = "CD8+ CD56+ T cells",
                                  "9_0" = "CD4+ CD103+ TRM",
                                  "9_1" = "CD8+ CD103+ TRM",
                                  "9_2" = "CD4+ CD103+ TRM", 
                                  "9_3" = "CD8+ CD103+ TEMRA", 
                                  "10" = "CD4+ Treg 1",
                                  "11" = "Proliferating",
                                  "12" = "Dead",
                                  "13" = "CD8+ CD103- TRM (nf)")
bri.integrated$CellType <- Idents(bri.integrated)

# Cell Types
Idents(bri.integrated) <- "seurat_clusters"
bri.integrated <- RenameIdents(object = bri.integrated,
                               "032_0" = "CD4",
                               "032_1" = "CD4",
                               "032_2" = "CD4",
                               "032_3" = "CD4",
                               "032_4" = "CD8",
                               "032_5" = "CD4",
                               "1_0" = "CD8",
                               "1_1_0" = "CD8", 
                               "1_1_1" = "CD8",
                               "1_1_2" = "CD8",
                               "1_2" = "CD8", 
                               "1_3" = "CD8",
                               "4" = "CD4+ Treg 2",   
                               "5" = "CD8",         
                               "6" = "CD4",
                               "7" = "CD8",            
                               "8_1" = "GD T cells",
                               "8_0" = "CD8",
                               "8_2" = "CD8",
                               "8_3" = "CD8",
                               "9_0" = "CD4",
                               "9_1" = "CD8",
                               "9_2" = "CD4", 
                               "9_3" = "CD8", 
                               "10" = "CD4+ Treg 1",
                               "11" = "Proliferating",
                               "12" = "Dead",
                               "13" = "CD8")
bri.integrated$CD4_CD8_status <- Idents(bri.integrated)

# save final RDS, which is annotated, and used for downstream analyses.
saveRDS(bri.integrated, "final-multimodal-sjs-clustered-tcr-annotated.rds")

