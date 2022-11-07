# libraries 
library(Seurat)
library(tidyverse)
library(anndata)

####
## import Data #####
#### 
bri.integrated <- readRDS("sjs_clustered_tcr.rds")
load("sjs_clustered_tcr_subclusters.Rdata")

# import subclusters
bri.integrated$seurat_clusters <- as.character(bri.integrated$wsnn_res.0.6)
bri.integrated$seurat_clusters[bri.integrated$seurat_clusters==1] <- paste(1,bri.integrated_cluster1$sub.cluster, sep = "_")
bri.integrated$seurat_clusters[bri.integrated$seurat_clusters==9] <- paste(9,bri.integrated_cluster9$wsnn_res.0.8, sep = "_")
bri.integrated$seurat_clusters[bri.integrated$seurat_clusters %in% c(0,3,2)] <- paste("032",bri.integrated_cluster0$wsnn_res.0.4, sep = "_")
bri.integrated$seurat_clusters[bri.integrated$seurat_clusters==8] <- paste(8,bri.integrated_cluster8$wsnn_res.0.4, sep = "_")

# cell types Wei and Artur
Idents(bri.integrated) <- "seurat_clusters"
bri.integrated <- RenameIdents(object = bri.integrated,
                                  #"032_0" = "CD4 Effector",
                                  "032_0" = "CD4 Effector Memory", # Pranali
                                  "032_1" = "CD4 Migratory",
                                  "032_2" = "CD4 Naive",
                                  "032_3" = "CD4 CD56+ NKT-like",
                                  "032_4" = "CD8 Central Memory",
                                  # "032_5" = "CD4 Effector Memory",
                                  "032_5" = "CD4 CD103- TRM",
                                  "1_0" = "CD8 SLEC",
                                  # "1_1" = "MPEC",
                                  "1_1_0" = "CD8 CD103- TRM", # Pranali
                                  "1_1_1" = "CD8 IL7Ra- Effector Memory (RO++)", # for now
                                  "1_1_2" = "CD8 MPEC (2)", # for now
                                  "1_2" = "CD8 IL7Ra- Effector Memory", # Pranali
                                  "1_3" = "CD8 MPEC",
                                  "10" = "LAMP-1+ CD335+",
                                  "11" = "Proliferating",
                                  "12" = "Dying Cells (MT++)",
                                  "13" = "CD8 NKT",
                                  "14" = "14?",
                                  "4" = "Treg",
                                  "5" = "CD8 Effector",
                                  "6" = "CD4 Cytotoxic",
                                  "7" = "CD8 Naive",
                                  # "8" = "Gamma Delta",
                                  "8_1" = "Gamma Delta",
                                  "8_0" = "CD8 CD56+ Effector",
                                  "8_2" = "CD8 CD56+ Effector",
                                  "8_3" = "CD8 CD56+ Effector",
                                  "9_0" = "CD4 CD103+ TRM",
                                  "9_1" = "CD8 CD45RO+ CD103+ TRM",
                                  "9_2" = "CD4 CD103+ TRM",
                                  # "9_3" = "CD8 MAIT",
                                  "9_3" = "CD8 CD45RA+ CD103+ TRM") #pranali

bri.integrated$CellType <- Idents(bri.integrated)

## Take out Cluster 14? ####
# bri.integrated <- bri.integrated[,bri.integrated$CellType!="14?"]
bri.integrated <- subset(bri.integrated, subset = CellType != "14?")
temp <- droplevels(bri.integrated$CellType)
bri.integrated$CellType <- temp

## Give New Names of Cell Types #####
levels(bri.integrated@meta.data$CellType)
bri.integrated@meta.data$CellTypeNumber <- bri.integrated@meta.data$CellType
temp <- bri.integrated@meta.data$CellTypeNumber
temp[temp == "CD8 IL7Ra- Effector Memory (RO++)"] <- "CD8 IL7Ra- Effector Memory"
temp[temp == "CD8 MPEC (2)"] <- "CD8 MPEC"
temp <- droplevels(temp)
bri.integrated@meta.data$CellTypeNumber <- temp
# levels(bri.integrated@meta.data$CellTypeNumber) <- paste("Cluster",1:length(levels(bri.integrated@meta.data$CellTypeNumber)), sep = "")

# Cell Type Numbers from 1-23
Idents(bri.integrated) <- "CellTypeNumber"
bri.integrated <- RenameIdents(object = bri.integrated,
                               "CD4 Effector Memory" = "18", 
                               "CD4 Migratory" = "12",
                               "CD4 Naive" = "13" ,
                               "CD4 CD56+ NKT-like" = "14",
                               "CD8 Central Memory" = "6",
                               "CD4 CD103- TRM" = "19",
                               "CD8 SLEC" = "1",
                               "CD8 CD45RO+ CD103+ TRM" = "8",
                               "CD8 IL7Ra- Effector Memory" = "10", 
                               "CD8 MPEC" = "7",
                               "LAMP-1+ CD335+" = "15",
                               "Proliferating" = "23",
                               "Dying Cells (MT++)" = "21",
                               "CD8 NKT" = "11",
                               "Treg" = "20",
                               "CD8 Effector" = "2",
                               "CD4 Cytotoxic" = "16",
                               "CD8 Naive" = "3",
                               "Gamma Delta" = "22",
                               "CD8 CD56+ Effector" = "5",
                               "CD4 CD103+ TRM" = "17",
                               "CD8 CD45RA+ CD103+ TRM" = "4",
                               "CD8 CD103- TRM" = "9")
bri.integrated$CellTypeNumber <- Idents(bri.integrated)

# Major Cell Types
Idents(bri.integrated) <- "CellType"
bri.integrated <- RenameIdents(object = bri.integrated,
                               "CD4 Effector Memory" = "CD4", 
                               "CD4 Migratory" = "CD4",
                               "CD4 Naive" = "CD4" ,
                               "CD4 CD56+ NKT-like" = "CD4",
                               "CD8 Central Memory" = "CD8",
                               "CD4 CD103- TRM" = "CD4",
                               "CD8 SLEC" = "CD8",
                               "CD8 CD45RO+ CD103+ TRM" = "CD8",
                               "CD8 IL7Ra- Effector Memory" = "CD8", 
                               "CD8 MPEC" = "CD8",
                               "LAMP-1+ CD335+" = "NA",
                               "Proliferating" = "NA",
                               "Dying Cells (MT++)" = "NA",
                               "CD8 NKT" = "CD8",
                               "Treg" = "CD4",
                               "CD8 Effector" = "CD8",
                               "CD4 Cytotoxic" = "CD4",
                               "CD8 Naive" = "CD8",
                               "Gamma Delta" = "NA",
                               "CD8 CD56+ Effector" = "CD8",
                               "CD4 CD103+ TRM" = "CD4",
                               "CD8 CD45RA+ CD103+ TRM" = "CD8",
                               "CD8 CD103- TRM" = "CD8",
                               "CD8 IL7Ra- Effector Memory (RO++)" = "CD8",
                               "CD8 MPEC (2)" = "CD8")
bri.integrated$CD4_CD8_status <- Idents(bri.integrated)

# Cell types Pranali and Sherrie
Idents(bri.integrated) <- "CellTypeNumber"
bri.integrated <- RenameIdents(object = bri.integrated,
                               "18" = "CD4 TCM", 
                               "12" = "CD4 CD45RA+ CD62L-low",
                               "13" = "CD4 CD62L-high Naive" ,
                               "14" = "CD4 CD45RA- CD62L-low",
                               "6" = "CD8 TCM",
                               "19" = "CD4 CD103- CD69+ TRM",
                               "1" = "CD8 Effector",
                               "8" = "CD8 CD103+ TRM",
                               "10" = "CD8 TEM", 
                               "7" = "CD8 Migratory Memory?",
                               "15" = "CD4 Treg 1",
                               "23" = "Proliferating",
                               "21" = "Dead",
                               "11" = "CD8 CD103- CD69+ TRM",
                               "20" = "CD4 Treg 2",
                               "2" = "CD8 CD103- TEMRA",
                               "16" = "CD4 CD45RA+ Effector",
                               "3" = "CD8 Naive",
                               "22" = "GD",
                               "5" = "CD8 CD56+",
                               "17" = "CD4 CD103+ CD69+ TRM",
                               "4" = "CD8 CD103+ TEMRA",
                               "9" = "CD8 CD103- TRM")
bri.integrated$CellTypePranali <- Idents(bri.integrated)


# ## OLD - separate out CD8 naive, should not be in with CD8 eff; and separate TEMRA. Complete the other names of cluster below.
# # Cell types SuperCluster (combined clusters)
# Idents(bri.integrated) <- "CellTypeNumber"
# bri.integrated <- RenameIdents(object = bri.integrated,
#                                '8' = 'CD8 TRM',
#                                '11' = 'CD8 TRM',
#                                '9' = 'CD8 TRM',
#                                '19' = 'CD4 TRM',
#                                '17' = 'CD4 TRM',
#                                '15' = 'Treg',
#                                '20' = 'Treg',
#                                '22' = 'GD',
#                                '1' = 'CD8 effector',
#                                '3' = 'CD8 effector',
#                                '5' = 'CD8 effector',
#                                '10' = 'CD8 effector memory',
#                                '2' = 'CD8 effector memory',
#                                '4' = 'CD8 effector memory',
#                                '12' = 'CD4 effector',
#                                '13' = 'CD4 effector',
#                                '14' = 'CD4 effector',
#                                '16' = 'CD4 effector')
# 
# bri.integrated$CellTypeSuperCluster_old <- Idents(bri.integrated)


# Cell types SuperCluster (combined clusters)
Idents(bri.integrated) <- "CellTypeNumber"
bri.integrated <- RenameIdents(object = bri.integrated,
                               '8' = 'CD8 TRM',
                               '11' = 'CD8 TRM',
                               '9' = 'CD8 TRM',
                               '19' = 'CD4 TRM',
                               '17' = 'CD4 TRM',
                               '15' = 'Treg',
                               '20' = 'Treg',
                               '22' = 'GD',
                               '1' = 'CD8 effector',
                               '3' = 'CD8 naive',
                               '5' = 'CD8 effector',
                               '10' = 'CD8 effector memory',
                               "7" = "CD8 effector memory",
                               '2' = 'CD8 TEMRA',
                               '4' = 'CD8 TEMRA',
                               "6" = "CD8 TCM",
                               "18" = "CD4 TCM",
                               '12' = 'CD4 naive',
                               '13' = 'CD4 effector',
                               '14' = 'CD4 effector',
                               '16' = 'CD4 effector',
                               "23" = "Proliferating",
                               "21" = "Dead")
              
bri.integrated$CellTypeSuperCluster <- Idents(bri.integrated)

# save
saveRDS(bri.integrated, "sjs_clustered_tcr_annotated.rds")

