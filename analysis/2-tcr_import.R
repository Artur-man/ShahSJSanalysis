# libraries
library(Seurat)
library(djvdj)
library(patchwork)
library(tidyverse)

# pull data
integrated_all <- readRDS("sjs_clustered.rds")
integrated_all_split <- SplitObject(integrated_all, split.by = "HTO_classification")

# tcr paths
path_795 <- c("../Data/Wei_Data/BRI-793-Wei/vdj/")
path_819 <- c("../Data/Wei_Data/BRI-817-Wei/vdj/")
path_822 <- c("../Data/Wei_Data/BRI-820-Wei/vdj/")

# names
BRI793_names <- c("Control-BCH-048-PBMC","MDE-006-PBMC","Control-BCH-001-PBMC","MDE-006-Skin","Control-039-Skin")
BRI817_names <- c("BCH026PBMC","MDE001PBMC","SJS003PBMC","Control031Skin","SJS003SKIN","MDE001SKIN")
BRI820_names <- c("SJS002PBMC","SJS001SKIN","MDE002SKIN","SJS001PBMC","SJS002SKIN","MDE002PBMC")
names_data <- rbind(cbind(BRI793_names,path_795),
                    cbind(BRI817_names,path_819),
                    cbind(BRI820_names,path_822))
names_data <- cbind(names_data,gsub("SKIN|Skin|PBMC","",names_data[,1]))
colnames(names_data) <- c("Sample","BioSample","Label")
names_data <- as.data.frame(names_data)
names_samples <- names(integrated_all_split)

# import vdj data for all demultiplexed samples
for(i in 1:length(integrated_all_split)){

  print(i)
  
  # exclude suffixes
  seu <- integrated_all_split[[i]]
  seu_cells <- sapply(Cells(seu), function(x){
    return(strsplit(x, split = "_")[[1]][1])
  })
  seu <- RenameCells(seu, new.names = seu_cells)
  
  # match tcrs and import to metadata
  path = names_data[names_data$Sample==names_samples[i],2]
  label = names_data[names_data$Sample==names_samples[i],3]
  seu <- import_vdj(
    sobj_in = seu,                       # Seurat object
    vdj_dir = path                       # Cellranger output directories
  )

  # attach sample names to clonotype IDs
  seu$clonotype_id[!is.na(seu$clonotype_id)] <- paste(seu$clonotype_id[!is.na(seu$clonotype_id)],label,sep = "_")
  seu <- RenameCells(seu, new.names = paste(Cells(seu),i, sep = "_"))
  integrated_all_split[[i]] <- seu
} 

# merge
integrated_all_merged <- merge(x = integrated_all_split[[1]], y = integrated_all_split[-1])
metadata <- integrated_all_merged@meta.data
integrated_all_merged_new <- integrated_all
integrated_all_merged_new@meta.data <- metadata

# save rds
saveRDS(integrated_all_merged_new, "sjs_clustered_tcr.rds")

