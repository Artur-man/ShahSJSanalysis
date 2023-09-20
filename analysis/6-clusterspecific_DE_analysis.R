# libraries
library(Seurat)
library(patchwork)
library(tidyverse)
library(anndata)
library(scales)
library(circlize)
library(reshape2)
library(speckle)
library(edgeR)
library(DESeq2)

####
## Import Data #####
####

bri.integrated <- readRDS("final-multimodal-sjs-clustered-tcr-annotated.rds")
bri.integrated$phenotype <- sapply(bri.integrated$phenotypeTissue, function(x) return(strsplit(x, split = "_")[[1]][1]))
bri.integrated$phenotype <- gsub("BCH","Control", bri.integrated$phenotype)

#### 
## Delete Dead cells from UMAP ###
####

bri.integrated <- bri.integrated[,!bri.integrated$CellType == "Dying Cells (MT++)"]

#### 
## DE analysis with edgeR Test on bulk level ####
####

allcelltypes <- unique(bri.integrated$CellType)
allcelltypes <- allcelltypes[!allcelltypes %in% c("GD T cells", "Dead", "Proliferating")]

# tissue
all_results <- NULL
for(tis in c("Skin","PBMC")){
  
  # cell type 
  for(celltype in allcelltypes){
    
    # select subset for tissue and cell type
    bri.integrated_subset <- bri.integrated[,bri.integrated$CellType == celltype & bri.integrated$Tissue == tis]
    bri.integrated_subset <- AggregateExpression(bri.integrated_subset, assays = "RNA", slot = "counts", group.by = "hash.ID")
    bri.integrated_subset_rna <- bri.integrated_subset$RNA
    
    # treatment
    treatment <- ifelse(!grepl("MDE|SJS", colnames(bri.integrated_subset_rna)), 
                        "Ctrl",
                        ifelse(grepl("SJS", colnames(bri.integrated_subset_rna)), "SJS","MDE"))
    combination <- combn(unique(treatment), 2)
    
    # print
    print(treatment)
    
    # pair of hash.ID
    for(kk in 1:ncol(combination)){
      
      # print 
      print(paste0(tis, "_", celltype, "_", paste(rev(combination[,kk]), collapse = "_vs_")))
      
      # select subset
      counts_data <- bri.integrated_subset_rna[,treatment %in% combination[,kk]]
      treat_cur <- treatment[treatment %in% combination[,kk]]

      # filtered
      max_count <- apply(counts_data,1,max)
      counts_data_filtered <- counts_data[which(max_count > 10),]
      
      # select variable
      rvar <- apply(counts_data_filtered, 1, var)
      idx.keep <- which(rvar > 0)
      counts_data_filtered <- counts_data_filtered[idx.keep, ]
      
      # test
      y <- edgeR::DGEList(counts = counts_data_filtered, group = factor(treat_cur))
      y <- edgeR::calcNormFactors(y)
      design <- model.matrix(~0 + treat_cur, data = y$samples)
      y <- edgeR::estimateDisp(y, design)
      fit <- edgeR::glmFit(y, design)
      lrt <- edgeR::glmLRT(fit, contrast = c(-1,1))
      res_table <- as.data.frame(edgeR::topTags(lrt, p.value = 1, n = Inf, sort.by = "logFC"))
      
      # get data
      res_table$ident1 <- combination[,kk][2]
      res_table$ident2 <- combination[,kk][1]
      res_table$Tissue <- tis
      res_table$CellType <- celltype
      res_table$gene <- rownames(res_table)
      all_results <- rbind(all_results, res_table)
    }
  }
}

write.table(all_results, file = "revision_files_17122022/DEResults_BulkedgeR.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

#### 
## DE analysis with edgeR Test on bulk level on ADT ####
####

allcelltypes <- unique(bri.integrated$CellType)
allcelltypes <- allcelltypes[!allcelltypes %in% c("GD T cells", "Dead", "Proliferating")]

# tissue
all_results <- NULL
for(tis in c("Skin","PBMC")){
  
  # cell type 
  for(celltype in allcelltypes){
    
    # select subset for tissue and cell type
    bri.integrated_subset <- bri.integrated[,bri.integrated$CellType == celltype & bri.integrated$Tissue == tis]
    bri.integrated_subset <- AggregateExpression(bri.integrated_subset, assays = "ADT", slot = "counts", group.by = "hash.ID")
    bri.integrated_subset_adt <- bri.integrated_subset$ADT
    
    # treatment
    treatment <- ifelse(!grepl("MDE|SJS", colnames(bri.integrated_subset_adt)), 
                        "Ctrl",
                        ifelse(grepl("SJS", colnames(bri.integrated_subset_adt)), "SJS","MDE"))
    combination <- combn(unique(treatment), 2)
    
    # print
    print(treatment)
    
    # pair of hash.ID
    for(kk in 1:ncol(combination)){
      
      # print 
      print(paste0(tis, "_", celltype, "_", paste(rev(combination[,kk]), collapse = "_vs_")))
      
      # select subset
      counts_data <- bri.integrated_subset_adt[,treatment %in% combination[,kk]]
      treat_cur <- treatment[treatment %in% combination[,kk]]
      
      # filtered
      max_count <- apply(counts_data,1,max)
      counts_data_filtered <- counts_data[which(max_count > 10),]
      
      # select variable
      rvar <- apply(counts_data_filtered, 1, var)
      idx.keep <- which(rvar > 0)
      counts_data_filtered <- counts_data_filtered[idx.keep, ]
      
      # test
      y <- edgeR::DGEList(counts = counts_data_filtered, group = factor(treat_cur))
      y <- edgeR::calcNormFactors(y)
      design <- model.matrix(~0 + treat_cur, data = y$samples)
      y <- edgeR::estimateDisp(y, design)
      fit <- edgeR::glmFit(y, design)
      lrt <- edgeR::glmLRT(fit, contrast = c(-1,1))
      res_table <- as.data.frame(edgeR::topTags(lrt, p.value = 1, n = Inf, sort.by = "logFC"))
      
      # get data
      # res_table$TestType <- paste0(tis, "_", celltype, "_", paste(rev(combination[,kk]), collapse = "_vs_"))
      res_table$ident1 <- combination[,kk][2]
      res_table$ident2 <- combination[,kk][1]
      res_table$Tissue <- tis
      res_table$CellType <- celltype
      res_table$gene <- rownames(res_table)
      all_results <- rbind(all_results, res_table)
    }
  }
}

write.table(all_results, file = "revision_files_17122022/DEResults_BulkedgeR_ADT.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

#### 
## DE analysis with edgeR Test on bulk level on combined Treg CD4 and CD8 ####
####

bri.integrated$CD4_CD8_Treg <- ifelse(grepl("Treg",bri.integrated$CellType), "Treg", 
                                      ifelse(grepl("CD4",bri.integrated$CellType), "CD4", 
                                             ifelse(grepl("CD8",bri.integrated$CellType), "CD8", as.character(bri.integrated$CellType))))
allcelltypes <- unique(bri.integrated$CD4_CD8_Treg)
allcelltypes <- allcelltypes[!allcelltypes %in% c("GD T cells", "Dead", "Proliferating")]

# tissue
all_results <- NULL
for(tis in c("Skin","PBMC")){
  
  # cell type 
  for(celltype in allcelltypes){
    
    # select subset for tissue and cell type
    bri.integrated_subset <- bri.integrated[,bri.integrated$CD4_CD8_Treg == celltype & bri.integrated$Tissue == tis]
    bri.integrated_subset <- AggregateExpression(bri.integrated_subset, assays = "RNA", slot = "counts", group.by = "hash.ID")
    bri.integrated_subset_rna <- bri.integrated_subset$RNA
    
    # treatment
    treatment <- ifelse(!grepl("MDE|SJS", colnames(bri.integrated_subset_rna)), 
                        "Ctrl",
                        ifelse(grepl("SJS", colnames(bri.integrated_subset_rna)), "SJS","MDE"))
    combination <- combn(unique(treatment), 2)
    
    # print
    print(treatment)
    
    # pair of hash.ID
    for(kk in 1:ncol(combination)){
      
      # print
      print(c(paste0(tis, "_", celltype, "_", paste(rev(combination[,kk]), collapse = "_vs_"))))
            
      # select subset
      counts_data <- bri.integrated_subset_rna[,treatment %in% combination[,kk]]
      treat_cur <- treatment[treatment %in% combination[,kk]]
      
      # filtered
      max_count <- apply(counts_data,1,max)
      counts_data_filtered <- counts_data[which(max_count > 10),]

      # # select variable
      rvar <- apply(counts_data_filtered, 1, var)
      idx.keep <- which(rvar > 0)
      counts_data_filtered <- counts_data_filtered[idx.keep, ]
      
      # test
      y <- edgeR::DGEList(counts = counts_data_filtered, group = factor(treat_cur))
      y <- edgeR::calcNormFactors(y)
      design <- model.matrix(~0 + treat_cur, data = y$samples)
      y <- edgeR::estimateDisp(y, design)
      fit <- edgeR::glmFit(y, design)
      lrt <- edgeR::glmLRT(fit, contrast = c(-1,1))
      res_table <- as.data.frame(edgeR::topTags(lrt, p.value = 1, n = Inf, sort.by = "logFC"))
      
      # get data
      res_table$ident1 <- combination[,kk][2]
      res_table$ident2 <- combination[,kk][1]
      res_table$Tissue <- tis
      res_table$CellType <- celltype
      res_table$gene <- rownames(res_table)
      all_results <- rbind(all_results, res_table)
    }
  }
}

write.table(all_results, file = "revision_files_17122022/DEResults_BulkedgeR_cd4cd8.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

#### 
## Merge DE Results ####
####

DEResults_BulkedgeR_cd4cd8 <- read.table("revision_files_17122022/DEResults_BulkedgeR_cd4cd8.tsv", sep = "\t", header = T)
DEResults_BulkedgeR <- read.table("revision_files_17122022/DEResults_BulkedgeR.tsv", sep = "\t", header = T)
DEResults_BulkedgeR_adt <- read.table("revision_files_17122022/DEResults_BulkedgeR_ADT.tsv", sep = "\t", header = T)
DEResults_BulkedgeR_cd4cd8$TestType <- "CD4-CD8-Treg"
DEResults_BulkedgeR$TestType <- "CellAnnotations"
DEResults_BulkedgeR_adt$TestType <- "CellAnnotations (ADT)"
DEResults_BulkedgeR_all <- rbind(DEResults_BulkedgeR_cd4cd8, DEResults_BulkedgeR, DEResults_BulkedgeR_adt)
common_colnames <- c("logFC", "PValue", "FDR", "Group1", "Group2", "Tissue", "CellType", "gene", "TestType")
DEResults_BulkedgeR_all <- DEResults_BulkedgeR_all[,c("logFC", "PValue", "FDR", "ident1", "ident2", "Tissue", "CellType", "gene", "TestType")]
colnames(DEResults_BulkedgeR_all) <- common_colnames
DEResults_BulkedgeR_all$Group2 <- gsub("Control","Ctrl",DEResults_BulkedgeR_all$Group2)
DEResults_BulkedgeR_all_sig <- DEResults_BulkedgeR_all[DEResults_BulkedgeR_all$FDR < 0.05,]
write.table(DEResults_BulkedgeR_all, file = "revision_files_17122022/DEResults_BulkedgeR_all.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

