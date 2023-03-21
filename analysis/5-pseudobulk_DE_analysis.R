# libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DESeq2)
library(debrowser)
library(vsn)

#### 
## Data Import ##############
####

# # read cellranger 10x filtered feature matrices 
bri793 <- Read10X(data.dir="BRI-793-Wei/cite/filtered_feature_bc_matrix/")
bri793hash <- Read10X(data.dir="BRI-793-Wei/hash/filtered_feature_bc_matrix/")
bri817 <- Read10X(data.dir="BRI-817-Wei/cite/filtered_feature_bc_matrix/")
bri817hash <- Read10X(data.dir="BRI-817-Wei/hash/filtered_feature_bc_matrix/")
bri820 <- Read10X(data.dir="BRI-820-Wei/cite/filtered_feature_bc_matrix/")
bri820hash <- Read10X(data.dir="BRI-820-Wei/hash/filtered_feature_bc_matrix/")

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

# Normalize HTO with central log transformation, use 99th quantile for cutoff
bri793.obj <- NormalizeData(bri793.obj, assay = "HTO", normalization.method = "CLR")
bri793.obj <- HTODemux(bri793.obj, assay = "HTO", positive.quantile = 0.99)
bri817.obj <- NormalizeData(bri817.obj, assay = "HTO", normalization.method = "CLR")
bri817.obj <- HTODemux(bri817.obj, assay = "HTO", positive.quantile = 0.99)
bri820.obj <- NormalizeData(bri820.obj, assay = "HTO", normalization.method = "CLR")
bri820.obj <- HTODemux(bri820.obj, assay = "HTO", positive.quantile = 0.99)

# filter out doublets and negatives
bri793.obj <- subset(bri793.obj, subset=HTO_classification.global == 'Singlet')
bri817.obj <- subset(bri817.obj, subset=HTO_classification.global == 'Singlet')
bri820.obj <- subset(bri820.obj, subset=HTO_classification.global == 'Singlet')

# define lists of Seurat Objects based on the the hash.id, 
bri793.list <- SplitObject(bri793.obj, split.by="hash.ID")
bri817.list <- SplitObject(bri817.obj, split.by="hash.ID")
bri820.list <- SplitObject(bri820.obj, split.by="hash.ID")
bri.list <- do.call(c, list(bri793.list,bri817.list,bri820.list))
bri.list.merged <- merge(bri.list[[1]], bri.list[-1])

#### 
## DESeq2 on pseudo-bulk data ###########
####

# Get Bulk RNA and ADT
counts <- AggregateExpression(bri.list.merged, group.by = "hash.ID", slot = "counts")
counts_RNA <- counts$RNA
counts_ADT <- counts$ADT
condition <- c(rep("Control",6),rep("MDE",6),rep("SJS",6))
tissue <- c("PBMC","Skin","Skin",
            "PBMC","PBMC","Skin",
            "PBMC","Skin","PBMC",
            "Skin","PBMC","Skin",
            "PBMC","Skin","PBMC",
            "Skin","PBMC","Skin")

# visualize genes
counts_normalized <- apply(counts$RNA, 2, function(x) x/sum(x) * 1000000)
gene <- c("GNLY","IFNG","GZMA","GZMB","PRF1")
counts_ggplot <- data.frame(sample = colnames(counts_normalized), t(counts_normalized[gene,]), tissue = tissue, condition = condition)
counts_ggplot <- melt(counts_ggplot, measure.vars = gene)
ggplot(counts_ggplot, aes(x = sample, y = value, fill = tissue)) + 
  geom_bar(stat = "identity") + 
  facet_grid(variable~condition, scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  labs(title = gene)

#### 
### DESeq2 on RNA: SJS vs MDE, tissue as covariate  ####
#### 

counts_RNA_SJSandMDE <- counts_RNA[,condition %in% c("SJS","MDE")]
tissue_SJSandMDE <- tissue[condition %in% c("SJS","MDE")]
condition_SJSandMDE <- condition[condition %in% c("SJS","MDE")]

# Filter low gene counts
max_count <- apply(counts_RNA_SJSandMDE,1,max)
counts_RNA_SJSandMDE <- counts_RNA_SJSandMDE[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_RNA_SJSandMDE,
                              DataFrame(condition_SJSandMDE, tissue_SJSandMDE),
                              design= ~ condition_SJSandMDE + tissue_SJSandMDE)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_RNA_SJSandMDE.txt")

#### 
### DESeq2 on RNA: SJS vs MDE, Skin and Blood separately ####
#### 

# Skin
counts_RNA_SJSandMDE_skin <- counts_RNA[,condition %in% c("SJS","MDE") & tissue %in% "Skin"]
condition_SJSandMDE_skin <- condition[condition %in% c("SJS","MDE") & tissue %in% "Skin"]

# Filter low gene counts
max_count <- apply(counts_RNA_SJSandMDE_skin,1,max)
counts_RNA_SJSandMDE_skin <- counts_RNA_SJSandMDE_skin[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_RNA_SJSandMDE_skin,
                              DataFrame(condition_SJSandMDE_skin),
                              design= ~ condition_SJSandMDE_skin)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_RNA_SJSandMDE_Skin.txt")

# Blood
counts_RNA_SJSandMDE_blood <- counts_RNA[,condition %in% c("SJS","MDE") & tissue %in% "PBMC"]
condition_SJSandMDE_blood <- condition[condition %in% c("SJS","MDE") & tissue %in% "PBMC"]

# Filter low gene counts
max_count <- apply(counts_RNA_SJSandMDE_blood,1,max)
counts_RNA_SJSandMDE_blood <- counts_RNA_SJSandMDE_blood[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_RNA_SJSandMDE_blood,
                              DataFrame(condition_SJSandMDE_blood),
                              design= ~ condition_SJSandMDE_blood)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_RNA_SJSandMDE_Blood.txt")

#### 
### DESeq2 on RNA: SJS vs Control, Skin and Blood separately ####
#### 

# Skin
counts_RNA_SJSandCon_skin <- counts_RNA[,condition %in% c("SJS","Control") & tissue %in% "Skin"]
condition_SJSandCon_skin <- condition[condition %in% c("SJS","Control") & tissue %in% "Skin"]

# Filter low gene counts
max_count <- apply(counts_RNA_SJSandCon_skin,1,max)
counts_RNA_SJSandCon_skin <- counts_RNA_SJSandCon_skin[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_RNA_SJSandCon_skin,
                              DataFrame(condition_SJSandCon_skin),
                              design= ~ condition_SJSandCon_skin)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_RNA_SJSandCon_Skin.txt")

# Blood
counts_RNA_SJSandCon_blood <- counts_RNA[,condition %in% c("SJS","Control") & tissue %in% "PBMC"]
condition_SJSandCon_blood <- condition[condition %in% c("SJS","Control") & tissue %in% "PBMC"]

# Filter low gene counts
max_count <- apply(counts_RNA_SJSandCon_blood,1,max)
counts_RNA_SJSandCon_blood <- counts_RNA_SJSandCon_blood[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_RNA_SJSandCon_blood,
                              DataFrame(condition_SJSandCon_blood),
                              design= ~ condition_SJSandCon_blood)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_RNA_SJSandCon_Blood.txt")


#### 
### DESeq2 on RNA: MDE vs Control, Skin and Blood separately ####
#### 

# Skin
counts_RNA_MDEandCon_skin <- counts_RNA[,condition %in% c("MDE","Control") & tissue %in% "Skin"]
condition_MDEandCon_skin <- condition[condition %in% c("MDE","Control") & tissue %in% "Skin"]

# Filter low gene counts
max_count <- apply(counts_RNA_MDEandCon_skin,1,max)
counts_RNA_MDEandCon_skin <- counts_RNA_MDEandCon_skin[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_RNA_MDEandCon_skin,
                              DataFrame(condition_MDEandCon_skin),
                              design= ~ condition_MDEandCon_skin)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_RNA_MDEandCon_Skin.txt")

# Blood
counts_RNA_MDEandCon_blood<- counts_RNA[,condition %in% c("MDE","Control") & tissue %in% "PBMC"]
condition_MDEandCon_blood <- condition[condition %in% c("MDE","Control") & tissue %in% "PBMC"]

# Filter low gene counts
max_count <- apply(counts_RNA_MDEandCon_blood,1,max)
counts_RNA_MDEandCon_blood <- counts_RNA_MDEandCon_blood[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_RNA_MDEandCon_blood,
                              DataFrame(condition_MDEandCon_blood),
                              design= ~ condition_MDEandCon_blood)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_RNA_MDEandCon_Blood.txt")

#### 
### DESeq2 on RNA: Blood vs Skin: SJS, MDE and healthy separately ####
#### 

# SJS
counts_RNA_BloodvsSkin_SJS <- counts_RNA[,condition %in% c("SJS")]
tissue_RNA_BloodvsSkin_SJS <- tissue[condition %in% c("SJS")]

# Filter low gene counts
max_count <- apply(counts_RNA_BloodvsSkin_SJS,1,max)
counts_RNA_BloodvsSkin_SJS <- counts_RNA_BloodvsSkin_SJS[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_RNA_BloodvsSkin_SJS,
                              DataFrame(tissue_RNA_BloodvsSkin_SJS),
                              design= ~ tissue_RNA_BloodvsSkin_SJS)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_RNA_BloodvsSkin_SJS.txt")

# MDE
counts_RNA_BloodvsSkin_MDE <- counts_RNA[,condition %in% c("MDE")]
tissue_RNA_BloodvsSkin_MDE <- tissue[condition %in% c("MDE")]

# Filter low gene counts
max_count <- apply(counts_RNA_BloodvsSkin_MDE,1,max)
counts_RNA_BloodvsSkin_MDE <- counts_RNA_BloodvsSkin_MDE[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_RNA_BloodvsSkin_MDE,
                              DataFrame(tissue_RNA_BloodvsSkin_MDE),
                              design= ~ tissue_RNA_BloodvsSkin_MDE)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_RNA_BloodvsSkin_MDE.txt")

# Healthy
counts_RNA_BloodvsSkin_MDE <- counts_RNA[,condition %in% c("Control")]
tissue_RNA_BloodvsSkin_MDE <- tissue[condition %in% c("Control")]

# Filter low gene counts
max_count <- apply(counts_RNA_BloodvsSkin_MDE,1,max)
counts_RNA_BloodvsSkin_MDE <- counts_RNA_BloodvsSkin_MDE[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_RNA_BloodvsSkin_MDE,
                              DataFrame(tissue_RNA_BloodvsSkin_MDE),
                              design= ~ tissue_RNA_BloodvsSkin_MDE)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_RNA_BloodvsSkin_MDE.txt")

#### 
### DESeq2 on ADT: SJS vs MDE, tissue as covariate  ####
#### 

counts_ADT_SJSandMDE <- counts_ADT[,condition %in% c("SJS","MDE")]
tissue_SJSandMDE <- tissue[condition %in% c("SJS","MDE")]
condition_SJSandMDE <- condition[condition %in% c("SJS","MDE")]

# Filter low gene counts
max_count <- apply(counts_ADT_SJSandMDE,1,max)
counts_ADT_SJSandMDE <- counts_ADT_SJSandMDE[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_ADT_SJSandMDE,
                              DataFrame(condition_SJSandMDE, tissue_SJSandMDE),
                              design= ~ condition_SJSandMDE + tissue_SJSandMDE)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_ADT_SJSandMDE.txt")

#### 
### DESeq2 on ADT: SJS vs MDE, Skin and Blood separately ####
#### 

# Skin
counts_ADT_SJSandMDE_skin <- counts_ADT[,condition %in% c("SJS","MDE") & tissue %in% "Skin"]
condition_SJSandMDE_skin <- condition[condition %in% c("SJS","MDE") & tissue %in% "Skin"]

# Filter low gene counts
max_count <- apply(counts_ADT_SJSandMDE_skin,1,max)
counts_ADT_SJSandMDE_skin <- counts_ADT_SJSandMDE_skin[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_ADT_SJSandMDE_skin,
                              DataFrame(condition_SJSandMDE_skin),
                              design= ~ condition_SJSandMDE_skin)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_ADT_SJSandMDE_Skin.txt")

# Blood
counts_ADT_SJSandMDE_blood <- counts_ADT[,condition %in% c("SJS","MDE") & tissue %in% "PBMC"]
condition_SJSandMDE_blood <- condition[condition %in% c("SJS","MDE") & tissue %in% "PBMC"]

# Filter low gene counts
max_count <- apply(counts_ADT_SJSandMDE_blood,1,max)
counts_ADT_SJSandMDE_blood <- counts_ADT_SJSandMDE_blood[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_ADT_SJSandMDE_blood,
                              DataFrame(condition_SJSandMDE_blood),
                              design= ~ condition_SJSandMDE_blood)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_ADT_SJSandMDE_Blood.txt")

#### 
### DESeq2 on ADT: SJS vs Control, Skin and Blood separately ####
#### 

# Skin
counts_ADT_SJSandCon_skin <- counts_ADT[,condition %in% c("SJS","Control") & tissue %in% "Skin"]
condition_SJSandCon_skin <- condition[condition %in% c("SJS","Control") & tissue %in% "Skin"]

# Filter low gene counts
max_count <- apply(counts_ADT_SJSandCon_skin,1,max)
counts_ADT_SJSandCon_skin <- counts_ADT_SJSandCon_skin[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_ADT_SJSandCon_skin,
                              DataFrame(condition_SJSandCon_skin),
                              design= ~ condition_SJSandCon_skin)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_ADT_SJSandCon_Skin.txt")

# Blood
counts_ADT_SJSandCon_blood <- counts_ADT[,condition %in% c("SJS","Control") & tissue %in% "PBMC"]
condition_SJSandCon_blood <- condition[condition %in% c("SJS","Control") & tissue %in% "PBMC"]

# Filter low gene counts
max_count <- apply(counts_ADT_SJSandCon_blood,1,max)
counts_ADT_SJSandCon_blood <- counts_ADT_SJSandCon_blood[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_ADT_SJSandCon_blood,
                              DataFrame(condition_SJSandCon_blood),
                              design= ~ condition_SJSandCon_blood)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_ADT_SJSandCon_Blood.txt")


#### 
### DESeq2 on ADT: MDE vs Control, Skin and Blood separately ####
#### 

# Skin
counts_ADT_MDEandCon_skin <- counts_ADT[,condition %in% c("MDE","Control") & tissue %in% "Skin"]
condition_MDEandCon_skin <- condition[condition %in% c("MDE","Control") & tissue %in% "Skin"]

# Filter low gene counts
max_count <- apply(counts_ADT_MDEandCon_skin,1,max)
counts_ADT_MDEandCon_skin <- counts_ADT_MDEandCon_skin[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_ADT_MDEandCon_skin,
                              DataFrame(condition_MDEandCon_skin),
                              design= ~ condition_MDEandCon_skin)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_ADT_MDEandCon_Skin.txt")

# Blood
counts_ADT_MDEandCon_blood<- counts_ADT[,condition %in% c("MDE","Control") & tissue %in% "PBMC"]
condition_MDEandCon_blood <- condition[condition %in% c("MDE","Control") & tissue %in% "PBMC"]

# Filter low gene counts
max_count <- apply(counts_ADT_MDEandCon_blood,1,max)
counts_ADT_MDEandCon_blood <- counts_ADT_MDEandCon_blood[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_ADT_MDEandCon_blood,
                              DataFrame(condition_MDEandCon_blood),
                              design= ~ condition_MDEandCon_blood)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_ADT_MDEandCon_Blood.txt")

#### 
### DESeq2 on ADT: Blood vs Skin: SJS and MDE separately ####
#### 

# SJS
counts_ADT_BloodvsSkin_SJS <- counts_ADT[,condition %in% c("SJS")]
tissue_ADT_BloodvsSkin_SJS <- tissue[condition %in% c("SJS")]

# Filter low gene counts
max_count <- apply(counts_ADT_BloodvsSkin_SJS,1,max)
counts_ADT_BloodvsSkin_SJS <- counts_ADT_BloodvsSkin_SJS[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_ADT_BloodvsSkin_SJS,
                              DataFrame(tissue_ADT_BloodvsSkin_SJS),
                              design= ~ tissue_ADT_BloodvsSkin_SJS)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_ADT_BloodvsSkin_SJS.txt")

# MDE
counts_ADT_BloodvsSkin_MDE <- counts_ADT[,condition %in% c("MDE")]
tissue_ADT_BloodvsSkin_MDE <- tissue[condition %in% c("MDE")]

# Filter low gene counts
max_count <- apply(counts_ADT_BloodvsSkin_MDE,1,max)
counts_ADT_BloodvsSkin_MDE <- counts_ADT_BloodvsSkin_MDE[which(max_count > 10),]

# test
dds <- DESeqDataSetFromMatrix(countData = counts_ADT_BloodvsSkin_MDE,
                              DataFrame(tissue_ADT_BloodvsSkin_MDE),
                              design= ~ tissue_ADT_BloodvsSkin_MDE)

# maximum might be better
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.table(as.data.frame(res), file = "counts_ADT_BloodvsSkin_MDE.txt")