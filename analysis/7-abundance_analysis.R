# libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(anndata)
library(scales)
library(ggsci) # palatte for NEJM
library(circlize)
library(reshape2)
library(xlsx)
library(speckle)
library(edgeR)
library(xlsx)

####
## Import and Configure Data #####
####

bri.integrated <- readRDS("sjs_clustered_tcr_annotated.rds")
bri.integrated$phenotype <- sapply(bri.integrated$phenotypeTissue, function(x) return(strsplit(x, split = "_")[[1]][1]))
bri.integrated$phenotype <- gsub("BCH","Control", bri.integrated$phenotype)
bri.integrated <- bri.integrated[,!bri.integrated$CellType == "Dying Cells (MT++)"]

#### 
## Get Cytotoxic Cells ####
####

# find cutvalues for 
cytotoxic_markers <- c("NKG7","GZMA", "GZMB", "GNLY", "PRF1")
hashes <- c("SJS001SKIN","SJS001PBMC","SJS002SKIN","SJS002PBMC","SJS003SKIN","SJS003PBMC",
            "MDE001SKIN","MDE001PBMC","MDE002SKIN","MDE002PBMC","MDE-006-Skin","MDE-006-PBMC","Control-BCH-048-PBMC","Control-BCH-001-PBMC","Control-039-Skin","Control031Skin","BCH026PBMC")

rna_counts <- bri.integrated@assays$integratedRNA@data
agg_cytotoxic <- colMeans(rna_counts[cytotoxic_markers,])
bri.integrated$cytotoxic_expression <- agg_cytotoxic
DefaultAssay(bri.integrated) <- "integratedRNA"
VlnPlot(bri.integrated, features = "NKG7", group.by = "CellType") + NoLegend() + geom_hline(yintercept = 2, linetype='dashed')
VlnPlot(bri.integrated, features = "GZMA", group.by = "CellType") + NoLegend() + geom_hline(yintercept = 1, linetype='dashed')
VlnPlot(bri.integrated, features = "GZMB", group.by = "CellType") + NoLegend() + geom_hline(yintercept = 0.5, linetype='dashed')
VlnPlot(bri.integrated, features = "GNLY", group.by = "CellType") + NoLegend() + geom_hline(yintercept = 2, linetype='dashed')
VlnPlot(bri.integrated, features = "PRF1", group.by = "CellType") + NoLegend() + geom_hline(yintercept = 0.5, linetype='dashed')
VlnPlot(bri.integrated, features = "cytotoxic_expression", group.by = "CellType") + NoLegend() + geom_hline(yintercept = 1.2, linetype='dashed', size = 1.5) + labs(title = "Average Cytotoxic Expression")
VlnPlot(bri.integrated, features = "cytotoxic_expression", group.by = "orig.ident") + NoLegend() + geom_hline(yintercept = 1.2, linetype='dashed') 

# label cytotoxic cells
bri.integrated$cytotoxic <- ifelse(bri.integrated$cytotoxic_expression > 1.2, "Cytotoxic", "Non-Cytotoxic")
bri.integrated$cytotoxic_CD4CD8 <- paste(bri.integrated$cytotoxic, bri.integrated$CD4_CD8_status, sep = "_")
bri.integrated$cytotoxic_clusters <- paste(bri.integrated$cytotoxic, bri.integrated$CellType, sep = "_")

table_cyto <- as.matrix(rbind(table(bri.integrated$cytotoxic, bri.integrated$hash.ID)))
table_cyto <- apply(table_cyto, 1, function(x){
  x/colSums(table_cyto)
})

#### 
## Propeller Abundance test with ANOVA ####
####

# all clusters
bri.integrated_test <- subset(bri.integrated, subset = Tissue == "Skin")
bri.integrated_test <- bri.integrated_test[,!bri.integrated_test$CellType %in% c("GD T cells", "CD8+ CD103+ TEMRA")]
prop.list <- getTransformedProps(clusters = droplevels(bri.integrated_test$CellType), 
                                 sample = bri.integrated_test$hash.ID, transform = "asin")
grp <- data.frame(cbind(bri.integrated_test$phenotype, bri.integrated_test$hash.ID)) %>% distinct()
rownames(grp) <- grp$X2
grp <- grp[colnames(prop.list$Proportions),"X1"]
design <- model.matrix(~0+grp)
results <- propeller.anova(prop.list = prop.list, design = design, coef = c(1,2,3), trend = FALSE, sort = TRUE, robust=TRUE)
write.xlsx(results, file = "differentialabundance_anova.xlsx", sheetName = "AllClusters_Skin")

# Skin
bri.integrated_test <- subset(bri.integrated, subset = Tissue == "Skin")
bri.integrated_test$cytotoxic_CD4CD8 <- ifelse(grepl("CD4|CD8", bri.integrated_test$CD4_CD8_status),
                                               paste(bri.integrated_test$cytotoxic, bri.integrated_test$CD4_CD8_status, sep = "_"),
                                               as.character(bri.integrated_test$CellType))
bri.integrated_test <- bri.integrated_test[,!bri.integrated_test$cytotoxic_CD4CD8 %in% c("GD T cells")]
prop.list <- getTransformedProps(clusters = bri.integrated_test$cytotoxic_CD4CD8, 
                                 sample = bri.integrated_test$hash.ID, transform = "asin")
grp <- data.frame(cbind(bri.integrated_test$phenotype, bri.integrated_test$hash.ID)) %>% distinct()
rownames(grp) <- grp$X2
grp <- grp[colnames(prop.list$Proportions),"X1"]
design <- model.matrix(~0+grp)
results <- propeller.anova(prop.list = prop.list, design = design, coef = c(1,2,3), trend = FALSE, sort = TRUE, robust=TRUE)
results <- results[grepl("Cytotoxic", rownames(results)),]
write.xlsx(results, file = "differentialabundance_anova.xlsx", sheetName = "Cytotoxic_Skin", append = TRUE)

# all clusters
bri.integrated_test <- subset(bri.integrated, subset = Tissue == "PBMC")
bri.integrated_test <- bri.integrated_test[,!bri.integrated_test$CellType %in% c("GD T cells", "CD8+ CD103+ TEMRA")]
prop.list <- getTransformedProps(clusters = droplevels(bri.integrated_test$CellType), 
                                 sample = bri.integrated_test$hash.ID, transform = "asin")
grp <- data.frame(cbind(bri.integrated_test$phenotype, bri.integrated_test$hash.ID)) %>% distinct()
rownames(grp) <- grp$X2
grp <- grp[colnames(prop.list$Proportions),"X1"]
design <- model.matrix(~0+grp)
results <- propeller.anova(prop.list = prop.list, design = design, coef = c(1,2,3), trend = FALSE, sort = TRUE, robust=TRUE)
write.xlsx(results, file = "differentialabundance_anova.xlsx", sheetName = "AllClusters_PBMC", append = TRUE)

# PBMC
bri.integrated_test <- subset(bri.integrated, subset = Tissue == "PBMC")
bri.integrated_test$cytotoxic_CD4CD8 <- ifelse(grepl("CD4|CD8", bri.integrated_test$CD4_CD8_status),
                                               paste(bri.integrated_test$cytotoxic, bri.integrated_test$CD4_CD8_status, sep = "_"),
                                               as.character(bri.integrated_test$CellType))
bri.integrated_test <- bri.integrated_test[,!bri.integrated_test$cytotoxic_CD4CD8 %in% c("GD T cells")]
prop.list <- getTransformedProps(clusters = bri.integrated_test$cytotoxic_CD4CD8, 
                                 sample = bri.integrated_test$hash.ID, transform = "asin")
grp <- data.frame(cbind(bri.integrated_test$phenotype, bri.integrated_test$hash.ID)) %>% distinct()
rownames(grp) <- grp$X2
grp <- grp[colnames(prop.list$Proportions),"X1"]
design <- model.matrix(~0+grp)
results <- propeller.anova(prop.list = prop.list, design = design, coef = c(1,2,3), trend = FALSE, sort = TRUE, robust=TRUE)
results <- results[grepl("Cytotoxic", rownames(results)),]
write.xlsx(results, file = "differentialabundance_anova.xlsx", sheetName = "Cytotoxic_PBMC", append = TRUE)

#### 
## Pairwise (Post hoc) Propeller t-test for CD4 Treg 2 and Cytotoxic CD8 T-cells ####
####

# comparisons
comparison <- c("MDE|SJS","Control|MDE","Control|SJS")

# Skin
alltest <- list()
for(i in 1:length(comparison)){
  bri.integrated_test <- subset(bri.integrated, subset = Tissue == "Skin")
  bri.integrated_test <- bri.integrated_test[,!bri.integrated_test$CellType %in% c("GD T cells","CD8+ CD103+ TEMRA")]
  bri.integrated_test <- bri.integrated_test[,grepl(comparison[i], bri.integrated_test$phenotype)]
  prop.list <- getTransformedProps(clusters = droplevels(bri.integrated_test$CellType), 
                                   sample = bri.integrated_test$hash.ID, transform = "asin")
  grp <- data.frame(cbind(bri.integrated_test$phenotype, bri.integrated_test$hash.ID)) %>% distinct()
  rownames(grp) <- grp$X2
  grp <- grp[colnames(prop.list$Proportions),"X1"]
  design <- model.matrix(~0+grp)
  contrasts <- c(-1,1)
  alltest[[i]] <- propeller.ttest(prop.list = prop.list, design = design, contrasts = contrasts, trend = FALSE, sort = TRUE, robust=TRUE) 
  alltest[[i]] <- data.frame(CellType = rownames(alltest[[i]]), Comparison = comparison[i], alltest[[i]], row.names = NULL)
  colnames(alltest[[i]])[grepl("PropMean", colnames(alltest[[i]]))] <- c("PropMean.1", "PropMean.2")
}
alltest <- do.call(rbind, alltest)
alltest <- alltest[alltest$CellType == "CD4+ Treg 2",]
alltest$FDR <- p.adjust(alltest$P.Value, method = "BH")
write.xlsx(alltest, file = "differentialabundance_pairwise.xlsx", sheetName = "CD4+ Treg 2 Skin", row.names = FALSE)

# Skin
alltest <- list()
for(i in 1:length(comparison)){
  bri.integrated_test <- subset(bri.integrated, subset = Tissue == "Skin")
  bri.integrated_test <- bri.integrated_test[,grepl(comparison[i], bri.integrated_test$phenotype)]
  bri.integrated_test$cytotoxic_CD4CD8 <- ifelse(grepl("CD4|CD8", bri.integrated_test$CD4_CD8_status),
                                                 paste(bri.integrated_test$cytotoxic, bri.integrated_test$CD4_CD8_status, sep = "_"),
                                                 as.character(bri.integrated_test$CellType))
  bri.integrated_test <- bri.integrated_test[,!bri.integrated_test$cytotoxic_CD4CD8 %in% "GD T cells"]
  prop.list <- getTransformedProps(clusters = bri.integrated_test$cytotoxic_CD4CD8, 
                                   sample = bri.integrated_test$hash.ID, transform = "asin")
  grp <- data.frame(cbind(bri.integrated_test$phenotype, bri.integrated_test$hash.ID)) %>% distinct()
  rownames(grp) <- grp$X2
  grp <- grp[colnames(prop.list$Proportions),"X1"]
  design <- model.matrix(~0+grp)
  contrasts <- c(-1,1)
  alltest[[i]] <- propeller.ttest(prop.list = prop.list, design = design, contrasts = contrasts, trend = FALSE, sort = TRUE, robust=TRUE) 
  alltest[[i]] <- data.frame(CellType = rownames(alltest[[i]]), Comparison = comparison[i], alltest[[i]], row.names = NULL)
  colnames(alltest[[i]])[grepl("PropMean", colnames(alltest[[i]]))] <- c("PropMean.1", "PropMean.2")
}
alltest <- do.call(rbind, alltest)
alltest_CD8 <- alltest[alltest$CellType == "Cytotoxic_CD8",]
alltest_CD8$FDR <- p.adjust(alltest_CD8$P.Value, method = "BH")
write.xlsx(alltest_CD8, file = "differentialabundance_pairwise.xlsx", sheetName = "Cytotoxic CD8 Skin", append = TRUE, row.names = FALSE)

# PBMC
alltest <- list()
for(i in 1:length(comparison)){
  bri.integrated_test <- subset(bri.integrated, subset = Tissue == "PBMC")
  bri.integrated_test <- bri.integrated_test[,grepl(comparison[i], bri.integrated_test$phenotype)]
  bri.integrated_test$cytotoxic_CD4CD8 <- ifelse(grepl("CD4|CD8", bri.integrated_test$CD4_CD8_status),
                                                 paste(bri.integrated_test$cytotoxic, bri.integrated_test$CD4_CD8_status, sep = "_"),
                                                 as.character(bri.integrated_test$CellType))
  bri.integrated_test <- bri.integrated_test[,!bri.integrated_test$cytotoxic_CD4CD8 %in% "GD T cells"]
  prop.list <- getTransformedProps(clusters = bri.integrated_test$cytotoxic_CD4CD8, 
                                   sample = bri.integrated_test$hash.ID, transform = "asin")
  grp <- data.frame(cbind(bri.integrated_test$phenotype, bri.integrated_test$hash.ID)) %>% distinct()
  rownames(grp) <- grp$X2
  grp <- grp[colnames(prop.list$Proportions),"X1"]
  design <- model.matrix(~0+grp)
  contrasts <- c(-1,1)
  alltest[[i]] <- propeller.ttest(prop.list = prop.list, design = design, contrasts = contrasts, trend = FALSE, sort = TRUE, robust=TRUE) 
  alltest[[i]] <- data.frame(CellType = rownames(alltest[[i]]), Comparison = comparison[i], alltest[[i]], row.names = NULL)
  colnames(alltest[[i]])[grepl("PropMean", colnames(alltest[[i]]))] <- c("PropMean.1", "PropMean.2")
}
alltest <- do.call(rbind, alltest)
alltest_CD8 <- alltest[alltest$CellType == "Cytotoxic_CD8",]
alltest_CD8$FDR <- p.adjust(alltest_CD8$P.Value, method = "BH")
write.xlsx(alltest_CD8, file = "differentialabundance_pairwise.xlsx", sheetName = "Cytotoxic CD8 PBMC", append = TRUE, row.names = FALSE)