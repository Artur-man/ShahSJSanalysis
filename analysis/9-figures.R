library(Seurat)
library(tidyverse)
source('plot_violin_seurat.R')

bri.integrated <- readRDS("final-multimodal-sjs-clustered-tcr-annotated.rds")
DefaultAssay(bri.integrated) <- "RNA"

## Configurations ####
####

pdfwidth = 8
pdfheight = 8

color_scheme <- list(`CD8+ Naive`="#c9a0dc",
                     `CD8+ TCM`="#87ff2a",
                     `CD8+ TMM`="#000f89",
                     `CD8+ CD103- CD69+ TRM (nf)`="#8b5a2b",
                     `CD8+ CD103- TEMRA`="#00f5ff",
                     `CD8+ CD103+ TEMRA`="#fa8072",
                     `CD8+ T Effector`="#ffdf00",
                     `CD8+ TEM`="#008b00",
                     `CD8+ CD103+ CD69+ TRM`="#ff0800",
                     `CD8+ CD103- CD69+ TRM`="#858585",
                     `CD8+ CD56+ T cells`="#ff1dce",
                     `GD T cells`="#03c03c",
                     `CD4+ Naive`="#ae2029",
                     `CD4+ CD45RA+ CD62L-low`="#9400d3",
                     `CD4+ CD45RA- CD62L-low`="#ff5800",
                     `CD4+ TCM`="#0000ff",
                     `CD4+ CD103- CD69+ TRM`="#5f9ea0",
                     `CD4+ CD103+ CD69+ TRM`="#fe28a2",
                     `CD4+ T Effector`="#9f00ff",
                     `CD4+ Treg 1`="#cd1076",
                     `CD4+ Treg 2`="#0abab5",
                     `Proliferating`="#cb4154")


#### Figure 2 ###

#CD4 / CD8
cd4cells <- WhichCells(bri.integrated,
                       cells = rownames(bri.integrated@meta.data[grep("CD4|Treg", bri.integrated$CellType),]))
cd8cells <- WhichCells(bri.integrated,
                       cells = rownames(bri.integrated@meta.data[grep("CD8", bri.integrated$CellType),]))

DimPlot(bri.integrated, label=F, reduction="wnn.umap", label.size = 5, sizes.highlight = 3, 
        cells.highlight = cd4cells, 
        cols.highlight = color_scheme[7]) + labs(title = "CD4 T Cells") + NoLegend()

DimPlot(bri.integrated, label=F, reduction="wnn.umap", label.size = 5, sizes.highlight = 3, 
        cells.highlight = cd8cells, 
        cols.highlight = color_scheme[8]) + labs(title = "CD8 T Cells") + NoLegend()

#2B
#heatmap of markers

#2C and D TCR Analysis




#### Figure 3 ####

bri.integrated.sub <- subset(bri.integrated, subset= CellType %in% c('CD8+ CD103+ TEMRA', 'CD8+ CD103- TEMRA', 'CD8+ T Effector', 
                                                                     'CD8+ TEM', 'CD8+ CD103+ CD69+ TRM', 'CD8+ CD103- CD69+ TRM', 'CD8+ CD56+ T cells'))
#Combine patients into disease groups
plotgene <- 'GNLY'
plotgene <- 'GZMB'
plotgene <- 'GZMA'
plotgene <- 'PRF1'
bri.integrated.sub.skin <- subset(bri.integrated.pranalisub, subset = Tissue=='Skin')
bri.integrated.sub.skin@meta.data$CellType <- 
  factor(bri.integrated.sub.skin@meta.data$CellType, 
         levels = c('CD8+ CD103- TEMRA','CD8+ CD103+ TEMRA','CD8+ T Effector',
                    'CD8+ TEM','CD8+ CD103+ CD69+ TRM',
                    'CD8+ CD103- CD69+ TRM', 'CD8+ CD56+ T cells'))
bri.integrated.sub.skin@meta.data$CellType <- bri.integrated.sub.skin$CellType
plot_violin_seurat(bri.integrated.sub.skin, gene=plotgene, color_by = 'CellType',
                   facet_by = 'phenotypeTissue',
                   facet_order = c('SJS_Skin','MDE_Skin','Control_Skin'),
                   colors = c('#00f5ff','#fa8072','#ffdf00','#008b00','#ff0800','#858585','#ff1dce'),
                   title_size = 20, axis_title_size = 15, axis_text_size = 5, legend_title_size = 10, 
                   legend_text_size = 8, facet_text_size = 10, number_label_text_size = 4)

#CCR8 and TNFSF4 in Tregs
bri.integrated.treg <- subset(bri.integrated, subset= CellType %in% c('CD4+ Treg 1', 'CD4+ Treg 2'))
bri.integrated.treg

# Assigng tissue type, specifically which type of PBMC or Skin
# care with the names - is case sensitive
diseaseType <- bri.integrated.treg@meta.data$phenotypeTissue
diseaseType[grepl('SJS',diseaseType)] <- 'SJS'
diseaseType[grepl('MDE',diseaseType)] <- 'MDE'
diseaseType[grepl('Control',diseaseType)] <- 'Control'
diseaseType[grepl('BCH',diseaseType)] <- 'Control'
bri.integrated.treg@meta.data$diseaseType <- diseaseType

bri.integrated.treg@meta.data$CellType <- bri.integrated.treg$CellTypePranali
#set order by setting factor levels
bri.integrated.treg$diseaseType <- 
  factor(bri.integrated.ptreg$diseaseType, levels = c('SJS','MDE','Control'))

plotgene <- 'CCR8'
plotgene <- 'TNFRSF4'
plot_violin_seurat(bri.integrated.treg, gene=plotgene, color_by = 'CellType',
                   facet_by = c('diseaseType','Tissue'),
                   color_order = c('CD4+ Treg 1', 'CD4+ Treg 2'),
                   colors = c('#cd1076','#0abab5'),
                   title_size = 20, axis_title_size = 15, axis_text_size = 10, legend_title_size = 10, 
                   legend_text_size = 8, facet_text_size = 10, number_label_text_size = 4)

#figure for IFNG violinplot across all clusters (except CD4 and CD8 naive and 2 Tregs) across skin phenotypes
bri.integrated.ifng <- 
  subset(bri.integrated, subset= CellType %in% 
           c('CD8+ TCM', 'CD8+ TMM', 'CD8+ CD103- CD69+ TRM (nf)',
             'CD8+ CD103- TEMRA', 'CD8+ CD103+ TEMRA', 'CD8+ T Effector', 'CD8+ TEM',
             'CD8+ CD103+ CD69+ TRM', 'CD8+ CD103- CD69+ TRM', 'CD8+ CD56+ T cells',
             'GD T cells', 'CD4+ CD45RA+ CD62L-low', 'CD4+ CD45RA- CD62L-low',
             'CD4+ TCM', 'CD4+ CD103- CD69+ TRM', 'CD4+ CD103+ CD69+ TRM',
             'CD4+ T Effector', 'Proliferating'))
bri.integrated.ifng@meta.data$CellType <- 
  factor(bri.integrated.ifng@meta.data$CellType, 
         levels = c('CD8+ TCM', 'CD8+ TMM', 'CD8+ CD103- CD69+ TRM (nf)',
                    'CD8+ CD103- TEMRA', 'CD8+ CD103+ TEMRA', 'CD8+ T Effector', 'CD8+ TEM',
                    'CD8+ CD103+ CD69+ TRM', 'CD8+ CD103- CD69+ TRM', 'CD8+ CD56+ T cells',
                    'GD T cells', 'CD4+ CD45RA+ CD62L-low', 'CD4+ CD45RA- CD62L-low',
                    'CD4+ TCM', 'CD4+ CD103- CD69+ TRM', 'CD4+ CD103+ CD69+ TRM',
                    'CD4+ T Effector', 'Proliferating'))
bri.integrated.ifng.skin <- subset(bri.integrated.ifng, subset= Tissue=='Skin')

plot_violin_seurat(bri.integrated.ifng.skin, gene='IFNG', color_by = 'CellType',
                   facet_by = 'phenotypeTissue', sig = 2,
                   facet_order = c('SJS_Skin','MDE_Skin','Control_Skin'),
                   colors = c('#87ff2a','#000f89','#8b5a2b','#00f5ff',
                              '#fa8072','#ffdf00','#008b00','#ff0800','#858585',
                              '#ff1dce','#03c03c','#9400d3','#ff5800',
                              '#0000ff','#5f9ea0','#fe28a2','#9f00ff','#cb4154'),
                   title_size = 20, axis_title_size = 15, axis_text_size = 12, legend_title_size = 10, 
                   legend_text_size = 8, facet_text_size = 10, number_label_text_size = 3)


#### Supplemental Files ####
#Supp Fig 3 - Final annotated object
DimPlot(bri.integrated, pt.size=1.5, cols=color_scheme, order = rev(names(color_scheme)), label = TRUE) + theme(legend.position = "bottom") + NoAxes() # For legend box and labels to write manually on NoLegend/Nolabel plot below
#Supp Fig 5 - skin/blood
DimPlot(bri.integrated, group.by = 'Tissue', cols=c('red','blue'), order = c('PBMC','Skin'), label = F) + NoAxes() + ggtitle('All samples') 

DimPlot(subset(bri.integrated, subset= hash.ID %in% c('SJS001SKIN','SJS001PBMC')), group.by = 'Tissue', cols=c('red','blue'), order = c('Skin','PBMC'), label = F) + NoAxes() + ggtitle("SJS001")
DimPlot(subset(bri.integrated, subset= hash.ID %in% c('SJS002SKIN','SJS002PBMC')), group.by = 'Tissue', cols=c('red','blue'), order = c('Skin','PBMC'), label = F) + NoAxes() + ggtitle("SJS002")
DimPlot(subset(bri.integrated, subset= hash.ID %in% c('SJS003SKIN','SJS003PBMC')), group.by = 'Tissue', cols=c('red','blue'), order = c('Skin','PBMC'), label = F) + NoAxes() + ggtitle("SJS003")
DimPlot(subset(bri.integrated, subset= hash.ID %in% c('MDE001SKIN','MDE001PBMC')), group.by = 'Tissue', cols=c('red','blue'), order = c('Skin','PBMC'), label = F) + NoAxes() + ggtitle("MDE001")
DimPlot(subset(bri.integrated, subset= hash.ID %in% c('MDE002SKIN','MDE002PBMC')), group.by = 'Tissue', cols=c('red','blue'), order = c('Skin','PBMC'), label = F) + NoAxes() + ggtitle("MDE002")
DimPlot(subset(bri.integrated, subset= hash.ID %in% c('MDE-006-Skin','MDE-006-PBMC')), group.by = 'Tissue', cols=c('red','blue'), order = c('Skin','PBMC'), label = F) + NoAxes() + ggtitle("MDE006")

DimPlot(subset(bri.integrated, subset= hash.ID %in% c('BCH026PBMC')), group.by = 'Tissue', cols=c('red'), label = F) + NoAxes() + ggtitle("BCH026PBMC")
DimPlot(subset(bri.integrated, subset= hash.ID %in% c('Control-BCH-001-PBMC')), group.by = 'Tissue', cols=c('red'), label = F) + NoAxes() + ggtitle("Control-BCH-001-PBMC")
DimPlot(subset(bri.integrated, subset= hash.ID %in% c('Control-BCH-048-PBMC')), group.by = 'Tissue', cols=c('red'), label = F) + NoAxes() + ggtitle("Control-BCH-048-PBMC")

DimPlot(subset(bri.integrated, subset= hash.ID %in% c('Control-039-Skin')), group.by = 'Tissue', cols=c('blue'), label = F) + NoAxes() + ggtitle("Control-039-Skin")
DimPlot(subset(bri.integrated, subset= hash.ID %in% c('Control031Skin')), group.by = 'Tissue', cols=c('blue'), label = F) + NoAxes() + ggtitle("Control031Skin")

#Individual plots of skin, blood, integrated
DimPlot(subset(bri.integrated, subset= hash.ID=="Control-039-Skin"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend() + NoAxes() + ggtitle("Control-039-Skin")
DimPlot(subset(bri.integrated, subset= hash.ID=="Control-BCH-001-PBMC"),pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend() + NoAxes() + ggtitle("Control-BCH-001-PBMC")
DimPlot(subset(bri.integrated, subset= hash.ID=="Control-BCH-048-PBMC"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("Control-BCH-048-PBMC")
DimPlot(subset(bri.integrated, subset= hash.ID=="BCH026PBMC"),pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend() + NoAxes() + ggtitle("BCH026PBMC")
DimPlot(subset(bri.integrated, subset= hash.ID=="Control031Skin"),pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("Control031Skin")
DimPlot(subset(bri.integrated, subset= hash.ID=="MDE-006-PBMC"),pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("MDE-006-PBMC")
DimPlot(subset(bri.integrated, subset= hash.ID=="MDE-006-Skin"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("MDE-006-Skin")
DimPlot(subset(bri.integrated, subset= hash.ID=="MDE001PBMC"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("MDE001PBMC")
DimPlot(subset(bri.integrated, subset= hash.ID=="MDE001SKIN"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("MDE001SKIN")
DimPlot(subset(bri.integrated, subset= hash.ID=="MDE002PBMC"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("MDE002PBMC")
DimPlot(subset(bri.integrated, subset= hash.ID=="MDE002SKIN"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("MDE002SKIN")
DimPlot(subset(bri.integrated, subset= hash.ID=="SJS001PBMC"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("SJS001PBMC")
DimPlot(subset(bri.integrated, subset= hash.ID=="SJS001SKIN"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("SJS001SKIN")
DimPlot(subset(bri.integrated, subset= hash.ID=="SJS002PBMC"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("SJS002PBMC")
DimPlot(subset(bri.integrated, subset= hash.ID=="SJS002SKIN"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("SJS002SKIN")
DimPlot(subset(bri.integrated, subset= hash.ID=="SJS003PBMC"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("SJS003PBMC")
DimPlot(subset(bri.integrated, subset= hash.ID=="SJS003SKIN"), pt.size = 2, cols=color_scheme, order = rev(names(color_scheme)), label = F) + NoLegend()+ NoAxes() + ggtitle("SJS003SKIN")


## Supp Fig 6
# use the independently clustered objects from script 8-pbmc-skin-independent-clustering.R
bri.integrated.skin <- readRDS("8-multimodal-sjs-skin-only.rds")
bri.integrated.pbmc_nonmt <- readRDS("8-multimodal-sjs-pbmc-only.rds")

bri.integrated.skin <- bri.integrated.skin[,bri.integrated.skin$CellTypePranali != 'NA']
bri.integrated.skin <- bri.integrated.skin[,bri.integrated.skin$CellTypePranali != 'Dead']
DimPlot(bri.integrated.skin, reduction = 'wnn.umap', label = TRUE, cols=color_scheme, order = rev(names(color_scheme)), label.size = 5, pt.size = 1, group.by = "CellType", repel = TRUE) + labs(title = "Cell Type Annotations (Integrated)")
DimPlot(bri.integrated.skin, reduction = 'wnn.umap', label = TRUE, label.size = 5, pt.size = 1, group.by = "wsnn_res.1.6", repel = TRUE) +NoAxes() + labs(title = "Cell Type Annotations (Independent clusters)")

bri.integrated.pbmc_nonmt <- bri.integrated.pbmc_nonmt[,bri.integrated.pbmc_nonmt$CellTypePranali != 'NA']
bri.integrated.pbmc_nonmt <- bri.integrated.pbmc_nonmt[,bri.integrated.pbmc_nonmt$CellTypePranali != 'Dead']
DimPlot(bri.integrated.pbmc_nonmt, reduction = 'wnn.umap', label = TRUE, cols=color_scheme, order = rev(names(color_scheme)), label.size = 5, pt.size = 1, group.by = "CellType", repel = TRUE) + labs(title = "Cell Type Annotations (Integrated)")
DimPlot(bri.integrated.pbmc_nonmt, reduction = 'wnn.umap', label = TRUE, label.size = 5, pt.size = 1, group.by = "wsnn_res.1.2", repel = TRUE) + NoAxes() + labs(title = "Cell Type Annotations (Independent)")


####
#### figure highlighting enrichment of cell types per blood type ####
#######

DimPlot(bri.integrated, pt.size=1, cols=color_scheme, order = rev(names(color_scheme)), label = TRUE) + theme(legend.position = "bottom") + NoAxes() # For legend box and labels to write manually on NoLegend/Nolabel plot below

length(color_scheme)
rep('gray',22)

DimPlot(subset(bri.integrated, subset = Tissue=='PBMC'), pt.size=1, cols=rep('gray',22), order = rev(names(color_scheme)), label = TRUE) + theme(legend.position = "bottom") + NoAxes() + ggtitle('Blood subset')
DimPlot(subset(bri.integrated, subset = Tissue=='Skin'), pt.size=1, cols=rep('gray',22), order = rev(names(color_scheme)), label = TRUE) + theme(legend.position = "bottom") + NoAxes() + ggtitle('Skin subset')

#### Heatmap for Supp Fig 6 ####

###


#### Supp Fig 7 ####
# clonal frequency and % of top clone in HC skin




### SUPP fig, which one for violin plots? Do we still supply it?

# GNLY and GZMB for SKIN for each individual patient SJS123, MDE123, HC1,2
#cd8+cd103-temra; cd8+cd103+temra;cd8+teff;cd8+tem;cd8+cd103+cd69+trm;cd8+cd103-cd69+trm;cd8+cd56+tcell
plotgene <- 'GNLY'
plot_violin_seurat(bri.integrated.sub.skin, gene=plotgene, color_by = 'CellType',
                   facet_by = 'hash.ID', sig = 2,
                   facet_order = c('SJS001SKIN','SJS002SKIN','SJS003SKIN','MDE001SKIN','MDE002SKIN','MDE-006-Skin','Control031Skin','Control-039-Skin'),
                   colors = c('#00f5ff','#fa8072','#ffdf00','#008b00','#ff0800','#858585','#ff1dce'),
                   title_size = 20, axis_title_size = 15, axis_text_size = 5, legend_title_size = 10, 
                   legend_text_size = 8, facet_text_size = 10, number_label_text_size = 4)







# RNA, ADT, WNN umaps
DimPlot(bri.integrated, group.by = "wsnn_res.0.6", reduction = "rnaintegrated.umap", label = TRUE, label.box = TRUE, repel = FALSE) + NoLegend()
DimPlot(bri.integrated, group.by = "wsnn_res.0.6", reduction = "adtintegrated.umap", label = TRUE, label.box = TRUE) + NoLegend()
DimPlot(bri.integrated, group.by = "CellType", reduction = "wnn.umap", label = TRUE, label.box = TRUE) + NoLegend()

DimPlot(bri.integrated, group.by = "CellType", reduction = "rnaintegrated.umap", label = TRUE, label.box = TRUE, repel = TRUE) + NoLegend()
DimPlot(bri.integrated, group.by = "CellType", reduction = "adtintegrated.umap", label = TRUE, label.box = TRUE) + NoLegend()
DimPlot(bri.integrated, group.by = "CellType", reduction = "wnn.umap", label = TRUE, label.box = TRUE) + NoLegend()

DimPlot(bri.integrated, group.by = "CellTypeNumber", reduction = "wnn.umap", label = TRUE, label.box = TRUE) + NoLegend()

#by hashID (all), and by individual samples
DimPlot(subset(bri.integrated, subset= Tissue=="PBMC"), group.by = "CellType", pt.size = .5, label = TRUE, label.size = 8, repel = TRUE) + ggtitle("All Blood") + NoLegend()
DimPlot(subset(bri.integrated, subset= Tissue=="Skin"), group.by = "CellType", pt.size = .5, label = TRUE, label.size = 8, repel = TRUE) + ggtitle("All Skin") + NoLegend()

DimPlot(subset(bri.integrated, subset= hash.ID=="BCH026PBMC"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("BCH026PBMC") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="Control-039-Skin"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("Control-039-Skin") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="Control-BCH-001-PBMC"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("Control-BCH-001-PBMC") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="Control-BCH-048-PBMC"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("Control-BCH-048-PBMC") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="Control031Skin"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("Control031Skin") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="MDE-006-PBMC"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("MDE-006-PBMC") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="MDE-006-Skin"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("MDE-006-Skin") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="MDE001PBMC"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("MDE001PBMC") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="MDE001SKIN"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("MDE001SKIN") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="MDE002PBMC"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("MDE002PBMC") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="MDE002SKIN"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("MDE002SKIN") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="SJS001PBMC"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("SJS001PBMC") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="SJS001SKIN"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("SJS001SKIN") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="SJS002PBMC"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("SJS002PBMC") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="SJS002SKIN"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("SJS002SKIN") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="SJS003PBMC"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE,) + ggtitle("SJS003PBMC") + NoLegend()
DimPlot(subset(bri.integrated, subset= hash.ID=="SJS003SKIN"), group.by = "CellType", pt.size = 1.5, label = TRUE, label.size = 5, repel = TRUE, label.box = TRUE) + ggtitle("SJS003SKIN") + NoLegend()



