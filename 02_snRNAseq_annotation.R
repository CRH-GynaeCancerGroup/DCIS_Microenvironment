########################
# snRNA-seq annotation #
########################

## This script provides the code used to annotate the snRNA-seq dataset, it follows the preprocessing script.


library(Seurat)
library(copykat)
set.seed(42)


#############################
# Annotate large cell types #
#############################

myeloidMarkers <- c("AIF1","CD68","CD163","CSF1R","CD14","ITGAM","ITGAX","FCGR3A","SIGLEC1")
lymphocyteMarkers <- c("CD3D","CD3G","CD3E","MS4A1","CD79A","CD19","NKG7","GZMB","CD8A")
endothelialMarkers <- c("PECAM1","CD34")
fibroblastMarkers <- c("ACTA2","FAP","PFGFRB","COL1A1","VIM","COL1A2","RGS5")
epithelialMarkers <- c("EPCAM","CDH1","SOX9","KRT19","KRT5","KRT6A")


VlnPlot(DCISSeurat,myeloidMarkers,ncol = 3, pt.size = 0)
VlnPlot(DCISSeurat,lymphocyteMarkers,ncol = 3, pt.size = 0) 
VlnPlot(DCISSeurat,endothelialMarkers,ncol = 3, pt.size = 0) 
VlnPlot(DCISSeurat,fibroblastMarkers,ncol = 3, pt.size = 0) 
VlnPlot(DCISSeurat,epithelialMarkers,ncol = 3, pt.size = 0) 
VlnPlot(DCISSeurat,c("SNAI1","SNAI2"),ncol = 3, pt.size = 0) 

DimPlot(DCISSeurat,label = TRUE, reduction = "umap.harmony")

largeCellType <- c("0" = "Fibroblast-pericyte", "1" = "Epithelial", "2" = "Epithelial", "3" = "Fibroblast-pericyte", "4" = "Epithelial",
                   "5" = "Epithelial", "6" = "Lymphocyte", "7" = "Lymphocyte", "8" = "Endothelial", "9" = "Lymphocyte",
                   "10" = "Myeloid", "11" = "Fibroblast-pericyte", "12" = "Fibroblast-pericyte", "13" = "Myeloid", "14" = "Epithelial",
                   "15" = "Epithelial", "16" = "Myeloid", "17" = "Lymphocyte", "18" = "Lymphocyte", "19" = "Epithelial",
                   "20" = "Fibroblast-pericyte", "21" = "Myeloid", "22" = "Endothelial", "23" = "Lymphocyte", "24" = "Lymphocyte")


newMeta <- data.frame(barcodes = colnames(DCISSeurat),large_cell_type = largeCellType[as.character(Idents(DCISSeurat))])
rownames(newMeta) = newMeta$barcodes

DCISSeurat <- AddMetaData(DCISSeurat,newMeta)

newMeta <- data.frame(barcodes = colnames(DCISSeurat),large_cell_type = DCISSeurat$large_cell_type)
rownames(newMeta) = newMeta$barcodes
newMeta$fine_cell_type <- NA
newMeta$subset <- NA

###############
# Lymphocytes #
###############

lymphocytes_DCIS <- subset(DCISSeurat, large_cell_type == "Lymphocyte")

lymphocytes_DCIS <- NormalizeData(lymphocytes_DCIS)
lymphocytes_DCIS <- FindVariableFeatures(lymphocytes_DCIS)
lymphocytes_DCIS <- ScaleData(lymphocytes_DCIS)
lymphocytes_DCIS <- RunPCA(lymphocytes_DCIS, verbose = TRUE)
lymphocytes_DCIS <- FindNeighbors(lymphocytes_DCIS)
lymphocytes_DCIS <- RunUMAP(lymphocytes_DCIS, dims = 1:30, reduction = "pca", reduction.name = "umap.lymphocyte")
lymphocytes_DCIS <- FindClusters(lymphocytes_DCIS,resolution=2.0)

fine_markers <- c("CD3D","CD3E","CD3G","CD4","CD8A","CD79A","MS4A1","JCHAIN","NCAM1")
VlnPlot(lymphocytes_DCIS,fine_markers,ncol = 3, pt.size = 0) ## clusters 
FeaturePlot(lymphocytes_DCIS,fine_markers,ncol = 3, reduction = "umap.lymphocyte") ## clusters 

BcellMarkers <- c("MS4A1","CD79A",sort(c("CD27","IGHD","IGHM","IGHG1","IGHG2","IGHG4","IGHA1","IGHA2","CD38","SDC1","MZB1","BCL6")))
VlnPlot(lymphocytes_DCIS,BcellMarkers,ncol = 4, pt.size = 0) ## clusters 

CD8TcellMarkers <- c("SELL", "CCR7", "TCF7", "LEF1", "IL7R", "GZMA", "GZMB", "GZMK", "PRF1", "IFNG", "TNF", "EOMES", "ITGAE", "CD69", "LAG3", "CTLA4")
VlnPlot(lymphocytes_DCIS,CD8TcellMarkers,ncol = 4, pt.size = 0) ## clusters 

CD4TcellMarkers <- c("CCR7", "IL7R", "SELL", "TCF7", "LEF1", "FOXP3", "MAF", "IKZF3", "CD40LG", "PRDM1", "RUNX3", "ITGB1", "CXCR5", "CXCL13", "GATA3", "FAS")
VlnPlot(lymphocytes_DCIS,CD4TcellMarkers,ncol = 4, pt.size = 0) ## clusters 

annotLymphocytes <- read.csv("lymphocyteReannot.csv",header = T,row.names = 1) ## table that has the clusters and the annotation that was decided based on the plots generated above, in two columns: fine_cell_type and subset

lymphocytes_DCIS$fine_cell_type <- annotLymphocytes[as.character(lymphocytes_DCIS$seurat_clusters),"Fine.cell.type"]
lymphocytes_DCIS$subset <- annotLymphocytes[as.character(lymphocytes_DCIS$seurat_clusters),"Subset"]

lymphoNewSubsets <- data.frame(Barcode = lymphocytes_DCIS$barcodes, fine_cell_type = lymphocytes_DCIS$fine_cell_type, subset = lymphocytes_DCIS$subset)

DimPlot(lymphocytes_DCIS, group.by = c("seurat_clusters","fine_cell_type","subset","grade_simplified"),reduction = "umap.lymphocyte",label = T, label.size = 6, ncol=2)


## modify in whole object
DCISSeurat$fine_cell_type[lymphoNewSubsets$Barcode] <- lymphoNewSubsets$fine_cell_type
DCISSeurat$subset[lymphoNewSubsets$Barcode] <- lymphoNewSubsets$subset

DCISSeurat$display_population = ifelse(is.na(DCISSeurat$subset),DCISSeurat$fine_cell_type,DCISSeurat$subset)

DCISSeurat <- subset(DCISSeurat, display_population != "DELETE") ## remove populations that were manually marked as to be deleted (contamination by other cell types)

table(DCISSeurat$fine_cell_type)
table(DCISSeurat$display_population)




###########
# Myeloid #
###########


myeloid_DCIS <- subset(DCISSeurat, large_cell_type == "Myeloid")

myeloid_DCIS <- NormalizeData(myeloid_DCIS)
myeloid_DCIS <- FindVariableFeatures(myeloid_DCIS)
myeloid_DCIS <- ScaleData(myeloid_DCIS)
myeloid_DCIS <- RunPCA(myeloid_DCIS, verbose = TRUE)
myeloid_DCIS <- FindNeighbors(myeloid_DCIS)
myeloid_DCIS <- RunUMAP(myeloid_DCIS, dims = 1:30, reduction = "pca", reduction.name = "umap.myeloid")
myeloid_DCIS <- FindClusters(myeloid_DCIS,resolution=3.0)

DimPlot(myeloid_DCIS,reduction = "umap.myeloid",label = T)

myeloidMarkers <- c("AIF1","CD68","CD163","CSF1R","CD14","ITGAM","ITGAX","FCGR3A","SIGLEC1","FOLR2","TREM2","SEPP1","SPP1","CXCL9","CD1C","XCR1","FCGR1A","CD123","CD15","LYVE1","KIT","C1QC","THBD")

VlnPlot(myeloid_DCIS,myeloidMarkers,ncol = 4, pt.size = 0)

markers <- FindAllMarkers(myeloid_DCIS)

finemyeloid <- c("0" = "Macrophage", "1" = "Macrophage", "2" = "Macrophage", "3" = "Macrophage", "4" = "Macrophage",
                 "5" = "Macrophage", "6" = "Macrophage", "7" = "Macrophage", "8" = "Macrophage", "9" = "Macrophage",
                 "10" = "Macrophage", "11" = "T cell", "12" = "Macrophage", "13" = "Macrophage", "14" = "Macrophage",
                 "15" = "Dendritic cell", "16" = "Macrophage", "17" = "Mast cell", "18" = "Macrophage", "19" = "Mast cell",
                 "20" = "Macrophage", "21" = "Macrophage", "22" = "Macrophage", "23" = "Dendritic cell", "24" = "Macrophage",
                 "25" = "Macrophage", "26" = "Macrophage", "27" = "Macrophage", "28" = "Macrophage", "29" = "Macrophage",
                 "30" = "Dendritic cell", "31" = "Macrophage", "32" = "Dendritic cell", "33" = "Dendritic cell", "34" = "Macrophage",
                 "35" = "Macrophage", "36" = "Macrophage", "37" = "Macrophage")

subsetmyeloid <- c("0" = "FOLR2 Macrophage", "1" = "FOLR2 Macrophage", "2" = "FOLR2 Macrophage", "3" = "FOLR2 Macrophage", "4" = "DCN Macrophage",
                   "5" = "FOLR2 Macrophage", "6" = "Unch. Macrophage", "7" = "CXCL9 Macrophage", "8" = "TREM2 Macrophage", "9" = "TREM2 Macrophage",
                   "10" = "FOLR2 Macrophage", "11" = "Contamination", "12" = "TREM2 Macrophage", "13" = "TREM2 Macrophage", "14" = "CLEC10A Macrophage",
                   "15" = "cDC2", "16" = "TREM2 Macrophage", "17" = NA, "18" = "TREM2 Macrophage", "19" = NA,
                   "20" = "FOLR2 Macrophage", "21" = "FOXA1 Macrophage", "22" = "TREM2 Macrophage", "23" = "cDC2", "24" = "CXCL9 Macrophage",
                   "25" = "TREM2 Macrophage", "26" = "S100A9 Macrophage", "27" = "Unch. Macrophage", "28" = "TREM2 Macrophage", "29" = "TREM2 Macrophage",
                   "30" = "cDC2", "31" = "FOLR2 Macrophage", "32" = "cDC1", "33" = "CCR7 DC", "34" = "S100A9 Macrophage",
                   "35" = "FOLR2 Macrophage", "36" = "TREM2 Macrophage", "37" = "TREM2 Macrophage")

newMeta[colnames(myeloid_DCIS),"fine_cell_type"] = finemyeloid[as.character(Idents(myeloid_DCIS))]
newMeta[colnames(myeloid_DCIS),"subset"] = subsetmyeloid[as.character(Idents(myeloid_DCIS))]
newMeta$display_population <- ifelse(is.na(newMeta$subset),newMeta$fine_cell_type,newMeta$subset)

myeloid_DCIS <- AddMetaData(myeloid_DCIS,newMeta)

DimPlot(myeloid_DCIS,group.by = c("fine_cell_type","display_population"),reduction = "umap.myeloid")

save(newMeta, file = "20250813_newMeta.RData")

table(newMeta$display_population,useNA="ifany")


#################################
# Revisit Macrophage annotation #
#################################


macrophages_DCIS <- subset(DCISSeurat, fine_cell_type == "Macrophage")

macrophages_DCIS <- NormalizeData(macrophages_DCIS)
macrophages_DCIS <- FindVariableFeatures(macrophages_DCIS)
macrophages_DCIS <- ScaleData(macrophages_DCIS)
macrophages_DCIS <- RunPCA(macrophages_DCIS, verbose = TRUE)
macrophages_DCIS <- FindNeighbors(macrophages_DCIS)
macrophages_DCIS <- RunUMAP(macrophages_DCIS, dims = 1:30, reduction = "pca", reduction.name = "umap.myeloid")
macrophages_DCIS <- FindClusters(macrophages_DCIS,resolution=3.0)

DimPlot(macrophages_DCIS,reduction = "umap.myeloid",label = T, label.size = 6)

myeloidMarkers <- c("FOLR2","TREM2","GPNMB","APOE","SIGLEC1","KRT19","SPP1","CXCL9","LYVE1","IL4I1","CD14","FOXA1","CCR2","VCAN","CD68","ACTA2")

VlnPlot(macrophages_DCIS,myeloidMarkers,ncol = 4, pt.size = 0)

markers <- FindAllMarkers(macrophages_DCIS)
markers$pct.diff <- markers$pct.1 - markers$pct.2


markers_28 <- markers[as.character(markers$cluster)=="28" & markers$avg_log2FC>0,]
markers_28$pct.diff <- markers_28$pct.1 - markers_28$pct.2

markers_13 <- markers[as.character(markers$cluster)=="13" & markers$avg_log2FC>0,]
markers[as.character(markers$cluster)=="9" & markers$avg_log2FC>0,] |> View()
markers[as.character(markers$cluster)=="11" & markers$avg_log2FC>0,] |> View()


FeaturePlot(macrophages_DCIS,c("FOLR2","TREM2","SPP1"),split.by = "grade_simplified", reduction = "umap.myeloid")
FeaturePlot(macrophages_DCIS,c("FOXA1","KRT19","GATA3","EPCAM"),split.by = "grade_simplified", reduction = "umap.myeloid")


macrophages_DCIS$old_annot <- macrophages_DCIS$display_population

annotMacroSubsets <- read.csv("macrophageReannot.csv",header = T,row.names = 1) ## table that has the clusters and the annotation that was decided based on the plots generated above, in two columns: fine_cell_type and subset

macrophages_DCIS$subset <- annotMacroSubsets[as.character(macrophages_DCIS$seurat_clusters),"Subset"]

DimPlot(macrophages_DCIS,group.by = c("seurat_clusters","subset","old_annot"), reduction = "umap.myeloid", label = TRUE)

macroNewSubsets <- data.frame(Barcode = macrophages_DCIS$barcodes, subset = macrophages_DCIS$subset)


## modify in whole object
DCISSeurat$subset[macroNewSubsets$Barcode] <- macroNewSubsets$subset

table(DCISSeurat$subset)

DCISSeurat$fine_cell_type[which(DCISSeurat$subset == "Classical monocyte")] <- "Monocyte"
DCISSeurat$display_population = ifelse(is.na(DCISSeurat$subset),DCISSeurat$fine_cell_type,DCISSeurat$subset)

DCISSeurat <- subset(DCISSeurat, display_population != "DELETE")



###############
# Endothelial #
###############


endothelials_DCIS <- subset(DCISSeurat, large_cell_type == "Endothelial")

endothelials_DCIS <- NormalizeData(endothelials_DCIS)
endothelials_DCIS <- FindVariableFeatures(endothelials_DCIS)
endothelials_DCIS <- ScaleData(endothelials_DCIS)
endothelials_DCIS <- RunPCA(endothelials_DCIS, verbose = TRUE)
endothelials_DCIS <- FindNeighbors(endothelials_DCIS)
endothelials_DCIS <- RunUMAP(endothelials_DCIS, dims = 1:30, reduction = "pca", reduction.name = "umap.endothelials")

endothelials_DCIS <- FindClusters(endothelials_DCIS,resolution = 1.5)

DimPlot(endothelials_DCIS,reduction = "umap.endothelials",label = T)

endothelialMarkers <- c("PECAM1","LYVE1","CD34","PDPN","PROX1","FLT4")

VlnPlot(endothelials_DCIS,endothelialMarkers,ncol = 3, pt.size = 0) ## clusters 

fineendothelials <- c("0" = "Blood Endothelial cell", "1" = "Blood Endothelial cell", "2" = "Blood Endothelial cell", "3" = "Blood Endothelial cell", "4" = "Blood Endothelial cell",
                      "5" = "Blood Endothelial cell", "6" = "Blood Endothelial cell", "7" = "Blood Endothelial cell", "8" = "Lymphatic Endothelial cell", "9" = "Blood Endothelial cell")

newMeta[colnames(endothelials_DCIS),"fine_cell_type"] = fineendothelials[as.character(Idents(endothelials_DCIS))]
newMeta$display_population <- ifelse(is.na(newMeta$subset),newMeta$fine_cell_type,newMeta$subset)

endothelials_DCIS <- AddMetaData(endothelials_DCIS,newMeta)

DimPlot(endothelials_DCIS,group.by = c("fine_cell_type","display_population"),reduction = "umap.endothelials")

#save(newMeta, file = "20250813_newMeta.RData")
rm(endothelials_DCIS)
gc()


#########################
# Fibroblasts-pericytes #
#########################


fibro_DCIS <- subset(DCISSeurat, large_cell_type == "Fibroblast-pericyte")

fibro_DCIS <- NormalizeData(fibro_DCIS)
fibro_DCIS <- FindVariableFeatures(fibro_DCIS)
fibro_DCIS <- ScaleData(fibro_DCIS)
fibro_DCIS <- RunPCA(fibro_DCIS, verbose = TRUE)
fibro_DCIS <- FindNeighbors(fibro_DCIS)
fibro_DCIS <- RunUMAP(fibro_DCIS, dims = 1:30, reduction = "pca", reduction.name = "umap.fibro")
fibro_DCIS <- FindClusters(fibro_DCIS,resolution=2.0)

DimPlot(fibro_DCIS, group.by = c("seurat_clusters","fine_cell_type","subset","grade_simplified"),reduction = "umap.fibro",label = T, label.size = 6, ncol=2)


fine_markers <- c("PDPN","PDGFRB","PDGFRA","S100A4","FAP","NES","CSPG4","MCAM","ITGB1","TNFS4","ACTA2","ADIPOQ","CXCL12","PTPRC","COL3A1","KRT5","COL5A1","THBS2","CD34","POSTN","FN1")
VlnPlot(fibro_DCIS,fine_markers,ncol = 4, pt.size = 0) ## clusters 

CAFMarkers <- c("CXCL12","IGF1", "APOD", "C7", "C3", "CXCL14", "IGBP6","ACTA2","COL3A1", "COL10A1", "THBS2", "VCAN", "POSTN", "CD34", "COL5A2", "COL1A1","FN1", "NES", "TNC","FABP4","ATXN1")
VlnPlot(fibro_DCIS,CAFMarkers,ncol = 4, pt.size = 0) ## clusters 



annotFibro <- read.csv("fibroReannot.csv",header = T,row.names = 1) ## table that has the clusters and the annotation that was decided based on the plots generated above, in two columns: fine_cell_type and subset

fibro_DCIS$fine_cell_type <- annotFibro[as.character(fibro_DCIS$seurat_clusters),"Large"]
fibro_DCIS$subset <- annotFibro[as.character(fibro_DCIS$seurat_clusters),"Fine"]


DimPlot(fibro_DCIS, group.by = c("seurat_clusters","fine_cell_type","subset","grade_simplified"),reduction = "umap.fibro",label = T, label.size = 6, ncol=2)



fibroNewSubsets <- data.frame(Barcode = fibro_DCIS$barcodes, fine_cell_type = fibro_DCIS$fine_cell_type, subset = fibro_DCIS$subset)

DimPlot(fibro_DCIS, group.by = c("seurat_clusters","fine_cell_type","subset","grade_simplified"),reduction = "umap.fibro",label = T, label.size = 6, ncol=2)



## modify in whole object
DCISSeurat$fine_cell_type[fibroNewSubsets$Barcode] <- fibroNewSubsets$fine_cell_type
DCISSeurat$subset[fibroNewSubsets$Barcode] <- fibroNewSubsets$subset

DCISSeurat$display_population = ifelse(is.na(DCISSeurat$subset),DCISSeurat$fine_cell_type,DCISSeurat$subset)

DCISSeurat <- subset(DCISSeurat, display_population != "DELETE")




#############################
# Annotate epithelial cells #
#############################


EpithelialFibro = subset(DCISSeurat,large_cell_type %in% c("Epithelial","Fibroblast-pericyte"))
copykat.EpiFibro <- copykat(rawmat=LayerData(EpithelialFibro, assay = "RNA", layer = "counts"),
                       id.type="S",# gene symbols
                       ngene.chr=5,# min number of genes/chr to enable discovery
                       win.size=25,
                       KS.cut=0.1,
                       sam.name = "Epithelial_Fibro",
                       distance="euclidean",
                       norm.cell.names="",
                       output.seg="FALSE",
                       plot.genes="FALSE",
                       genome="hg20",
                       n.cores=32)


copycatPred = copykat.EpiFibro$prediction

copykatPred <- read.table("./Copykat/Epithelial_Fibro_copykat_prediction.txt",sep="\t",header=1)
rownames(copykatPred) = copykatPred$cell.names
DCISSeurat$ploidy = copykatPred[DCISSeurat$barcodes,"copykat.pred"]


addmargins(table(ploidy = DCISSeurat$ploidy, cellType = DCISSeurat$large_cell_type, useNA = "ifany"))
addmargins(table(ploidy = DCISSeurat$ploidy, grade = DCISSeurat$grade))

identical(rownames(newMeta),colnames(DCISSeurat))

table(DCISSeurat$ploidy[rownames(newMeta)],useNA="ifany")

newMeta$fine_cell_type = ifelse(newMeta$large_cell_type == "Epithelial" & DCISSeurat$ploidy[rownames(newMeta)] == "aneuploid", "Malignant Epithelial cell", newMeta$fine_cell_type)
newMeta$fine_cell_type = ifelse(newMeta$large_cell_type == "Epithelial" & newMeta$fine_cell_type == "Epithelial", "Normal Epithelial cell", newMeta$fine_cell_type)

table(newMeta$fine_cell_type)

newMeta$display_population = ifelse(newMeta$large_cell_type == "Epithelial",newMeta$fine_cell_type,newMeta$display_population)

DCISSeurat <- AddMetaData(DCISSeurat, newMeta)