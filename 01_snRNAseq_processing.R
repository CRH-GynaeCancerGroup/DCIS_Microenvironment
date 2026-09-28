########################
# snRNA-seq processing #
########################


library(Seurat)
library(SeuratWrappers)
library(parallel)
source("~/R/Rutils/FlorentUtils.R")

sampleNames = list.dirs("./Preprocessed",recursive = FALSE) |> fp.strsplit(sep = "/", index = 3) ## read in folders containing the output from cellranger

sampleAnnot <- fp.quickRead("sampleAnnot.csv",sep=",") #table containing the grade annotation for each sample
rownames(sampleAnnot) <- sampleAnnot$Sample

countsList <- mclapply(sampleNames, function(s){
  
  h5File <- NULL
  fileList <- list.files(paste0("./Preprocessed/",s),recursive = TRUE)
  h5File <- fileList[grep("filtered_feature_bc_matrix.h5",fileList)]
  dgCMatrix <- Read10X_h5(paste0("./Preprocessed/",s,"/",h5File))
  colnames(dgCMatrix) <- paste(s,colnames(dgCMatrix),sep="_")
  
  return(dgCMatrix)
})


metaList <- mclapply(sampleNames, function(s){
  
  h5File <- NULL
  fileList <- list.files(paste0("./Preprocessed/",s),recursive = TRUE)
  h5File <- fileList[grep("filtered_feature_bc_matrix.h5",fileList)]
  dgCMatrix <- Read10X_h5(paste0("./Preprocessed/",s,"/",h5File))
  
  sampleID = s
  
  meta <- data.frame(barcodes = paste(s,colnames(dgCMatrix),sep="_"),
                     sampleID = rep(sampleID,ncol(dgCMatrix)),
                     grade = rep(sampleAnnot[sampleID,"Grade"],ncol(dgCMatrix)))
  rownames(meta) <- meta$barcodes
  
  return(meta)
})
names(countsList) <- names(metaList) <- sampleNames


DCISSeurat <- CreateSeuratObject(counts = countsList,meta.data = dplyr::bind_rows(metaList))
Layers(DCISSeurat)

DCISSeurat[["percent.mt"]] <- PercentageFeatureSet(DCISSeurat, pattern = "^MT-")

## remove low-quality cells
VlnPlot(DCISSeurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,)
print(paste("Before filtering:",length(DCISSeurat$orig.ident),"cells",sep=" "))
DCISSeurat <- subset(DCISSeurat, nFeature_RNA > 600 & percent.mt < 25)
print(paste("After filtering:",length(DCISSeurat$orig.ident),"cells",sep=" "))

#pre-processing
DCISSeurat <- NormalizeData(DCISSeurat)
DCISSeurat <- FindVariableFeatures(DCISSeurat)
DCISSeurat <- ScaleData(DCISSeurat)
DCISSeurat <- RunPCA(DCISSeurat, verbose = TRUE)
DCISSeurat <- FindNeighbors(DCISSeurat)
DCISSeurat <- FindClusters(DCISSeurat)
DCISSeurat <- RunUMAP(DCISSeurat, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

## Integration
DCISSeurat <- IntegrateLayers(
  object = DCISSeurat, method = HarmonyIntegration,
  new.reduction = "integrated.harmony"
)

DCISSeurat <- FindNeighbors(DCISSeurat, reduction = "integrated.harmony", dims = 1:30)
DCISSeurat <- FindClusters(DCISSeurat, cluster.name = "harmony_clusters")

DCISSeurat <- RunUMAP(DCISSeurat, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")

DimPlot(DCISSeurat, reduction = "umap.unintegrated", shuffle = TRUE, group.by = "sampleID") +
  DimPlot(DCISSeurat, reduction = "umap.harmony", shuffle = TRUE, group.by = "sampleID") +
  DimPlot(DCISSeurat, reduction = "umap.harmony", shuffle = TRUE, group.by = "harmony_clusters") +
  DimPlot(DCISSeurat, reduction = "umap.harmony", shuffle = TRUE, group.by = "grade")

DCISSeurat <- JoinLayers(DCISSeurat)
