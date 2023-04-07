library(Seurat)
library(dplyr)
library(patchwork)


#####Read files#####

setwd("/Volumes/LaCie SSD/All Raw Data Cardiotoxicity/scRNAseq")
CTRL1 <- readRDS("8_20200903_Mu_heart_captisol_1h_Balbc_total_cells_v3_2020_10_05_seurat.rds")
CTRL2 <- readRDS("9_20200903_Mu_heart_captisol_1h_Balbc_totoal_cells_v3_2020_10_05_seurat.rds")
CFZ1 <- readRDS("10_20200903_Mu_heart_CFZ_1h_Balbc_total_cells_v3_2020_10_05_seurat.rds")
CFZ2 <- readRDS("11_20200903_Mu_heart_CFZ_1h_Balbc_total_cells_v3_2020_10_05_seurat.rds")
CFZATRA1 <-readRDS("12_20200903_Mu_heart_CFZATRA_1h_Balbc_total_cells_v3_2020_10_05_seurat.rds")
CFZATRA2 <- readRDS("13_20200903_Mu_heart_CFZATRA_1h_Balbc_total_cells_v3_2020_10_05_seurat.rds")

CTRLmerged <- merge(CTRL1, CTRL2)
CTRLmerged$dataset
CTRLmerged$grp <- "Vehicle"
CFZmerged <- merge(CFZ1, CFZ2)
CFZmerged$dataset
CFZmerged$grp <- "CFZ"

CFZATRAmerged <- merge(CFZATRA1, CFZATRA2)
CFZATRAmerged$dataset
CFZATRAmerged$grp <- "CFZATRA"

#####set resolution#####
res = c(0.2,0.4,0.6,0.8)

##########################################################################################################################################################################
##########################################################################################################################################################################

# Identify the 10 most highly variable genes
CFZvsCTRL <- NormalizeData(CFZvsCFZATRA)
CFZvsCTRL <- merge(CFZmerged, y = CTRLmerged, add.cell.ids = c("CFZ", "CTRL"), project = "Cardio",
                   merge.data = TRUE)

CFZvsCTRL <- FindVariableFeatures(CFZvsCTRL, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(CFZvsCTRL), 10)
VlnPlot(CFZvsCTRL, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

###Scale, PCA, UMAP###
all.genes <- rownames(CFZvsCTRL)
CFZvsCTRL <- ScaleData(CFZvsCTRL, features = all.genes)

CFZvsCTRL <- RunPCA(object = CFZvsCTRL, npcs = 20, verbose = FALSE)
CFZvsCTRL <- RunTSNE(object = CFZvsCTRL, reduction = "pca", dims = 1:20)
CFZvsCTRL <- RunUMAP(object = CFZvsCTRL, reduction = "pca", dims = 1:20)
CFZvsCTRL <- FindNeighbors(object = CFZvsCTRL, reduction = "pca", 
                           dims = 1:20)

for (i in 1:length(res)) {
  CFZvsCTRL <- FindClusters(object = CFZvsCTRL, resolution = 0.6, #res[i], 
                            random.seed = 1234)
}


DimPlot(CFZvsCTRL, label = TRUE)
ElbowPlot(CFZvsCTRL)

CFZvsCTRL$seurat_clusters <-CFZvsCTRL$RNA_snn_res.0.6
IDCluster_CFZvsCTRL <- data.frame(ID=colnames(CFZvsCTRL)) %>% mutate(cluster=CFZvsCTRL$seurat_clusters)
VlnPlot(CFZvsCTRL, features=unlist(TopFeatures(CFZvsCTRL[["pca"]])), ncol=3)

VlnPlot(CFZvsCTRL, features= c("ENSMUSG00000091898.Tnnc1", "ENSMUSG00000026414.Tnnt2", 
                               "ENSMUSG00000035458.Tnni3"), ncol=3)

cardiom<-FeaturePlot(object = CFZvsCTRL, features= c("ENSMUSG00000091898.Tnnc1", "ENSMUSG00000026414.Tnnt2", 
                                                     "ENSMUSG00000035458.Tnni3"), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZvsCTRL, features= c("ENSMUSG00000091898.Tnnc1", "ENSMUSG00000026414.Tnnt2", 
                               "ENSMUSG00000035458.Tnni3"), ncol=3)

fibr<-FeaturePlot(object = CFZvsCTRL, features= c("ENSMUSG00000045680.Tcf21", "ENSMUSG00000001506.Col1a1", 
                                                  "ENSMUSG00000029661.Col1a2"), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZvsCTRL, features= c("ENSMUSG00000045680.Tcf21", "ENSMUSG00000001506.Col1a1", 
                               "ENSMUSG00000029661.Col1a2"), ncol=3)

macrop<-FeaturePlot(object = CFZvsCTRL, features= c("ENSMUSG00000024610.Cd74", "ENSMUSG00000051439.Cd14", 
                                                    "ENSMUSG00000052336.Cx3cr1"), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZvsCTRL, features= c("ENSMUSG00000024610.Cd74", "ENSMUSG00000051439.Cd14", 
                               "ENSMUSG00000052336.Cx3cr1"), ncol=3)

VEC<-FeaturePlot(object = CFZvsCTRL, features= c("ENSMUSG00000002944.Cd36", "ENSMUSG00000031871.Cdh5", 
                                                 "ENSMUSG00000020717.Pecam1"), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZvsCTRL, features= c("ENSMUSG00000002944.Cd36", "ENSMUSG00000031871.Cdh5", 
                               "ENSMUSG00000020717.Pecam1"), ncol=3)

Pericy<-FeaturePlot(object = CFZvsCTRL, features= c("ENSMUSG00000032135.Mcam", "ENSMUSG00000025348.Itga7" 
), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZvsCTRL, features= c("ENSMUSG00000032135.Mcam", "ENSMUSG00000025348.Itga7"), ncol=3)


table(Idents(CFZvsCTRL))

p1<-DimPlot(CFZvsCTRL, label = TRUE, reduction = "umap", pt.size = 0.5, group.by= "grp") #+ NoLegend()
p2<-DimPlot(CFZvsCTRL, label = TRUE, reduction = "umap", pt.size = 0.5, group.by= "seurat_clusters")

library(cowplot)
plot_grid(p1, p2)

#####Find Marker genes#####
CFZvsCaptisol.All.rna.markers <- FindAllMarkers(CFZvsCTRL, max.cells.per.ident = 200, min.diff.pct = 0.3, test.use = "wilcox", only.pos = FALSE)

All.hearts.rna.markers.cardiomyvsfibroblasts<- FindMarkers(CFZvsCTRL, ident.1 = 13, ident.2 = c(0,1,2,3,5) , min.pct = 0.5, grouping.var = "grp")
All.hearts.rna.markers.cardiomyvsfibroblasts<- cbind(Row.Names = rownames(All.hearts.rna.markers.cardiomyvsfibroblasts), All.hearts.rna.markers.cardiomyvsfibroblasts)

All.hearts.rna.markers.cardiomyvsVEC<- FindMarkers(CFZvsCTRL, ident.1 = 13, ident.2 = c(4,9) , min.pct = 0.5, grouping.var = "grp")
All.hearts.rna.markers.cardiomyvsVEC<- cbind(Row.Names = rownames(All.hearts.rna.markers.cardiomyvsVEC), All.hearts.rna.markers.cardiomyvsVEC)

All.hearts.rna.markers.cardiomyvsMacroph<- FindMarkers(CFZvsCTRL, ident.1 = 13, ident.2 = c(7,8,11) , min.pct = 0.5, grouping.var = "grp")
All.hearts.rna.markers.cardiomyvsMacroph<- cbind(Row.Names = rownames(All.hearts.rna.markers.cardiomyvsMacroph), All.hearts.rna.markers.cardiomyvsMacroph)

All.hearts.rna.markers.cardiomyvsPericyt<- FindMarkers(CFZvsCTRL, ident.1 = 13, ident.2 = c(9,10) , min.pct = 0.5, grouping.var = "grp")
All.hearts.rna.markers.cardiomyvsPericyt<- cbind(Row.Names = rownames(All.hearts.rna.markers.cardiomyvsPericyt), All.hearts.rna.markers.cardiomyvsPericyt)

All.hearts.rna.markers.VEC<- FindMarkers(CFZvsCTRL, ident.1 = c(4,9), min.pct = 0.5, grouping.var = "grp")
All.hearts.rna.markers.VEC<- cbind(Row.Names = rownames(All.hearts.rna.markers.VEC), All.hearts.rna.markers.VEC)

All.hearts.rna.markers.macroph<- FindMarkers(CFZvsCTRL, ident.1 = c(7,8,11), min.pct = 0.5, grouping.var = "grp")
All.hearts.rna.markers.macroph<- cbind(Row.Names = rownames(All.hearts.rna.markers.macroph), All.hearts.rna.markers.macroph)

##########################################################################################################################################################################
###CFZ vs CFZATRA###
CFZvsCFZATRA <- merge(CFZmerged, y = CFZATRAmerged, add.cell.ids = c("CFZ", "CFZATRA"), project = "Cardio",
                      merge.data = TRUE)

# visualization of clusters ALL groups
CFZvsCFZATRA <- NormalizeData(CFZvsCFZATRA)
CFZvsCFZATRA <- FindVariableFeatures(CFZvsCFZATRA, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(CFZvsCFZATRA), 10)
VlnPlot(CFZvsCFZATRA, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

plot1 <- VariableFeaturePlot(CFZvsCFZATRA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

###Scale, PCA, UMAP###
all.genes <- rownames(CFZvsCFZATRA)
CFZvsCFZATRA <- ScaleData(CFZvsCFZATRA, features = all.genes)

CFZvsCFZATRA <- RunPCA(object = CFZvsCFZATRA, npcs = 20, verbose = FALSE)
CFZvsCFZATRA <- RunTSNE(object = CFZvsCFZATRA, reduction = "pca", dims = 1:15)
CFZvsCFZATRA <- RunUMAP(object = CFZvsCFZATRA, reduction = "pca", dims = 1:15)
CFZvsCFZATRA <- FindNeighbors(object = CFZvsCFZATRA, reduction = "pca", 
                              dims = 1:20)
for (i in 1:length(res)) {
  CFZvsCFZATRA <- FindClusters(object = CFZvsCFZATRA, resolution = 0.4, #res[i], 
                               random.seed = 1234)
}

DimPlot(CFZvsCFZATRA, label = TRUE)
ElbowPlot(CFZvsCFZATRA)

CFZvsCFZATRA$seurat_clusters <-CFZvsCFZATRA$RNA_snn_res.0.6
IDCluster_CFZvsCFZATRA <- data.frame(ID=colnames(CFZvsCFZATRA)) %>% mutate(cluster=CFZvsCFZATRA$seurat_clusters)
IDCluster_CFZvsCFZATRA <- write_tsv(IDCluster_CFZvsCFZATRA, "IDCluster_CFZvsCFZATRA.tsv")



VlnPlot(CFZvsCFZATRA, features=unlist(TopFeatures(CFZvsCFZATRA[["pca"]])), ncol=3)

VlnPlot(CFZvsCFZATRA, features=unlist(TopFeatures(CFZvsCFZATRA[["pca"]])), ncol=3)
VlnPlot(CFZvsCFZATRA, features= c("ENSMUSG00000091898.Tnnc1", "ENSMUSG00000026414.Tnnt2", 
                                  "ENSMUSG00000035458.Tnni3"), ncol=3)

cardiom<-FeaturePlot(object = CFZvsCFZATRA, features= c("ENSMUSG00000091898.Tnnc1", "ENSMUSG00000026414.Tnnt2", 
                                                        "ENSMUSG00000035458.Tnni3"), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZvsCFZATRA, features= c("ENSMUSG00000091898.Tnnc1", "ENSMUSG00000026414.Tnnt2", 
                                  "ENSMUSG00000035458.Tnni3"), ncol=3)

fibr<-FeaturePlot(object = CFZvsCFZATRA, features= c("ENSMUSG00000045680.Tcf21", "ENSMUSG00000001506.Col1a1", 
                                                     "ENSMUSG00000029661.Col1a2"), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZvsCFZATRA, features= c("ENSMUSG00000045680.Tcf21", "ENSMUSG00000001506.Col1a1", 
                                  "ENSMUSG00000029661.Col1a2"), ncol=3)

macrop<-FeaturePlot(object = CFZvsCFZATRA, features= c("ENSMUSG00000024610.Cd74", "ENSMUSG00000051439.Cd14", 
                                                       "ENSMUSG00000052336.Cx3cr1"), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZvsCFZATRA, features= c("ENSMUSG00000024610.Cd74", "ENSMUSG00000051439.Cd14", 
                                  "ENSMUSG00000052336.Cx3cr1"), ncol=3)

VEC<-FeaturePlot(object = CFZvsCFZATRA, features= c("ENSMUSG00000002944.Cd36", "ENSMUSG00000031871.Cdh5", 
                                                    "ENSMUSG00000020717.Pecam1"), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZvsCFZATRA, features= c("ENSMUSG00000002944.Cd36", "ENSMUSG00000031871.Cdh5", 
                                  "ENSMUSG00000020717.Pecam1"), ncol=3)

Pericy<-FeaturePlot(object = CFZvsCFZATRA, features= c("ENSMUSG00000032135.Mcam", "ENSMUSG00000025348.Itga7" 
), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZvsCFZATRA, features= c("ENSMUSG00000032135.Mcam", "ENSMUSG00000025348.Itga7"), ncol=3)


table(Idents(CFZvsCFZATRA))
p1<-DimPlot(CFZvsCFZATRA, label = FALSE, reduction = "umap", pt.size = 0.5, group.by= "dataset") + NoLegend()
p2<-DimPlot(CFZvsCFZATRA, label = TRUE, reduction = "umap", pt.size = 0.5, group.by= "seurat_clusters")

library(cowplot)
plot_grid(p1, p2)

#####get marker genes###
CFZvsCFZATRA.All.rna.markers <- FindAllMarkers(CFZvsCFZATRA, max.cells.per.ident = 200, min.diff.pct = 0.3, test.use = "wilcox", only.pos = FALSE) #cardiomyocytes=9

All.hearts.rna.markers.cardiomyvsfibroblasts<- FindMarkers(CFZvsCFZATRA, ident.1 = 15, ident.2 = c(0,1,2,3,6) , min.pct = 0.5, grouping.var = "dataset")
All.hearts.rna.markers.cardiomyvsfibroblasts<- cbind(Row.Names = rownames(All.hearts.rna.markers.cardiomyvsfibroblasts), All.hearts.rna.markers.cardiomyvsfibroblasts)

All.hearts.rna.markers.cardiomyvsVEC<- FindMarkers(CFZvsCTRL, ident.1 = 15, ident.2 = c(4,8) , min.pct = 0.5, grouping.var = "dataset")
All.hearts.rna.markers.cardiomyvsVEC<- cbind(Row.Names = rownames(All.hearts.rna.markers.cardiomyvsVEC), All.hearts.rna.markers.cardiomyvsVEC)

All.hearts.rna.markers.cardiomyvsMacroph<- FindMarkers(CFZvsCTRL, ident.1 = 15, ident.2 = c(5,13,14) , min.pct = 0.5, grouping.var = "dataset")
All.hearts.rna.markers.cardiomyvsMacroph<- cbind(Row.Names = rownames(All.hearts.rna.markers.cardiomyvsMacroph), All.hearts.rna.markers.cardiomyvsMacroph)

All.hearts.rna.markers.cardiomyvsPericyt<- FindMarkers(CFZvsCTRL, ident.1 = 15, ident.2 = c(9,10) , min.pct = 0.5, grouping.var = "dataset")
All.hearts.rna.markers.cardiomyvsPericyt<- cbind(Row.Names = rownames(All.hearts.rna.markers.cardiomyvsPericyt), All.hearts.rna.markers.cardiomyvsPericyt)

All.hearts.rna.markers.VEC<- FindMarkers(CFZvsCFZATRA, ident.1 = c(4,8), min.pct = 0.5, grouping.var = "dataset")
All.hearts.rna.markers.VEC<- cbind(Row.Names = rownames(All.hearts.rna.markers.VEC), All.hearts.rna.markers.VEC)

All.hearts.rna.markers.macroph<- FindMarkers(CFZvsCFZATRA, ident.1 = c(5,13,14), min.pct = 0.5, grouping.var = "dataset")
All.hearts.rna.markers.macroph<- cbind(Row.Names = rownames(All.hearts.rna.markers.macroph), All.hearts.rna.markers.VEC)


##########################################################################################################################################################################
###CFZATRAvsCTRL###

CFZATRAvsCTRL <- merge(CFZATRAmerged, y = CTRLmerged, add.cell.ids = c("CFZATRA", "CTRL"), project = "Cardio",
                       merge.data = TRUE)

# visualization of clusters ALL groups
CFZATRAvsCTRL <- NormalizeData(CFZATRAvsCTRL)
CFZATRAvsCTRL <- FindVariableFeatures(CFZATRAvsCTRL, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(CFZATRAvsCTRL), 10)
VlnPlot(CFZATRAvsCTRL, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

plot1 <- VariableFeaturePlot(CFZATRAvsCTRL)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


###Scale, PCA, UMAP###
all.genes <- rownames(CFZATRAvsCTRL)
CFZATRAvsCTRL <- ScaleData(CFZATRAvsCTRL, features = all.genes)

CFZATRAvsCTRL <- RunPCA(object = CFZATRAvsCTRL, npcs = 20, verbose = FALSE)
ElbowPlot(CFZATRAvsCTRL)

CFZATRAvsCTRL <- RunTSNE(object = CFZATRAvsCTRL, reduction = "pca", dims = 1:10)
CFZATRAvsCTRL <- RunUMAP(object = CFZATRAvsCTRL, reduction = "pca", dims = 1:10)
CFZATRAvsCTRL <- FindNeighbors(object = CFZATRAvsCTRL, reduction = "pca", 
                               dims = 1:10)

for (i in 1:length(res)) {
  CFZATRAvsCTRL <- FindClusters(object = CFZATRAvsCTRL, resolution = 0.6, #res[i], 
                                random.seed = 1234)
}

DimPlot(CFZATRAvsCTRL, label = TRUE)

CFZATRAvsCTRL$seurat_clusters <-CFZATRAvsCTRL$RNA_snn_res.0.6
IDCluster_CFZATRAvsCTRL <- data.frame(ID=colnames(CFZATRAvsCTRL)) %>% mutate(cluster=CFZATRAvsCTRL$seurat_clusters)
IDCluster_CFZATRAvsCTRL <- write_tsv(IDCluster_CFZATRAvsCTRL, "IDCluster_CFZATRAvsCTRL.tsv")


VlnPlot(CFZATRAvsCTRL, features=unlist(TopFeatures(CFZATRAvsCTRL[["pca"]])), ncol=3)
VlnPlot(CFZATRAvsCTRL, features= c("ENSMUSG00000091898.Tnnc1", "ENSMUSG00000026414.Tnnt2", 
                                   "ENSMUSG00000035458.Tnni3"), ncol=3)

cardiom<-FeaturePlot(object = CFZATRAvsCTRL, features= c("ENSMUSG00000091898.Tnnc1", "ENSMUSG00000026414.Tnnt2", 
                                                         "ENSMUSG00000035458.Tnni3"), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZATRAvsCTRL, features= c("ENSMUSG00000091898.Tnnc1", "ENSMUSG00000026414.Tnnt2", 
                                   "ENSMUSG00000035458.Tnni3"), ncol=3)

fibr<-FeaturePlot(object = CFZATRAvsCTRL, features= c("ENSMUSG00000045680.Tcf21", "ENSMUSG00000001506.Col1a1", 
                                                      "ENSMUSG00000029661.Col1a2"), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZATRAvsCTRL, features= c("ENSMUSG00000045680.Tcf21", "ENSMUSG00000001506.Col1a1", 
                                   "ENSMUSG00000029661.Col1a2"), ncol=3)

macrop<-FeaturePlot(object = CFZATRAvsCTRL, features= c("ENSMUSG00000024610.Cd74", "ENSMUSG00000051439.Cd14", 
                                                        "ENSMUSG00000052336.Cx3cr1"), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZATRAvsCTRL, features= c("ENSMUSG00000024610.Cd74", "ENSMUSG00000051439.Cd14", 
                                   "ENSMUSG00000052336.Cx3cr1"), ncol=3)

VEC<-FeaturePlot(object = CFZATRAvsCTRL, features= c("ENSMUSG00000002944.Cd36", "ENSMUSG00000031871.Cdh5", 
                                                     "ENSMUSG00000020717.Pecam1"), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZATRAvsCTRL, features= c("ENSMUSG00000002944.Cd36", "ENSMUSG00000031871.Cdh5", 
                                   "ENSMUSG00000020717.Pecam1"), ncol=3)

LEC<-FeaturePlot(object = CFZATRAvsCTRL, features= c("ENSMUSG00000002944.Cd36", "ENSMUSG00000031871.Cdh5", 
                                                     "ENSMUSG00000020717.Pecam1"), min.cutoff = "q9", pt.size = 0.5)

VlnPlot(CFZATRAvsCTRL, features= c("ENSMUSG00000002944.Cd36", "ENSMUSG00000031871.Cdh5", 
                                   "ENSMUSG00000020717.Pecam1"), ncol=3)


Pericy<-FeaturePlot(object = CFZATRAvsCTRL, features= c("ENSMUSG00000032135.Mcam", "ENSMUSG00000025348.Itga7" 
), min.cutoff = "q9", pt.size = 0.5)
VlnPlot(CFZATRAvsCTRL, features= c("ENSMUSG00000032135.Mcam", "ENSMUSG00000025348.Itga7"), ncol=3)


table(Idents(CFZATRAvsCTRL))

p1<-DimPlot(CFZATRAvsCTRL, label = TRUE, reduction = "umap", pt.size = 0.5, group.by= "grp") #+ NoLegend()
p2<-DimPlot(CFZATRAvsCTRL, label = TRUE, reduction = "umap", pt.size = 0.5, group.by= "seurat_clusters")

#####get marker genes###
CFZATRAvsCTRL.All.rna.markers <- FindAllMarkers(CFZATRAvsCTRL, max.cells.per.ident = 200, min.diff.pct = 0.3, test.use = "wilcox", only.pos = FALSE) #cardiomyocytes=9

All.hearts.rna.markers.cardiomyvsfibroblasts<- FindMarkers(CFZvsCFZATRA, ident.1 = 15, ident.2 = c(0,1,2,3,6) , min.pct = 0.5, grouping.var = "dataset")
All.hearts.rna.markers.cardiomyvsfibroblasts<- cbind(Row.Names = rownames(All.hearts.rna.markers.cardiomyvsfibroblasts), All.hearts.rna.markers.cardiomyvsfibroblasts)

All.hearts.rna.markers.cardiomyvsVEC<- FindMarkers(CFZvsCTRL, ident.1 = 15, ident.2 = c(4,8) , min.pct = 0.5, grouping.var = "dataset")
All.hearts.rna.markers.cardiomyvsVEC<- cbind(Row.Names = rownames(All.hearts.rna.markers.cardiomyvsVEC), All.hearts.rna.markers.cardiomyvsVEC)

All.hearts.rna.markers.cardiomyvsMacroph<- FindMarkers(CFZvsCTRL, ident.1 = 15, ident.2 = c(5,13,14) , min.pct = 0.5, grouping.var = "dataset")
All.hearts.rna.markers.cardiomyvsMacroph<- cbind(Row.Names = rownames(All.hearts.rna.markers.cardiomyvsMacroph), All.hearts.rna.markers.cardiomyvsMacroph)

All.hearts.rna.markers.cardiomyvsPericyt<- FindMarkers(CFZvsCTRL, ident.1 = 15, ident.2 = c(9,10) , min.pct = 0.5, grouping.var = "dataset")
All.hearts.rna.markers.cardiomyvsPericyt<- cbind(Row.Names = rownames(All.hearts.rna.markers.cardiomyvsPericyt), All.hearts.rna.markers.cardiomyvsPericyt)

All.hearts.rna.markers.VEC<- FindMarkers(CFZvsCFZATRA, ident.1 = c(4,8), min.pct = 0.5, grouping.var = "dataset")
All.hearts.rna.markers.VEC<- cbind(Row.Names = rownames(All.hearts.rna.markers.VEC), All.hearts.rna.markers.VEC)

All.hearts.rna.markers.macroph<- FindMarkers(CFZvsCFZATRA, ident.1 = c(5,13,14), min.pct = 0.5, grouping.var = "dataset")
All.hearts.rna.markers.macroph<- cbind(Row.Names = rownames(All.hearts.rna.markers.macroph), All.hearts.rna.markers.VEC)


p2<- FeaturePlot(object = CFZvsCTRL, features= c("ENSMUSG00000091898.Tnnc1"), min.cutoff = "q9", pt.size = 0.5)
p2

##########################################################################################################################################################################
##########################################################################################################################################################################

#####RUN :
##Fibroblasts, Endothelial and macrophages --> CFZ vs Captisol; CFZATRA vs Captisol ; CFZ vs CFZATRA ###

##################################################CFZvsCaptisol
CFZvsCaptisol.All.rna.markers <- FindAllMarkers(CFZvsCTRL, max.cells.per.ident = 200, min.diff.pct = 0.3, test.use = "wilcox", only.pos = FALSE)

CFZvsCaptisol.rna.markers.macroph<- FindMarkers(CFZvsCTRL, ident.1 = c(7,8), min.pct = 0.5, grouping.var = "grp")
CFZvsCaptisol.rna.markers.macroph<- cbind(Row.Names = rownames(CFZvsCaptisol.rna.markers.macroph), CFZvsCaptisol.rna.markers.macroph)

CFZvsCaptisol.rna.markers.Vascularendothelial<- FindMarkers(CFZvsCTRL, ident.1 = c(4), min.pct = 0.5, grouping.var = "grp")
CFZvsCaptisol.rna.markers.Vascularendothelial<- cbind(Row.Names = rownames(CFZvsCaptisol.rna.markers.Vascularendothelial), CFZvsCaptisol.rna.markers.Vascularendothelial)

CFZvsCaptisol.rna.markers.fibrobl<- FindMarkers(CFZvsCTRL, ident.1 = c(0,1,2,3,5,6), min.pct = 0.5, grouping.var = "grp")
CFZvsCaptisol.rna.markers.fibrobl<- cbind(Row.Names = rownames(CFZvsCaptisol.rna.markers.fibrobl), CFZvsCaptisol.rna.markers.fibrobl)

##################################################CFZvsCFZATRA
CFZvsCFZATRA.All.rna.markers <- FindAllMarkers(CFZvsCFZATRA, max.cells.per.ident = 200, min.diff.pct = 0.3, test.use = "wilcox", only.pos = FALSE)

CFZvsCFZATRA.rna.markers.macroph<- FindMarkers(CFZvsCFZATRA, ident.1 = c(5,7), min.pct = 0.5, grouping.var = "grp")
CFZvsCFZATRA.rna.markers.macroph<- cbind(Row.Names = rownames(CFZvsCFZATRA.rna.markers.macroph), CFZvsCFZATRA.rna.markers.macroph)

CFZvsCFZATRA.rna.markers.Vascularendothelial<- FindMarkers(CFZvsCFZATRA, ident.1 = c(3), min.pct = 0.5, grouping.var = "grp")
CFZvsCFZATRA.rna.markers.Vascularendothelial<- cbind(Row.Names = rownames(CFZvsCFZATRA.rna.markers.Vascularendothelial), CFZvsCFZATRA.rna.markers.Vascularendothelial)

CFZvsCFZATRA.rna.markers.fibrobl<- FindMarkers(CFZvsCFZATRA, ident.1 = c(0,1,2,4,6), min.pct = 0.5, grouping.var = "grp")
CFZvsCFZATRA.rna.markers.fibrobl<- cbind(Row.Names = rownames(CFZvsCFZATRA.rna.markers.fibrobl), CFZvsCFZATRA.rna.markers.fibrobl)

##################################################CFZATRAvsCTRL
CFZATRAvsCTRL.All.rna.markers <- FindAllMarkers(CFZATRAvsCTRL, max.cells.per.ident = 200, min.diff.pct = 0.3, test.use = "wilcox", only.pos = FALSE) #cardiomyocytes=9

CFZATRAvsCTRL.rna.markers.macroph<- FindMarkers(CFZATRAvsCTRL, ident.1 = c(6,8), min.pct = 0.5, grouping.var = "grp")
CFZATRAvsCTRL.rna.markers.macroph<- cbind(Row.Names = rownames(CFZATRAvsCTRL.rna.markers.macroph), CFZATRAvsCTRL.rna.markers.macroph)

CFZATRAvsCTRL.rna.markers.Vascularendothelial<- FindMarkers(CFZATRAvsCTRL, ident.1 = c(4), min.pct = 0.5, grouping.var = "grp")
CFZATRAvsCTRL.rna.markers.Vascularendothelial<- cbind(Row.Names = rownames(CFZATRAvsCTRL.rna.markers.Vascularendothelial), CFZATRAvsCTRL.rna.markers.Vascularendothelial)

CFZATRAvsCTRL.rna.markers.fibrobl<- FindMarkers(CFZATRAvsCTRL, ident.1 = c(0,1,2,3,5), min.pct = 0.5, grouping.var = "grp")
CFZATRAvsCTRL.rna.markers.fibrobl<- cbind(Row.Names = rownames(CFZATRAvsCTRL.rna.markers.fibrobl), CFZATRAvsCTRL.rna.markers.fibrobl)

##########################################################################################################################################################################
##########################################################################################################################################################################

macrophages_all <- CFZvsCaptisol.rna.markers.macroph %>%
  left_join(CFZvsCFZATRA.rna.markers.macroph, by='Row.Names') %>%  left_join(CFZATRAvsCTRL.rna.markers.macroph, by='Row.Names')

colnames(macrophages_all) <- c("Gene", "p.val_CFZvsCaptisol", "av_log2FC_CFZvsCaptisol", "pct.1_CFZ", "pct.2_Captisol", "p_val_adj_CFZvsCaptisol", 
                               "p.val_CFZvsCFZATRA", "av_log2FC_CFZvsCFZATRA", "pct.1_CFZ", "pct.2_CFZATRA", "p_val_adj_CFZvsCFZATRA", 
                               "p.val_CFZATRAvsCaptisol", "av_log2FC_CFZATRAvsCaptisol", "pct.1_CFZATRA", "pct.2_Captisol", "p_val_adj_CFZATRAvsCaptisol")


VEC_all <- CFZvsCaptisol.rna.markers.Vascularendothelial %>%
  left_join(CFZvsCFZATRA.rna.markers.Vascularendothelial, by='Row.Names') %>%  left_join(CFZATRAvsCTRL.rna.markers.Vascularendothelial, by='Row.Names')

colnames(VEC_all) <- c("Gene", "p.val_CFZvsCaptisol", "av_log2FC_CFZvsCaptisol", "pct.1_CFZ", "pct.2_Captisol", "p_val_adj_CFZvsCaptisol", 
                       "p.val_CFZvsCFZATRA", "av_log2FC_CFZvsCFZATRA", "pct.1_CFZ", "pct.2_CFZATRA", "p_val_adj_CFZvsCFZATRA", 
                       "p.val_CFZATRAvsCaptisol", "av_log2FC_CFZATRAvsCaptisol", "pct.1_CFZATRA", "pct.2_Captisol", "p_val_adj_CFZATRAvsCaptisol")


Fibroblasts_all <- CFZvsCaptisol.rna.markers.fibrobl %>%
  left_join(CFZvsCFZATRA.rna.markers.fibrobl, by='Row.Names') %>%  left_join(CFZATRAvsCTRL.rna.markers.fibrobl, by='Row.Names')

colnames(Fibroblasts_all) <- c("Gene", "p.val_CFZvsCaptisol", "av_log2FC_CFZvsCaptisol", "pct.1_CFZ", "pct.2_Captisol", "p_val_adj_CFZvsCaptisol", 
                               "p.val_CFZvsCFZATRA", "av_log2FC_CFZvsCFZATRA", "pct.1_CFZ", "pct.2_CFZATRA", "p_val_adj_CFZvsCFZATRA", 
                               "p.val_CFZATRAvsCaptisol", "av_log2FC_CFZATRAvsCaptisol", "pct.1_CFZATRA", "pct.2_Captisol", "p_val_adj_CFZATRAvsCaptisol")
write.csv(Fibroblasts_all, file= "DE.Fibroblasts_all.csv")