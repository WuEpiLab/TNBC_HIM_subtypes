#----- title: "HIM_single_cell_scRNA" ------------------------------------------
# dir.create("/public/home/scRNA/R/HIM_env", recursive = TRUE)
# renv::init(project = "/public/home/scRNA/R/HIM_env", bare = TRUE, restart = TRUE)
# renv::deactivate()

#----- loading R packages ------------------------------------------------------
library(dplyr)
library(Seurat)
library(SeuratObject)
library(lattice)
library(patchwork)
library(ggplot2)
library(magrittr)
library(pheatmap)
library(tidyverse)
library(cowplot)
library(gridExtra)

#----- set environment ---------------------------------------------------------
Sys.setenv(LANGUAGE = "en")
setwd("/public/home/scRNA/HIM-omic")
set.seed(123)

# # 使用多线程计算
# library(future)
# future::plan("multicore", workers = 2) # do parallel
# future::plan("sequential") # 切换回顺序计算
# # 设置并行计算的线程数
# num_cores <- 2  # 设置你想要的线程数
# registerDoParallel(cores = num_cores)


#----- loading data ------------------------------------------------------------
dir1 <- c("G:/HIM-omic/HIM_single cell/LC-X20220510001-TNBC001/filtered_feature_bc_matrix", # HIM253
          "G:/HIM-omic/HIM_single cell/LC-X20220530012-TNBC002/filtered_feature_bc_matrix", # HIM252
          "G:/HIM-omic/HIM_single cell/LC-X20220605001-TNBC003/filtered_feature_bc_matrix", # HIM251
          "G:/HIM-omic/HIM_single cell/LC-X20220613009-TNBC004/filtered_feature_bc_matrix", # HIM250
          "G:/HIM-omic/HIM_single cell/LC-X20220621008-TNBC005/filtered_feature_bc_matrix", # HIM248
          "G:/HIM-omic/HIM_single cell/LC-X20221114013-TNBC006/filtered_feature_bc_matrix",
          "G:/HIM-omic/HIM_single cell/LC-X20221116012-TNBC007/filtered_feature_bc_matrix",
          "G:/HIM-omic/HIM_single cell/LC-X20221111007-10XVDJ-HIM-LEN003/mRNA/HIM_LEN003_T/filtered_feature_bc_matrix") # HIM135
dir2 <- c("G:/HIM-omic/HIM_single cell/rawdata/CID3946",
          "G:/HIM-omic/HIM_single cell/rawdata/CID3963", 
          "G:/HIM-omic/HIM_single cell/rawdata/CID4465", 
          "G:/HIM-omic/HIM_single cell/rawdata/CID4495",
          "G:/HIM-omic/HIM_single cell/rawdata/CID4513",
          "G:/HIM-omic/HIM_single cell/rawdata/CID4515",
          "G:/HIM-omic/HIM_single cell/rawdata/CID4523",
          "G:/HIM-omic/HIM_single cell/rawdata/CID44041",
          "G:/HIM-omic/HIM_single cell/rawdata/CID44971",
          "G:/HIM-omic/HIM_single cell/rawdata/CID44991")

dir3 <- c("G:/HIM-omic/HIM_single cell/rawdata/TN-MH0126",
          "G:/HIM-omic/HIM_single cell/rawdata/TN-MH0135",
          "G:/HIM-omic/HIM_single cell/rawdata/TN-SH0106",
          "G:/HIM-omic/HIM_single cell/rawdata/TN-MH0114-T2",
          "G:/HIM-omic/HIM_single cell/rawdata/TN-B1-Tum0554",
          "G:/HIM-omic/HIM_single cell/rawdata/TN-B1-MH4031",
          "G:/HIM-omic/HIM_single cell/rawdata/TN-B1-MH0177",
          "G:/HIM-omic/HIM_single cell/rawdata/TN-B1-MH0131")

names(dir1) <- c("TNBC001","TNBC002","TNBC003","TNBC004","TNBC005","TNBC006","TNBC007","LEN003")
names(dir2) <- c("CID3946","CID3963","CID4465","CID4495","CID4513","CID4515","CID4523","CID44041","CID44971","CID44991")
names(dir3) <- c("TN-MH0126","TN-MH0135","TN-SH0106","TN-MH0114-T2","TN-B1-Tum0554","TN-B1-MH4031","TN-B1-MH0177","TN-B1-MH0131")

scRNAlist1 <- list()
scRNAlist2 <- list()
scRNAlist3 <- list()

#----- CreateSeuratObject
#----- dir1
# nFeature_RNA代表每个细胞测到的基因数目，nCount代表每个细胞测到所有基因的表达量之和
# min_cells指基因至少在1个细胞中表达
for (i in 1:length(dir1)) {
  counts1 <- Read10X(data.dir = dir1[i])
  scRNAlist1[[i]] = CreateSeuratObject(counts1,project = names(dir1[i]),min_cells = 1, min_ngene = 500, max_ngene = inf)
  scRNAlist1[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist1[[i]],pattern = "^MT-")}
#----- dir2
for (i in 1:length(dir2)) {
  counts2 <- Read10X(data.dir = dir2[i],gene.column = 1)
  scRNAlist2[[i]] = CreateSeuratObject(counts2,project = names(dir2[i]),min_cells = 1, min_ngene = 500, max_ngene = inf)
  scRNAlist2[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist2[[i]],pattern = "^MT-")}
#----- dir3
for (i in 1:length(dir3)) {
  counts3 <- Read10X(data.dir = dir3[i],gene.column = 2)
  scRNAlist3[[i]] = CreateSeuratObject(counts3,project = names(dir3[i]),min_cells = 1, min_ngene = 500, max_ngene = inf)
  scRNAlist3[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist3[[i]],pattern = "^MT-")}

#----- remove aptmer ----------------------------------------------------------#
tail(rownames(scRNAlist1[[1]]),3)
tail(rownames(scRNAlist1[[2]]),3)
tail(rownames(scRNAlist1[[3]]),3)
tail(rownames(scRNAlist1[[4]]),75)
tail(rownames(scRNAlist1[[5]]),75)
tail(rownames(scRNAlist1[[6]]),3)
tail(rownames(scRNAlist1[[7]]),3)
tail(rownames(scRNAlist1[[8]]),3)

# remove aptmer
scRNAlist1[[4]] <- subset(scRNAlist1[[4]], features = rownames(scRNAlist1[[4]])[1:(length(rownames(scRNAlist1[[4]] ))-73)])
scRNAlist1[[5]] <- subset(scRNAlist1[[5]], features = rownames(scRNAlist1[[5]])[1:(length(rownames(scRNAlist1[[5]] ))-73)])
tail(rownames(scRNAlist1[[3]]),3)
tail(rownames(scRNAlist1[[4]]),3)
tail(rownames(scRNAlist1[[5]]),3)


#----- before quality control -------------------------------------------------#
p1 <- VlnPlot(scRNAlist1[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 <- VlnPlot(scRNAlist1[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p3 <- VlnPlot(scRNAlist1[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p4 <- VlnPlot(scRNAlist1[[4]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p5 <- VlnPlot(scRNAlist1[[5]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p6 <- VlnPlot(scRNAlist1[[6]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p7 <- VlnPlot(scRNAlist1[[7]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p8 <- VlnPlot(scRNAlist1[[8]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plotc = (p1|p2)/(p3|p4)/(p5|p6)/(p7|p8) #(p9|plot_spacer())
ggsave(plotc,filename = "figure/inhouse_pre_control.png",height = 14,width = 18)

#----- scRNAlist2
p1 <- VlnPlot(scRNAlist2[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 <- VlnPlot(scRNAlist2[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p3 <- VlnPlot(scRNAlist2[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p4 <- VlnPlot(scRNAlist2[[4]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p5 <- VlnPlot(scRNAlist2[[5]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p6 <- VlnPlot(scRNAlist2[[6]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p7 <- VlnPlot(scRNAlist2[[7]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p8 <- VlnPlot(scRNAlist2[[8]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p9 <- VlnPlot(scRNAlist2[[9]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p10 <- VlnPlot(scRNAlist2[[10]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plotc = (p1|p2)/(p3|p4)/(p5|p6)/(p7|p8)/(p9|p10)
ggsave(plotc,filename = "figure/CID_pre_control.png",height = 14,width = 18)

#----- scRNAlist3
p1 <- VlnPlot(scRNAlist3[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 <- VlnPlot(scRNAlist3[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p3 <- VlnPlot(scRNAlist3[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p4 <- VlnPlot(scRNAlist3[[4]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p5 <- VlnPlot(scRNAlist3[[5]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p6 <- VlnPlot(scRNAlist3[[6]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p7 <- VlnPlot(scRNAlist3[[7]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p8 <- VlnPlot(scRNAlist3[[8]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plotc = (p1|p2)/(p3|p4)/(p5|p6)/(p7|p8)
ggsave(plotc,filename = "figure/TN_pre_control.png",height = 14,width = 18)


#----- quality control --------------------------------------------------------#
for (i in 1:length(dir1)) {scRNAlist1[[i]] <- subset(scRNAlist1[[i]],subset = percent.mt < 25)}
for (i in 1:length(dir2)) {scRNAlist2[[i]] <- subset(scRNAlist2[[i]],subset = percent.mt < 25)}
for (i in 1:length(dir3)) {scRNAlist3[[i]] <- subset(scRNAlist3[[i]],subset = percent.mt < 25)}


#----- agter quality control -------------------------------------------------#
p1 <- VlnPlot(scRNAlist1[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 <- VlnPlot(scRNAlist1[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p3 <- VlnPlot(scRNAlist1[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p4 <- VlnPlot(scRNAlist1[[4]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p5 <- VlnPlot(scRNAlist1[[5]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p6 <- VlnPlot(scRNAlist1[[6]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p7 <- VlnPlot(scRNAlist1[[7]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p8 <- VlnPlot(scRNAlist1[[8]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plotc = (p1|p2)/(p3|p4)/(p5|p6)/(p7|p8)#/(p9|plot_spacer())
ggsave(plotc,filename = "figure/inhouse_after_control.png",height = 14,width = 18)

#----- scRNAlist2
p1 <- VlnPlot(scRNAlist2[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 <- VlnPlot(scRNAlist2[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p3 <- VlnPlot(scRNAlist2[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p4 <- VlnPlot(scRNAlist2[[4]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p5 <- VlnPlot(scRNAlist2[[5]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p6 <- VlnPlot(scRNAlist2[[6]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p7 <- VlnPlot(scRNAlist2[[7]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p8 <- VlnPlot(scRNAlist2[[8]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p9 <- VlnPlot(scRNAlist2[[9]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p10 <- VlnPlot(scRNAlist2[[10]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plotc = (p1|p2)/(p3|p4)/(p5|p6)/(p7|p8)/(p9|p10)
ggsave(plotc,filename = "figure/CID_after_control.png",height = 14,width = 18)

#----- scRNAlist3
p1 <- VlnPlot(scRNAlist3[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 <- VlnPlot(scRNAlist3[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p3 <- VlnPlot(scRNAlist3[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p4 <- VlnPlot(scRNAlist3[[4]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p5 <- VlnPlot(scRNAlist3[[5]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p6 <- VlnPlot(scRNAlist3[[6]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p7 <- VlnPlot(scRNAlist3[[7]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p8 <- VlnPlot(scRNAlist3[[8]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plotc = (p1|p2)/(p3|p4)/(p5|p6)/(p7|p8)
ggsave(plotc,filename = "figure/TN_after_control.png",height = 14,width = 18)






#----- integration data --------------------------------------------------------
# using merge 
list1 <- scRNAlist1[2:length(scRNAlist1)]
integrate1 <- merge(scRNAlist1[[1]],list1)

list2 <- scRNAlist2[2:length(scRNAlist2)]
integrate2 <- merge(scRNAlist2[[1]],list2)

list3 <- scRNAlist3[2:length(scRNAlist3)]
integrate3 <- merge(scRNAlist3[[1]],list3)

tmp <- merge(integrate1,integrate2)
HIM_integrate <- merge(tmp,integrate3)

saveRDS(HIM_integrate,file = "HIM_integrate.rds")



HIM_integrate <-  readRDS(file = "HIM_integrate.rds")

#----- remove LEN002
HIM_integrate <- HIM_integrate[,HIM_integrate@meta.data$orig.ident != "LEN002"]

#----- HIM_intergrate workflow ---------------------------------
# https://zhuanlan.zhihu.com/p/493477082
# parallelly::availableCores()
# future::plan("multicore", workers = 2) # use two core

# add data batch 
batch1 <- c("TNBC001","TNBC002","TNBC003","TNBC004","TNBC005","TNBC006","TNBC007","LEN003")
batch2 <- c("CID3946","CID3963","CID4465","CID4495","CID4513","CID4515","CID4523","CID44041","CID44971","CID44991")
batch3 <- c("TN-MH0126","TN-MH0135","TN-SH0106","TN-MH0114-T2","TN-B1-Tum0554","TN-B1-MH4031","TN-B1-MH0177","TN-B1-MH0131")

batch <- ifelse(HIM_integrate$orig.ident %in% batch1,"batch1",
                ifelse(HIM_integrate$orig.ident %in% batch2,"batch2",
                       ifelse(HIM_integrate$orig.ident %in% batch3,"batch3","reset")))
table(batch)
HIM_integrate$batch <- batch


HIM_integrate <- NormalizeData(HIM_integrate)
HIM_integrate <- FindVariableFeatures(HIM_integrate, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(pbmc)
HIM_integrate <- ScaleData(HIM_integrate, features = VariableFeatures(HIM_integrate))
HIM_integrate <- RunPCA(HIM_integrate, verbose = T)
library(harmony)
HIM_integrate <- RunHarmony(HIM_integrate, group.by.vars = 'batch') # "orig.ident"
ElbowPlot(object = HIM_integrate, ndims = 50)
HIM_integrate <- FindNeighbors(HIM_integrate, reduction = "harmony", dims = 1:30)
library(clustree)
HIM_integrate <- FindClusters(object = HIM_integrate,resolution = c(seq(.1,1.6,.2)))
p1 <- clustree(HIM_integrate@meta.data, prefix = "RNA_snn_res.")
ggsave(p1, filename = "figure/clustree.png",width = 12,height = 12)

HIM_integrate <- RunUMAP(HIM_integrate, reduction = "harmony", dims = 1:30)
HIM_integrate <- RunTSNE(HIM_integrate, reduction = "harmony", dims = 1:30)


p1 <- DimPlot(HIM_integrate, reduction = "tsne", group.by = "batch",label = F) + 
  scale_color_manual(values=c("batch1" = "#B3C254",
                              "batch2" = "#F46328",
                              "batch3" = "#FA8721")) + NoLegend()+ 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))

p2 <- DimPlot(HIM_integrate, reduction = "tsne", group.by = "orig.ident",label = F) + NoLegend()+ 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))

ggsave(p1, filename = "HIM-integrate/figure/HIM_integrate_batch.png",width = 8, height = 8)
ggsave(p2, filename = "HIM-integrate/figure/HIM_integrate_orig_ident.pdf",width = 6, height = 6)


saveRDS(HIM_integrate,file = "HIM_integrate.rds")



HIM_integrate <-  readRDS(file = "HIM_integrate.rds")




#----- pesudo bulk prediction HIM subtype --------------------------------------
HIM_pesudo <- AverageExpression(HIM_integrate, group.by = "orig.ident")
HIM_pesudo <- as.data.frame(HIM_pesudo)

write.csv(HIM_pesudo,file = "HIM_pesudo.csv")

HIM_pesudo <- read.csv(file = "HIM_pesudo.csv")




#----- add HIM subtype 
table(HIM_integrate$orig.ident)
subtype <- read.csv(file = "PamrRes_scRNAseq_pesudo.csv")

CS1 <- subtype[subtype$clust == "CS1",]$X
H <- colnames(HIM_integrate[,HIM_integrate@meta.data$orig.ident %in% CS1])
# "CID3963","CID4495","CID44971","CID4513","TN-B1-MH0177","TN-B1-Tum0554","TN-MH0114-T2","TN-SH0106",
# "TNBC001","TNBC004","TNBC005" 
# 
CS2 <- subtype[subtype$clust == "CS2",]$X
M <- colnames(HIM_integrate[,HIM_integrate@meta.data$orig.ident %in% CS2])
# "CID3946","CID44041","CID4465","LEN003","TNBC002","TNBC003","TNBC006"
CS3 <- subtype[subtype$clust == "CS3",]$X
I <- colnames(HIM_integrate[,HIM_integrate@meta.data$orig.ident %in% CS3])
# "CID44991","CID4515","CID4523","TN-B1-MH0131","TN-B1-MH4031","TN-MH0126","TN-MH0135","TNBC007"
group <- ifelse(colnames(HIM_integrate) %in% H,"H",
                ifelse(colnames(HIM_integrate) %in% M,"M",
                       ifelse(colnames(HIM_integrate) %in% I,"I","rest")))

cluster <- ifelse(colnames(HIM_integrate) %in% H,"CS1",
                  ifelse(colnames(HIM_integrate) %in% M,"CS2",
                         ifelse(colnames(HIM_integrate) %in% I,"CS3","rest")))
table(group)
table(cluster)
HIM_integrate@meta.data$group <- group
HIM_integrate@meta.data$cluster <- cluster


saveRDS(HIM_integrate,file = "HIM_integrate.rds")



HIM_integrate <-  readRDS(file = "HIM_integrate.rds")







#----- singleR cluster --------------------------------------------------
library(SingleR)
library(celldex)
# BlueprintEncodeData (Obtain human bulk RNA-seq data from Blueprint and ENCODE)
# DatabaseImmuneCellExpressionData (Obtain human bulk RNA-seq data from DICE)
# HumanPrimaryCellAtlasData (Obtain the HPCA data)
# MonacoImmuneData (Obtain bulk RNA-seq data of sorted human immune cells)
# NovershternHematopoieticData (Obtain bulk microarray expression for sorted hematopoietic cells)

# ImmGenData (Obtain mouse bulk expression data from the Immunologic Genome Project)
# MouseRNAseqData (Obtain mouse bulk expression data of sorted cell populations (RNA-seq)

ref <- readRDS(file = "H:/referencedata/celldex/BlueprintEncodeData.rds")
tmp <- GetAssayData(HIM_integrate,slot = "data")
cluster <- HIM_integrate@meta.data$RNA_snn_res.0.3
singler <- SingleR(test = tmp, ref = ref, labels = ref$label.main, clusters = cluster, fine.tune = T)
new.cluster.ids <- singler$first.labels
HIM_integrate@active.ident <- HIM_integrate$seurat_clusters
names(new.cluster.ids) <- levels(HIM_integrate)
HIM_integrate <- RenameIdents(HIM_integrate, new.cluster.ids)
HIM_integrate@meta.data$singleR_cluster <- HIM_integrate@active.ident # add singleR celltype to meta.data

p1 <- DimPlot(HIM_integrate, reduction = "umap",label = T)
ggsave(p1, filename = "figure/singleR.png",width = 8, height = 6)





#----- identify cell type ------------------------------------------------------
# 对HIM_integrate重命名
table(Tcells$cell_type)
for (i in c(as.character(levels(Tcells$cell_type)))){
  tmp <- rownames(Tcells@meta.data[Tcells$cell_type == i,])
  HIM_integrate$re_celltype[rownames(HIM_integrate@meta.data) %in% tmp] <- i
}
table(Tcells$cell_type)
table(HIM_integrate$re_celltype)

for (i in c(as.character(levels(Mac$cell_type)))){
  tmp <- rownames(Mac@meta.data[Mac$cell_type == i,])
  HIM_integrate$re_celltype[rownames(HIM_integrate@meta.data) %in% tmp] <- i
}
table(Mac$cell_type)
table(HIM_integrate$re_celltype)






table(stromal$cell_type)
for (i in c(as.character(levels(stromal$cell_type)))){
  tmp <- rownames(stromal@meta.data[stromal$cell_type == i,])
  HIM_integrate$re_celltype[rownames(HIM_integrate@meta.data) %in% tmp] <- i
}
table(stromal$cell_type)
table(HIM_integrate$re_celltype)

HIM_integrate$cell_type[HIM_integrate@meta.data$cell_type == 'cDCs'] <- 'Macrophages'


# cancer_cell <- c("CD24","KRT19","SCGB2A2","EPCAM", "PROM1", "ALDH1A1")
# basal ("KRT14","MYLK")
# T_cell <- c("CD2", "CD3D", "CD3E", "CD3G", "CD8A")
# NK_cell <- c("GNLY", "KLRD1", "NCR1","XCL1","NCAM1","FGFBP2", "FCG3RA", "CX3CR1","NKG7")
# B_cell <- c("CD79A","MZB1", "MS4A1","CD19", "CD79B")
# Plasma <- c("JCHAIN", "IGHG1", "SDC1")
# myeloid_cell <- c("LYZ", "CD68", "TYROBP")
# macrophages <- c("CD163", "CSF1R", "CD14")
# Granulocyte <- c("ENPP3","FUT4","CSF3R")
# mast_cell <- c("CPA3","TPSAB1","TPSB2")
# cDCs <- c("CLEC10A", "CD1C", "CD1E")
# pDCs<- c("IL3RA","LILRA4","CXCR3","IRF7")
# endothelial_cells <- c("PECAM1", "S100A4","CLDN5","FLT1","RAMP2", "VWF")
# fibroblasts <- c("COL1A1", "COL1A2", "COL3A1", "PDGFRB", "THY1", "ITGB1","DCN","C1R","FGF7", "MME", 'LUM', 'GSN')

# cluster_markers <- c("CD24","KRT19","SCGB2A2","EPCAM", "PROM1", "ALDH1A1",
#                      "PTPRC","CD2", "CD3D", "CD3E", "CD3G", "CD8A",
#                      "GNLY", "KLRD1", "NCR1","XCL1","NCAM1","FGFBP2", "CX3CR1","NKG7",
#                      "CD79A","MZB1", "MS4A1","CD19", "CD79B",
#                      "JCHAIN", "IGHG1", "SDC1",
#                      "LYZ", "CD68", "TYROBP",
#                      "CD163", "CSF1R", "CD14",
#                      "ENPP3","FUT4","CSF3R",
#                      "CPA3","TPSAB1","TPSB2",
#                      "CLEC10A", "CD1C", "CD1E",
#                      "IL3RA","LILRA4","CXCR3","IRF7",
#                      "PECAM1", "S100A4","CLDN5","FLT1","RAMP2", "VWF",
#                      "COL1A1", "COL1A2", "COL3A1", "PDGFRB", "THY1", "ITGB1","DCN","C1R","FGF7", "MME", 'LUM', 'GSN',
#                      'C1QA','C1QB','C1QC')

# HIM_integrate@active.ident <- HIM_integrate$RNA_snn_res.1.5
# marker_plot <- DotPlot(HIM_integrate,features = cluster_markers,group.by = "RNA_snn_res.0.5") + coord_flip()
# marker_plot
# ggsave(marker_plot, filename = 'figure/cluster marker.png',  width = 14, height = 10,bg = "white")

p1 <- DimPlot(HIM_integrate, group.by = "RNA_snn_res.0.5",raster = FALSE,label = T, reduction = "tsne") +
  scale_color_manual(values=c("Cancer cells" = '#90c9f8', "Endothelial cells" = '#ffcf6f',"CAFs" = '#ffab6f', "Mast cells" = '#f9869e',
                              "pDCs" ='#f9869e', "cDCs" ='#f9869e', "Macrophages" = '#f9869e', "B cells" = '#f9869e',
                              "NK cells" = '#f9869e', "T cells" = '#f9869e'))+ NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))

p2 <- DimPlot(HIM_integrate, group.by = "cell_type",raster = FALSE,label = F, reduction = "tsne")  +
  scale_color_manual(values=c("Cancer cells" = '#2EC4B6', "Endothelial cells" = '#FF9F1C',"CAFs" = '#ffab6f', "Mast cells" = '#fad5da',
                              "pDCs" ='#f5a7b1', "cDCs" ='#ee7362', "Macrophages" = '#ee6297', "B cells" = '#ee6274',
                              "NK cells" = '#ec4b60', "T cells" = '#E71D36'))+ NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))
ggsave(p1, filename = "figure/HIM_integrate_tsne.pdf",width = 8, height = 8)
ggsave(p1, filename = "figure/HIM_integrate_tsne.png",width = 8, height = 8)
ggsave(p2, filename = "figure/HIM_integrate_cell_type.pdf",width = 8, height = 6)
ggsave(p2, filename = "figure/HIM_integrate_cell_type.png",width = 8, height = 8)

saveRDS(HIM_integrate,file = "HIM_integrate.rds")

HIM_integrate <- readRDS(file = "HIM-integrate/rds_file/HIM_integrate.rds")







#----- marker gene heat-map ---------------------------------------------------
cluster_markers <- c("CD24","KRT19","SCGB2A2","EPCAM", "PROM1",
                     "PTPRC","CD2", "CD3D", "CD3E", "CD3G", "CD8A","CXCL13",
                     "GNLY", "KLRD1", "NCR1","XCL1","NCAM1","FGFBP2", "CX3CR1","NKG7",
                     "CD79A","MZB1", "MS4A1","CD19", "CD79B","IGHG1",
                     "LYZ", "CD68", "TYROBP","CD163", "CSF1R", "CD14","CSF3R","FUT4",
                     "CLEC10A", "CD1C", "CD1E","C1QA", "C1QB", "C1QC",
                     "IL3RA","LILRA4","CXCR3","IRF7",
                     "ENPP3","CPA3","TPSAB1","TPSB2",
                     "LEPR","PECAM1","CD34","VWF","VEGFC","MCAM", "CAV1","ACKR1","SELE", "SELP","VCAM1",
                     "PDPN","ACTA2", "TAGLN", "MYH11", "MYLK","PDGFRB", "RGS5",
                     "FAP","COL1A1", "COL1A2","POSTN","DCN", "FIGF", "C3", "IGF1")

library(pheatmap)
for (i in 1:length(cluster_markers)) {
  ids = HIM_integrate[cluster_markers[i],]
  ids = ids@assays$RNA@data %>% as.numeric()
  assign(cluster_markers[i],tapply(ids, HIM_integrate@meta.data$cell_type,mean))
}

heatmap_matrix <- NULL
for (i in 1:length(cluster_markers)) {
  heatmap_matrix = rbind(heatmap_matrix,get(cluster_markers[i]))
}

row.names(heatmap_matrix) = cluster_markers
ids = c("Cancer cells","T cells","NK cells","B cells","Macrophages","pDCs","Mast cells",
        'Endothelial cells','CAFs')
heatmap_matrix = heatmap_matrix[,ids]

col <- c(colorRampPalette(c('#3ab2e8','white'))(36),
         colorRampPalette(c('white','#da3446'))(36))

p1 <- pheatmap(heatmap_matrix,color = col,scale = 'row',cluster_cols = F,cluster_rows = F) #breaks = unique(c(seq(-2.5,2.5, length=100)))
ggsave(p1, filename = '/public/home/scRNA/HIM-omic/HIM_single_cell/result/HIM-integrate/figure/HIM_integrate_markers.pdf',height = 12,width = 4)




#----- HIM subtypes Cellratio --------------------------------------------------
#----- bar plot of each patients
Cellratio <- as.data.frame(prop.table(table(HIM_integrate$cell_type,HIM_integrate$orig.ident),margin = 2))
Cellratio$Var1 <- factor(Cellratio$Var1, levels = c("Cancer cells","Endothelial cells","CAFs","Mast cells",
                                                    "pDCs","Macrophages","B cells","NK cells","T cells"))
Cellratio$Var2 <- factor(Cellratio$Var2, levels = c("CID3963","CID4495","CID44971","CID4513","LEN002","TN-B1-MH0177","TN-B1-Tum0554",
                                                    "TN-MH0114-T2","TN-SH0106","TNBC001","TNBC004","TNBC005","TNBC003",
                                                    "CID3946","CID44041","CID4465","LEN003","TNBC002","TNBC006",
                                                    "CID44991","CID4515","CID4523","TN-B1-MH0131","TN-B1-MH4031","TN-MH0126","TN-MH0135","TNBC007"))
library(ggsci)
colors <- c('#2EC4B6','#FF9F1C','#ffab6f','#fad5da', '#f5a7b1',
            '#ee6297','#ee6274','#ec4b60',
            '#E71D36')
p1 <- ggplot(Cellratio,aes(x = Var2, y= Freq,fill = Var1)) +
  geom_bar(stat ="identity", position = "stack")+
  labs(x='sample',y = 'ratio') +
  RotatedAxis()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

plotc <- p1 + scale_fill_manual(values=colors)
ggsave(plotc, filename = "figure/HIM ident cell ratio.pdf",width = 6, height = 5)



#----- top and right bar plot 
Cellratio <- as.data.frame(prop.table(table(HIM_integrate$cell_type,HIM_integrate$cluster),margin = 2))

cluster_number <- as.data.frame(table(HIM_integrate$cell_type,HIM_integrate$cluster))
cell_number <- as.data.frame(table(HIM_integrate$cluster, HIM_integrate$cell_type))

colourCount = length(unique(cluster_number$Var1))
cluster_number$Var1 <- factor(cluster_number$Var1, levels = c("Cancer cells","Endothelial cells","CAFs","Mast cells", 
                                                              "pDCs","Macrophages","B cells","NK cells","T cells"))
cluster_number$Var2 <- factor(cluster_number$Var2, levels = c("CS1","CS2","CS3"))
library(ggsci)
colors <- c('#2EC4B6','#FF9F1C','#ffab6f','#fad5da', '#f5a7b1',
            '#ee7362', '#ee6297','#ee6274','#ec4b60',
            '#E71D36')
top <- ggplot(cluster_number,aes(x = Var2, y= Freq,fill = Var1)) +
  geom_bar(stat ="identity", position = "stack")+
  labs(x='',y = '') +
  RotatedAxis()+
  theme_classic()+
  theme(axis.text.x = NULL) # +coord_flip() element_text(angle = 45,hjust = 1)

# top + coord_polar(theta = 'y') + scale_fill_manual(values=colors)+ facet_wrap(~Var2,ncol=3) # 绘制环形图

plotc <- top + scale_fill_manual(values=colors)
ggsave(plotc, filename = "figure/HIM ident cell counts top.pdf",width = 5, height = 5)


right <- ggplot(cluster_number,aes(x= Freq,y = Var1,fill = Var2)) +
  geom_bar(stat ="identity", position = "stack")+
  labs(x='',y = '') +
  RotatedAxis()+
  theme_classic()+
  theme(axis.text.x = NULL) # +coord_flip() element_text(angle = 45,hjust = 1)

plotc <- right + scale_fill_manual(values=c('CS1'='#E71D36','CS2'='#FF9F1C','CS3'='#2EC4B6'))
ggsave(plotc, filename = "figure/HIM ident cell counts right.pdf",width = 5, height = 5)


#----- bubble plot
Cellratio$Var1 <- factor(Cellratio$Var1, levels = c("Cancer cells", "Endothelial cells" ,"CAFs", "Mast cells",
                                                    "pDCs", "Macrophages", "B cells",
                                                    "NK cells", "T cells"))


# 创建基础绘图对象
color_values <- c('#2EC4B6','#FF9F1C','#ffab6f','#fad5da', '#f5a7b1',
                  '#ee7362', '#ee6297','#ee6274','#ec4b60',
                  '#E71D36')

plot <- ggplot(Cellratio, aes()) + 
  geom_point(aes(x = Var2, y = Var1, size = Freq,color = Var1)) +
  scale_color_manual(values = color_values)+
  scale_size_continuous(range = c(2, 15)) +
  theme_classic() 

ggsave(plot, filename = "figure/HIM ident cell bubble.pdf",width = 4, height = 5)





#----- bar plot of each cluster
Cellratio <- as.data.frame(prop.table(table(HIM_integrate$cell_type,HIM_integrate$cluster),margin = 2))
Cellratio$Var1 <- factor(Cellratio$Var1, levels = c("Cancer cells","Endothelial cells","CAFs","Mast cells",
                                                    "pDCs","Macrophages","B cells","NK cells","T cells"))
Cellratio$Var2 <- factor(Cellratio$Var2, levels = c("CS1","CS2","CS3"))

library(ggsci)
colors <- c('#2EC4B6','#FF9F1C','#ffab6f','#fad5da', '#f5a7b1',
            '#ee6297','#ee6274','#ec4b60',
            '#E71D36')
p1 <- ggplot(Cellratio,aes(x = Var2, y= Freq,fill = Var1)) +
  geom_bar(stat ="identity", position = "stack")+
  labs(x='sample',y = 'ratio') +
  RotatedAxis()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

plotc <- p1 + scale_fill_manual(values = colors)
ggsave(plotc, filename = "figure/HIM cluster cell ratio.pdf",width = 6, height = 5)




#----- 饼图
library(tidyverse)
library(PieGlyph)
library(dplyr)
library(tidyr)


data <- Cellratio %>% pivot_wider(names_from = Var2,
                                  values_from = Freq)# 用pivot_wider将长格式转换为宽格式
desired_order <- c("Cancer cells","Endothelial cells","CAFs","Mast cells",
                   "pDCs","Macrophages","B cells","NK cells","T cells")
data <- data %>%
  filter(Var1 %in% desired_order) %>%
  arrange(factor(Var1, levels = desired_order))

pdf(file = "/public/home/scRNA/HIM-omic/HIM_single_cell/result/HIM-integrate/figure/CS3-pie.pdf",width = 6, height = 5)
pie(data$CS3, labels = data$Var1, col = colors)
dev.off()










#----- HIM_integrate membrane protein score --------------------------------------------------
load('/public/home/scRNA/HIM-omic/Proteogenomic/in_house/result/fuclust_mem_3sub.rda') #膜蛋白 一致性聚类 的三分型 的fuclust # 分类结果
load('/public/home/scRNA/HIM-omic/Proteogenomic/in_house/result/HIM_proteomic_data_v2.rda') # pt蛋白表达数据
gs_tcsa <- read.csv("/public/home/scRNA/HIM-omic/Proteogenomic/in_house/result/HIM_subtype/TCSA_membrane_gene_list.CSV",check.names = F) #3567
CS1_deg <- read.csv("/public/home/scRNA/HIM-omic/Proteogenomic/in_house/result/Specific_gene_of_CS1.csv",row.names = 1)
CS2_deg <- read.csv("/public/home/scRNA/HIM-omic/Proteogenomic/in_house/result/Specific_gene_of_CS2.csv",row.names = 1)
CS3_deg <- read.csv("/public/home/scRNA/HIM-omic/Proteogenomic/in_house/result/Specific_gene_of_CS3.csv",row.names = 1)

CS1_surface <- CS1_deg[rownames(CS1_deg) %in% gs_tcsa$`HGNC Symbol`,]
CS2_surface <- CS2_deg[rownames(CS2_deg) %in% gs_tcsa$`HGNC Symbol`,]
CS3_surface <- CS3_deg[rownames(CS3_deg) %in% gs_tcsa$`HGNC Symbol`,]

CS1_surface <- CS1_surface[CS1_surface$logFC > 0.5 & CS1_surface$P.Value < 0.05,]
CS2_surface <- CS2_surface[CS2_surface$logFC > 0.5 & CS2_surface$P.Value < 0.05,]
CS3_surface <- CS3_surface[CS3_surface$logFC > 0.5 & CS3_surface$P.Value < 0.05,]

CS1_surface <-  CS1_surface[order(CS1_surface$logFC, decreasing = T),][1:30,]
CS2_surface <-  CS2_surface[order(CS2_surface$logFC, decreasing = T),][1:30,]
CS3_surface <-  CS3_surface[order(CS3_surface$logFC, decreasing = T),][1:30,]


CS1_UP <- list(row.names(CS1_surface))
CS2_UP <- list(row.names(CS2_surface))
CS3_UP <- list(row.names(CS3_surface))

save(CS1_UP,CS2_UP,CS3_UP,file = '/public/home/scRNA/HIM-omic/Proteogenomic/in_house/result/cluster_surface')

load(file = '/public/home/scRNA/HIM-omic/Proteogenomic/in_house/result/cluster_surface')

HIM_integrate <- AddModuleScore(object = HIM_integrate, features = CS1_UP,name = "CS1_UP")
HIM_integrate <- AddModuleScore(object = HIM_integrate, features = CS2_UP,name = "CS2_UP")
HIM_integrate <- AddModuleScore(object = HIM_integrate, features = CS3_UP,name = "CS3_UP")

p1 <- FeaturePlot(object = HIM_integrate,reduction = 'tsne', features = "CS1_UP1",raster=FALSE, cols = c('grey','#ee6274','#E71D36')) + NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))
p2 <- FeaturePlot(object = HIM_integrate,reduction = 'tsne', features = "CS2_UP1",raster=FALSE, cols = c('grey','#ffbf69','#FF9F1C')) + NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))
p3 <- FeaturePlot(object = HIM_integrate,reduction = 'tsne', features = "CS3_UP1",raster=FALSE, cols = c('grey','#64dbd0','#2EC4B6')) + NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))

ggsave(p1,filename = "HIM-integrate/figure/CS1 up AddModuleScore in HIM_integrate.pdf",width = 6, height = 6) 
ggsave(p2,filename = "HIM-integrate/figure/CS2 up AddModuleScore in HIM_integrate.pdf",width = 6, height = 6)
ggsave(p3,filename = "HIM-integrate/figure/CS3 up AddModuleScore in HIM_integrate.pdf",width = 6, height = 6) 


#----- AUCells
library(AUCell)
# https://zhuanlan.zhihu.com/p/482523999
names(CS1_UP) <- 'CS1_up'
names(CS2_UP) <- 'CS2_up'
names(CS3_UP )<- 'CS3_up'
cells_rankings <- AUCell_buildRankings(HIM_integrate@assays$RNA@data,splitByBlocks=TRUE) 

CS1_AUC <- AUCell_calcAUC(CS1_UP, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
CS2_AUC <- AUCell_calcAUC(CS2_UP, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
CS3_AUC <- AUCell_calcAUC(CS3_UP, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

HIM_integrate$CS1_AUC  <- as.numeric(getAUC(CS1_AUC)['CS1_up',])
HIM_integrate$CS2_AUC  <- as.numeric(getAUC(CS2_AUC)['CS2_up',])
HIM_integrate$CS3_AUC  <- as.numeric(getAUC(CS3_AUC)['CS3_up',])


p1 <- FeaturePlot(object = HIM_integrate, features = "CS1_AUC",reduction = "tsne",raster = F) +
  scale_color_gradientn(colours = c('grey','#ffb3d2','red')) +
  theme(legend.position = "right") +
  labs(x= NULL,y = NULL)+
  ggtitle('cs1 up score AUCell')

p1 <- FeaturePlot(object = HIM_integrate, features = "CS1_AUC",reduction = "tsne",raster = F) +
  scale_color_gradientn(colours = c('grey','#ffb3d2','red')) + NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))

p2 <- FeaturePlot(object = HIM_integrate, features = "CS2_AUC",reduction = "tsne",raster = F) +
  scale_color_gradientn(colours = c('grey','#ffd59c','#d07600')) +
  theme(legend.position = "right") +
  labs(x= NULL,y = NULL)+
  ggtitle('cs2 up score AUCell')

p2 <- FeaturePlot(object = HIM_integrate, features = "CS2_AUC",reduction = "tsne",raster = F)  +
  scale_color_gradientn(colours = c('grey','#ffd59c','#d07600')) + NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))

p3 <- FeaturePlot(object = HIM_integrate, features = "CS3_AUC",reduction = "tsne",raster = F) +
  scale_color_gradientn(colours = c('grey','#64dbd0','#2EC4B6')) +
  theme(legend.position = "right") +
  labs(x= NULL,y = NULL)+
  ggtitle('cs3 up score AUCell')

p3 <- FeaturePlot(object = HIM_integrate, features = "CS3_AUC",reduction = "tsne",raster = F) +
  scale_color_gradientn(colours = c('grey','#64dbd0','#2EC4B6')) + NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))


ggsave(p1,filename = "HIM-integrate/figure/CS1 up AUC score in HIM_integrate.png",width = 6, height = 6) 
ggsave(p2,filename = "HIM-integrate/figure/CS2 up AUC score in HIM_integrate.pdf",width = 6, height = 6) 
ggsave(p3,filename = "HIM-integrate/figure/CS3 up AUC score in HIM_integrate.pdf",width = 6, height = 6) 





#----- cancer cells ------------------------------------------------------------
cancer <- HIM_integrate[,HIM_integrate@meta.data$cell_type == "Cancer cells"]

cancer <- NormalizeData(cancer)
cancer <- FindVariableFeatures(cancer, selection.method = "dispersion",
                               nfeatures = 2000)
cancer <- ScaleData(cancer, features = VariableFeatures(cancer))
cancer <- RunPCA(cancer, verbose = T)

ElbowPlot(object = cancer, ndims = 50)
cancer <- FindNeighbors(cancer, dims = 1:20)
cancer <- FindClusters(object = cancer,resolution = c(seq(.1,1.6,.2)))
library(clustree)
p1 <- clustree(cancer@meta.data, prefix = "RNA_snn_res.")
ggsave(p1, filename = "figure/cancer_clustree.png",width = 8, height = 6)
#-------resolution decide using the number of cluster
#-------dims decide using the number of PCA
cancer <- RunUMAP(cancer, dims = 1:20)
cancer <- RunTSNE(cancer, dims = 1:20)


p1 <- DimPlot(cancer, reduction = "umap", group.by = "orig.ident",label = T)
p2 <- DimPlot(cancer, reduction = "tsne", group.by = "RNA_snn_res.0.9", label = F)
p3 <- DimPlot(cancer, reduction = "umap", group.by = "group") +
  scale_color_manual(values=c('H'='#E71D36','M'='#FF9F1C','I'='#2EC4B6'))


ggsave(p1, filename = "figure/cancer_umap.png",width = 6, height = 4)
ggsave(p2, filename = "figure/cancer_tsne.png",width = 7, height = 6)
ggsave(p3, filename = "figure/cancer_umap_group.png",width = 6, height = 4)


saveRDS(cancer,file = "HIM_cancer/rds_file/cancer.rds")

cancer <- readRDS(file = "HIM_cancer/rds_file/cancer.rds")

DotPlot(cancer,features = 'CD274',group.by = 'cluster')


#----- remove T cell -----------------------------------------------------#
DotPlot(cancer,features = cluster_markers,group.by = "RNA_snn_res.0.7") + coord_flip()
cancer@meta.data[which(cancer@meta.data$RNA_snn_res.0.7 %in% c(7,8,10,13,36)),'cell_type'] = "T cells"
redefine_cancer <- cancer[,cancer@meta.data$cell_type == "T cells"]
#----- add Cancer cells to HIM_integrate 
HIM_integrate$cell_type <- as.character(HIM_integrate$cell_type)
HIM_integrate$cell_type[match(colnames(redefine_cancer),colnames(HIM_integrate))] <- as.character(redefine_cancer$cell_type)
HIM_integrate$cell_type <- as.factor(HIM_integrate$cell_type)

table(HIM_integrate$cell_type)
table(cancer$cell_type)
table(redefine_cancer$cell_type)

saveRDS(HIM_integrate,file = "HIM_integrate.rds")

HIM_integrate <- readRDS(file = "HIM_integrate.rds")


HIM_integrate$cell_type <- as.character(HIM_integrate$cell_type)
HIM_integrate$re_celltype <- as.character(HIM_integrate$re_celltype)

tmp <-  rownames(HIM_integrate[,HIM_integrate$cell_type == 'Cancer cells']@meta.data)
tmp1 <- rownames(HIM_integrate[,HIM_integrate$re_celltype == 'Cancer cells']@meta.data)
rename_cells <- setdiff(tmp1,tmp) # 顺序对取差集有影响

for (i in rename_cells){
  HIM_integrate$re_celltype[rownames(HIM_integrate@meta.data) == i] <- "T cells"
}

table(HIM_integrate$cell_type) 
table(HIM_integrate$re_celltype)

saveRDS(cancer,file = "cancer.rds")

cancer <- readRDS(file = "cancer.rds")




#----- cancer cell infercnv -------------------------------------------------------
# https://www.jianshu.com/p/72005719ed7c
# https://zhuanlan.zhihu.com/p/433234064
# https://zhuanlan.zhihu.com/p/376438151
# https://zhuanlan.zhihu.com/p/376438151
# https://www.jianshu.com/p/b12d165e5085

# https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV
# BiocManager::install("infercnv")
# install JAGS-4.3.0.exe in windows
# BiocManager::install("rjags",force = TRUE)

#---- inferCNV requires
# https://github.com/broadinstitute/inferCNV/wiki
# 1. a raw counts matrix of single-cell RNA-Seq expression
# 2. an annotations file which indicates which cells are tumor vs. normal.
# 3. a gene/chromosome positions file

# 基因定位文件相对不好准备，有三个途径获得：
# 1. 如果用自己测序的数据分析，可以找测序公司要，这个最准确;
# 2. 使用broad准备好的文件，问题是版本太老了，会有一些数据丢失；
# 3. 下载基因组对应的GTF文件，自己提取基因坐标信息，这个对初学者可能有点困难。
#     download gtf file from https://www.gencodegenes.org/human/
#     zcat  gencode.v42.annotation.gtf.gz |perl -alne '{next unless $F[2] eq "gene" ;/gene_name /"(.*?)/";/; print "$F[0]/t$F[3]/t$F[4]/t$1" }' > gencode.annotation.gtf

library(tidyverse)
library(infercnv)
library(rjags)
library(AnnoProbe) # if cannot install AnnoProbe then change to use python
library(Seurat)
library(Matrix)

HIM_integrate <- readRDS(file = "HIM_integrate.rds")
unique(HIM_integrate$orig.ident)
dir.create("infercnv")


#----- 使用循环计算肿瘤CNV ----------------------------------------------------#
for(i in unique(HIM_integrate$orig.ident)){
  setwd('G:/HIM-omic/HIM_single cell/result/infercnv')
  if (!dir.exists(i)) {
    dir.create(i)
  }
  setwd(i)
  
  tmp <- HIM_integrate[,HIM_integrate@meta.data$orig.ident == i]
  tmp_infercnv <- tmp[,tmp@meta.data$cell_type %in% c("T cells","Cancer cells")]
  
  # 1. prepare raw counts matrix and annotations file
  cellAnnota <- subset(tmp_infercnv@meta.data, select='cell_type')
  cellAnnota$cell_type <- as.character(cellAnnota$cell_type)
  
  exprMatrix <- GetAssayData(tmp_infercnv,assay = "RNA", slot='counts')
  exprMatrix[1:5,1:5]
  exprMatrix <- as.matrix(exprMatrix)
  
  # 2. using AnnoProbe change symbol ID
  geneInfor=annoGene(rownames(exprMatrix),"SYMBOL",'human')
  geneInfor_kmeans <- annoGene(rownames(exprMatrix),"SYMBOL",'human')
  colnames(geneInfor)
  geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]      
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  length(unique(geneInfor[,1]))
  head(geneInfor)
  
  # 检查矩阵的一致性
  exprMatrix = exprMatrix[match(geneInfor[,1], rownames(exprMatrix)),] 
  
  rownames(geneInfor) <- geneInfor$SYMBOL
  geneInfor <- geneInfor[,-1]
  
  identical(colnames(exprMatrix),rownames(cellAnnota))  
  identical(rownames(exprMatrix),rownames(geneInfor))
  # write.table(exprMatrix, 'inferCNV/exprMatrix.txt', col.names=NA, sep='\t')
  # write.table(cellAnnota, 'inferCNV/cellAnnota.txt', col.names=F, sep='\t')
  
  # run infercnv
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = exprMatrix,
                                      annotations_file = cellAnnota,
                                      delim="\t",
                                      gene_order_file = geneInfor,
                                      ref_group_names = c("T cells"))
  
  # https://github.com/broadinstitute/infercnv/issues/533
  # options(bitmapType="Xlib")
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir = "cnv1/", # dir is auto-created for storing outputs
                               cluster_by_groups = FALSE,  # cluster
                               denoise = F, 
                               HMM= FALSE, # HMM set as F for save time
                               num_threads= 8)
  # infercnv_obj_medianfiltered = infercnv::apply_median_filtering(infercnv_obj)
  # infercnv::plot_cnv(infercnv_obj_medianfiltered,
  #                    output_filename='infercnv.median_filtered',
  #                    x.range=c(0.9,1.1),
  #                    x.center=1,
  #                    title = "infercnv",
  #                    color_safe_pal = FALSE)
  # save(infercnv_obj, file = 'infer_cnv.RData')
}


#----- 通过聚类的方式判断肿瘤细胞 ---------------------------------------------#
# 定义函数，计算CNV score并绘制热图
calculate_CNV_score <- function(infercnv_obj, folder) {
  expr <- infercnv_obj@expr.data
  normal_loc <- infercnv_obj@reference_grouped_cell_indices$`T cells`
  cancer_loc <- infercnv_obj@observation_grouped_cell_indices$`Cancer cells`
  
  anno.df <- data.frame(
    CB = c(colnames(expr)[normal_loc], colnames(expr)[cancer_loc]),
    class = c(rep("normal", length(normal_loc)), rep("cancer", length(cancer_loc)))
  )
  
  # 聚类，7类，提取结果
  kmeans.result <- kmeans(t(expr), 7)
  kmeans_df <- data.frame(kmeans_class = kmeans.result$cluster)
  kmeans_df$CB <- rownames(kmeans_df)
  kmeans_df <- inner_join(kmeans_df, anno.df, by = "CB") # 合并
  kmeans_df_s <- arrange(kmeans_df, kmeans_class) # 排序
  rownames(kmeans_df_s) <- kmeans_df_s$CB
  
  # 将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:7排序
  kmeans_df_s$kmeans_class <- as.factor(kmeans_df_s$kmeans_class) 
  # 计算CNV score
  expr2 <- (expr - 1) ^ 2 # 放大CNV score差异
  CNV_score <- as.data.frame(colMeans(expr2))
  colnames(CNV_score) <- "CNV_score"
  CNV_score$CB <- rownames(CNV_score)
  CNV_score <- inner_join(CNV_score, kmeans_df_s, by = "CB")
  
  # 定义热图的注释，及配色
  color_v <- RColorBrewer::brewer.pal(8, "Dark2")[1:7] # 类别数
  CNV_score %>% ggplot(aes(kmeans_class, CNV_score)) +
    geom_violin(aes(fill = kmeans_class), color = "NA") +
    scale_fill_manual(values = color_v) +
    theme_bw()
  
  # 保存结果
  ggsave(file = paste0(folder,"","/CNV level.pdf"), width = 10, height = 6, units = "cm")
  write.csv(CNV_score,file = paste0(folder,"","/CNV level.csv"))
}


folders <- list.dirs("G:/HIM-omic/HIM_single cell/result/infercnv",recursive = FALSE, full.names = TRUE)
# recursive = FALSE表示只列出当前目录的直接子文件夹，不递归列出子文件夹。
# full.names = FALSE表示只返回文件夹的名称，不包括完整的路径。

# 循环处理每个文件夹
for (i in folders) {
  tmp <- paste0(i,"", "/cnv1/run.final.infercnv_obj")
  infercnv_obj <- readRDS(tmp)
  
  calculate_CNV_score(infercnv_obj, i)
}




#----- 针对每个肿瘤样本定义 ---------------------------------------------------#
setwd(dir = "G:/HIM-omic/HIM_single cell/result/infercnv")
folders <- list.dirs("G:/HIM-omic/HIM_single cell/result/infercnv",recursive = FALSE, full.names = FALSE)

CNV_list <- list()
for (i in folders) {
  file_path <- paste0(i,'','/CNV level.csv')
  CNV_list[[i]] <- read.csv(file_path,row.names = 1)
}

# table(CNV_list[['CID3946']]$kmeans_class,CNV_list[['CID3946']]$class) # 4个肿瘤细胞，不做比较
table(CNV_list[['CID3963']]$kmeans_class,CNV_list[['CID3963']]$class) # 全部是肿瘤细胞
table(CNV_list[['CID44041']]$kmeans_class,CNV_list[['CID44041']]$class) # 151个细胞，聚类完全分开
table(CNV_list[['CID4465']]$kmeans_class,CNV_list[['CID4465']]$class) # 3,4,7 肿瘤细胞
table(CNV_list[['CID4495']]$kmeans_class,CNV_list[['CID4495']]$class) ### 1,2,3,5 肿瘤细胞
table(CNV_list[['CID44971']]$kmeans_class,CNV_list[['CID44971']]$class) # 6,7 肿瘤细胞
table(CNV_list[['CID44991']]$kmeans_class,CNV_list[['CID44991']]$class) ### 1,2,3,5,6,7 肿瘤细胞 
table(CNV_list[['CID4513']]$kmeans_class,CNV_list[['CID4513']]$class) # 3,7 肿瘤细胞
table(CNV_list[['CID4515']]$kmeans_class,CNV_list[['CID4515']]$class) ### 1,2,4,5,6,7 肿瘤细胞
table(CNV_list[['CID4523']]$kmeans_class,CNV_list[['CID4523']]$class) ### 2,3,4,5,6,7 肿瘤细胞

table(CNV_list[['TN-B1-MH0131']]$kmeans_class,CNV_list[['TN-B1-MH0131']]$class) ### 拷贝数低，没有与其他样本聚集在一起，全部为肿瘤细胞
table(CNV_list[['TN-B1-MH0177']]$kmeans_class,CNV_list[['TN-B1-MH0177']]$class) # 1 肿瘤细胞
table(CNV_list[['TN-B1-MH4031']]$kmeans_class,CNV_list[['TN-B1-MH4031']]$class) ### 1,2,3,5,6,7 肿瘤细胞
table(CNV_list[['TN-B1-Tum0554']]$kmeans_class,CNV_list[['TN-B1-Tum0554']]$class) # 1,2,4,5 肿瘤细胞
table(CNV_list[['TN-MH0114-T2']]$kmeans_class,CNV_list[['TN-MH0114-T2']]$class) # 1,2,7 肿瘤细胞
table(CNV_list[['TN-MH0126']]$kmeans_class,CNV_list[['TN-MH0126']]$class) # 1,2,4,7 肿瘤细胞
table(CNV_list[['TN-MH0135']]$kmeans_class,CNV_list[['TN-MH0135']]$class) ### 拷贝数低，没有与其他样本聚集在一起，全部为肿瘤细胞
table(CNV_list[['TN-SH0106']]$kmeans_class,CNV_list[['TN-SH0106']]$class) ### 1,3,7 肿瘤细胞

table(CNV_list[['TNBC001']]$kmeans_class,CNV_list[['TNBC001']]$class) # 2,4 肿瘤细胞
table(CNV_list[['TNBC002']]$kmeans_class,CNV_list[['TNBC002']]$class) # 6,7 肿瘤细胞
table(CNV_list[['TNBC003']]$kmeans_class,CNV_list[['TNBC003']]$class) # 1,4,5,6 肿瘤细胞
table(CNV_list[['TNBC004']]$kmeans_class,CNV_list[['TNBC004']]$class) # 2,3,4,7 肿瘤细胞
table(CNV_list[['TNBC005']]$kmeans_class,CNV_list[['TNBC005']]$class) # 1,5,6 肿瘤细胞
table(CNV_list[['TNBC006']]$kmeans_class,CNV_list[['TNBC006']]$class) # 6,7 肿瘤细胞
table(CNV_list[['TNBC007']]$kmeans_class,CNV_list[['TNBC007']]$class) ### 1,2,3,4,5,7 肿瘤细胞
table(CNV_list[['LEN003']]$kmeans_class,CNV_list[['LEN003']]$class) ### 1,2,5,6 肿瘤细胞


#----- 标记肿瘤和非肿瘤细胞 ---------------------------------------------------#
condition_list <- list()
condition_list[['CID4465']] <- c(3,4,7)
condition_list[['CID4495']] <- c(1,2,3,5)
condition_list[['CID44971']] <- c(6,7)
condition_list[['CID44991']] <- c(1,2,3,5,6,7)
condition_list[['CID4513']] <- c(3,7)
condition_list[['CID4515']] <- c(1,2,4,5,6,7)
condition_list[['CID4523']] <- c(2,3,4,5,6,7)

condition_list[['TN-B1-MH0177']] <- c(1)
condition_list[['TN-B1-MH4031']] <- c(1,2,3,5,6,7)
condition_list[['TN-B1-Tum0554']] <- c(1,2,4,5)
condition_list[['TN-MH0114-T2']] <- c(1,2,7)
condition_list[['TN-MH0126']] <- c(1,2,4,7)
condition_list[['TN-SH0106']] <- c(1,3,7)

condition_list[['TNBC001']] <- c(2,4)
condition_list[['TNBC002']] <- c(6,7)
condition_list[['TNBC003']] <- c(1,4,5,6)
condition_list[['TNBC004']] <- c(2,3,4,7)
condition_list[['TNBC005']] <- c(1,5,6)
condition_list[['TNBC006']] <- c(6,7)
condition_list[['TNBC007']] <- c(1,2,3,4,5,7)
condition_list[['LEN003']] <- c(1,2,5,6)

# 提取normal cells
normal_list <- list()
for (i in names(condition_list)) {
  tmp <- CNV_list[[i]][CNV_list[[i]]$class == 'cancer',]
  normal_list[[i]] <- tmp[!(tmp$kmeans_class %in% condition_list[[i]]),]
}

dim(normal_list[['TNBC007']])

cancer$cell_type <- as.character(cancer$cell_type)
# 对肿瘤细胞重命名
for (i in names(normal_list)){
  tmp <- normal_list[[i]]$CB
  cancer$cell_type[rownames(cancer@meta.data) %in% tmp] <- "Normal cells"
}
table(cancer$orig.ident,cancer$cell_type) # 49768 4031

p1 <- DimPlot(cancer, reduction = "umap", group.by = "cell_type",label = F,cols = c('black','green')) + NoLegend()+ 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))
ggsave(p1,file = "figure/cancer_normal.pdf", width = 7, height = 6)
ggsave(p1,file = "figure/cancer_normal.png", width = 6, height = 6)

saveRDS(cancer,file = "cancer.rds")

cancer <- readRDS(file = "cancer.rds")




#----- cancer cell remove normal cells -------------------------------------------------
cancer_remor <- cancer[,cancer@meta.data$cell_type == "Cancer cells"]

batch1 <- c("TNBC001","TNBC002","TNBC003","TNBC004","TNBC005","TNBC006","TNBC007","LEN003")
batch2 <- c("CID3946","CID3963","CID4465","CID4495","CID4513","CID4515","CID4523","CID44041","CID44971","CID44991")
batch3 <- c("TN-MH0126","TN-MH0135","TN-SH0106","TN-MH0114-T2","TN-B1-Tum0554","TN-B1-MH4031","TN-B1-MH0177","TN-B1-MH0131")
batch <- ifelse(cancer_remor$orig.ident %in% batch1,"batch1",
                ifelse(cancer_remor$orig.ident %in% batch2,"batch2",
                       ifelse(cancer_remor$orig.ident %in% batch3,"batch3","reset")))
table(batch)
cancer_remor$batch <- batch

cancer_remor <- NormalizeData(cancer_remor)
cancer_remor <- FindVariableFeatures(cancer_remor, selection.method = "dispersion", nfeatures = 2000)
cancer_remor <- ScaleData(cancer_remor, features = VariableFeatures(cancer_remor))
cancer_remor <- RunPCA(cancer_remor, verbose = T)

library(harmony)
cancer_remor <- RunHarmony(cancer_remor, group.by.vars = "batch")
ElbowPlot(object = cancer_remor, ndims = 50)

cancer_remor <- FindNeighbors(cancer_remor, reduction = "harmony", dims = 1:20)
cancer_remor <- FindClusters(object = cancer_remor,resolution = c(seq(.1,1.6,.2)))
library(clustree)
p1 <- clustree(cancer_remor@meta.data, prefix = "RNA_snn_res.")
ggsave(p1,filename = "figure/cancer_remor_harmony.png",width = 8, height = 6)

cancer_remor <- RunUMAP(cancer_remor, reduction = "harmony", dims = 1:20)
cancer_remor <- RunTSNE(cancer_remor, reduction = "harmony", dims = 1:20)

p1 <- DimPlot(cancer_remor, reduction = "umap", group.by = "cluster",raster = F) +
  scale_color_manual(values=c('CS1'='#E71D36','CS2'='#FF9F1C','CS3'='#2EC4B6')) + NoLegend()+ 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))

p2 <- DimPlot(cancer_remor, reduction = "umap", group.by = "RNA_snn_res.0.5",raster = F,label = T)
# p3 <- DimPlot(cancer_remor, reduction = "umap", group.by = "RNA_snn_res.0.7",label = T)

ggsave(p1, filename = "figure/cancer_remor_umap.pdf",width = 6, height = 6)
ggsave(p1, filename = "figure/cancer_remor_umap.png",width = 6, height = 6)

ggsave(p2, filename = "figure/cancer_remor/umap0.5.pdf",width = 6, height = 6)

saveRDS(cancer_remor,file = "HIM_cancer/rds_file/cancer_remor.rds")

cancer_remor <- readRDS(file = "cancer_remor.rds")

prop.table(table(cancer_remor$cluster,cancer_remor$RNA_snn_res.0.5),margin = 2)

#----- cancer cell addmodule score ---------------------------------------------
load('G:/HIM-omic/Proteogenomic/in_house/result/fuclust_mem_3sub.rda') #膜蛋白 一致性聚类 的三分型 的fuclust # 分类结果
load('G:/HIM-omic/Proteogenomic/in_house/result/HIM_proteomic_data_v2.rda') # pt蛋白表达数据
gs_tcsa <- read.csv("G:/HIM-omic/Proteogenomic/in_house/result/HIM_subtype/TCSA_membrane_gene_list.CSV",check.names = F) #3567

CS3_deg <- read.csv("G:/HIM-omic/Proteogenomic/in_house/result/Specific_gene_of_CS3.csv",row.names = 1)
CS3_UP <- intersect(rownames(CS3_deg[CS3_deg$logFC > 0.2 & CS3_deg$P.Value < 0.05,]),gs_tcsa$`HGNC Symbol`)

CS3_UP <- intersect(rownames(cancer_remor),CS3_UP)

DotPlot(object = cancer_remor,features = CS3_UP,group.by = 'RNA_snn_res.0.7') +
  coord_flip()+ theme(axis.text.x = element_text(angle = 45,hjust = 1))



#----- T_cells -----------------------------------------------------------------
Tcells <- HIM_integrate[,HIM_integrate@meta.data$cell_type %in% c('T cells')]

noiseGenes <- c("MALAT1", "^HSPA", "^MT-", "^MT\\.", "^RP[SL]", "^LOC(0-9)", "^TR(A|B|G|D)V",
                "^MTRNR",'^AC(0-9)','^CTB-','^RP11-','^CTD-','^ENSG','^LINC(0-9)')
Tcells <- Tcells[!grepl(paste0(noiseGenes, collapse = "|"), rownames(Tcells)),] # 54402 50031
Tcells <- subset(Tcells, subset = nFeature_RNA > 200) # 54402 49837


Tcells <- NormalizeData(Tcells)
Tcells <- FindVariableFeatures(Tcells, selection.method = "dispersion",
                               nfeatures = 2000)
Tcells <- ScaleData(Tcells, features = VariableFeatures(Tcells))
Tcells <- RunPCA(Tcells, verbose = T)
library(harmony)
Tcells <- RunHarmony(Tcells, group.by.vars = "batch")
ElbowPlot(object = Tcells, ndims = 50)
Tcells <- FindNeighbors(Tcells, reduction = "harmony", dims = 1:20)
Tcells <- FindClusters(object = Tcells,resolution = c(seq(.1,1.6,.2)))
library(clustree)
p1 <- clustree(Tcells@meta.data, prefix = "RNA_snn_res.")
ggsave(p1, filename = "figure/Tcell_clustree.png",width = 8, height = 6)
#-------resolution decide using the number of cluster
#-------dims decide using the number of PCA
Tcells <- RunUMAP(Tcells, reduction = "harmony", dims = 1:20)
Tcells <- RunTSNE(Tcells, reduction = "harmony", dims = 1:20)


p1 <- DimPlot(Tcells, reduction = "umap", group.by = "RNA_snn_res.1.5",label = T) 
# +scale_color_manual(values=c('11'='#E71D36','NKT'='black')) + NoLegend()

p2 <- DimPlot(Tcells, reduction = "umap", group.by = "RNA_snn_res.1.5", label = T) + NoLegend()
p3 <- DimPlot(Tcells, reduction = "tsne", group.by = "cluster",label = T)

ggsave(p1, filename = "HIM_immune/figure/Tcell_umap.png",width = 5, height = 4)
ggsave(p2, filename = "HIM_immune/figure/Tcell_tsne.png",width = 8, height = 8)
ggsave(p3, filename = "HIM_immune/figure/Tcell_tsne_cluster.png",width = 8, height = 6)


# saveRDS(Tcells,file = "HIM_immune/rds_file/Tcells.rds")
# Tcells <- readRDS(file = "HIM_immune/rds_file/Tcells.rds")

saveRDS(Tcells,file = "HIM_immune/rds_file/Tcells_1219.rds")


Tcells <- readRDS(file = "HIM_immune/rds_file/Tcells_1219.rds")



#----- T_cell type -------------------------------------------------------------
# Error: useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE.
# remotes::install_version("matrixStats", version="1.1.0")

#----- singleR
library(SingleR)
library(celldex)
# BlueprintEncodeData (Obtain human bulk RNA-seq data from Blueprint and ENCODE)
# DatabaseImmuneCellExpressionData (Obtain human bulk RNA-seq data from DICE)
# HumanPrimaryCellAtlasData (Obtain the HPCA data)
# MonacoImmuneData (Obtain bulk RNA-seq data of sorted human immune cells)
# NovershternHematopoieticData (Obtain bulk microarray expression for sorted hematopoietic cells)

# ImmGenData (Obtain mouse bulk expression data from the Immunologic Genome Project)
# MouseRNAseqData (Obtain mouse bulk expression data of sorted cell populations (RNA-seq)
ref <- readRDS(file = "/public/home/scRNA/referencedata/celldex/NovershternHematopoieticData.rds")

tmp <- GetAssayData(Tcells, slot = "data")
cluster <- Tcells@meta.data$RNA_snn_res.1.5
singler <- SingleR(test = tmp, ref = ref, labels = ref$label.main, clusters = cluster,
                   fine.tune = T)
new.cluster.ids <- singler$labels
Tcells@active.ident <- Tcells$RNA_snn_res.1.5
names(new.cluster.ids) <- levels(Tcells)
Tcells <- RenameIdents(Tcells, new.cluster.ids)
Tcells@meta.data$singleR_cluster <- Tcells@active.ident # add singleR celltype to meta.data


p1 <- DimPlot(Tcells, reduction = "tsne", group.by = "singleR_cluster",label = T) #+
# scale_color_manual(values=c('Dendritic cells'='#E71D36','B cells'='black')) + NoLegend()
p2 <- DimPlot(Tcells, reduction = "umap", group.by = "singleR_cluster",label = T)

ggsave(p1, filename = "figure/Tcells_singleR.png",width = 8, height = 6)
ggsave(p2, filename = "figure/Tcells_singleR_umap.png",width = 8, height = 6)



#----- celltypist
library(Seurat)
library(SeuratDisk)

# Seu2scan <- DietSeurat(Tcells,counts = TRUE,data = TRUE,scale.data = FALSE,features = NULL,
#                        assays = NULL,dimreducs = NULL,graphs = NULL,misc = TRUE) # 目前还不能删除data

# 格式转换需要两步完成
SaveH5Seurat(Seu2scan,filename = "HIM_immune/data/Tcells.h5Seurat",overwrite = T)
Convert(source = "HIM_immune/data/Tcells.h5Seurat",dest = "h5ad",overwrite = T)

# 读取分析完成的数据
# 使用 Convert 和 LoadH5Seurat 读取数据会报错，将meta data 输出为CSV 然后再读取
meta_data <- read.csv(file = 'HIM_immune/data/Tcells.csv',row.names = 1)
meta_data[1:5,1:5]
Tcells@meta.data$majority_voting <- meta_data[,"majority_voting"]
head(Tcells)
table(Tcells$majority_voting)

p1 <- DimPlot(Tcells, reduction = "tsne", group.by = "majority_voting",label = T)+ NoLegend()
# + scale_color_manual(values=c('Tem/Trm cytotoxic T cells'='#E71D36')) 
ggsave(p1,filename = "HIM_immune/figure/Tcells_majority_voting.png",width = 6, height = 5) 


#----- 参考singleR，将Tcells 分为 CD4 and CD8
table(Tcells$singleR_cluster)

tmp <- rownames(Tcells@meta.data[Tcells$singleR_cluster == 'Granulocytes',])
Tcells$singleR_cluster[rownames(Tcells@meta.data) %in% tmp] <- 'CD8+ T cells'

table(Tcells$singleR_cluster)





#----- 手动注释
#  https://www.biocompare.com/Editorial-Articles/569888-A-Guide-to-T-Cell-Markers/
# 'PTPRC','CD3D','CD3E','CD3G','CD4','CD8A','CD8B','CD69','CD38',
# 'CCR7','IL7R','SELL','TCF7','LEF1','CCR5', 'FOXP1','PRDM1','CD127' ----------(naïve, Tscm, Tcm)
# 'GZMA','GZMK','GZMH','GZMB' --------------------------------------------------(Teff)
# 'HLA-DRB1','HLA-DQB1','HLA-DRA','HLA-DQA1','HLA-DRB5' ------------------------(activated T cells, Tem, Teff)
# 'PDCD1','CXCL13','CTLA4' -----------------------------------------------------(Tex)
# 'ISG15','IFIT1','CCR1', ------------------------------------------------------(T-ISG)
# 'FOXP3', 'TGFB1' -------------------------------------------------------------(Tregs)
# 'CXCR3','IFNG','IL2','TNF','STAT4','TBX21' -----------------------------------(Th1, Tc1)
# 'PTGDR2','CCR4','IL13','IL4','IL5','GATA3'  ----------------------------------(Th2, Tc2)
# 'IL9','IRF4' -----------------------------------------------------------------(Th9)
# 'CCR6','RORA','RORC','IL17A','IL17F','KLRB1','IL21','IL25','IL26' ------------(Th17, Tc17)
# 'CCR10','IL22', 'AHR','FOXO4' ------------------------------------------------(Th22)
# 'NCAM1' ----------------------------------------------------------------------(NKT)

DefaultAssay(Tcells) <- 'RNA'

p1 <- FeaturePlot(Tcells,features = "CD8A",reduction = "tsne")+ NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))
p2 <- FeaturePlot(Tcells,features = "CD4",reduction = "tsne")+ NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))
ggsave(p1, filename = "HIM_immune/figure/Tcell_CD8A.pdf",width = 6, height = 6)
ggsave(p2, filename = "HIM_immune/figure/Tcell_CD4.pdf",width = 6, height = 6)



p1 <- DotPlot(Tcells,features = c('PTPRC','CD3D','CD3E','CD3G','CD4','CD8A','CD8B','CD69','CD38',
                                  'CCR7','IL7R','SELL','TCF7','LEF1','CCR5', 'FOXP1','PRDM1','CD127',
                                  'GZMA','GZMK','GZMH','GZMB',
                                  'HLA-DRB1','HLA-DQB1','HLA-DRA','HLA-DQA1','HLA-DRB5',
                                  'PDCD1','CXCL13','CTLA4',
                                  'ISG15','IFIT1','CCR1',
                                  'FOXP3', 'TGFB1',
                                  'CXCR3','IFNG','IL2','TNF','STAT4','TBX21',
                                  'PTGDR2','CCR4','IL13','IL4','IL5','GATA3',
                                  'IL9','IRF4',
                                  'CCR6','RORA','RORC','IL17A','IL17F','KLRB1','IL21','IL25','IL26',
                                  'CCR10','IL22', 'AHR','FOXO4',
                                  'NCAM1',
                                  'SEPTIN1','SEPTIN6','SEPTIN7','SEPTIN9',
                                  'MKI67','TOP2A'),
              group.by = 'RNA_snn_res.1.5') + coord_flip()

ggsave(p1, filename = "HIM_immune/figure/Tcell_dotplot.png",width = 12, height = 12)



cluster.ids <- c('0'= 'CD4-Tn','1'= 'CD4-Tn','15'= 'CD4-Tn','19'= 'CD4-Tn',
                 '2'= 'CD8-Tn','24'= 'CD8-Tn','26'= 'CD8-Tn',
                 '3'= 'CD8-Teff','5'= 'CD8-Teff','7'= 'CD8-Teff','10'= 'CD8-Teff','14'= 'CD8-Teff',
                 '17'= 'CD8-Teff','22'= 'CD8-Teff','23'= 'CD8-Teff','25'= 'CD8-Teff','32'= 'CD8-Teff','33'= 'CD8-Teff',
                 '9'= 'CD8-Tex',
                 '11'= 'CD4-Th','20'= 'CD4-Th',
                 '12'= 'CD8-ISG',
                 '13'= 'CD4-ISG',
                 '4'= 'CD4-Treg','8'= 'CD4-Treg',
                 '6'= 'T-SEPTIN','21'= 'T-SEPTIN','27'= 'T-SEPTIN','28'= 'T-SEPTIN',
                 '29'= 'T-SEPTIN','30'= 'T-SEPTIN',
                 '16'= 'CD8-Tpro','31'= 'CD8-Tpro',
                 '18'= 'CD4-Tpro')
Tcells@active.ident <- Tcells$RNA_snn_res.1.5
Tcells <- RenameIdents(Tcells, cluster.ids)                        
Tcells$cell_type <- Tcells@active.ident

p1 <- DotPlot(Tcells,features = c('CD4','CD8A','CD8B',
                                  'CCR7','IL7R','SELL','TCF7','LEF1',
                                  'GZMA','GZMK','GZMH','GZMB',
                                  'PDCD1','CXCL13',
                                  'ISG15','IFIT1',
                                  'FOXP3',
                                  'SEPTIN1','SEPTIN6','SEPTIN7','SEPTIN9',
                                  'MKI67','TOP2A'),
              group.by = 'cell_type') + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))


p2 <- DimPlot(Tcells, reduction = "umap", group.by = "cell_type",label = T) + NoLegend()
p3 <- DimPlot(Tcells, reduction = "tsne", group.by = "cell_type",label = T) + NoLegend()
p4 <- DimPlot(Tcells, reduction = "tsne", group.by = "cell_type",label = T)
p5 <- DimPlot(Tcells, reduction = "tsne", group.by = "cell_type",label = F) + NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),
        axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))

ggsave(p1,filename = "HIM_immune/figure/Tcells_type_marker.pdf",width = 6, height = 6)
ggsave(p2,filename = "HIM_immune/figure/Tcells_type_umap.png",width = 6, height = 6)
ggsave(p3,filename = "HIM_immune/figure/Tcells_type_tsne.png",width = 6, height = 6)
ggsave(p4,filename = "HIM_immune/figure/Tcells_type_tsne_legend.pdf",width = 6, height = 6)
ggsave(p5,filename = "HIM_immune/figure/Tcells_type_tsne.pdf",width = 6, height = 6)



saveRDS(Tcells,file = "HIM_immune/rds_file/Tcells_1219.rds")


Tcells <- readRDS(file = "HIM_immune/rds_file/Tcells_1219.rds")








#----- T_cells cell ratio ------------------------------------------------------
Cellratio <- prop.table(table(Tcells$cell_type, Tcells$cluster),margin = 1)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
Cellratio$Var2 <- factor(Cellratio$Var2, levels = c('CS1','CS2','CS3'))

# 柱状图
library(ggsci)
p1 <- ggplot(Cellratio,aes(x = Var1, y= Freq,fill = Var2)) +
  geom_bar(stat = "identity", position = "stack")+ # dodge
  theme_classic() +
  labs(x='Sample',y = 'Ratio') +
  RotatedAxis() 
plotc <- p1 + scale_fill_manual(values=c('CS1'='#ee6274','CS2'='#ffbf69','CS3'='#64dbd0'))

ggsave(plotc, filename = "HIM_immune/figure/T cell ratio.pdf",width = 6, height = 4)

# 饼图
library(tidyverse)
library(PieGlyph)

data <- data.frame(table(Tcells$cell_type))
data$CS1 <- Cellratio[Cellratio$Var2 ==  "CS1",]$Freq
data$CS2 <- Cellratio[Cellratio$Var2 ==  "CS2",]$Freq
data$CS3 <- Cellratio[Cellratio$Var2 ==  "CS3",]$Freq

p1 <- ggplot(data = data, aes(x = Var1, y = Freq))+
  geom_pie_glyph(aes(radius = Freq), 
                 slices = c('CS1', 'CS2', 'CS3'), 
                 colour = 'black')+
  scale_fill_manual(values=c('CS1'='#f5a7b1','CS2'='#ffe0b6','CS3'='#a2e9e3'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) 
ggsave(p1,filename = "HIM_immune/figure/T cell ratio pie.pdf",width = 8, height = 6)





#----- T cells function annotation -------------------------------------------
Tcells@active.ident <- Tcells$cell_type

T_DEGs <- FindAllMarkers(Tcells, only.pos = T, test.use = "wilcox")

write.csv(T_DEGs,file = 'HIM_immune/data/T_RNA15_DEGs.csv')

T_DEGs <- read.csv(file = 'HIM_immune/data/T_RNA15_DEGs.csv',row.names = 1)

top10 <- T_DEGs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p1 <- DoHeatmap(Tcells, # input data is scale.data
                features = unique(top10$gene),
                group.by = "cell_type",
                slot = 'scale.data',
                group.bar = T, size = 4)
ggsave(p1,filename = "HIM_immune/figure/Tcell_DEG15_top10.png",width = 12,height = 16)


#----- T_celltype_DEGs
Tcells@active.ident <- Tcells$cell_type
T_celltype_DEGs <- FindAllMarkers(Tcells, only.pos = T, test.use = "wilcox")
write.csv(T_celltype_DEGs,file = 'HIM_immune/data/T_celltype_DEGs.csv')

T_celltype_DEGs <- read.csv(file = 'HIM_immune/data/T_celltype_DEGs.csv',row.names = 1)

top10 <- T_celltype_DEGs %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

#----- enrichment analysis 
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(enrichplot)
library(ggstar)
library(GOplot)
library(dplyr)
GO <- read.gmt("/public/home/scRNA/referencedata/GSEA/c5 Ontology gene sets/c5.go.v7.5.1.entrez.gmt") #GO
HM <- read.gmt("/public/home/scRNA/referencedata/GSEA/hallmark gene sets/h.all.v7.5.1.entrez.gmt")
KEGG <- read.gmt("/public/home/scRNA/referencedata/GSEA/c2 curated gene sets/c2.cp.kegg.v7.5.1.entrez.gmt")
Reactome <- read.gmt("/public/home/scRNA/referencedata/GSEA/c2 curated gene sets/c2.cp.reactome.v7.5.1.entrez.gmt")


# gene ID转换
ids = bitr(T_celltype_DEGs$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
T_celltype_DEGs = merge(T_celltype_DEGs,ids,by.x='gene',by.y='SYMBOL')

# T_celltype_DEGs$cluster <- factor(T_celltype_DEGs$cluster,levels = c(0:29))
# make group
gcSample = split(T_celltype_DEGs$ENTREZID, T_celltype_DEGs$cluster)


# GO
go <- compareCluster(gcSample,fun = "enricher", TERM2GENE = GO, pvalueCutoff = 0.05)
p1 <- dotplot(go, showCategory = 5,shape = F) + theme(axis.text.y = element_text(size = 7),
                                                      axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(p1,filename = "HIM_immune/figure/T_celltype_DEGs_GO.pdf",height = 10,width = 8)

saveRDS(go,file = 'HIM_immune/rds_file/T_celltype_DEGs_GO.rds')


# HM
hm <- compareCluster(gcSample,fun = "enricher", TERM2GENE = HM, pvalueCutoff = 0.05)
p1 <- dotplot(hm, showCategory = 5,shape = F) + theme(axis.text.y = element_text(size = 7),
                                                      axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(p1,filename = "HIM_immune/figure/T_celltype_DEGs_HM.pdf",height = 10,width = 8)

saveRDS(go,file = 'HIM_immune/rds_file/T_celltype_DEGs_HM.rds')


# KEGG
kegg <- compareCluster(gcSample,fun = "enricher", TERM2GENE = KEGG, pvalueCutoff = 0.05)
p1 <- dotplot(kegg, showCategory = 10,shape = F) + theme(axis.text.y = element_text(size = 7),
                                                         axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(p1,filename = "HIM_immune/figure/T_celltype_DEGs_KEGG.pdf",height = 10,width = 8)

saveRDS(kegg,file = 'HIM_immune/rds_file/T_celltype_DEGs_KEGG.rds')


# reactome
reactome <- compareCluster(gcSample,fun = "enricher", TERM2GENE = Reactome, pvalueCutoff = 0.05)
p1 <- dotplot(reactome, showCategory = 10,shape = F) + theme(axis.text.y = element_text(size = 7),
                                                             axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(p1,filename = "HIM_immune/figure/T_celltype_DEGs_reactome.pdf",height = 20,width = 8)

saveRDS(reactome,file = 'HIM_immune/rds_file/T_celltype_DEGs_reactome.rds')

reactome <- readRDS(file = 'HIM_immune/rds_file/T_celltype_DEGs_reactome.rds')




#----- T-septiin senescence signature ------------------------------------------
# 参考文献 A new gene set identifies senescent cells and predicts senescence-associated pathways across tissues
senescence_geneset <- read.table(file = 'senescene_seneset.txt')
senescence_geneset <- as.character(senescence_geneset$V1)
senescence_geneset <- list(senescence_geneset)

Tcells <- AddModuleScore(object = Tcells, features = senescence_geneset,name = "senescence_geneset")


VlnPlot(Tcells,features = 'senescence_geneset1',group.by = 'cell_type',pt.size = 0)

DotPlot(Tcells,features = c('TP53','CDKN1A','CDKN2A','ATM','H2AX','IL6','CXCL8','IL1A','IL1B','TNF','VEGFA',
                            'IFNG','B3GAT1','KLRG1','HAVCR2','TIGIT','GLB1','TERT'),group.by = 'cell_type')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 12))



p1 <- VlnPlot(Tcells,features = c('H2AX'),group.by = 'cell_type') + NoLegend()
ggsave(p1,filename = 'HIM_immune/figure/Tcell_H2AX.pdf',height = 3,width = 4)


#----- T cells immune checkpoint --------------------------------------------
check_markers <- c('ICOSLG','CD40','PDCD1LG2','CD80','CD274','HAVCR2','TNFSF18','CD276','PVR','LGALS9','PVRL2',
                   'CD86','TNFRSF14','C10orf54','TNFSF9','KLRC2','CD226','CD40LG','ADORA2A','TNFSF4','LAG3','PDCD1',
                   'BTLA','TNFRSF18','TNFRSF4','FAS','ICOS','TNFRSF9','TIGIT','CTLA4','CD70','CD27','CXCL13','CDKN1A','TGFB1')
# check_markers <- c('CD40LG','BTLA','GITR','PDCD1','HAVCR2','LAG3','CD28','CTLA4','TIGIT',
#                    'TNFRSF4','TNFRSF9','ICOS','CD27','SELPLG') # TNFRSF4(OX40),CD40LG(CD40L),TNFRSF9(4-1BB),HAVCR2(TIM3)
p1 <- DotPlot(Tcells,features = check_markers,group.by = "cell_type") + coord_flip() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = NULL))
ggsave(p1,file = "HIM_immune/figure/T ICIs.pdf")



#----- Mac cells -------------------------------------------------------------
Mac <- HIM_integrate[,HIM_integrate@meta.data$cell_type %in% c('Macrophages')]

noiseGenes <- c("MALAT1", "^HSPA", "^MT-", "^MT\\.", "^RP[SL]", "^LOC(0-9)", "^TR(A|B|G|D)V",
                "^MTRNR",'^AC(0-9)','^CTB-','^RP11-','^CTD-','^ENSG','^LINC(0-9)')
Mac <- Mac[!grepl(paste0(noiseGenes, collapse = "|"), rownames(Mac)),] # 64646 24021
Mac <- subset(Mac, subset = nFeature_RNA > 200) # 64646 23965


Mac <- NormalizeData(Mac)
Mac <- FindVariableFeatures(Mac, selection.method = "dispersion",
                            nfeatures = 2000)
Mac <- ScaleData(Mac, features = VariableFeatures(Mac))
Mac <- RunPCA(Mac, verbose = T)
library(harmony)
Mac <- RunHarmony(Mac, group.by.vars = "batch")
ElbowPlot(object = Mac, ndims = 50)
Mac <- FindNeighbors(Mac, reduction = "harmony", dims = 1:30)
Mac <- FindClusters(object = Mac,resolution = c(seq(.1,1.6,.2)))
library(clustree)
p1 <- clustree(Mac@meta.data, prefix = "RNA_snn_res.")
ggsave(p1, filename = "HIM_immune/figure/Mac_clustree.png",width = 8, height = 6)
#------ resolution decide using the number of cluster
#------ dims decide using the number of PCA
Mac <- RunUMAP(Mac, reduction = "harmony", dims = 1:30)
Mac <- RunTSNE(Mac, reduction = "harmony", dims = 1:30)


p1 <- DimPlot(Mac, reduction = "tsne", group.by = "RNA_snn_res.1.5",label = T) 
p2 <- DimPlot(Mac, reduction = "umap", group.by = "RNA_snn_res.1.5",label = T)

ggsave(p1, filename = "HIM_immune/figure/Mac_tsne.png",width = 5, height = 4)
ggsave(p2, filename = "HIM_immune/figure/Mac_umap.png",width = 5, height = 4)





saveRDS(Mac,file = "HIM_immune/rds_file/Mac_1226.rds")

Mac <- readRDS(file = "HIM_immune/rds_file/Mac_1226.rds")


#----- Mac cell type -----------------------------------------------------------
Mac@active.ident <- Mac$RNA_snn_res.1.5
Mac_DEGs <- FindAllMarkers(Mac,
                           only.pos = T,
                           test.use = "wilcox")

write.csv(Mac_DEGs,file = 'HIM_immune/data/Mac_DEGs.csv')

Mac_DEGs <- read.csv(file = 'HIM_immune/data/Mac_DEGs.csv',row.names = 1)
top10 <- Mac_DEGs %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)


#----- singleR
# Error: useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE.
# remotes::install_version("matrixStats", version="1.1.0")

library(SingleR)
library(celldex)
# BlueprintEncodeData (Obtain human bulk RNA-seq data from Blueprint and ENCODE)
# DatabaseImmuneCellExpressionData (Obtain human bulk RNA-seq data from DICE)
# HumanPrimaryCellAtlasData (Obtain the HPCA data)
# MonacoImmuneData (Obtain bulk RNA-seq data of sorted human immune cells)
# NovershternHematopoieticData (Obtain bulk microarray expression for sorted hematopoietic cells)

# ImmGenData (Obtain mouse bulk expression data from the Immunologic Genome Project)
# MouseRNAseqData (Obtain mouse bulk expression data of sorted cell populations (RNA-seq)
ref <- readRDS(file = "/public/home/scRNA/referencedata/celldex/NovershternHematopoieticData.rds")

tmp <- GetAssayData(Mac, slot = "data")
cluster <- Mac@meta.data$RNA_snn_res.1.5
singler <- SingleR(test = tmp, ref = ref, labels = ref$label.main, clusters = cluster,
                   fine.tune = T)
new.cluster.ids <- singler$labels
Mac@active.ident <- Mac$RNA_snn_res.1.5
names(new.cluster.ids) <- levels(Mac)
Mac <- RenameIdents(Mac, new.cluster.ids)
Mac@meta.data$singleR_cluster <- Mac@active.ident # add singleR celltype to meta.data


p1 <- DimPlot(Mac, reduction = "tsne", group.by = "singleR_cluster",label = T) #+
# scale_color_manual(values=c('Dendritic cells'='#E71D36','B cells'='black')) + NoLegend()
p2 <- DimPlot(Mac, reduction = "umap", group.by = "singleR_cluster",label = T)

ggsave(p1, filename = "HIM_immune/figure/Mac_singleR.png",width = 8, height = 6)
ggsave(p2, filename = "HIM_immune/figure/Mac_singleR_umap.png",width = 8, height = 6)



#----- celltypist
library(Seurat)
library(SeuratDisk)

Seu2scan <- DietSeurat(Mac,counts = TRUE,data = TRUE,scale.data = FALSE,features = NULL,
                       assays = NULL,dimreducs = TRUE,graphs = TRUE,misc = TRUE) # 目前还不能删除data

# 格式转换需要两步完成
SaveH5Seurat(Seu2scan,filename = "HIM_immune/data/Mac.h5Seurat",overwrite = T)
Convert(source = "HIM_immune/data/Mac.h5Seurat",dest = "h5ad",overwrite = T)

# 读取分析完成的数据
# 使用 Convert 和 LoadH5Seurat 读取数据会报错，将meta data 输出为CSV 然后再读取
meta_data <- read.csv(file = 'HIM_immune/data/Mac.csv',row.names = 1)
meta_data[1:5,1:5]
Mac@meta.data$majority_voting <- meta_data[,"majority_voting"]
head(Mac)
table(Mac$majority_voting)

p1 <- DimPlot(Mac, reduction = "tsne", group.by = "majority_voting",label = T)+ NoLegend()
# + scale_color_manual(values=c('Tem/Trm cytotoxic T cells'='#E71D36'))
p2 <- DimPlot(Mac, reduction = "umap", group.by = "majority_voting",label = T)+ NoLegend()

ggsave(p1,filename = "HIM_immune/figure/Mac_majority_voting_tsne.png",width = 6, height = 5) 
ggsave(p2,filename = "HIM_immune/figure/Mac_majority_voting_umap.png",width = 6, height = 5) 




#----- 手动注释
# A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells
# cDC1_CLEC9A <- c("CLEC9A","FLT3", "IDO1")
# cDC2_CD1C <- c("CD1C", "FCER1A", "HLA-DQA1")
# cDC3_LAMP3 <- c("LAMP3", "CCR7", "FSCN1")
# Mono_CD14 <- c("FCN1", "S100A9", "S100A8")
# Mono_CD16 <- c("FCGR3A","LST1",'LILRB2') #FCGR3A:CD16
# Macro_INHBA <- c("INHBA","IL1RN","CCL4")
# Macro_NLRP3 <- c("NLRP3", "EREG", "IL1B")
# Macro_LYVE1 <- c("LYVE1","PLTP","SEPP1")
# Macro_C1QC <- c("C1QC", "C1QA","APOE")
# Neutrophil <- c('S100A9','S100A8','ENPP3','FUT4','CSF3R')
# ISG15+ TAMs: 'ISG15','IFIT1','CCR1',

cluster_markers <- c('CD14',"FCN1",'THBS1',"S100A9", "S100A8",
                     "FCGR3A","LST1",'LILRB2',
                     'CD163','CD68',"C1QC", "C1QA","APOE",'SPP1','TMEM2',
                     "CCL4","IL1B",'APOC1',
                     'CSF3R','ENPP3','FUT4',
                     'ISG15','IFIT1','CCR1',
                     "CLEC9A","FLT3","IDO1",
                     "CD1C", "FCER1A", "HLA-DQA1",
                     "LAMP3", "CCR7", "FSCN1",
                     'CD33','CCR3',
                     'MKI67','TOP2A')
DotPlot(Mac,features = cluster_markers,group.by = 'cell_type') + coord_flip()

FeaturePlot(Mac,features = c('CD14',"FCGR3A",'SPP1','TMEM2',"C1QC",'ITGAM'))




# tumor-infiltrating monocytes expressed higher levels of inflammatory cytokines and chemokines (IL1B, CCL4, CXCL2,and CXCR4),
cluster.ids <- c("0" = "TAM_C1QC", "1" = "TAM_C1QC", "3" = "TAM_C1QC", '5' = 'TAM_C1QC','6' = 'TAM_C1QC',
                 "9" = "TAM_C1QC","10" = "TAM_C1QC","11" = "TAM_C1QC","13" = "TAM_C1QC","14" = "TAM_C1QC",
                 "16" = "TAM_C1QC","17" = "TAM_C1QC","26" = "TAM_C1QC","36" = "TAM_C1QC",
                 
                 "8" = "TAM_C1QC","21" = "TAM_C1QC","37" = "TAM_C1QC",
                 
                 "29" = "Neutrophil","35" = "Neutrophil", # 'S100A9','S100A8', CSF3R
                 
                 "4" = "TAM_SPP1","7" = "TAM_SPP1","15" = "TAM_SPP1","19" = "TAM_SPP1","22" = "TAM_SPP1",
                 "24" = "TAM_SPP1","27" = "TAM_SPP1","34" = "TAM_SPP1", # spp1
                 
                 "2" = "Mono_CD14","25" = "Mono_CD14", "31" = "Mono_CD14",
                 "33" = "TAM_ISG15",
                 "32" = "cDC1",
                 "18" = "cDC2","23" = "cDC2","38" = "cDC2",
                 "28" = "cDC3",
                 '12' = 'TAM-pro','20' = 'TAM-pro','30' = 'TAM-pro')
Mac@active.ident <- Mac$RNA_snn_res.1.5
Mac <- RenameIdents(Mac, cluster.ids)                        
Mac$cell_type <- Mac@active.ident


cluster_markers <- c('CD14',"FCGR3A",'CD163','CD68',
                     "FCN1",'THBS1',
                     "C1QC", "C1QA","APOE",
                     'SPP1',
                     'ISG15','IFIT1','CCR1',
                     "CLEC9A","FLT3","IDO1",
                     "CD1C", "FCER1A", "HLA-DQA1",
                     "LAMP3", "CCR7", "FSCN1",
                     "S100A9", "S100A8",'FCGR3B','CSF3R',
                     'MKI67','TOP2A')
DotPlot(Mac,features = cluster_markers,group.by = 'cell_type') + coord_flip()

p1 <- DimPlot(Mac, reduction = "tsne", group.by = "cell_type",label = T)  + NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))

ggsave(p1, filename = "HIM_immune/figure/Mac_tsne.pdf",width = 7.5, height = 6)
ggsave(p1, filename = "HIM_immune/figure/Mac_tsne.png",width = 6.5, height = 6)


library(pheatmap)
for (i in 1:length(cluster_markers)) {
  ids = Mac[cluster_markers[i],]
  ids = ids@assays$RNA@data %>% as.numeric()
  assign(cluster_markers[i],tapply(ids, Mac@meta.data$cell_type,mean))
}

heatmap_matrix <- NULL
for (i in 1:length(cluster_markers)) {
  heatmap_matrix = rbind(heatmap_matrix,get(cluster_markers[i]))
}

row.names(heatmap_matrix) = cluster_markers
ids = c('Mono_CD14','TAM_C1QC','TAM_SPP1','TAM_ISG15','cDC1','cDC2','cDC3','Neutrophil','TAM-pro')
heatmap_matrix = heatmap_matrix[,ids]

col <- c(colorRampPalette(c('#3ab2e8','white'))(36),
         colorRampPalette(c('white','#da3446'))(36))

p1 <- pheatmap(heatmap_matrix,color = col, scale = 'row',cluster_cols = F,cluster_rows = F) #breaks = unique(c(seq(-2.5,2.5, length=100)))
ggsave(p1, filename = '/public/home/scRNA/HIM-omic/HIM_single_cell/result/HIM_immune/figure/Mac_markers.pdf',height = 6,width = 4)




saveRDS(Mac,file = "HIM_immune/rds_file/Mac_1226.rds")

Mac <- readRDS(file = "/public/home/scRNA/HIM-omic/HIM_single_cell/result/HIM_immune/rds_file/Mac_1226.rds")





#----- Mac cells cell ratio -----------------------------------------------
Cellratio <- prop.table(table(Mac$cell_type, Mac$cluster),margin = 1)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
Cellratio$Var2 <- factor(Cellratio$Var2, levels = c('CS1','CS2','CS3'))

library(ggsci)
p1 <- ggplot(Cellratio,aes(x = Var1, y= Freq,fill = Var2)) +
  geom_bar(stat = "identity", position = "stack")+ # dodge
  theme_classic() +
  labs(x='Sample',y = 'Ratio') +
  RotatedAxis() 
plotc <- p1 + scale_fill_manual(values=c('CS1'='#ee6274','CS2'='#ffbf69','CS3'='#64dbd0'))

ggsave(plotc, filename = "HIM_immune/figure/HIM_immune cell ratio.pdf",width = 6, height = 4)


library(tidyverse)
library(PieGlyph)
data <- data.frame(table(Mac$cell_type))
data$CS1 <- Cellratio[Cellratio$Var2 ==  "CS1",]$Freq
data$CS2 <- Cellratio[Cellratio$Var2 ==  "CS2",]$Freq
data$CS3 <- Cellratio[Cellratio$Var2 ==  "CS3",]$Freq

p1 <- ggplot(data = data, aes(x = Var1, y = Freq))+
  geom_pie_glyph(aes(radius = Freq), 
                 slices = c('CS1', 'CS2', 'CS3'), 
                 colour = 'black')+
  scale_fill_manual(values=c('CS1'='#f5a7b1','CS2'='#ffe0b6','CS3'='#a2e9e3'))+
  theme_classic()
ggsave(p1,filename = "HIM_immune/figure/Mac cell ratio pie.pdf",width = 8, height = 6)





#----- B cells -----------------------------------------------------------------
Bcells <- HIM_integrate[,HIM_integrate@meta.data$cell_type == "B cells"]

noiseGenes <- c("MALAT1", "^HSPA", "^MT-", "^MT\\.", "^RP[SL]", "^LOC(0-9)", "^TR(A|B|G|D)V",
                "^MTRNR",'^AC(0-9)','^CTB-','^RP11-','^CTD-','^ENSG','^LINC(0-9)')
Bcells <- Bcells[!grepl(paste0(noiseGenes, collapse = "|"), rownames(Bcells)),] # 64646 14050
Bcells <- subset(Bcells, subset = nFeature_RNA > 200) # 64646 13580

# Bcells <- Bcells[,Bcells@meta.data$RNA_snn_res.0.1 != "4"]

Bcells <- NormalizeData(Bcells)
Bcells <- FindVariableFeatures(Bcells, selection.method = "dispersion",
                               nfeatures = 2000)
Bcells <- ScaleData(Bcells, features = VariableFeatures(Bcells))
Bcells <- RunPCA(Bcells, verbose = T)
library(harmony)
Bcells <- RunHarmony(Bcells, group.by.vars = "orig.ident")
ElbowPlot(object = Bcells, ndims = 50)
Bcells <- FindNeighbors(Bcells, reduction = "harmony", dims = 1:10)
Bcells <- FindClusters(object = Bcells,resolution = c(seq(.1,1.6,.2)))
library(clustree)
p1 <- clustree(Bcells@meta.data, prefix = "RNA_snn_res.")
ggsave(p1, filename = "HIM_immune/figure/Bcell_clustree.png",width = 8, height = 6)
#-------resolution decide using the number of cluster
#-------dims decide using the number of PCA
Bcells <- RunUMAP(Bcells, reduction = "harmony", dims = 1:10)
Bcells <- RunTSNE(Bcells, reduction = "harmony", dims = 1:10)

p1 <- DimPlot(Bcells, reduction = "tsne", group.by = "RNA_snn_res.0.1",label = T) 
p2 <- DimPlot(Bcells, reduction = "umap", group.by = "RNA_snn_res.0.1",label = T)

ggsave(p1, filename = "HIM_immune/figure/Bcells_tsne.png",width = 5, height = 4)
ggsave(p2, filename = "HIM_immune/figure/Bcells_umap.png",width = 5, height = 4)



saveRDS(Bcells,file = "HIM_immune/rds_file/Bcells.rds")

Bcells <- readRDS(file = "HIM_immune/rds_file/Bcells.rds")



#----- B cells cell type -----------------------------------------------------------
Bcells@active.ident <- Bcells$cell_type
Bcells_DEGs <- FindAllMarkers(Bcells,
                              only.pos = T,
                              test.use = "wilcox")

write.csv(Bcells_DEGs,file = 'HIM_immune/data/Bcells_DEGs.csv')

Bcells_DEGs <- read.csv(file = 'HIM_immune/data/Bcells_DEGs.csv',row.names = 1)
top10 <- Bcells_DEGs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


#----- singleR
# Error: useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE.
# remotes::install_version("matrixStats", version="1.1.0")

library(SingleR)
library(celldex)
# BlueprintEncodeData (Obtain human bulk RNA-seq data from Blueprint and ENCODE)
# DatabaseImmuneCellExpressionData (Obtain human bulk RNA-seq data from DICE)
# HumanPrimaryCellAtlasData (Obtain the HPCA data)
# MonacoImmuneData (Obtain bulk RNA-seq data of sorted human immune cells)
# NovershternHematopoieticData (Obtain bulk microarray expression for sorted hematopoietic cells)

# ImmGenData (Obtain mouse bulk expression data from the Immunologic Genome Project)
# MouseRNAseqData (Obtain mouse bulk expression data of sorted cell populations (RNA-seq)
ref <- readRDS(file = "/public/home/scRNA/referencedata/celldex/NovershternHematopoieticData.rds")

tmp <- GetAssayData(Bcells, slot = "data")
cluster <- Bcells@meta.data$RNA_snn_res.0.1
singler <- SingleR(test = tmp, ref = ref, labels = ref$label.main, clusters = cluster,
                   fine.tune = T)
new.cluster.ids <- singler$labels
Bcells@active.ident <- Bcells$RNA_snn_res.0.1
names(new.cluster.ids) <- levels(Bcells)
Bcells <- RenameIdents(Bcells, new.cluster.ids)
Bcells@meta.data$singleR_cluster <- Bcells@active.ident # add singleR celltype to meta.data


p1 <- DimPlot(Bcells, reduction = "tsne", group.by = "singleR_cluster",label = T) #+
# scale_color_manual(values=c('Dendritic cells'='#E71D36','B cells'='black')) + NoLegend()
# p2 <- DimPlot(Bcells, reduction = "umap", group.by = "singleR_cluster",label = T)

ggsave(p1, filename = "HIM_immune/figure/Bcells_singleR.png",width = 8, height = 6)
# ggsave(p2, filename = "HIM_immune/figure/Bcells_singleR_umap.png",width = 8, height = 6)



#----- celltypist
library(Seurat)
library(SeuratDisk)

Seu2scan <- DietSeurat(Bcells,counts = TRUE,data = TRUE,scale.data = FALSE,features = NULL,
                       assays = NULL,dimreducs = TRUE,graphs = TRUE,misc = TRUE) # 目前还不能删除data

# 格式转换需要两步完成
SaveH5Seurat(Seu2scan,filename = "HIM_immune/data/Bcells.h5Seurat",overwrite = T)
Convert(source = "HIM_immune/data/Bcells.h5Seurat",dest = "h5ad",overwrite = T)

# 读取分析完成的数据
# 使用 Convert 和 LoadH5Seurat 读取数据会报错，将meta data 输出为CSV 然后再读取
meta_data <- read.csv(file = 'HIM_immune/data/Bcells.csv',row.names = 1)
meta_data[1:5,1:5]
Bcells@meta.data$majority_voting <- meta_data[,"majority_voting"]
head(Bcells)
table(Bcells$majority_voting)

p1 <- DimPlot(Bcells, reduction = "tsne", group.by = "majority_voting",label = T)+ NoLegend()
# + scale_color_manual(values=c('Tem/Trm cytotoxic T cells'='#E71D36'))
# p2 <- DimPlot(Bcells, reduction = "umap", group.by = "majority_voting",label = T)+ NoLegend()

ggsave(p1,filename = "HIM_immune/figure/Bcells_majority_voting_tsne.png",width = 6, height = 5) 
# ggsave(p2,filename = "HIM_immune/figure/Bcells_majority_voting_umap.png",width = 6, height = 5) 




# 手动注释
#B cells (form marker gene: CD79A) 
#CD20+B:MS4A1; CD138+ Plasma:SDC1
#Naive B: (CD20+, CD27−, and CD38−)主要的基因是 IGHD, FCER2, TCL1A, and IL4R
#memory B cells (CD20+, CD27+, and CD38–), 主要的基因是 CD27, AIM2, TNFRSF13B
#germinal center (GC) B cells (CD20+, CD27+, CD38+, and CD138−),主要的基因是S1PI2, LRMP, SUGCT, MME, MKI67, and AICDA
cluster.ids <- c("0" = "B cells", "1" = "plasma cells", "2" = "plasma cells", '3' = 'plasma cells')
Bcells@active.ident <- Bcells$RNA_snn_res.0.1
Bcells <- RenameIdents(Bcells, cluster.ids)                        
Bcells$cell_type <- Bcells@active.ident



cluster_markers <- c('CD79A',"MS4A1", "IGHD", "FCER2", "TCL1A", "IL4R","CD27", 'IGHG1','CD38',"SDC1")
p1 <- DotPlot(Bcells,features = cluster_markers,group.by = 'cell_type') + coord_flip()
ggsave(p1, filename = 'HIM_immune/figure/Bcells_markers.pdf',height = 4,width = 6)

p1 <- DimPlot(Bcells, reduction = "tsne", group.by = "cell_type",label = T)
p2 <- DimPlot(Bcells, reduction = "tsne", group.by = "cell_type",label = F) + NoLegend() + 
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))
ggsave(p1,filename = "HIM_immune/figure/Bcells_cell_type.pdf",width = 5, height = 4)
ggsave(p2,filename = "HIM_immune/figure/Bcells_cell_type.png",width = 4, height = 4)



saveRDS(Bcells,file = "HIM_immune/rds_file/Bcells.rds")

Bcells <- readRDS(file = "HIM_immune/rds_file/Bcells.rds")



#----- B cells cells cell ratio -----------------------------------------------
Cellratio <- prop.table(table(Bcells$cell_type, Bcells$cluster),margin = 1)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
Cellratio$Var2 <- factor(Cellratio$Var2, levels = c('CS1','CS2','CS3'))

library(ggsci)
p1 <- ggplot(Cellratio,aes(x = Var1, y= Freq,fill = Var2)) +
  geom_bar(stat = "identity", position = "stack")+ # dodge
  theme_classic() +
  labs(x='Sample',y = 'Ratio') +
  RotatedAxis() 
plotc <- p1 + scale_fill_manual(values=c('CS1'='#ee6274','CS2'='#ffbf69','CS3'='#64dbd0'))

ggsave(plotc, filename = "HIM_immune/figure/Bcells ratio.pdf",width = 6, height = 4)


library(tidyverse)
library(PieGlyph)
data <- data.frame(table(Bcells$cell_type))
data$CS1 <- Cellratio[Cellratio$Var2 ==  "CS1",]$Freq
data$CS2 <- Cellratio[Cellratio$Var2 ==  "CS2",]$Freq
data$CS3 <- Cellratio[Cellratio$Var2 ==  "CS3",]$Freq

p1 <- ggplot(data = data, aes(x = Var1, y = Freq))+
  geom_pie_glyph(aes(radius = Freq), 
                 slices = c('CS1', 'CS2', 'CS3'), 
                 colour = 'black')+
  scale_fill_manual(values=c('CS1'='#f5a7b1','CS2'='#ffe0b6','CS3'='#a2e9e3'))+
  theme_classic()
ggsave(p1,filename = "HIM_immune/figure/Bcells ratio pie.pdf",width = 8, height = 6)







#----- CAFs cells -----------------------------------------------------------
CAFs <- HIM_integrate[,HIM_integrate@meta.data$cor_type %in% c('myCAFs','apCAFs','Cycling CAFs',
                                                               'iCAFs','imPVL', 'Perivascular-like cells')]

# add data batch 
# batch1 <- c("TNBC001","TNBC002","TNBC003","TNBC004","TNBC005","TNBC006","TNBC007","LEN003")
# batch2 <- c("CID3946","CID3963","CID4465","CID4495","CID4513","CID4515","CID4523","CID44041","CID44971","CID44991")
# batch3 <- c("TN-MH0126","TN-MH0135","TN-SH0106","TN-MH0114-T2","TN-B1-Tum0554","TN-B1-MH4031","TN-B1-MH0177","TN-B1-MH0131")
# 
# batch <- ifelse(CAFs$orig.ident %in% batch1,"batch1",
#                 ifelse(CAFs$orig.ident %in% batch2,"batch2",
#                        ifelse(CAFs$orig.ident %in% batch3,"batch3","reset")))
# table(batch)
# CAFs$batch <- batch


# noiseGenes <- c("MALAT1", "^HSPA", "^MT-", "^MT\\.", "^RP[SL]", "^LOC(0-9)", "^TR(A|B|G|D)V",
#                 "^MTRNR",'^AC(0-9)','^CTB-','^RP11-','^CTD-','^ENSG','^LINC(0-9)')
# CAFs <- CAFs[!grepl(paste0(noiseGenes, collapse = "|"), rownames(CAFs)),] # 64646 23394
# CAFs <- subset(CAFs, subset = nFeature_RNA > 200) # 64646 23269


CAFs <- NormalizeData(CAFs)
CAFs <- FindVariableFeatures(CAFs, selection.method = "vst",
                             nfeatures = 2000)
CAFs <- ScaleData(CAFs, features = VariableFeatures(CAFs))
CAFs <- RunPCA(CAFs, verbose = T)
library(harmony)
CAFs <- RunHarmony(CAFs, group.by.vars = "batch")
ElbowPlot(object = CAFs, ndims = 50)
CAFs <- FindNeighbors(CAFs, reduction = "harmony", dims = 1:20)
library(clustree)
CAFs <- FindClusters(object = CAFs,resolution = c(seq(.1,1.6,.2)))
p1 <- clustree(CAFs@meta.data, prefix = "RNA_snn_res.")
ggsave(p1, filename = "HIM_stromal//figure/CAFs_clustree.png",width = 8, height = 6)
#-------resolution decide using the number of cluster
#-------dims decide using the number of PCA
CAFs <- RunUMAP(CAFs, reduction = "harmony", dims = 1:20)
CAFs <- RunTSNE(CAFs, reduction = "harmony", dims = 1:20)


p1 <- DimPlot(CAFs, reduction = "tsne", group.by = "RNA_snn_res.1.5", label = T)

ggsave(p1, filename = "HIM_stromal/figure/CAFs_tsne.png",width = 12, height = 5)


saveRDS(CAFs,file = "HIM_stromal/rds_file/CAFs.rds")

CAFs <- readRDS(file = "/public/home/scRNA/HIM-omic/HIM_single_cell/result/HIM_stromal/rds_file/CAFs.rds")





#----- CAFs cells type ------------------------------------------------------
CAFs@active.ident <- CAFs$cell_type
CAFs_DEGs <- FindAllMarkers(CAFs,
                            only.pos = T,
                            test.use = "wilcox")

write.csv(CAFs_DEGs,file = 'HIM_CAFs/data/CAFs_DEGs.csv')

CAFs_DEGs <- read.csv(file = 'HIM_CAFs/data/CAFs_DEGs.csv',row.names = 1)
top10 <- CAFs_DEGs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


#----- singleR
# Error: useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE.
# remotes::install_version("matrixStats", version="1.1.0")

library(SingleR)
library(celldex)
# BlueprintEncodeData (Obtain human bulk RNA-seq data from Blueprint and ENCODE)
# DatabaseImmuneCellExpressionData (Obtain human bulk RNA-seq data from DICE)
# HumanPrimaryCellAtlasData (Obtain the HPCA data)
# MonacoImmuneData (Obtain bulk RNA-seq data of sorted human immune cells)
# NovershternHematopoieticData (Obtain bulk microarray expression for sorted hematopoietic cells)

# ImmGenData (Obtain mouse bulk expression data from the Immunologic Genome Project)
# MouseRNAseqData (Obtain mouse bulk expression data of sorted cell populations (RNA-seq)
ref <- readRDS(file = "/public/home/scRNA/referencedata/celldex/BlueprintEncodeData.rds")

tmp <- GetAssayData(CAFs, slot = "data")
cluster <- CAFs@meta.data$RNA_snn_res.1.5
singler <- SingleR(test = tmp, ref = ref, labels = ref$label.main, clusters = cluster,
                   fine.tune = T)
new.cluster.ids <- singler$labels
CAFs@active.ident <- CAFs$RNA_snn_res.1.5
names(new.cluster.ids) <- levels(CAFs)
CAFs <- RenameIdents(CAFs, new.cluster.ids)
CAFs@meta.data$singleR_cluster <- CAFs@active.ident # add singleR celltype to meta.data


p1 <- DimPlot(CAFs, reduction = "tsne", group.by = "singleR_cluster",label = T) #+
# scale_color_manual(values=c('Dendritic cells'='#E71D36','B cells'='black')) + NoLegend()
# p2 <- DimPlot(CAFs, reduction = "umap", group.by = "singleR_cluster",label = T)

ggsave(p1, filename = "HIM_CAFs/figure/CAFs_singleR.png",width = 8, height = 6)
# ggsave(p2, filename = "HIM_CAFs/figure/CAFs_singleR_umap.png",width = 8, height = 6)


#----- marker gene annoation
# CAFs: "PDGFRA", "COL1A1"
# perivascular-like (PVL) cells: "MCAM"/CD146, "ACTA2" and "PDGFRB"
# endothelial cells: "PECAM1"/CD31 and "CD34";
# lymphatic endothelial cells: "LYVE1"
# cycling PVL cells: "MKI67","AURKA"

# myCAFs: "FAP","PDPN", "COL1A1", "COL1A2","POSTN","CFD","DCN", "FBLN2","SEPP1",
# antigen-presenting CAFs: "HLA-DRA", "HLA-DPA1", "HLA-DQA1", "CD74", "C1QA", "C1QB", "C1QC"
# inflammatory CAF: "CXCL12", "C3", "IGF1", "FIGF", "PDGFD", "ALDH1A1","ID2", "EGFR"
# endothelial cells: PECAM1/CD31, CD34, VWF,"VEGFC"
# venular endothelial cells: ACKR1,SELE, SELP, ICAM1, VCAM1,HLA-DRA
# perivascular markers: "ACTA2", "MCAM", "CAV1", "TAGLN", "MYH11", "MYLK",
# immature-PVL (imPVL): "PDGFRB", "CD36", "RGS5"
# lymphatic endothelial cells: LYVE1

# "KDR":VEGFR,"FGFR1","FGFR2","FGFR3","FGFR4","PDGFRB","KIT","RET"

# DotPlot(CAFs, features = c("ACTA2","FAP","S100A4","PDPN", "COL1A1", "COL1A2","POSTN","CFD","DCN", "FBLN2","SEPP1",
#                               "HLA-DRA", "HLA-DPA1", "HLA-DQA1", "CD74", "C1QA", "C1QB", "C1QC",
#                               "CXCL12","CXCL13","C3", "IGF1", "FIGF", "PDGFD", "ALDH1A1","ID2", "EGFR","KLF4","LEPR",
#                               "PECAM1","CD34","VWF","VEGFC","ACKR1","SELE", "SELP","ICAM1","VCAM1",
#                               "MCAM", "CAV1", "TAGLN", "MYH11", "MYLK","PDGFRB", "RGS5","MKI67","AURKA","LYVE1"),
#         group.by = "RNA_snn_res.1.5") +
#   coord_flip()
# CAFs$cell_type <- as.character(CAFs$cell_type)
# CAFs$cell_type[CAFs@meta.data$RNA_snn_res.1.5 %in% c(2,7,11,16,17,18,21,22,24,27,28,32,34,35)] <- 'myCAFs'
# CAFs$cell_type[CAFs@meta.data$RNA_snn_res.1.5 %in% c(8,13,14,31)] <- 'imPVL'
# CAFs$cell_type[CAFs@meta.data$RNA_snn_res.1.5 %in% c(1,9,26,38)] <- 'iCAFs'
# CAFs$cell_type[CAFs@meta.data$RNA_snn_res.1.5 %in% c(19)] <- 'Cycling CAFs'
# CAFs$cell_type[CAFs@meta.data$RNA_snn_res.1.5 %in% c(5,29,37,39)] <- 'apCAFs'
# CAFs$cell_type[CAFs@meta.data$RNA_snn_res.1.5 %in% c(0,10,20,23)] <- 'Perivascular-like cells'
# CAFs$cell_type[CAFs@meta.data$RNA_snn_res.1.5 %in% c(36)] <- 'Cycling Endothelial cells'
# CAFs$cell_type[CAFs@meta.data$RNA_snn_res.1.5 %in% c(6,12,15)] <- 'Endothelial cells'
# CAFs$cell_type[CAFs@meta.data$RNA_snn_res.1.5 %in% c(3,4,25,30)] <- 'Venular endothelial cells'
# CAFs$cell_type[CAFs@meta.data$RNA_snn_res.1.5 %in% c(33)] <- 'Lymphatic endothelial cells'
# CAFs$cell_type <- as.factor(CAFs$cell_type)

p1 <- DimPlot(CAFs, group.by = "cor_type",label = T,reduction = "tsne")
p2 <- DimPlot(CAFs, group.by = "cell_type",label = F,reduction = "tsne",raster = F) + NoLegend() +
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))
ggsave(p1, filename = "HIM_stromal/figure/CAFs_cell_type.pdf",width = 8.3, height = 6)
ggsave(p2, filename = "HIM_stromal/figure/CAFs_cell_type.png",width = 4, height = 4)

# saveRDS(CAFs,file = "HIM_CAFs/rds_file/CAFs.rds")
# 
# CAFs <- readRDS(file = "HIM_CAFs/rds_file/CAFs.rds")





#----- CAFs cells cell ratio -----------------------------------------------
Cellratio <- prop.table(table(CAFs$cor_type, CAFs$cluster),margin = 1)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
Cellratio$Var2 <- factor(Cellratio$Var2, levels = c('CS1','CS2','CS3'))

library(ggsci)
p1 <- ggplot(Cellratio,aes(x = Var1, y= Freq,fill = Var2)) +
  geom_bar(stat = "identity", position = "stack")+ # dodge
  theme_classic() +
  labs(x='Sample',y = 'Ratio') +
  RotatedAxis() 
plotc <- p1 + scale_fill_manual(values=c('CS1'='#ee6274','CS2'='#ffbf69','CS3'='#64dbd0'))

ggsave(plotc, filename = "HIM_stromal/figure/HIM_immune cell ratio.pdf",width = 6, height = 4)


library(tidyverse)
library(PieGlyph)
data <- data.frame(table(CAFs$cor_type))
data$CS1 <- Cellratio[Cellratio$Var2 ==  "CS1",]$Freq
data$CS2 <- Cellratio[Cellratio$Var2 ==  "CS2",]$Freq
data$CS3 <- Cellratio[Cellratio$Var2 ==  "CS3",]$Freq

p1 <- ggplot(data = data, aes(x = Var1, y = Freq))+
  geom_pie_glyph(aes(radius = Freq), 
                 slices = c('CS1', 'CS2', 'CS3'), 
                 colour = 'black')+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values=c('CS1'='#f5a7b1','CS2'='#ffe0b6','CS3'='#a2e9e3'))
ggsave(p1,filename = "HIM_stromal/figure/CAFs cell ratio pie.pdf",width = 6, height = 6)






#----- CAFs marker gene heatmap ---------------------------------------------
#----- cell type heatmap
cluster_markers <- c("ACTA2","FAP","S100A4", "COL1A1", "COL1A2","POSTN",
                     "HLA-DRA", "HLA-DPA1", "HLA-DQA1", "CD74", "C1QA", "C1QB", "C1QC","CXCL13","CFD","DCN", "FBLN2","SEPP1",
                     "CXCL12","C3", "IGF1", "FIGF", "PDGFD", "ID2", "EGFR","MKI67","AURKA",
                     "MCAM", "CAV1", "TAGLN", "MYH11", "MYLK","PDGFRB", "RGS5")
library(pheatmap)
for (i in 1:length(cluster_markers)) {
  ids = CAFs[cluster_markers[i],]
  ids = ids@assays$RNA@data %>% as.numeric()
  assign(cluster_markers[i],tapply(ids, CAFs@meta.data$cor_type,mean))
}

heatmap_matrix <- NULL
for (i in 1:length(cluster_markers)) {
  heatmap_matrix = rbind(heatmap_matrix,get(cluster_markers[i]))
}

row.names(heatmap_matrix) = cluster_markers
ids = c('myCAFs','Cycling CAFs','apCAFs','iCAFs',
        'Perivascular-like cells','imPVL')
heatmap_matrix = heatmap_matrix[,ids]

col <- c(colorRampPalette(c('#3ab2e8','white'))(36),
         colorRampPalette(c('white','#da3446'))(36))

p1 <- pheatmap(heatmap_matrix,color = col,scale = 'row',cluster_cols = F,cluster_rows = F)
ggsave(p1, filename = '/public/home/scRNA/HIM-omic/HIM_single_cell/result/HIM_stromal/figure/CAFscell_annotion.pdf',height = 8,width = 4)


# #----- recell type ------------------------------------------------------------#
# cluster_markers <- c("PDGFRA", "FAP","PDPN","COL1A1","COL1A2","COL3A1","S100A4",
#                      "ACTA2", "MCAM", "CAV1", "TAGLN", "MYH11", "MYLK","MKI67",
#                      "PECAM1", "CD34",
#                      "LYVE1")
# library(pheatmap)
# for (i in 1:length(cluster_markers)) {
#   ids = CAFs[cluster_markers[i],]
#   ids = ids@assays$RNA@data %>% as.numeric()
#   assign(cluster_markers[i],tapply(ids, CAFs@meta.data$recell_type,mean))
# }
# 
# heatmap_matrix <- NULL
# for (i in 1:length(cluster_markers)) {
#   heatmap_matrix = rbind(heatmap_matrix,get(cluster_markers[i]))
# }
# 
# row.names(heatmap_matrix) = cluster_markers
# ids = c("CAFs","PVL cells","Endothelial cells")
# heatmap_matrix = heatmap_matrix[,ids]
# p1 <- pheatmap(heatmap_matrix,scale = 'row',cluster_cols = F,cluster_rows = F)
# ggsave(p1, filename = 'figure/CAFscell_annotion-1.pdf',height = 8,width = 4)




#----- CAFs DEGs ------------------------------------------------------------
CAFs_degs <- FindMarkers(HIM_integrate,
                         ident.1 = "CAFs cells",
                         group.by = "cell_type",
                         only.pos = F,
                         logfc.threshold = 0.25,
                         test.use = "wilcox")
write.csv(CAFs_degs,file = 'CAFs_degs.csv')

CAFs_hm  <- CAFs_degs$avg_log2FC
names(CAFs_hm) <- rownames(CAFs_degs)
CAFs_hm <-  sort(CAFs_hm,decreasing = T)
head(CAFs_hm)

# GO
gsea_CAFs <- GSEA(CAFs_hm,TERM2GENE = GO, pvalueCutoff = 1)
save(gsea_CAFs,file = "CAFs_GOBP.rda")
View(gsea_CAFs@result)

# mapping 
CAFs <- gsea_CAFs@result[order(gsea_CAFs@result$NES,decreasing = T),]
up_pathway <- CAFs[c(1:10),]
up_pathway$Direction <- 'UP'

CAFs <- CAFs[order(CAFs$NES,decreasing = F),]
dow_pathway <- gsea_CAFs[c(1:10),]
dow_pathway$Direction <- 'Down'

merge_pathway <- rbind(up_pathway,dow_pathway)
merge_pathway$ID <- tolower(merge_pathway$ID) 
merge_pathway$IDgsub('gobp_','',merge_pathway$ID)


# merge_pathway <- merge_pathway[order(merge_pathway$NES,decreasing = T),]
# 按照富集值的绝对值大小重新排序术语
merge_pathway$ID <- reorder(merge_pathway$ID, merge_pathway$NES,decreasing = FALSE)


# 绘制双向条形图
p1 <- ggplot(merge_pathway, aes(x=ID, y=NES)) + 
  geom_bar(stat='identity', aes(fill=Direction), width=.7)+ 
  scale_fill_manual(values = c("UP" = "#D14237", "Down" = "#3763B1")) +
  coord_flip()+
  theme_light()
ggsave(p1,filename = "figure/CAFs_GOBP.pdf", width = 8, height = 4)






#----- Endothelial cells -----------------------------------------------------------
Endothelial <- HIM_integrate[,HIM_integrate@meta.data$cor_type %in% c('Cycling Endothelial cells','Endothelial cells',
                                                                      'Lymphatic endothelial cells',
                                                                      'Venular endothelial cells')]

# # add data batch 
# batch1 <- c("TNBC001","TNBC002","TNBC003","TNBC004","TNBC005","TNBC006","TNBC007","LEN003")
# batch2 <- c("CID3946","CID3963","CID4465","CID4495","CID4513","CID4515","CID4523","CID44041","CID44971","CID44991")
# batch3 <- c("TN-MH0126","TN-MH0135","TN-SH0106","TN-MH0114-T2","TN-B1-Tum0554","TN-B1-MH4031","TN-B1-MH0177","TN-B1-MH0131")
# 
# batch <- ifelse(Endothelial$orig.ident %in% batch1,"batch1",
#                 ifelse(Endothelial$orig.ident %in% batch2,"batch2",
#                        ifelse(Endothelial$orig.ident %in% batch3,"batch3","reset")))
# table(batch)
# Endothelial$batch <- batch


# noiseGenes <- c("MALAT1", "^HSPA", "^MT-", "^MT\\.", "^RP[SL]", "^LOC(0-9)", "^TR(A|B|G|D)V",
#                 "^MTRNR",'^AC(0-9)','^CTB-','^RP11-','^CTD-','^ENSG','^LINC(0-9)')
# Endothelial <- Endothelial[!grepl(paste0(noiseGenes, collapse = "|"), rownames(Endothelial)),] # 64646 23394
# Endothelial <- subset(Endothelial, subset = nFeature_RNA > 200) # 64646 23269


Endothelial <- NormalizeData(Endothelial)
Endothelial <- FindVariableFeatures(Endothelial, selection.method = "vst",
                                    nfeatures = 2000)
Endothelial <- ScaleData(Endothelial, features = VariableFeatures(Endothelial))
Endothelial <- RunPCA(Endothelial, verbose = T)
library(harmony)
Endothelial <- RunHarmony(Endothelial, group.by.vars = "batch")
ElbowPlot(object = Endothelial, ndims = 50)
Endothelial <- FindNeighbors(Endothelial, reduction = "harmony", dims = 1:20)
library(clustree)
Endothelial <- FindClusters(object = Endothelial,resolution = c(seq(.1,1.6,.2)))
p1 <- clustree(Endothelial@meta.data, prefix = "RNA_snn_res.")
ggsave(p1, filename = "HIM_stromal/figure/Endothelial_clustree.png",width = 8, height = 6)
#-------resolution decide using the number of cluster
#-------dims decide using the number of PCA
Endothelial <- RunUMAP(Endothelial, reduction = "harmony", dims = 1:20)
Endothelial <- RunTSNE(Endothelial, reduction = "harmony", dims = 1:20)


p1 <- DimPlot(Endothelial, reduction = "tsne", group.by = "RNA_snn_res.1.5", label = T)

ggsave(p1, filename = "HIM_stromal/figure/Endothelial_tsne.png",width = 12, height = 5)


saveRDS(Endothelial,file = "HIM_stromal/rds_file/Endothelial.rds")

Endothelial <- readRDS(file = "HIM_stromal/rds_file/Endothelial.rds")





#----- Endothelial cells type ------------------------------------------------------
Endothelial@active.ident <- Endothelial$RNA_snn_res.1.5
Endothelial_DEGs <- FindAllMarkers(Endothelial,
                                   only.pos = T,
                                   test.use = "wilcox")

write.csv(Endothelial_DEGs,file = 'HIM_stromal/data/Endothelial_DEGs.csv')

Endothelial_DEGs <- read.csv(file = 'HIM_stromal/data/Endothelial_DEGs.csv',row.names = 1)

top10 <- Endothelial_DEGs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
plot = DoHeatmap(Endothelial, # input data is scale.data
                 features = unique(top10$gene),
                 group.by = "RNA_snn_res.1.5",
                 group.bar = T)


#----- singleR
# Error: useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE.
# remotes::install_version("matrixStats", version="1.1.0")

library(SingleR)
library(celldex)
# BlueprintEncodeData (Obtain human bulk RNA-seq data from Blueprint and ENCODE)
# DatabaseImmuneCellExpressionData (Obtain human bulk RNA-seq data from DICE)
# HumanPrimaryCellAtlasData (Obtain the HPCA data)
# MonacoImmuneData (Obtain bulk RNA-seq data of sorted human immune cells)
# NovershternHematopoieticData (Obtain bulk microarray expression for sorted hematopoietic cells)

# ImmGenData (Obtain mouse bulk expression data from the Immunologic Genome Project)
# MouseRNAseqData (Obtain mouse bulk expression data of sorted cell populations (RNA-seq)
ref <- readRDS(file = "/public/home/scRNA/referencedata/celldex/BlueprintEncodeData.rds")

tmp <- GetAssayData(Endothelial, slot = "data")
cluster <- Endothelial@meta.data$RNA_snn_res.1.5
singler <- SingleR(test = tmp, ref = ref, labels = ref$label.main, clusters = cluster,
                   fine.tune = T)
new.cluster.ids <- singler$labels
Endothelial@active.ident <- Endothelial$RNA_snn_res.1.5
names(new.cluster.ids) <- levels(Endothelial)
Endothelial <- RenameIdents(Endothelial, new.cluster.ids)
Endothelial@meta.data$singleR_cluster <- Endothelial@active.ident # add singleR celltype to meta.data


p1 <- DimPlot(Endothelial, reduction = "tsne", group.by = "singleR_cluster",label = T) #+
# scale_color_manual(values=c('Dendritic cells'='#E71D36','B cells'='black')) + NoLegend()
# p2 <- DimPlot(Endothelial, reduction = "umap", group.by = "singleR_cluster",label = T)

ggsave(p1, filename = "HIM_stromal/figure/Endothelial_singleR.png",width = 8, height = 6)
# ggsave(p2, filename = "HIM_stromal/figure/Endothelial_singleR_umap.png",width = 8, height = 6)








#----- marker gene annoation 
# CAFs: "PDGFRA", "COL1A1"
# perivascular-like (PVL) cells: "MCAM"/CD146, "ACTA2" and "PDGFRB"
# endothelial cells: "PECAM1"/CD31 and "CD34";
# lymphatic endothelial cells: "LYVE1"
# cycling PVL cells: "MKI67","AURKA"

# myCAFs: "FAP","PDPN", "COL1A1", "COL1A2","POSTN","CFD","DCN", "FBLN2","SEPP1",
# antigen-presenting CAFs: "HLA-DRA", "HLA-DPA1", "HLA-DQA1", "CD74", "C1QA", "C1QB", "C1QC"
# inflammatory CAF: "CXCL12", "C3", "IGF1", "FIGF", "PDGFD", "ALDH1A1","ID2", "EGFR"
# endothelial cells: PECAM1/CD31, CD34, VWF,"VEGFC"
# venular endothelial cells: ACKR1,SELE, SELP, ICAM1, VCAM1,HLA-DRA
# perivascular markers: "ACTA2", "MCAM", "CAV1", "TAGLN", "MYH11", "MYLK",
# immature-PVL (imPVL): "PDGFRB", "CD36", "RGS5"
# lymphatic endothelial cells: LYVE1

# "KDR":VEGFR,"FGFR1","FGFR2","FGFR3","FGFR4","PDGFRB","KIT","RET"

# DotPlot(Endothelial, features = c("KLF4","LEPR","PECAM1","CD34","VWF","VEGFC","ACKR1","SELE", "SELP","ICAM1","VCAM1",
#                               "MCAM", "CAV1", "TAGLN", "MYH11", "MYLK","PDGFRB", "RGS5","MKI67","AURKA","LYVE1"),
#         group.by = "RNA_snn_res.1.5") +
#   coord_flip()
# Endothelial$cell_type <- as.character(Endothelial$cell_type)
# Endothelial$cell_type[Endothelial@meta.data$RNA_snn_res.1.5 %in% c(2,7,11,16,17,18,21,22,24,27,28,32,34,35)] <- 'myCAFs'
# Endothelial$cell_type[Endothelial@meta.data$RNA_snn_res.1.5 %in% c(8,13,14,31)] <- 'imPVL'
# Endothelial$cell_type[Endothelial@meta.data$RNA_snn_res.1.5 %in% c(1,9,26,38)] <- 'iCAFs'
# Endothelial$cell_type[Endothelial@meta.data$RNA_snn_res.1.5 %in% c(19)] <- 'Cycling CAFs'
# Endothelial$cell_type[Endothelial@meta.data$RNA_snn_res.1.5 %in% c(5,29,37,39)] <- 'apCAFs'
# Endothelial$cell_type[Endothelial@meta.data$RNA_snn_res.1.5 %in% c(0,10,20,23)] <- 'Perivascular-like cells'
# Endothelial$cell_type[Endothelial@meta.data$RNA_snn_res.1.5 %in% c(36)] <- 'Cycling Endothelial cells'
# Endothelial$cell_type[Endothelial@meta.data$RNA_snn_res.1.5 %in% c(6,12,15)] <- 'Endothelial cells'
# Endothelial$cell_type[Endothelial@meta.data$RNA_snn_res.1.5 %in% c(3,4,25,30)] <- 'Venular endothelial cells'
# Endothelial$cell_type[Endothelial@meta.data$RNA_snn_res.1.5 %in% c(33)] <- 'Lymphatic endothelial cells'
# Endothelial$cell_type <- as.factor(Endothelial$cell_type)

p1 <- DimPlot(Endothelial, group.by = "cor_type",label = T,reduction = "tsne")
p2 <- DimPlot(Endothelial, group.by = "cor_type",label = F,reduction = "tsne",raster = F) + NoLegend() +
  theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_blank(),axis.line = element_blank(),axis.ticks=element_blank(),
        plot.margin =  margin(t = -1, r = -1, b = -1, l = -1))
ggsave(p1, filename = "HIM_stromal/figure/Endothelial_cell_type-1.pdf",width = 8.3, height = 6)
ggsave(p2, filename = "HIM_stromal/figure/Endothelial_cell_type.pdf",width = 4, height = 4)

# saveRDS(Endothelial,file = "HIM_stromal/rds_file/Endothelial.rds")

Endothelial <- readRDS(file = "/public/home/scRNA/HIM-omic/HIM_single_cell/result/HIM_stromal/rds_file/Endothelial.rds")


#----- Endothelial cells cell ratio -----------------------------------------------
Cellratio <- prop.table(table(Endothelial$cor_type, Endothelial$cluster),margin = 1)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
Cellratio$Var2 <- factor(Cellratio$Var2, levels = c('CS1','CS2','CS3'))

library(ggsci)
p1 <- ggplot(Cellratio,aes(x = Var1, y= Freq,fill = Var2)) +
  geom_bar(stat = "identity", position = "stack")+ # dodge
  theme_classic() +
  labs(x='Sample',y = 'Ratio') +
  RotatedAxis() 
plotc <- p1 + scale_fill_manual(values=c('CS1'='#ee6274','CS2'='#ffbf69','CS3'='#64dbd0'))

ggsave(plotc, filename = "HIM_stromal/figure/HIM_immune cell ratio.pdf",width = 6, height = 4)


library(tidyverse)
library(PieGlyph)
data <- data.frame(table(Endothelial$cor_type))
data$CS1 <- Cellratio[Cellratio$Var2 ==  "CS1",]$Freq
data$CS2 <- Cellratio[Cellratio$Var2 ==  "CS2",]$Freq
data$CS3 <- Cellratio[Cellratio$Var2 ==  "CS3",]$Freq

p1 <- ggplot(data = data, aes(x = Var1, y = Freq))+
  geom_pie_glyph(aes(radius = Freq), 
                 slices = c('CS1', 'CS2', 'CS3'), 
                 colour = 'black')+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values=c('CS1'='#f5a7b1','CS2'='#ffe0b6','CS3'='#a2e9e3'))
ggsave(p1,filename = "HIM_stromal/figure/Endothelial cell ratio pie.pdf",width = 6, height = 6)






#----- Endothelial marker gene heatmap ---------------------------------------------
#----- cell type heatmap
cluster_markers <- c("KLF4","LEPR","PECAM1","CD34","VWF","VEGFC","ACKR1","SELE", 
                     "SELP","ICAM1","MKI67","AURKA","LYVE1")
library(pheatmap)
for (i in 1:length(cluster_markers)) {
  ids = Endothelial[cluster_markers[i],]
  ids = ids@assays$RNA@data %>% as.numeric()
  assign(cluster_markers[i],tapply(ids, Endothelial@meta.data$cor_type,mean))
}

heatmap_matrix <- NULL
for (i in 1:length(cluster_markers)) {
  heatmap_matrix = rbind(heatmap_matrix,get(cluster_markers[i]))
}

row.names(heatmap_matrix) = cluster_markers
ids = c('Endothelial cells','Venular endothelial cells','Cycling Endothelial cells',
        'Lymphatic endothelial cells')
heatmap_matrix = heatmap_matrix[,ids]

col <- c(colorRampPalette(c('#3ab2e8','white'))(36),
         colorRampPalette(c('white','#da3446'))(36))

p1 <- pheatmap(heatmap_matrix, color = col,scale = 'row',cluster_cols = F,cluster_rows = F)
ggsave(p1, filename = '/public/home/scRNA/HIM-omic/HIM_single_cell/result/HIM_stromal/figure/Endothelial_annotion.pdf',height = 8,width = 4)


#----- recell type ------------------------------------------------------------#
cluster_markers <- c("PDGFRA", "FAP","PDPN","COL1A1","COL1A2","COL3A1","S100A4",
                     "ACTA2", "MCAM", "CAV1", "TAGLN", "MYH11", "MYLK","MKI67",
                     "PECAM1", "CD34",
                     "LYVE1")
library(pheatmap)
for (i in 1:length(cluster_markers)) {
  ids = Endothelial[cluster_markers[i],]
  ids = ids@assays$RNA@data %>% as.numeric()
  assign(cluster_markers[i],tapply(ids, Endothelial@meta.data$recell_type,mean))
}

heatmap_matrix <- NULL
for (i in 1:length(cluster_markers)) {
  heatmap_matrix = rbind(heatmap_matrix,get(cluster_markers[i]))
}

row.names(heatmap_matrix) = cluster_markers
ids = c("CAFs","PVL cells","Endothelial cells")
heatmap_matrix = heatmap_matrix[,ids]


p1 <- pheatmap(heatmap_matrix,scale = 'row',cluster_cols = F,cluster_rows = F)
ggsave(p1, filename = 'figure/Endothelialcell_annotion-1.pdf',height = 8,width = 4)





#----- correlation across cell subsets ----------------------------------------------------
# 对HIM_integrate重命名
HIM_integrate$cor_type <- 'normal'

HIM_integrate$cor_type[HIM_integrate@meta.data$cell_type == 'Mast cells'] <- 'Mast cells'
HIM_integrate$cor_type[HIM_integrate@meta.data$cell_type == 'pDCs'] <- 'pDCs'
HIM_integrate$cor_type[HIM_integrate@meta.data$cell_type == 'NK cells'] <- 'NK cells'


class(Tcells$cell_type)
for (i in c(as.character(levels(Tcells$cell_type)))){
  tmp <- rownames(Tcells@meta.data[Tcells$cell_type == i,])
  HIM_integrate$cor_type[rownames(HIM_integrate@meta.data) %in% tmp] <- i 
}

class(Bcells$cell_type)
for (i in c(as.character(levels(Bcells$cell_type)))){
  tmp <- rownames(Bcells@meta.data[Bcells$cell_type == i,])
  HIM_integrate$cor_type[rownames(HIM_integrate@meta.data) %in% tmp] <- i
}

class(Mac$cell_type)
for (i in c(as.character(levels(Mac$cell_type)))){
  tmp <- rownames(Mac@meta.data[Mac$cell_type == i,])
  HIM_integrate$cor_type[rownames(HIM_integrate@meta.data) %in% tmp] <- i
}

class(Endothelial$cell_type)
for (i in c(as.character(levels(Endothelial$cell_type)))){
  tmp <- rownames(Endothelial@meta.data[Endothelial$cell_type == i,])
  HIM_integrate$cor_type[rownames(HIM_integrate@meta.data) %in% tmp] <- i
}

class(cancer_remor$cluster)
cancer_remor$cluster <- as.factor(cancer_remor$cluster)
for (i in c(as.character(levels(cancer_remor$cluster)))){
  tmp <- rownames(cancer_remor@meta.data[cancer_remor$cluster == i,])
  HIM_integrate$cor_type[rownames(HIM_integrate@meta.data) %in% tmp] <- i
}





# verify
levels(HIM-integrate$cor_type)
levels(Tcells$cell_type)
levels(Bcells$cell_type)
levels(Mac$cell_type)
levels(stromal$cell_type)
levels(cancer_remor$cluster)




saveRDS(HIM_integrate,file = "HIM-integrate/rds_file/HIM_integrate_326.rds") 


HIM_integrate <- readRDS(file = "HIM_single_cell/result/HIM-integrate/rds_file/HIM_integrate_326.rds") 




#----- 绘制相关性
# https://zhuanlan.zhihu.com/p/28076189
library(corrplot)

# remove normal levels cells
dim(HIM_integrate) # 92782 173137
tmp <- HIM_integrate[,HIM_integrate@meta.data$cor_type != 'normal']
dim(tmp)  # 92782 167966

class(tmp$cor_type)
tmp$cor_type <- as.factor(tmp$cor_type)

# 计算相关性
cor_table <- prop.table(table(tmp$orig.ident,tmp$cor_type),margin = 2)
cor_table <- as.data.frame(unclass(cor_table))
sum(cor_table$apCAFs)

col=c(rev(brewer.pal(9,"Blues")),"White",brewer.pal(9,"Reds"))
# 计算相关性矩阵
M <- cor(cor_table,method="spearman")



order_vector <- c("CD4-ISG","CD4-Th","CD4-Tn","CD4-Tpro","CD4-Treg","CD8-ISG",
                  "CD8-Teff","CD8-Tex","CD8-Tn","CD8-Tpro","T-SEPTIN","cDC1","NK cells","B cells",
                  "cDC2","cDC3","pDCs","TAM_C1QC","TAM_ISG15","TAM_SPP1","TAM-pro",
                  "Mono_CD14","Mast cells","plasma cells","Neutrophil",
                  "apCAFs","Cycling CAFs","myCAFs","Cycling Endothelial cells","imPVL","iCAFs","Endothelial cells",
                  "Lymphatic endothelial cells","Perivascular-like cells","Venular endothelial cells",
                  "CS1","CS2","CS3") # "normal",

M <- M[order_vector, order_vector]  # 重新排序相关性矩阵
# testRes$p_adjusted <- testRes$p_adjusted[order_vector,order_vector]

library(pheatmap)
library(RColorBrewer)
# col = c(rev(brewer.pal(9,'Blues')),"White",brewer.pal(9,'qual'))

col <- c(colorRampPalette(c('#3ab2e8', 'white'))(30),
         rep("white", 12),  # 插入更多白色
         colorRampPalette(c('white', '#da3446'))(30))

# 检查颜色数量
length(col)  # 应该是72个颜色

# breaks <- seq(-1, 1, length.out = length(col) + 1) # 调整 breaks

pdf('HIM_single_cell/result/HIM-integrate/figure/cell cluster correlation-6.pdf',width = 10,height = 10)
pheatmap(M,color = col,cluster_rows = F, cluster_cols = F,)
dev.off()

pdf('HIM_single_cell/result/HIM-integrate/figure/cell cluster correlation-5.pdf',width = 10,height = 10)
corrplot(as.matrix(M),
         type = "full",
         # p.mat = testRes$p_adjusted,
         method = "square",
         tl.cex = 0.85,
         cl.cex = 0.85,
         col = col,
         tl.col = "black",
         insig='blank',
         # sig.level = 0.05,
         cl.pos = "b")# 图例放在底部
dev.off()





# # 计算p值
# testRes = cor.mtest(M, conf.level = 0.95) # 一次性检验多个相关性的函数
# 
# # 扁平化p值
# p_values <- as.vector(testRes$p)
# 
# # Benjamini-Hochberg校正
# p_adjusted <- p.adjust(p_values, method = "BH")
# 
# # 将校正后的p值重新组织成与相关性矩阵相同的形状
# p_adjusted_matrix <- matrix(p_adjusted, nrow = nrow(testRes$p))
# rownames(p_adjusted_matrix) <- rownames(testRes$p)
# colnames(p_adjusted_matrix) <- colnames(testRes$p)
# 
# # 过滤显著性相关性
# testRes$p_adjusted <- p_adjusted_matrix


# color_vector <- c(colorRampPalette(c('#313695','white'))(36),
#                   colorRampPalette(c('white','#ff7f00'))(36))

# pdf('HIM-integrate/figure/cell cluster correlation.pdf',width = 10,height = 10)
# corrplot(M, p.mat = testRes$p_adjusted, type = "full",col = color_vector,insig='blank',sig.level = 0.05)
# dev.off()
# 
# 
# pdf('HIM-integrate/figure/cell cluster correlation-1.pdf',width = 10,height = 10)
# corrplot(M, p.mat = testRes$p, type = "upper",col = color_vector,insig='blank',sig.level = 0.01)
# dev.off()












#----- HIM_integrate seurat2scanpy ---------------------------------------------
library(Seurat)
library(SeuratDisk)
HIM_integrate <- readRDS(file = "HIM-integrate/rds_file/HIM_integrate_1227.rds") # 92782 173137

HIM_integrate@active.ident <- as.factor(HIM_integrate$cell_type) 

Seurat::SetDefaultAssay(HIM_integrate, assay = "RNA")

options(Seurat.object.assay.version = "v3")
Seu2scan <- CreateSeuratObject(counts = GetAssayData(HIM_integrate, assay="RNA", layer = "counts"),
                               #data = GetAssayData(HIM_integrate, assay="RNA", layer = "data"),
                               meta.data = HIM_integrate@meta.data)

# 格式转换需要两步完成
SaveH5Seurat(Seu2scan,filename = "HIM-integrate/data/HIM_integrate.h5Seurat",overwrite = T)
Convert(source = "HIM-integrate/data/HIM_integrate.h5Seurat",dest = "h5ad",overwrite = T)













#----- HIM_integrate cellchat ---------------------------------------------------
# 细胞通讯分析的工具主要有: cellphoneDB、Celltalker、iTALK 、NicheNet、CellChat
# https://github.com/sqjin/CellChat
# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
# https://cloud.tencent.com/developer/article/1865540
# https://www.jianshu.com/p/58a2ac0e0448
# https://cloud.tencent.com/developer/article/1853994

#----- cell chat
library(CellChat)
# CellChat需要用户输入两个东西：细胞的基因表达矩阵和细胞的标签。
# 对于表达矩阵，行是基因（行名），列是细胞（列名），矩阵需要是标准化以后的，若输入的是原始的counts，可以用CellChat提供的normalizeData进行标准化。
# 细胞的标签信息需要做成一个dataframe，要求行名是细胞的名字。
# loading data
HIM_integrate <- readRDS(file = "HIM-integrate/rds_file/HIM_integrate_326.rds") 

tmp <- HIM_integrate[,HIM_integrate@meta.data$cor_type != 'normal']

data.input = tmp@assays$RNA@data # normalized data matrix
meta = subset(tmp@meta.data,select = c(orig.ident,cor_type)) # meta data 行名为细胞名
unique(meta$cor_type)
head(meta)

# creat cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cor_type")
# add cell information and labels to meta
cellchat <- addMeta(cellchat, meta = meta)# Add the cell information into meta slot
cellchat <- setIdent(cellchat, ident.use = "cor_type") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# loading ligand-receptor database
CellChatDB <- CellChatDB.human #mouse use CellChatDB.mouse
showDatabaseCategory(CellChatDB) # Show the description of CellChatDB databse
dplyr::glimpse(CellChatDB$interaction)    # Show the structure of the database

CellChatDB.use <- subsetDB(CellChatDB, search =  'Cell-Cell Contact')    # search =  Secreted Signaling, Cell-Cell Contact, ECM-Receptor
cellchat@DB <- CellChatDB.use    # set the used database in the object

# The expression data was preprocessed
# 首先识别出在一类细胞中过表达的配体或受体，然后把基因表达数据投射到蛋白质互作网络中。
# 只要配体或受体有一个过表达，则该配体-受体互作对就被识别出来。
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost

# # 使用多线程计算
# library(future)
# future::plan("multicore", workers = 2) # do parallel
# future::plan("sequential") # 切换回顺序计算


# Identify over-expressed signaling genes associated with each cell group
cellchat <- identifyOverExpressedGenes(cellchat)
# Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
cellchat <- identifyOverExpressedInteractions(cellchat) 
# Project gene expression data onto a protein-protein interaction network
cellchat <- projectData(cellchat, PPI.human)


# To infer cell-to-cell communication
# CellChat通过给每个相互作用分配一个概率值并进行置换检验来推断具有生物学意义的细胞间通信。
# 推测的每个配受体对的细胞通信网络和信号通路分别存储在“net”和“netP”中。
# future = TRUE
cellchat <- computeCommunProb(cellchat) # 较慢
cellchat <- filterCommunication(cellchat, min.cells = 10) # 过滤掉小于10个细胞的胞间通讯网络
cellchat <- computeCommunProbPathway(cellchat)
# Calculate the aggregated network by counting the number of links or summarizing the communication probability
cellchat <- aggregateNet(cellchat) 



saveRDS(cellchat, file = "HIM-integrate/data/Secreted.rds")
saveRDS(cellchat, file = "HIM-integrate/data/Contact.rds")
saveRDS(cellchat, file = "HIM-integrate/data/ECM.rds")


#----- Load CellChat object of tumor Secreted ---------------------------------#
cellchat_Secreted <- readRDS(file = "HIM-integrate/data/Secreted.rds")
cellchat_Contact <- readRDS(file = "HIM-integrate/data/Contact.rds")
cellchat_ECM <- readRDS(file = "HIM-integrate/data/ECM.rds")

cellchat_Secreted <- netAnalysis_computeCentrality(cellchat_Secreted, slot.name = "netP")
cellchat_Contact <- netAnalysis_computeCentrality(cellchat_Contact, slot.name = "netP")
cellchat_ECM <- netAnalysis_computeCentrality(cellchat_ECM, slot.name = "netP")

cellchat_Secreted <- updateCellChat(cellchat_Secreted)
cellchat_Contact <- updateCellChat(cellchat_Contact)
cellchat_ECM <- updateCellChat(cellchat_ECM)

# Chord diagram
# CS1
levels(cellchat_Secreted@idents)
pdf("HIM-integrate/figure/cellchat_CS1_Secreted.pdf",height = 8,width = 8)
netVisual_chord_gene(cellchat_Secreted, sources.use = 16, targets.use = c(1,2,3:15,19,23,25,26,29,30,32,34:37), slot.name = "netP", lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_Secreted, sources.use = c(2,3:10,12:15,19,29,33,34), targets.use = 16, slot.name = "netP", lab.cex = 0.5,legend.pos.y = 30)
dev.off()


pdf("HIM-integrate/figure/cellchat_CS1_Contact.pdf",height = 8,width = 8)
netVisual_chord_gene(cellchat_Contact, sources.use = 16, targets.use = c(1,2,3:15,19,23,25,26,29,30,32,34:37), slot.name = "netP", lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_Contact, sources.use = c(2,3:10,12:15,19,29,33,34), targets.use = 16, slot.name = "netP", lab.cex = 0.5,legend.pos.y = 30)
dev.off()


# CS2
levels(cellchat_Secreted@idents)
pdf("HIM-integrate/figure/cellchat_CS2_Secreted.pdf",height = 8,width = 8)
netVisual_chord_gene(cellchat_Secreted, sources.use = 17, targets.use = c(20:22,24,31,38), slot.name = "netP", lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_Secreted, sources.use = c(20:22,24,31,38), targets.use = 17, slot.name = "netP", lab.cex = 0.5,legend.pos.y = 30)
dev.off()

levels(cellchat_Contact@idents)
pdf("HIM-integrate/figure/cellchat_CS2_Contact.pdf",height = 8,width = 8)
netVisual_chord_gene(cellchat_Contact, sources.use = 17, targets.use = c(20:22,24,31,38), slot.name = "netP", lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_Contact, sources.use = c(20:22,24,31,38), targets.use = 17, slot.name = "netP", lab.cex = 0.5,legend.pos.y = 30)
dev.off()

pdf("HIM-integrate/figure/cellchat_cellchat_ECM.pdf",height = 8,width = 8)
netVisual_chord_gene(cellchat_ECM, sources.use = 17, targets.use = c(1,20,21,23,24,26,28,30,35), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_ECM, sources.use = c(1,20,21,23,24,26,28,30,35), targets.use = 17, lab.cex = 0.5,legend.pos.y = 30)
dev.off()


# T-SEPTIN
levels(cellchat_Secreted@idents)
pdf("HIM-integrate/figure/cellchat_T-SEPTIN_Secreted.pdf",height = 8,width = 8)
netVisual_chord_gene(cellchat_Secreted, sources.use = 33, targets.use = c(21,22,31), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_Secreted, sources.use = c(21,22,31), targets.use = 33, lab.cex = 0.5,legend.pos.y = 30)
dev.off()

levels(cellchat_Contact@idents)
pdf("HIM-integrate/figure/cellchat_T-SEPTIN_Contact.pdf",height = 8,width = 8)
netVisual_chord_gene(cellchat_Contact, sources.use = 33, targets.use = c(21,22,31), slot.name = "netP", lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_Contact, sources.use = c(21,22,31), targets.use = 33, slot.name = "netP", lab.cex = 0.5,legend.pos.y = 30)
dev.off()


#----- all cell type Chord diagram
levels(cellchat_Secreted@idents)
groupSize <- as.numeric(table(cellchat_Secreted@idents)) # number of cells in each cell group

pdf("HIM-integrate/figure/cellchat_Secreted.pdf",height = 20,width = 20)
par(mfrow = c(1,1), xpd=TRUE)
par(mar = c(0, 0, 0, 0)) 
mat <- cellchat_Secreted@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


levels(cellchat_Contact@idents)
groupSize <- as.numeric(table(cellchat_Contact@idents)) # number of cells in each cell group

pdf("HIM-integrate/figure/cellchat_Contact.pdf",height = 20,width = 20)
par(mfrow = c(1,1), xpd=TRUE)
par(mar = c(0, 0, 0, 0)) 
mat <- cellchat_Contact@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

levels(cellchat_ECM@idents)
groupSize <- as.numeric(table(cellchat_ECM@idents)) # number of cells in each cell group

pdf("HIM-integrate/figure/cellchat_ECM.pdf",height = 20,width = 20)
par(mfrow = c(1,1), xpd=TRUE)
par(mar = c(0, 0, 0, 0)) 
mat <- cellchat_ECM@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


#----- all cell type Heatmap
pdf("HIM-integrate/figure/cellchat_Heatmap.pdf",height = 10,width = 20)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_Contact, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_Contact, pattern = "incoming")
ht1 + ht2
dev.off()


# Bubble plot
# CS1
levels(cellchat_Secreted@idents)
netVisual_bubble(cellchat_Secreted, sources.use = 16, targets.use = c(2,3:10,12:15,19,29,33,34), remove.isolate = FALSE)

levels(cellchat_Contact@idents)
netVisual_bubble(cellchat_Contact, sources.use = 16, targets.use = c(2,3:10,12:15,19,29,33,34), remove.isolate = FALSE)

# CS2
levels(cellchat_Secreted@idents)
netVisual_bubble(cellchat_Secreted, sources.use = 17, targets.use = c(1,20,21,23,24,26,28,30,35), remove.isolate = FALSE)

levels(cellchat_Contact@idents)
netVisual_bubble(cellchat_Contact, sources.use = 17, targets.use = c(1,20,21,23,24,26,28,30,35), remove.isolate = FALSE)

levels(cellchat_ECM@idents)
netVisual_bubble(cellchat_ECM, sources.use = 17, targets.use = c(1,20,21,23,24,26,28,30,35), remove.isolate = FALSE)





















#----- HIM_integrate cellphoneDB -----------------------------------------------
HIM_integrate <- readRDS(file = "HIM-integrate/rds_file/HIM_integrate_1227.rds") 


# 导出metadata行名和细胞类型
meta_data <- cbind(rownames(HIM_integrate@meta.data),
                   HIM_integrate@meta.data[,'cor_type', drop=F]) 
meta_data <- as.matrix(meta_data)
table(is.na(meta_data)) #细胞类型不能为空
colnames(meta_data) <- c('barcode_sample','cell_type')
head(meta_data)
write.csv(meta_data,'HIM-integrate/data/cellphonedb/metadata.csv',row.names=F)


# 将counts导出为h5ad文件
# https://cloud.tencent.com/developer/article/2193912?areaSource=104001.13&traceId=H74biujC4GYFLR5uInZha
library(Seurat)
library(SeuratDisk)

# Loom文件以矩阵的形式存储单细胞数据，其中行表示细胞，列表示基因或特征。导致使用scanpy读取数据后有缺失
# h5ad文件也以类似的方式存储数据，但同时还可以包含其他附加的注释信息，如样本信息、细胞类型等


# 将scale.data删除
annadata <- DietSeurat(HIM_integrate,counts = TRUE,data = TRUE,scale.data = FALSE,features = NULL,
                       assays = NULL,dimreducs = NULL,graphs = NULL,misc = TRUE) # 目前还不能删除data
annadata@reductions <- HIM_integrate@reductions


# 格式转换需要两步完成
SaveH5Seurat(annadata,filename = "HIM-integrate/data/HIM_integrate.h5Seurat",overwrite = T)
Convert(source = "HIM-integrate/data/HIM_integrate.h5Seurat",dest = "h5ad",overwrite = T)




# 使用R对结果进行可视化
# devtools::install_github('zktuong/ktplots', dependencies = TRUE)
library(Seurat)
library(dplyr)
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
library(ktplots)
library(SingleCellExperiment)
library(reticulate)

pvals <- read.delim("HIM-integrate/data/cellphonedb/statistical_analysis_pvalues_12_28_2023_082602.txt", check.names = FALSE)
means <- read.delim("HIM-integrate/data/cellphonedb/statistical_analysis_means_12_28_2023_082602.txt", check.names = FALSE)

pdf(file = 'HIM-integrate/figure/cellphonedb_plot.pdf',height = 80,width = 8)
plot_cpdb(scdata = HIM_integrate, cell_type1 = '', cell_type2 = '', #这里的cell_type1、2指定互作细胞，如果不指定则默认所有细胞,
          celltype_key = 'cor_type',
          means = means, pvals = pvals,
          keep_significant_only = TRUE)  + coord_flip() + small_axis(fontsize = 5) 
dev.off()

pdf(file = 'HIM-integrate/figure/cellphonedb_T−SEPTIN.pdf',height = 80,width = 8)
plot_cpdb(scdata = HIM_integrate, cell_type1 = 'T−SEPTIN', cell_type2 = 'Endothelial cells', #这里的cell_type1、2指定互作细胞，如果不指定则默认所有细胞,
          celltype_key = 'cor_type',
          means = means, pvals = pvals,
          keep_significant_only = TRUE)  + coord_flip() + small_axis(fontsize = 5) 
dev.off()












#----- CS2 drug2cell -----------------------------------------------------------
HIM_integrate <- readRDS(file = "/public/home/scRNA/HIM-omic/HIM_single_cell/result/HIM-integrate/rds_file/HIM_integrate_326.rds")


library(SeuratDisk)
# 格式转换需要两步完成
SaveH5Seurat(HIM_integrate,filename = "/public/home/scRNA/HIM-omic/HIM_single_cell/result/HIM-integrate/data/HIM_integrate.h5Seurat",overwrite = T)
Convert(source = "/public/home/scRNA/HIM-omic/HIM_single_cell/result/HIM-integrate/data/HIM_integrate.h5Seurat",dest = "h5ad",overwrite = T)
# 最终得到integrate_pre.h5Seurat和integrate_pre.h5ad两个文件


#----- Lenvatinib targets ---------------------------------------------------------
CS2_drug2cell <- HIM_integrate[,HIM_integrate@meta.data$cluster %in% c('CS2')]

lenva_target <- c('FLT1','KDR','FLT4','FGFR1','FGFR2','FGFR3','FGFR4',
                  'PDGFRA','KIT','RET') # FLT1(VEGFR1),KDR(VEGFR2),FLT4(VEGFR3)
RTK_family <- c('ALK', 'ACVR1C', 'NTRK3', 'KIT', 'EGFR', 'PDGFRA', 'TGFBR2', 'EPHB1',
                'MET', 'NTRK2', 'TIE1', 'DDR2', 'TEK', 'ACVRL1', 'EPHA5', 'ROS1',
                'FLT4', 'EPHB6', 'TYRO3', 'PDGFRB', 'EPHA3', 'EPHA4', 'LTK', 'EPHA2',
                'MERTK', 'CSF1R', 'ACVR2A', 'INSRR', 'KDR', 'AXL', 'FLT3', 'FGFR1',
                'BMPR1A', 'ACVR1', 'NTRK1', 'FLT1', 'MUSK', 'BMPR2', 'ROR2', 'EPHA7',
                'EPHB2', 'TGFBR1', 'ACVR2B', 'AMHR2', 'IGF1R', 'INSR', 'ACVR1B', 'EPHB4',
                'EPHA6', 'EPHA1', 'ERBB4', 'EPHB3', 'DDR1', 'FGFR2', 'ERBB3', 'ERBB2',
                'BMPR1B', 'MST1R', 'FGFR3', 'RET', 'EPHA10', 'FGFR4', 'EPHA8')
library(pheatmap)
for (i in 1:length(RTK_family)) {
  ids = CS2_drug2cell[RTK_family[i],]
  ids = ids@assays$RNA@data %>% as.numeric()
  assign(RTK_family[i],tapply(ids, CS2_drug2cell@meta.data$cor_type, mean))
}

matrix <- NULL
for (i in 1:length(RTK_family)) {
  matrix = rbind(matrix,get(RTK_family[i]))
}

row.names(matrix) = RTK_family

heatmap_matrix <- matrix[,colnames(matrix) != "normal"]


col <- c(colorRampPalette(c('#3ab2e8', 'white'))(36),
         colorRampPalette(c('white', '#da3446'))(36))

p1 <- pheatmap(t(heatmap_matrix), color = col,scale = 'column',cluster_cols = T,cluster_rows = T)
# breaks = unique(c(seq(-2,2, length=100))))
ggsave(p1, filename = '/public/home/scRNA/HIM-omic/HIM_single_cell/result/HIM-integrate/figure/RTK_family.pdf',
       height = 8,width = 14)















