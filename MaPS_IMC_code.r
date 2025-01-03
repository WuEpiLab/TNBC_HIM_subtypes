'''
steinbock包做细胞分割、获取每张图片的细胞的蛋白表达矩阵：参考https://bodenmillergroup.github.io/steinbock/latest/
cycombine进行批次校正：参考https://biosurf.org/cyCombine_ref_manual.html
'''
library(imcRtools)
library(cytomapper)
library(tidyverse)
library(RColorBrewer)
library(scuttle)
library(patchwork)
library(dittoSeq)
library(viridis)
library(data.table)
library(kohonen)
library(ConsensusClusterPlus)
library(bluster)
library(dplyr)
library(paletteer)
library(dittoSeq)
library(forcats)
library(ggplot2)
library(scales)
library(ComplexHeatmap)

set.seed(220225)
mypath<-"./work_data/steinbock/intensities/"
spe <- read_steinbock("work_data/steinbock/")
spe


#get pateint_id
spe$patient_id <- str_extract(spe$sample_id, "CS[1-3]_HIM[1-9][1-9][1-9]")
spe$ROI <- str_extract(spe$sample_id, "00[1-8]")

#arcsin transform counts with cofactor 1
assay(spe, "exprs") <- asinh(counts(spe)/1)

#select marker
rowData(spe)$use_channel <- !grepl("DNA|Histone|PDL1", rownames(spe))

color_vectors <- list()

ROI <- setNames(brewer.pal(length(unique(spe$ROI)), name = "BrBG"), 
                unique(spe$ROI))
patient_id <- setNames(brewer.pal(length(unique(spe$patient_id)), name = "Paired"), 
                       unique(spe$patient_id))
sample_id <- setNames(c(brewer.pal(6, "YlOrRd")[3:6],
                        brewer.pal(6, "PuBu")[3:6],
                        brewer.pal(6, "YlGn")[3:6],
                        brewer.pal(6, "BuPu")[3:6],
                        brewer.pal(6,"Reds")[3:6],
                        brewer.pal(6,"Set1")[3:6],
                        brewer.pal(6,"PuOr")[3:6],
                        brewer.pal(6,"PiYG")[3:6],
                        brewer.pal(6,"BrBG")[3:6],
                        brewer.pal(6,"PRGn")[3:6],
                        brewer.pal(6,"RdGy")[3:6],
                        brewer.pal(6,"RdBu")[3:6]
                        ),unique(spe$sample_id))
color_vectors$ROI <- ROI
color_vectors$patient_id <- patient_id
color_vectors$sample_id <- sample_id
metadata(spe)$color_vectors <- color_vectors

#load Images
images <- loadImages("work_data/steinbock/img/")
#load Image mask
masks <- loadImages("work_data/steinbock/masks/", as.is = TRUE)
channelNames(images) <- rownames(spe)
images

patient_id <- rep(unique(spe$patient_id),each=4)
mcols(images) <- mcols(masks) <- DataFrame(sample_id = unique(spe$sample_id),
                                           patient_id = patient_id)
										   

clip99 <- function(x){ 
  p99 <- quantile(x, probs = 0.99, na.rm = T)
  print(p99)
  x[x > p99 & !is.na(x)] <- p99
  return(x)
}
exprs <- assay(spe,'exprs')
exprs1 <- as.data.table(t(exprs))

exprs1[, eval(rownames(spe)) := 
         lapply(.SD, clip99), .SDcols = rownames(spe)]
assay(spe,'exprs') <- as.matrix(t(exprs1))

#umap、tsne before batch correction
spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs") 
spe <- runTSNE(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs")

#umap、tsne after cycombine correction
spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs_corrected") 
spe <- runTSNE(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs_corrected")

#marker_cluster heatmap
#get exprs
exprs <- assay(spe,'exprs')
library(data.table)
exprs1 <- as.data.table(t(exprs))
#增加一列
exprs1[,cluster:=spe$celltype]
clusters <- 'cluster'
get_hmap_data <- function(zscore = F) {
  channels <- rownames(spe)[rowData(spe)$use_channel]
  bygrp <- clusters
  toplot <- exprs1[,lapply(.SD, median), .SDcols = channels, by = bygrp]
  newcols <- toplot[, get(bygrp)]
  toplot <- t(toplot[, -c(bygrp), with = F])
  if(zscore) toplot <- t(apply(toplot, 1, scale_clip))
  colnames(toplot) <- newcols
  return(toplot)
}
zscored_dat <- get_hmap_data(zscore=T)
zscored_dat <- as.data.frame(zscored_dat)
row_order <- c("Cancer cells","CD11b+ Cancer cells","CD73+ Cancer cells","CD57+ Cancer cells","Ki-67+/CD73+ Cancer cells","αSMA+ Cancer cells","EGFR+ Cancer cells","Ki-67+ Cancer cells","Ki-67+ Epi","CD38+ Epi","EGFR+ Epi",
               "Stromal cells","Endothelial cells","CAFs","myCAFs","vCAFs","iCAFs","M1-like Mø","M2a-like Mø","M2c-like Mø","Ki-67+ M2-like Mø","B cells",
               "Ki-67+ B cells","NK cells","CD4+ T cells","Memory CD4+ T cells","CD8+ T cells","Memory CD8+ T cells","Treg cells",
              "CD38+ immune suppressive lymphocytes","Undefined T cells"
                           )
row_levels <- factor(row_order,levels=c("Cancer cells","CD11b+ Cancer cells","CD73+ Cancer cells","CD57+ Cancer cells","Ki-67+/CD73+ Cancer cells","αSMA+ Cancer cells","EGFR+ Cancer cells","Ki-67+ Cancer cells","Ki-67+ Epi","CD38+ Epi","EGFR+ Epi",
                                        "Stromal cells","Endothelial cells","CAFs","myCAFs","vCAFs","iCAFs","M1-like Mø","M2a-like Mø","M2c-like Mø","Ki-67+ M2-like Mø","B cells",
                                        "Ki-67+ B cells","NK cells","CD4+ T cells","Memory CD4+ T cells","CD8+ T cells","Memory CD8+ T cells","Treg cells",
                                        "CD38+ immune suppressive lymphocytes","Undefined T cells"
))

cell_count <- read.table('./work_data/cell_count.csv',header = TRUE,sep=',')
split_anno = factor(c(rep('Tumor', each = 11), rep('Stroma', each = 6), rep('Immune', each = 14)),levels = c("Tumor","Stroma","Immune"))
cp<-rowAnnotation(type = anno_block(
  gp = gpar(fill = rev(metadata(spe)$color_vectors$HIM)),
  labels = c("n = 140,500 (35.42%)", "n = 77,806 (19.61%)", "n = 178,366 (44.97%)"), 
  labels_gp = gpar(col = "white", fontsize = 10), height = unit(5,"cm")),
  CellTypes=row_levels, col=list(CellTypes=metadata(spe)$color_vectors$celltype),
  "Cell Count"=anno_barplot(cell_count$Freq,gp = gpar(fill=metadata(spe)$color_vectors$celltype),border = F,width = unit(8,"cm"))
)
zscored_dat <- as.data.frame(zscored_dat)
pdf(file="./work_data/分型文章图/cell_marker_heatmap.pdf", width=100, height=100)
Heatmap(t(zscored_dat[row_order]),name="z_score",cluster_rows = F,cluster_columns = T,row_split = split_anno,row_title_gp=gpar(fontsize=15),right_annotation = cp,show_row_names = F)
dev.off()

#spe$cluster:clustering 
spe$celltype <- recode(spe$cluster,
                           "1" = "Ki-67+ B cells",
                           "2" = "B cells",
                           "3" = "EGFR+ Epi",
                           "4" = "NK cells",
                           "5" = "NK cells",
                           "6" = "CD11b+ Cancer cells",
                           "7" = "CD73+ Cancer cells",
                           "8" = "αSMA+ Cancer cells",
                           "9" = "CD4+ T cells",
                           "10" = "Undefined T cells",
                           "11" = "B cells",
                           "12" = "M2a-like Mø",
                           "13" = "EGFR+ Cancer cells",
                           "14" = "Ki-67+/CD73+ Cancer cells",
                           "15" = "CD57+ Cancer cells",
                           "16" = "CD73+ Cancer cells",
                           "17" = "CD4+ T cells",
                           "18" = "CD8+ T cells",
                           "19" = "B cells",
                           "20" = "Ki-67+ M2-like Mø",
                           "21" = "Ki-67+ Cancer cells",
                           "22" = "Ki-67+ Cancer cells",
                           "23" = "αSMA+ Cancer cells",
                           "24" = "Cancer cells",
                           "25" = "CD8+ T cells",
                           "26" = "Treg cells",
                           "27" = "CD8+ T cells",
                           "28" = "CD57+ Cancer cells",
                           "29" = "Ki-67+ M2-like Mø",
                           "30" = "Ki-67+ Cancer cells",
                           "31" = "Ki-67+ Cancer cells",
                           "32" = "Cancer cells",
                           "33" = "M2c-like Mø",
                           "34" = "Memory CD8+ T cells",
                           "35" = "Memory CD4+ T cells",
                           "36" = "CD8+ T cells",
                           "37" = "CD38+ Epi",
                           "38" = "Ki-67+ Epi",
                           "39" = "Cancer cells",
                           "40" = "Cancer cells",
                           "41" = "M2a-like Mø",
                           "42" = "M1-like Mø",
                           "43" = "M1-like Mø",
                           "44" = "CD38+ immune suppressive lymphocytes",
                           "45" = "CD38+ Epi",
                           "46" = "CD11b+ Cancer cells",
                           "47" = "CD38+ immune suppressive lymphocytes",
                           "48" = "Cancer cells",
                           "49" = "M2c-like Mø",
                           "50" = "Endothelial cells",
                           "51" = "Memory CD4+ T cells",
                           "52" = "iCAFs",
                           "53" = "B cells",
                           "54" = "myCAFs",
                           "55" = "myCAFs",
                           "56" = "Cancer cells",
                           "57" = "CAFs",
                           "58" = "Endothelial cells",
                           "59" = "vCAFs",
                           "60" = "myCAFs",
                           "61" = "Stromal cells",
                           "62" = "Stromal cells",
                           "63" = "M2c-like Mø",
                           "64" = "iCAFs"
)

metadata(spe)$color_vectors$celltype <- setNames(c("#95e6df","#85e2da","#75dfd5","#65dbd0","#56d7cb","#46d4c6","#36d0c1","#2ec4b6","#2ab4a7","#27a498","#23948a",
                                                   "#ffb857","#ffb043","#ffa730","#ff9f1c","#ff9708","#f48d00",
                                                   "#f28795","#f07685","#ef6475","#ed5265","#eb4056","#e92f46","#e71d36","#dd273d","#e71d36","#f1132f","#fa0a28",
                                                   "#d9172f","#c8152b","#b61327"),row_order)


												   
#celltype frequency barplot												   
spe$celltype <- fct_relevel(spe$celltype,row_order)
spe$HIM <- as.factor(spe$HIM)
metadata(spe)$color_vectors$HIM <- setNames(c('#E71D36','#FF9F1C','#2EC4B6'),unique(spe$HIM))
dittoBarPlot(spe, var = "celltype",var.labels.reorder = c(3,4,9,8,19,31,11,16,17,5,12,27,13,2,25,30,14,20,21,22,18,1,15,26,7,23,10,24,28,6,29), group.by = "HIM") +
  scale_fill_manual(values = metadata(spe)$color_vectors$celltype)


#celltype interaction analysis
out <- testInteractions(spe, 
                        group_by = "sample_id",
                        label = "celltype", 
                        method= "classic",
                        colPairName = "neighborhood",
                        BPPARAM = BiocParallel::SerialParam(RNGseed = 123)
                        )
out$from_label <- as.factor(out$from_label)
out$to_label <- as.factor(out$to_label)
out$from_label <- fct_relevel(out$from_label,row_order)
out$to_label <- fct_relevel(out$to_label,row_order)
out %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval))+ggtitle("CS1")+
  scale_fill_gradient2(limits=c(-16,16),breaks=c(-15,-10,-5,0,5,10,15),low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title=element_text(hjust=0.5))
 
#叠加三个分型的互作信息
heatmap_HIM <- read.csv("./work_data/Heatmap_classic1.csv",header=TRUE,check.names = FALSE)
rownames(heatmap_HIM) = heatmap_HIM[,1]
heatmap_HIM=heatmap_HIM[,-1]
colnames(heatmap_HIM) <- row_order
celltype_to_color <- setNames(c("#95e6df","#85e2da","#75dfd5","#65dbd0","#56d7cb","#46d4c6","#36d0c1","#2ec4b6","#2ab4a7","#27a498","#23948a",
           "#ffb857","#ffb043","#ffa730","#ff9f1c","#ff9708","#f48d00",
           "#f28795","#f07685","#ef6475","#ed5265","#eb4056","#e92f46","#e71d36","#dd273d","#e71d36","#f1132f","#fa0a28",
           "#d9172f","#c8152b","#b61327"),row_order)

ha=HeatmapAnnotation(celltype = row_order, col = list(celltype=celltype_to_color),simple_anno_size=unit(0.3,"cm"),show_annotation_name = FALSE,show_legend = FALSE)
f1 = colorRamp2(seq(min(heatmap_HIM), max(heatmap_HIM), length = 3), c("#2166AC", "white", "#B2182B"), space = "RGB")
Heatmap(as.matrix(heatmap_HIM),col=f1,cluster_rows = F,cluster_columns = F,show_column_names = FALSE,bottom_annotation=ha,
        row_split = factor(c(rep(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"),each=3)),levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31")),row_title_rot = 0,row_title = row_order,row_gap = unit(0, "mm"), column_gap = unit(0.28, "mm"), border = TRUE,show_row_names = F,
        column_split = factor(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"),levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31")),column_title = NULL,width = unit(10, "cm"), height = unit(12, "cm"),row_title_gp = gpar(col=celltype_to_color,fontsize=10),column_names_gp = gpar(fontsize=10),heatmap_legend_param = list(
  title = NULL, at=c(-16,16),
  labels = c("avoidance", "interaction")
))


#neighborhood analysis
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "knn", k = 10)
spe <- aggregateNeighbors(spe, 
                          colPairName = "knn_interaction_graph", 
                          aggregate_by = "metadata", 
                          count_by = "celltype")
cn_1 <- kmeans(spe$aggregatedNeighbors, centers = 20)
spe$cn_celltypes <- as.factor(cn_1$cluster)
for_plot <- prop.table(table(spe$cn_celltypes, spe$celltype33), 
                       margin = 1)
#cn_ct heatmap
pheatmap(for_plot, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column",cluster_rows = FALSE,cluster_cols = FALSE,name = "score")
		 
cn_names<- recode(spe$cn_celltypes,
                         "1"="CD8+T/Mø/CAFs cells",
                         "2"="Mø/Epi/T cells",
                         "3"="Epi/Mø/T cells",
                         "4"="Ki67+/CD73+ Cancer cells",
                         "5"="B/Ki67+B/Ki-67+ M2-like Mø",
                         "6"="CD4+T/Treg/Mø cells",
                         "7"="Stromal/myCAFs/M2-like Mø cells",
                         "8"="Endothelial/CAFs cells",
                         "9"="Epi/Ki-67+ M2-like Mø",
                         "10"="T cells zone",
                         "11"="Ki67+ Cancer/Ki67+ Epi cells",
                         "12"="CD57+/Ki67+/CD73+ Cancer cells",
                         "13"="B cells zone",
                         "14"="Cancer cells zone I",
                         "15"="CAFs/Epi/CD38+ lymphocyte cells",
                         "16"="Cancer/M2a-like Mø cells",
                         "17"="CD38+ lymphocyte/iCAFs/Epi cells",
                         "18"="NK/Epi cells",
                         "19"="Cancer cells zone II",
                         "20"="Cancer/CD38+ lymphocyte cells"
)
spe$cn_names <- cn_names

#plot cn frequency barplot
metadata(spe)$color_vectors$cn_names <- setNames(c("#FFEDA0","#78C679","#e92f46","#DE77AE","#F1B6DA","#eb4056","#238443","#ed5265","#23948a","#f48d00","#41AB5D","#238443","#FED976","#FC4E2A","#88419D","#27a498","#ed5265","#2ab4a7","#2ec4b6","#36d0c1"),unique(spe$cn_names))

metadata(spe)$color_vectors$cn_names <- setNames(c("#36d0c1","#2ec4b6","#2ab4a7","#27a498","#23948a","#78C679","#41AB5D","#28984d","#238443","#f48d00","#FFEDA0","#FED976","#e92f46","#eb4056","#ed5265","#ef6475","#F1B6DA","#FF69B4","#DE77AE","#88419D"),levels(spe$cn_names))
dittoBarPlot(spe, var = "cn_names",var.labels.reorder = c(4,5,10,15,16,1,2,9,20,12,6,7,13,14,19,3,8,17,11,18),group.by = "HIM") +
  scale_fill_manual(values = metadata(spe)$color_vectors$cn_names)
  



