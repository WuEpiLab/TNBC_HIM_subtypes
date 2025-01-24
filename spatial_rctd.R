library(spacexr)
library(Seurat)
library(ggplot2)

sprna_134 <- read.VisiumSpatialRNA('data/spatial/scanpy_format/HIM134/')
sprna_135 <- read.VisiumSpatialRNA('data/spatial/scanpy_format/HIM-LEN003-LXY/')
sprna_263 <- read.VisiumSpatialRNA('data/spatial/scanpy_format/HIM263/')

sprna_list <- list(sprna_134, sprna_135, sprna_263)
name_list <- c('134', '135', '263')

him_seurat <- readRDS('data/HIM_integrate_326.rds')

counts <- GetAssayData(him_seurat, slot='counts')
cell_types <- him_seurat$cor_type
reference <- Reference(counts, as.factor(cell_types))


rctd_list <- lapply(sprna_list, function(x, reference) {
  return(create.RCTD(x, reference, max_cores = 4))
}, reference=reference)

rctd_res_list <- lapply(rctd_list, function(x) {
  return(run.RCTD(x, doublet_mode = 'doublet'))
})

# save(rctd_res_list, file='saves/rctd_doublet.rds')
load('saves/rctd_doublet.rds')

dir.create('export/RCTD_Plots/doublet/')

# Running RCTD and processing results

for (i in c(1:length(rctd_res_list))) {
  results <- rctd_res_list[[i]]@results
  norm_weights = normalize_weights(results$weights) 
  cell_type_names <- rctd_res_list[[i]]@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- rctd_res_list[[i]]@spatialRNA
  
  resultsdir <- paste0('export/RCTD_Plots/doublet/', name_list[i]) ## you may change this to a more accessible directory on your computer.
  dir.create(resultsdir)
  
  plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
  plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
  plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                       results$results_df) 
  plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA) + theme(axis.text.x = element_text(angle=45, hjust=1))
  plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir) 
  
}


results <- myRCTD@results
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA

resultsdir <- 'RCTD_Plots' ## you may change this to a more accessible directory on your computer.
dir.create(resultsdir)


plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir) + theme_classic() + scale_colour_brewer(palette = "Set1")

# Integrating with STdeconvolve

library(STdeconvolve)

i=2

results <- rctd_res_list[[i]]@results
norm_weights = normalize_weights(results$weights) 
cell_type_names <- rctd_res_list[[i]]@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- rctd_res_list[[i]]@spatialRNA

# resultsdir <- paste0('export/RCTD_Plots/doublet/', name_list[i])

# plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir) + theme_classic() + scale_colour_discrete() + ggplot2::geom_point(ggplot2::aes(shape=19,color=class))


m <- as.matrix(norm_weights)

coords <- sprna_list[[i]]@coords
coords$x <- coords$x * 1.7

cell_order <- c('normal','CS1','CS2','CS3','Endothelial cells','Lymphatic endothelial cells','Venular endothelial cells',
                'Cycling Endothelial cells','Perivascular-like cells','imPVL','apCAFs','myCAFs','iCAFs','Cycling CAFs',
                'TAM_C1QC','TAM_ISG15','TAM_SPP1','TAM-pro','B cells','plasma cells','cDC1','cDC2','cDC3','pDCs','Mast cells',
                'Mono_CD14','Neutrophil','NK cells','CD4-ISG','CD4-Th','CD4-Tn','CD4-Tpro','CD4-Treg','CD8-ISG','CD8-Teff',
                'CD8-Tex','CD8-Tn','CD8-Tpro','T-SEPTIN')

m <- m[,cell_order]

m[,2] <- apply(m[,c(2:4)], 1, sum)
colnames(m)[2] <- 'Cancer cells'
m <- m[,-c(3,4)]

# saveRDS(m, 'him263_matrix.rds')
saveRDS(coords, 'him135_coords.rds')


plt <- vizAllTopics(theta = m,
                    pos = coords,
                    topicOrder=seq(ncol(m)),
                    topicCols=rainbow(ncol(m)),
                    groups = NA,
                    group_cols = NA,
                    r = 0.95, # size of scatterpies; adjust depending on the coordinates of the pixels
                    lwd = 0.005,
                    showLegend = TRUE,
                    plotTitle = "scatterpies")


#plt + ggplot2::guides(fill=ggplot2::guide_legend(ncol=2))


#plt
# save 20 x 12
plt + scale_fill_manual(values = c('Topic.normal' = '#cccccc',
                                   'Topic.Cancer.cells' = '#46d4c6',
                                   'Topic.Endothelial.cells' = '#ffeda8',
                                   'Topic.Lymphatic.endothelial.cells' = '#ffe689',
                                   'Topic.Venular.endothelial.cells' = '#f5d14b',
                                   'Topic.Cycling.Endothelial.cells' = '#eac129',
                                   'Topic.Perivascular.like.cells' = '#f5c47f',
                                   'Topic.imPVL' = '#efac4e',
                                   'Topic.apCAFs' = '#ffb347',
                                   'Topic.myCAFs' = '#ffa82d',
                                   'Topic.iCAFs' = '#e79d34',
                                   'Topic.Cycling.CAFs' = '#d89434',
                                   'Topic.TAM_C1QC' = '#e4c0f4',
                                   'Topic.TAM_ISG15' = '#c573e7',
                                   'Topic.TAM_SPP1' = '#b441e4',
                                   'Topic.TAM.pro' = '#a927df',
                                   'Topic.B.cells' = '#fdbadc',
                                   'Topic.plasma.cells' = '#f8a4cf',
                                   'Topic.cDC1' = '#f591c4',
                                   'Topic.cDC2' = '#f284bc',
                                   'Topic.cDC3' = '#ee76b3',
                                   'Topic.pDCs' = '#ea6bac',
                                   'Topic.Mast.cells' = '#e760a5',
                                   'Topic.Mono_CD14' = '#e2559d',
                                   'Topic.Neutrophil' = '#d84590',
                                   'Topic.NK.cells' = '#ef6475',
                                   'Topic.CD4.ISG' = '#ed5265',
                                   'Topic.CD4.Th' = '#eb4056',
                                   'Topic.CD4.Tn' = '#e92f46',
                                   'Topic.CD4.Tpro' = '#e71d36',
                                   'Topic.CD4.Treg' = '#dd273d',
                                   'Topic.CD8.ISG' = '#e71d36',
                                   'Topic.CD8.Teff' = '#f1132f',
                                   'Topic.CD8.Tex' = '#fa0a28',
                                   'Topic.CD8.Tn' = '#d9172f',
                                   'Topic.CD8.Tpro' = '#c8152b',
                                   'Topic.T.SEPTIN' = '#b61327'))


# "#95e6df","#85e2da","#75dfd5","#65dbd0","#56d7cb","#46d4c6","#36d0c1","#2ec4b6","#2ab4a7","#27a498","#23948a",
# "#ffb857","#ffb043","#ffa730","#ff9f1c","#ff9708","#f48d00",
# "#f28795","#f07685","#ef6475","#ed5265","#eb4056",
# "#e92f46","#e71d36","#dd273d","#e71d36","#f1132f","#fa0a28","#d9172f","#c8152b","#b61327"
# 
# "Cancer cells","αSMA+ Cancer cells","EGFR+ Cancer cells","Ki67+ Cancer cells","CD11b+ Cancer cells","CD73+ Cancer cells","CD57+ Cancer cells","Ki67+/CD73+ Cancer cells","ki67+ Epi","CD38+ Epi","EGFR+ Epi",
# "Stromal cells","Endothelial cells","CAFs","myCAFs","vCAFs","iCAFs",
# "M1-like Mø","M2a-like Mø","M2c-like Mø","B cells","ki67+ B cells","NK cells","CD4+ T cells",
# "Memory CD4+ T cells","CD8+ T cells","Memory CD8+ T cells","Treg cells","Activated T cells","CD38+ immune suppressive lymphocytes","Undefined T cells"

pdf('rctd_263_separated.pdf', height=7, width=9)

for (i in cell_order) {
  print(vizTopic(theta = m, pos = coords, topic = i, plotTitle = i,
           size = 2.5, stroke = 0.5, alpha = 1,
           low = "white",
           high = "red"))
}

dev.off()

vizTopic(theta = m, pos = coords, topic = "CS2", plotTitle = "X5",
         size = 2.5, stroke = 1, alpha = 0.5,
         low = "white",
         high = "red")

vizTopic(theta = m, pos = coords, topic = "CS3", plotTitle = "X5",
         size = 2.5, stroke = 1, alpha = 0.5,
         low = "white",
         high = "red")

# T-septin CAF spatial co-expression -------------------------------------------

library(PieGlyph)

celltype_tseptin <- c("T-SEPTIN")

celltype_caf <- c("apCAFs", "Cycling CAFs", "Cycling Endothelial cells", "Endothelial cells", "iCAFs", "imPVL", "Lymphatic endothelial cells", "myCAFs", "Perivascular-like cells", "Venular endothelial cells" ) 

celltype_cancer <- c("CS1", "CS2", "CS3")

celltype_other <- c("B cells", "CD4-ISG", "CD4-Th", "CD4-Tn", "CD4-Tpro", "CD4-Treg", "CD8-ISG", "CD8-Teff", "CD8-Tex", "CD8-Tn", "CD8-Tpro", "cDC1", "cDC2", "cDC3", 
                    "Mast cells", "Mono_CD14", "Neutrophil", "NK cells", "normal", "pDCs", "plasma cells", "TAM_C1QC", "TAM_ISG15", "TAM_SPP1", "TAM-pro")

m2 <- as.data.frame(cbind(m[, colnames(m) %in% celltype_tseptin],
                          apply(m[, colnames(m) %in% celltype_caf], 1, sum),
                          apply(m[, colnames(m) %in% celltype_cancer], 1, sum),
                          apply(m[, colnames(m) %in% celltype_other], 1, sum)))

colnames(m2) <- c('T-septin', 'CAFs', 'Tumor', 'Others')

m2 <- merge(coords, m2, by=0)

# landscape 12*16
ggplot(m2, aes(x=x, y=y)) +
  geom_pie_glyph(slices = c('T-septin', 'CAFs', 'Tumor', 'Others')) +
  theme_classic()

plt <- vizAllTopics(theta = m2,
                    pos = coords,
                    topicOrder=seq(ncol(m2)),
                    topicCols=rainbow(ncol(m2)),
                    groups = NA,
                    group_cols = NA,
                    r = 0.95, # size of scatterpies; adjust depending on the coordinates of the pixels
                    lwd = 0.005,
                    showLegend = TRUE,
                    plotTitle = "scatterpies")

# Visualizing distribution per category of cells ------------------------------------

library(reshape2)

zone_135 <- readRDS('him135_label.rds')
zone_263 <- readRDS('him263_label.rds')

align_rctd_zones <- function(m, rds) {
  zone <- readRDS(rds)
  
  zone <- zone[rownames(m),]
  zone_dist <- data.frame(row.names=colnames(m))
  
  for (i in unique(zone$label)) {
    zone_dist[,i] <- apply(m[zone$barcode[zone$label == i],], 2, mean)
  }
  
  return(zone_dist)
}

# save heatmaps

zone_dist_134 <- align_rctd_zones(m, 'him134_label.rds')
zone_dist_134 <- zone_dist_134[,c('Immune','Epithelial','Stromal','Boundary')]
pheatmap(zone_dist_134, cluster_rows = F, cluster_cols = F, scale = 'row') # save portrait 10 x 4.5

zone_dist_135 <- align_rctd_zones(m, 'him135_label.rds')
zone_dist_135 <- zone_dist_135[,c('Immune','Epithelial','Stromal')]
pheatmap(zone_dist_135, cluster_rows = F, cluster_cols = F, scale = 'row') # save portrait 10 x 4

zone_dist_263 <- align_rctd_zones(m, 'him263_label.rds')
zone_dist_263 <- zone_dist_263[,c('Immune','Epithelial','Stromal')]
pheatmap(zone_dist_263, cluster_rows = F, cluster_cols = F, scale = 'row')

# Bar plot

zone_dist_135$celltype <- rownames(zone_dist_135)
zone_dist_135 <- melt(zone_dist_135)

ggplot(zone_dist_135, aes(x=variable, y=value, fill=celltype)) + 
  geom_bar(position="stack", stat="identity")

pheatmap(zone_dist_135, cluster_rows = F, cluster_cols = F, scale = 'row')

# Linear plot 

library(ggchromatic)
library(ggnewscale)

ggplot(m, aes(x=`T-SEPTIN`, y=`Endothelial cells`, color=`T-SEPTIN`)) +
  geom_point() + theme_classic()

df <- merge(coords, m, by='row.names')

ggplot(df, aes(x=x, y=y)) +
  geom_point(aes(colour = rgb_spec(`T-SEPTIN`, iCAFs, myCAFs)))

ggplot(df, aes(x=x, y=y)) +
  geom_point(aes(color=`T-SEPTIN`), alpha=0.7, size=4) +
  scale_colour_gradient(low = "#eeeeee", high = "red") +
  new_scale_color() +
  geom_point(aes(color=iCAFs), alpha=0.5, size=4) +
  scale_colour_gradient(low = "#eeeeee", high = "blue") +
  theme_void()

ggplot(df, aes(x=x, y=y)) +
  geom_point(aes(color=iCAFs), alpha=1, size=4) +
  scale_colour_gradient(low = "#eeeeee", high = "blue") +
  new_scale_color() +
  geom_point(aes(color=`T-SEPTIN`), alpha=1, size=2) +
  scale_colour_gradient(low = "#eeeeee", high = "red") +
  theme_void()

ggplot(df, aes(x=x, y=y)) +
  geom_point(aes(color=`Endothelial cells`), alpha=1, size=4) +
  scale_colour_gradient(low = "#eeeeee", high = "darkgreen") +
  new_scale_color() +
  geom_point(aes(color=`T-SEPTIN`), alpha=1, size=2) +
  scale_colour_gradient(low = "#eeeeee", high = "red") +
  theme_void()

ggplot(df, aes(x=x, y=y)) +
  geom_point(aes(fill=`T-SEPTIN`, color=iCAFs), stroke=1) +
  scale_colour_gradient(low = "lightgrey", high = "red") +
  scale_fill_gradient(low = "lightgrey", high = "blue")

# Clustering with RCTD result -----------------------------------------------

weights_cs1 <- rctd_res_list[[1]]@results$weights
rownames(weights_cs1) <- paste0(rownames(weights_cs1), '-cs1')

weights_cs2 <- rctd_res_list[[2]]@results$weights
rownames(weights_cs2) <- paste0(rownames(weights_cs2), '-cs2')

weights_cs3 <- rctd_res_list[[3]]@results$weights
rownames(weights_cs3) <- paste0(rownames(weights_cs3), '-cs3')

rctd_res_matrix <- rbind(weights_cs1,
                         weights_cs2,
                         weights_cs3)


pheatmap(rctd_res_matrix, cluster_rows = F, show_rownames = F)



rctd_seurat <- CreateSeuratObject(t(rctd_res_matrix))

rctd_seurat$sample <- factor(c(rep('HIM134', nrow(weights_cs1)),
                               rep('HIM135', nrow(weights_cs2)),
                               rep('HIM263', nrow(weights_cs3))), 
                             levels=c('HIM134', 'HIM135', 'HIM263'))

rctd_seurat <- ScaleData(rctd_seurat, features = rownames(rctd_seurat))
rctd_seurat <- RunPCA(rctd_seurat, features = rownames(rctd_seurat))

rctd_seurat <- FindNeighbors(rctd_seurat, dims = 1:10)
rctd_seurat <- FindClusters(rctd_seurat, resolution = 0.05)

rctd_seurat <- RunUMAP(rctd_seurat, dims = 1:10)
DimPlot(rctd_seurat, reduction = "umap")

DimPlot(rctd_seurat, reduction = "umap", group.by='sample')


