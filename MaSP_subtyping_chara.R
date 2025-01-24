library(limma)
library(dplyr)
library(pamr)
library(fgsea)
library(tibble)
library(reshape2)

source('R/functions.R')


# Protein expression tally -----------------------------------------------------

pt <- read.csv("data/Protein_Matrix.csv",row.names = 1)#9888 360
pt$symbol <- gsub('^.*_','',rownames(pt))
pt <- pt[!duplicated(pt$symbol),]
rownames(pt) <- pt$symbol
pt <- pt[-ncol(pt)]#9534 330

de <- read.csv("data/要去掉的重复的数据V_2.csv")
pt <- pt[setdiff(colnames(pt),de$Barcode)]#9875 326
#save(pt,file='HIM_proteomic_data_v2.rda')

tcsa <- read.csv('data/TCSA_membrane_gene_list.csv')
pt_mem <- pt %>% filter(row.names(pt) %in% tcsa$HGNC.Symbol)

count = apply(pt, 2, function(x) {
  return(length(na.omit(x)))
})

count <- data.frame(count = count, label = 'count')

ggplot(count, aes(x=label, y=count)) +
  geom_violin() +
  geom_jitter() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

count <- count[order(count$count),]
count$rank <- c(1:nrow(count))

p1 <- ggplot(count, aes(x=rank, y=count)) +
  geom_point() +
  ylim(0, max(count$count)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


count_mem = apply(pt_mem, 2, function(x) {
  return(length(na.omit(x)))
})

count_mem <- data.frame(count = count_mem, label = 'count')

ggplot(count_mem, aes(x=label, y=count)) +
  geom_violin() +
  geom_jitter() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

count_mem <- count_mem[order(count_mem$count),]
count_mem$rank <- c(1:nrow(count_mem))

p2 <- ggplot(count_mem, aes(x=rank, y=count)) +
  geom_point() +
  ylim(0, max(count_mem$count)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p1 + p2

pt_count <- data.frame(row.names = rownames(pt))

pt_count$count = apply(pt, 1, function(x) {
  return(length(na.omit(x))/312)
})

pt_count$exp = apply(pt, 1, function(x) {
  return(mean(na.omit(x)))
})

pt_count$label <- rownames(pt_count)

ggplot(pt_count, aes(x=exp, y=count)) +
  geom_point() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text_repel(aes(label = label), 
                  max.overlaps = 50,
                  min.segment.length = 1,
                  force = 3, force_pull = 3)

#-------------------------------------------------------------------------------

# Making heatmap of proteins and pathways --------------------------------------

library(ConsensusClusterPlus)

load('data/HIM_proteomic_data_v2.rda')

pam.metb.data <- pamr.from.excel('data/PamAlg_traindata_metabric_for_newispy2.txt', 304, sample.labels=T)
pam.metb.data <- pamr.from.excel('data/PamAlg_traindata_metabric_new3clust.txt', 304, sample.labels=T)

tcsa <- read.csv('data/TCSA_membrane_gene_list.csv')

fill <- scale(pt[which(rownames(pt) %in% tcsa$`HGNC.Symbol`),])#1128

workDir <- 'results/clust_mem'

results <- ConsensusClusterPlus(as.matrix(fill),
                                maxK = 10,
                                reps = 1000,
                                pItem = 0.8,
                                pFeature = 1,
                                title = workDir,
                                clusterAlg ="pam",
                                distance = "euclidean",
                                seed = 123456,
                                plot = "pdf")

calcICL(results, title="consensusScore", plot="pdf")

result_3c <- results[[3]][[3]]
result_3c <- result_3c[order(result_3c, decreasing=F)]

him_orig <- list(expression = pt,
                 subtype = paste0('CS', result_3c),
                 samplelabels = names(result_3c))

saveRDS(him_orig, file='HIM_orig.rds')
him_orig <- readRDS(file='HIM_orig.rds')

# Checking if results matches with previous clustering

library(caret)

df_a <- data.frame(row.names = him_orig[[3]], type_new = him_orig[[2]])
df_b <- data.frame(row.names = pam.metb.data[['samplelabels']], type_old = pam.metb.data[['y']])


df_a <- df_a %>% filter(row.names(df_a) %in% rownames(df_b))
df_a <- df_a[rownames(df_b),]
df_b <- df_b[rownames(df_b),]

confusionMatrix(as.factor(df_a), as.factor(df_b))

# Mapping clustering quality

library(factoextra)
library(NbClust)

fviz_nbclust(t(fill), pam, method = 'wss') +
  labs(subtitle = 'Elbow method')

fviz_nbclust(t(fill), pam, method = 'silhouette') +
  labs(subtitle = 'Silhoouette mothod')

fviz_nbclust(t(fill), pam, method = 'gap_stat', nstart=25, nboot=50) +
  labs(subtitle = 'Gap statistic method')

fviz_nbclust(t(fill), kmeans, method = 'wss') +
  labs(subtitle = 'Elbow method')

fviz_nbclust(t(fill), kmeans, method = 'silhouette') +
  labs(subtitle = 'Silhoouette mothod')

# Hallmark for 4 subtypes

library(GSVA)
library(pheatmap)

hallmark_label <- read.csv('data/hallmark_labelling.csv', row.name = 1)

pathways <- getGmt("D:/Data/GSEA/h.all.v7.4.symbols.gmt")

gsva_pt <- gsva(expr=as.matrix(pt), gset.idx.list=pathways, 
                method="gsva", kcdf="Gaussian")

gsva_pt <- gsva_pt[,names(results[[4]][[3]][order(results[[4]][[3]])])]
gsva_pt <- as.data.frame(gsva_pt)

rownames(gsva_pt) <- gsub('HALLMARK_', '', rownames(gsva_pt))

gsva_pt <- gsva_pt[rownames(hallmark_label),]

annotation_col <- data.frame(row.names = names(results[[4]][[3]]), subtype = results[[4]][[3]])

ann_colors <- list(
  subtype = c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6')
)

pheatmap(gsva_pt, scale = 'row',
         cluster_cols = F, cluster_rows = F,
         annotation_row = hallmark_label,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         #breaks = seq(-1.5, 1.5, length.out = 100),
         breaks = seq(-2, 2, length.out = 100),
         color=colorRampPalette(c('#b6e8ff','white','#da3446'))(100),
         gaps_row = c(8,13,19,22,32,37),
         gaps_col = c(74,225),
         show_colnames = F)

pheatmap(gsva_pt, scale = 'row',
         cluster_cols = F, cluster_rows = F,
         annotation_row = hallmark_label,
         annotation_col = annotation_col,
         #breaks = seq(-1.5, 1.5, length.out = 100),
         breaks = seq(-2, 2, length.out = 100),
         color=colorRampPalette(c('#b6e8ff','white','#da3446'))(100),
         gaps_row = c(8,13,19,22,32,37),
         gaps_col = c(71,114,244),
         show_colnames = F)

sum(results[[3]][[3]] == 2)

# Sankey plot between 3 and 4 subtypes

library(ggplot2)
library(ggsankey)

cluster_sankey_3v4 <- data.frame(cluster_3 = paste0('3-', results[[3]][[3]]),
                                 cluster_4 = paste0('4-', results[[4]][[3]]))

cluster_sankey_3v4 <- cluster_sankey_3v4 %>% make_long(cluster_3, cluster_4)

cluster_3v4_colors <- c(`3-1`='#E71D36', `3-2`='#FF9F1C', `3-3`='#2EC4B6',
                        `4-1`='#E71D36', `4-2`='#FF9F1C', `4-3`='#2EC4B6', `4-4`='grey')

# landscape 8*6
ggplot(cluster_sankey_3v4, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = .6) +
  geom_alluvial_text(size = 3, color = "white") +
  scale_fill_manual(values = cluster_3v4_colors) +
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("3 cluster vs 4 cluster")




# Finding DE proteins and corresponding pathways

library(pheatmap)
library(RColorBrewer)

pt <- pt[,him_orig[['samplelabels']]]

de_cs1 = limma_de_HIM(pt, him_orig, 'CS1')
de_cs2 = limma_de_HIM(pt, him_orig, 'CS2')
de_cs3 = limma_de_HIM(pt, him_orig, 'CS3')

de_cs1_mem = de_cs1 %>% filter(row.names(de_cs1) %in% tcsa$HGNC.Symbol)
de_cs2_mem = de_cs2 %>% filter(row.names(de_cs2) %in% tcsa$HGNC.Symbol)
de_cs3_mem = de_cs3 %>% filter(row.names(de_cs3) %in% tcsa$HGNC.Symbol)

diff_prot <- c(rownames(de_cs1 %>% filter(adj.P.Val < 0.05)),
               rownames(de_cs2 %>% filter(adj.P.Val < 0.05)),
               rownames(de_cs3 %>% filter(adj.P.Val < 0.05)))

diff_prot <- diff_prot[!duplicated(diff_prot)]

#pathways <- gmtPathways("files/h.all.v7.4.symbols.gmt")
pathways <- gmtPathways("files/c5.go.v7.5.1.symbols.gmt")
pathways <- gmtPathways("files/c2.cp.kegg.v7.5.1.symbols.gmt")

gsea_cs1 = gsea_limma(de_cs1, pathways)
gsea_cs2 = gsea_limma(de_cs2, pathways)
gsea_cs3 = gsea_limma(de_cs3, pathways)

gsea_cs1 = gsea_cs1[order(gsea_cs1$NES, decreasing=T),]
gsea_cs2 = gsea_cs2[order(gsea_cs2$NES, decreasing=T),]
gsea_cs3 = gsea_cs3[order(gsea_cs3$NES, decreasing=T),]

gsea_cs1_filtered <- gsea_cs1 %>% filter(pathway %in% (gsea_cs2 %>% filter(NES < 0))$pathway &
                                           pathway %in% (gsea_cs3 %>% filter(NES < 0))$pathway)
gsea_cs2_filtered <- gsea_cs2 %>% filter(pathway %in% (gsea_cs1 %>% filter(NES < 0))$pathway &
                                           pathway %in% (gsea_cs3 %>% filter(NES < 0))$pathway)
gsea_cs3_filtered <- gsea_cs3 %>% filter(pathway %in% (gsea_cs1 %>% filter(NES < 0))$pathway &
                                           pathway %in% (gsea_cs2 %>% filter(NES < 0))$pathway)


save(gsea_cs1_filtered,
     gsea_cs2_filtered,
     gsea_cs3_filtered,
     file='gsea_filtered.Rdata')

load('gsea_filtered.Rdata')

# fill <- as.matrix(pt %>% filter(row.names(pt) %in% diff_prot))
# 
# pheatmap(fill, cluster_cols = F, scale = 'row',
#          breaks = seq(-5, 5, length.out = 100),
#          color=colorRampPalette(c('#2255aa', 'white', '#aa5555'))(100))



#-------------------------------------------------------------------------------

# Constructing Major Heatmap ---------------------------------------------------


top_genes <- function(de_res, gsea_res, n, m=c(1,2,3)) {
  de_res <- de_res[order(de_res$adj.P.Val),]
  gsea_res <- gsea_res[order(gsea_res$NES, decreasing=T),]
  
  top <- rownames(de_res %>% filter(row.names(de_res) %in% gsea_res$leadingEdge[[m[1]]]))[1:n]
  de_res <- de_res %>% filter(!row.names(de_res) %in% top)
  
  top <- c(top, rownames(de_res %>% filter(row.names(de_res) %in% gsea_res$leadingEdge[[m[2]]]))[1:n])
  de_res <- de_res %>% filter(!row.names(de_res) %in% top)
  
  top <- c(top, rownames(de_res %>% filter(row.names(de_res) %in% gsea_res$leadingEdge[[m[3]]]))[1:n])
  
  return(top)
}

n = 5
m3 = c(1,10, 17)

top_cs1 <- top_genes(de_cs1 %>% filter(row.names(de_cs1) %in% tcsa$HGNC.Symbol), gsea_cs1_filtered, n)
top_cs2 <- top_genes(de_cs2 %>% filter(row.names(de_cs2) %in% tcsa$HGNC.Symbol), gsea_cs2_filtered, n)
top_cs3 <- top_genes(de_cs3 %>% filter(row.names(de_cs3) %in% tcsa$HGNC.Symbol), gsea_cs3_filtered, n, m3)

top_cs1 <- top_genes(de_cs1, gsea_cs1_filtered, n)
top_cs2 <- top_genes(de_cs2, gsea_cs2_filtered, n)
top_cs3 <- top_genes(de_cs3, gsea_cs3_filtered, n, m3)

top_prots <- c(top_cs1, top_cs2, top_cs3)
top_paths <- c(gsea_cs1_filtered$pathway[1:3],
               gsea_cs2_filtered$pathway[1:3],
               gsea_cs3_filtered$pathway[m3])
  
sum(duplicated(c(top_cs1, top_cs2, top_cs3)))

# Constructing Fig 1 heatmap #### 

# Distinctly expressed membrane protein part

# pick_sig_proteins = function(de_res) {
#   de_res <- de_res %>% filter(logFC > 0)
#   df <- rbind(de_res[order(de_res$logFC, decreasing=T),][1:6,],
#               de_res[order(de_res$adj.P.Val),][1:6,])
#   df <- df[!duplicated(df),]
#   return(df)
# }

pick_sig_proteins = function(de_res) {
  de_res <- de_res %>% filter(adj.P.Val < 0.05)
  de_res <- de_res %>% filter(logFC > 0.35)
  df <- de_res[order(de_res$t, decreasing=T),][1:10,]
  return(df)
}

hm_exp <- rbind(pt[rownames(pt) %in% rownames(pick_sig_proteins(de_cs1_mem)),], #10
                pt[rownames(pt) %in% rownames(pick_sig_proteins(de_cs2_mem)),], #10
                pt[rownames(pt) %in% rownames(pick_sig_proteins(de_cs3_mem)),]) #9

pheatmap(hm_exp, cluster_cols = F, cluster_rows = F, scale = 'row',
         breaks = seq(-3, 3, length.out = 100),
         color=colorRampPalette(c('#2255aa', 'white', '#aa5555'))(100),
         gaps_col = c(74,225),
         show_colnames = F)

# GSVA GO GSEA Construction

library(GSVA)
library(GSEABase)


{
  # 12
  list_diff_cs1 <- list('T Cell Related' = c(2,5,14),
                        'Antigen Presentation' = c(1,8,9,10,12,15),
                        'Viral Response' = c(6,11,53))
  # 18
  list_diff_cs2 <- list('ECM' = c(1,2,3,4,5,6,9,10,14,15,17),
                        'Tissue Development' = c(64,95,371,328,232),
                        'Clotting' = c(7,8))
  # 16
  list_diff_cs3 <- list('Cell Development' = c(10,27),
                        'Transcription Related' = c(7,8,12,28,29),
                        'Mitochondrial Related' = c(1,3,9,11,14),
                        'Metabolic Processes' = c(2,5,6,15))
} # Defining GO pathways for each subtype

gsea_diff_list <- c(gsea_cs1_filtered$pathway[unlist(list_diff_cs1)],
                    gsea_cs2_filtered$pathway[unlist(list_diff_cs2)],
                    gsea_cs3_filtered$pathway[unlist(list_diff_cs3)])

{
  # 21
  list_diff_cs1 <- list('MHC Antigen Presentation' = c('GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_ANTIGEN',
                                                       'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN',
                                                       'GOBP_POSITIVE_REGULATION_OF_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY',
                                                       'GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_PROTEIN_COMPLEX',
                                                       'GOMF_MHC_CLASS_II_PROTEIN_COMPLEX_BINDING',
                                                       'GOBP_REGULATION_OF_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY',
                                                       'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I'),
                        'T-Cell Activation' = c('GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY',
                                                'GOCC_T_CELL_RECEPTOR_COMPLEX',
                                                'GOBP_POSITIVE_REGULATION_OF_T_CELL_RECEPTOR_SIGNALING_PATHWAY',
                                                'GOBP_POSITIVE_REGULATION_OF_ALPHA_BETA_T_CELL_DIFFERENTIATION',
                                                'GOBP_T_CELL_DIFFERENTIATION_IN_THYMUS',
                                                'GOBP_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL',
                                                'GOBP_T_CELL_LINEAGE_COMMITMENT',
                                                'GOBP_GAMMA_DELTA_T_CELL_DIFFERENTIATION'),
                        'Adaptive Immune Response' = c('GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY',
                                                       'GOBP_B_CELL_DIFFERENTIATION',
                                                       'GOBP_DENDRITIC_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION',
                                                       'GOBP_POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY',
                                                       'GOBP_POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY',
                                                       'GOBP_B_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE'))
  # 30
  list_diff_cs2 <- list('ECM' = c('GOBP_ELASTIC_FIBER_ASSEMBLY',
                                  'GOBP_COLLAGEN_FIBRIL_ORGANIZATION',
                                  'GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY',
                                  'GOCC_COMPLEX_OF_COLLAGEN_TRIMERS',
                                  'GOCC_BASEMENT_MEMBRANE',
                                  'GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT_CONFERRING_TENSILE_STRENGTH',
                                  'GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION',
                                  'GOCC_COLLAGEN_TRIMER',
                                  'GOCC_EXTERNAL_ENCAPSULATING_STRUCTURE',
                                  'GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX',
                                  'GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT'),
                        'Humoral Immune' = c('GOBP_COMPLEMENT_ACTIVATION',
                                             'GOCC_IMMUNOGLOBULIN_COMPLEX',
                                             'GOBP_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_CIRCULATING_IMMUNOGLOBULIN',
                                             'GOBP_HUMORAL_IMMUNE_RESPONSE',
                                             'GOCC_IMMUNOGLOBULIN_COMPLEX_CIRCULATING',
                                             'GOMF_IMMUNOGLOBULIN_RECEPTOR_BINDING',
                                             'GOBP_ACUTE_INFLAMMATORY_RESPONSE'),
                        'Angiogenesis' = c('GOCC_BLOOD_MICROPARTICLE',
                                           'GOBP_BLOOD_VESSEL_MORPHOGENESIS',
                                           'GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION',
                                           'GOBP_POSITIVE_REGULATION_OF_VASCULATURE_DEVELOPMENT',
                                           'GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT',
                                           'GOBP_VASCULATURE_DEVELOPMENT',
                                           'GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_STIMULUS',
                                           'GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_GROWTH_FACTOR_STIMULUS',
                                           'GOBP_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_SIGNALING_PATHWAY',
                                           'GOBP_VASOCONSTRICTION',
                                           'GOBP_REGULATION_OF_VASOCONSTRICTION',
                                           'GOBP_VASCULOGENESIS'))
  # 30
  list_diff_cs3 <- list('Metabolic Processes' = c('GOBP_MITOCHONDRIAL_MEMBRANE_ORGANIZATION',
                                                  'GOCC_MITOCHONDRIAL_ENVELOPE',
                                                  'GOBP_INNER_MITOCHONDRIAL_MEMBRANE_ORGANIZATION',
                                                  'GOCC_MITOCHONDRIAL_MATRIX',
                                                  'GOBP_REGULATION_OF_MITOCHONDRIAL_GENE_EXPRESSION',
                                                  'GOBP_PYRIMIDINE_NUCLEOSIDE_TRIPHOSPHATE_METABOLIC_PROCESS',
                                                  'GOBP_PYRIMIDINE_RIBONUCLEOTIDE_BIOSYNTHETIC_PROCESS',
                                                  'GOBP_PYRIMIDINE_RIBONUCLEOSIDE_TRIPHOSPHATE_BIOSYNTHETIC_PROCESS',
                                                  'GOBP_CELLULAR_AMINO_ACID_METABOLIC_PROCESS'),
                        'Cell Proliferation' = c('GOBP_MAMMARY_GLAND_ALVEOLUS_DEVELOPMENT',
                                                 'GOBP_REGULATION_OF_NEUROBLAST_PROLIFERATION',
                                                 'GOBP_POSITIVE_REGULATION_OF_NEUROBLAST_PROLIFERATION',
                                                 'GOBP_MEIOSIS_I_CELL_CYCLE_PROCESS',
                                                 'GOBP_INNER_CELL_MASS_CELL_PROLIFERATION',
                                                 'GOBP_POSITIVE_REGULATION_OF_MAMMARY_GLAND_EPITHELIAL_CELL_PROLIFERATION',
                                                 'GOBP_POSITIVE_REGULATION_OF_G0_TO_G1_TRANSITION',
                                                 'GOBP_CARDIOBLAST_PROLIFERATION',
                                                 'GOBP_MESENCHYMAL_CELL_PROLIFERATION_INVOLVED_IN_LUNG_DEVELOPMENT',
                                                 'GOBP_POSITIVE_REGULATION_OF_NEURAL_PRECURSOR_CELL_PROLIFERATION',
                                                 'GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION'),
                        'Cell Adhesion' = c('GOBP_CARDIAC_MUSCLE_CELL_CARDIAC_MUSCLE_CELL_ADHESION',
                                            'GOMF_CADHERIN_BINDING_INVOLVED_IN_CELL_CELL_ADHESION',
                                            'GOMF_PROTEIN_BINDING_INVOLVED_IN_HETEROTYPIC_CELL_CELL_ADHESION'),
                        'Stemness' = c('GOMF_WNT_ACTIVATED_RECEPTOR_ACTIVITY',
                                       'GOMF_CORECEPTOR_ACTIVITY_INVOLVED_IN_WNT_SIGNALING_PATHWAY_PLANAR_CELL_POLARITY_PATHWAY',
                                       'GOBP_CANONICAL_WNT_SIGNALING_PATHWAY_INVOLVED_IN_REGULATION_OF_CELL_PROLIFERATION',
                                       'GOBP_WNT_PROTEIN_SECRETION',
                                       'GOBP_REGULATION_OF_STEM_CELL_DIFFERENTIATION',
                                       'GOBP_POSITIVE_REGULATION_OF_STEM_CELL_DIFFERENTIATION',
                                       'GOBP_NEGATIVE_REGULATION_OF_STEM_CELL_PROLIFERATION'))
} # Defining GO pathways for each subtype 2

gsea_diff_list <- c(unlist(list_diff_cs1), unlist(list_diff_cs2), unlist(list_diff_cs3))

# gsea_diff_list <- c(gsea_cs1_filtered$pathway[1:10],
#                     gsea_cs2_filtered$pathway[1:10],
#                     gsea_cs3_filtered$pathway[1:10])


pathways <- getGmt("files/c5.go.v7.5.1.symbols.gmt")
# pathways <- getGmt("files/c8.all.v2023.2.Hs.symbols.gmt")


gsva_pt <- gsva(expr=as.matrix(pt), gset.idx.list=pathways, 
                method="gsva", kcdf="Gaussian")

hm_gsva_go <- gsva_pt[rownames(gsva_pt) %in% gsea_diff_list,]
gsea_diff_list <- gsea_diff_list[gsea_diff_list %in% rownames(hm_gsva_go)]
hm_gsva_go <- hm_gsva_go[gsea_diff_list,]

pheatmap(hm_gsva_go, scale = 'row',
         gaps_row = c(21,51), 
         gaps_col = c(sum(him_orig[['subtype']] %in% 'CS1'), 
                      sum(him_orig[['subtype']] %in% c('CS1', 'CS2'))),
         cluster_rows = F, cluster_cols = F,
         show_colnames = F,
         breaks = seq(-3, 3, length.out = 100),
         color=colorRampPalette(c('#b6e8ff','white','#da3446'))(100))

# Label proteins as pathways

prot_sig_cs1 <- rownames(pick_sig_proteins(de_cs1_mem))
prot_go_label_cs1 <- list()
for (i in prot_sig_cs1) {
  prot_go_label_cs1[[i]] <- gsub('0|1|2|3|4|5|6|7|8|9', '', names(unlist(pathways)[which(unlist(pathways) %in% i)]))
  prot_go_label_cs1[[i]] <- gsub('0|1|2|3|4|5|6|7|8|9', '', names(unlist(list_diff_cs1)[which(unlist(list_diff_cs1) %in% prot_go_label_cs1[[i]])]))
  prot_go_label_cs1[[i]] <- prot_go_label_cs1[[i]][!duplicated(prot_go_label_cs1[[i]])]
}

prot_sig_cs2 <- rownames(pick_sig_proteins(de_cs2_mem))
prot_go_label_cs2 <- list()
for (i in prot_sig_cs2) {
  prot_go_label_cs2[[i]] <- gsub('0|1|2|3|4|5|6|7|8|9', '', names(unlist(pathways)[which(unlist(pathways) %in% i)]))
  prot_go_label_cs2[[i]] <- gsub('0|1|2|3|4|5|6|7|8|9', '', names(unlist(list_diff_cs2)[which(unlist(list_diff_cs2) %in% prot_go_label_cs2[[i]])]))
  prot_go_label_cs2[[i]] <- prot_go_label_cs2[[i]][!duplicated(prot_go_label_cs2[[i]])]
}

prot_sig_cs3 <- rownames(pick_sig_proteins(de_cs3_mem))
prot_go_label_cs3 <- list()
for (i in prot_sig_cs3) {
  prot_go_label_cs3[[i]] <- gsub('0|1|2|3|4|5|6|7|8|9', '', names(unlist(pathways)[which(unlist(pathways) %in% i)]))
  prot_go_label_cs3[[i]] <- gsub('0|1|2|3|4|5|6|7|8|9', '', names(unlist(list_diff_cs3)[which(unlist(list_diff_cs3) %in% prot_go_label_cs3[[i]])]))
  prot_go_label_cs3[[i]] <- prot_go_label_cs3[[i]][!duplicated(prot_go_label_cs3[[i]])]
}

# Rearranging protein orders
hm_exp <- hm_exp[c("HLA-DRA","PTPRC","B2M","CYBB","LCP1","ITGAL","CD38","ICAM3","LSP1","CLCNKB",
                   "CAV1","LRP1","CD248","COL6A2","PDGFRB","PALM","CD99","ADD1","ANK2","EHD2",
                   "EPCAM","TACSTD2","ABCC5","ECEL1","RAC3","SLC7A5","PTK6","ABCB5","PLXNB3"),]


# Adding estimate score

hm_esti <- read.csv('files/estimate_sore.csv', row.names = 1)

hm_esti <- as.data.frame(t(hm_esti))

hm_esti <- hm_esti[,him_orig$samplelabels]
hm_esti <- hm_esti[c('ImmuneScore', 'StromalScore', 'TumorPurity'),]

# Adding 4 subtypes score

result_4c <- results[[4]][[3]]
result_4c <- result_4c[him_orig[['samplelabels']]]


# Adding clinical labels


info <- read.csv("data/addmessage_PamrRes_3clust_proteo.csv", row.names = 1)

info$lymph.node.status <- as.character(info$lymph.node.status)
info$KI67 <- as.character(info$KI67)
info$tumor.grade <- as.character(info$tumor.grade)

info$lymph.node.status[is.na(info$lymph.node.status)] <- 'na'

info$tumor.grade[is.na(info$tumor.grade)] <- 'na'

info$KI67[is.na(info$KI67)] <- 'na'
info$KI67[info$KI67 == '1'] <- '<30%'
info$KI67[info$KI67 == '2'] <- '>=30%'




library(pheatmap)
library(enrichplot)
library(fgsea)


hm_full <- rbind(hm_esti, hm_gsva_go, hm_exp)


annotations <- data.frame(row.names = him_orig[['samplelabels']],
                          subtype = him_orig[['subtype']])

ann_colors <- list(
  subtype = c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6')
)

annotations <- merge(annotations, info[,-c(2:5)], by='row.names')
rownames(annotations) <- annotations$Row.names
annotations <- annotations[,-c(1,4)]
annotations$age <- ifelse(annotations$age <= 40,"<= 40",
                             ifelse(annotations$age > 40 & annotations$age < 55,"40~55",
                                    ifelse(annotations$age >= 55,">= 55","reset")))


ann_colors <- list(
  subtype = c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6'),
  subtype4 = c(CS1='#E71D36', CS2='#2c750b', CS3='#FF9F1C', CS4='#2EC4B6'),
  age = c("<= 40" = "#d6a2e0","40~55" = "#c76cd9",">= 55" = "#b008d1"),
  tumor.grade = c("na" = "#ffffff","1" = "#d7ab7a","2" = "#d1ad4a","3" = "#cd7229"),
  KI67 = c("na" =  "#ffffff","<30%" = "#d18899",">=30%" = "#b02145"),
  lymph.node.status = c("na" =  "#ffffff","0" = "#888888","1" = "black")
)


pheatmap(hm_full, scale = 'row',
         gaps_row = c(3,13,24,36,46,56), 
         gaps_col = c(sum(him_orig[['subtype']] %in% 'CS1'), 
                      sum(him_orig[['subtype']] %in% c('CS1', 'CS2'))),
         cluster_rows = F, cluster_cols = F,
         show_colnames = F,
         breaks = seq(-3, 3, length.out = 100),
         color=colorRampPalette(c('#2255aa', 'white', '#aa5555'))(100))


pheatmap(hm_full, scale = 'row',
         gaps_row = c(3,24,54,84,94,104), 
         gaps_col = c(sum(him_orig[['subtype']] %in% 'CS1'), 
                      sum(him_orig[['subtype']] %in% c('CS1', 'CS2'))),
         annotation_col = annotations,
         annotation_colors = ann_colors,
         cluster_rows = F, cluster_cols = F,
         show_colnames = F,
         breaks = seq(-3, 3, length.out = 100),
         color=colorRampPalette(c('#b6e8ff','white','#da3446'))(100))

#ggsave('export/fig1_hm.pdf', width=17, height=15)

#####

