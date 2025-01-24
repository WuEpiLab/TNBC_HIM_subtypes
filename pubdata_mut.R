library(maftools)
library(data.table)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggpubr)

source('R/functions.R')

# Data Pre-processing ----------------------------------------------------------

him_orig = readRDS('data/HIM_orig.rds')

# TCGA processing

exp_tcga = read.table('data/TCGA/HiSeqV2', sep='\t', header=T)
colnames(exp_tcga) = gsub('.', '', colnames(exp_tcga), fixed=T)

pam_tcga = tcga_pam50(exp_tcga)

write.csv(pam_tcga, 'data/TCGA_PAM50_label.csv')
pam_tcga <- read.csv('data/TCGA_PAM50_label.csv', row.names = 1)

rownames(exp_tcga) = exp_tcga[,1]
exp_tcga = exp_tcga[,-1]
exp_tcga = exp_tcga[,colnames(exp_tcga) %in% filter(pam_tcga, subtype=='Basal')$id]
exp_tcga = exp_tcga[,as.integer(substr(colnames(exp_tcga), 11, 12)) < 10]
#exp_tcga = exp_tcga[,colnames(exp_tcga) %in% unlist(lapply(maf_tcga_tnbc, function(x) x[1,16]))]

him_tcga <- subtype_casting(na.omit(exp_tcga), him_orig)

nrow(him_tcga %>% filter(clust == 'CS1')) # 37
nrow(him_tcga %>% filter(clust == 'CS2')) # 95
nrow(him_tcga %>% filter(clust == 'CS3')) # 72

write.csv(him_tcga, 'data/TCGA_HIM_label.csv')
him_tcga <- read.csv('data/TCGA_HIM_label.csv', row.names=1)

# Metabric processing

exp_meta = read.table("data/META/data_expression_median.txt", sep = '\t', header=T)
exp_meta = exp_meta[,-2]
exp_meta = exp_meta[!duplicated(exp_meta$Hugo_Symbol),]
colnames(exp_meta) <- gsub('.', '', colnames(exp_meta), fixed=T)

pam_meta = tcga_pam50(exp_meta)

write.csv(pam_tcga, 'data/META_PAM50_label.csv')
pam_tcga <- read.csv('data/META_PAM50_label.csv', row.names = 1)

rownames(exp_meta) = exp_meta$Hugo_Symbol
exp_meta = exp_meta[,-1]
exp_meta = exp_meta[,colnames(exp_meta) %in% filter(pam_meta, subtype=='Basal')$id]
#exp_meta = exp_meta[,colnames(exp_meta) %in% unique(maf_meta@data$Tumor_Sample_Barcode)]

him_meta = subtype_casting(na.omit(exp_meta), him_orig)

nrow(him_meta %>% filter(clust == 'CS1')) # 44
nrow(him_meta %>% filter(clust == 'CS2')) # 99
nrow(him_meta %>% filter(clust == 'CS3')) # 101

write.csv(him_meta, 'data/META_HIM_label.csv')
him_meta <- read.csv('data/META_HIM_label.csv', row.names = 1)

# CPTAC processing

exp_cptac = read.table("data/CPTAC/data_mrna_seq_fpkm.txt", header=T)
exp_cptac = exp_cptac[!duplicated(exp_cptac$Hugo_Symbol),]
rownames(exp_cptac) = exp_cptac$Hugo_Symbol
exp_cptac = exp_cptac[,-1]

info_cptac = read.table("data/CPTAC/data_clinical_sample.txt", sep='\t', header=T)
exp_cptac = exp_cptac[,colnames(exp_cptac) %in% filter(info_cptac, PAM50=='Basal')$SAMPLE_ID]

him_cptac = subtype_casting(na.omit(exp_cptac), him_orig)

nrow(him_cptac %>% filter(clust == 'CS1')) # 5
nrow(him_cptac %>% filter(clust == 'CS2')) # 12
nrow(him_cptac %>% filter(clust == 'CS3')) # 12

write.csv(him_cptac, 'data/CPTAC_HIM_label.csv')

# Mutation processing METABRIC

mut_meta = read.table('data/META/data_mutations_extended.txt', sep='\t', header=T)
mut_meta$Tumor_Sample_Barcode = gsub('-', '', mut_meta$Tumor_Sample_Barcode, fixed=T)
mut_meta = filter(mut_meta, Tumor_Sample_Barcode %in% rownames(him_meta))

write.table(mut_meta, file='data/maf_META.txt', row.names=F, sep='\t', quote=F)

# Mutation processing CPTAC

mut_cptac = read.table('data/CPTAC/data_mutations.txt', sep='\t', header=T)
mut_cptac = filter(mut_cptac, Tumor_Sample_Barcode %in% rownames(him_cptac))

write.table(mut_cptac, file='data/maf_CPTAC.txt', row.names=F, sep='\t', quote=F)

#####

# Maftools mutation signature --------------------------------------------------

# Metabric

maf_meta <- read.maf("data/maf_META.txt")

maf_meta@data$clust <- him_meta[match(maf_meta@data$Tumor_Sample_Barcode,rownames(him_meta)),]$clust
maf_meta@clinical.data$clust <- him_meta[match(maf_meta@clinical.data$Tumor_Sample_Barcode,rownames(him_meta)),]$clust
maf_meta@clinical.data <- maf_meta@clinical.data[order(maf_meta@clinical.data$clust,decreasing = F),]

maf_cptac <- read.maf("data/maf_CPTAC.txt")

maf_cptac@data$clust <- him_cptac[match(maf_cptac@data$Tumor_Sample_Barcode,rownames(him_cptac)),]$clust
maf_cptac@clinical.data$clust <- him_cptac[match(maf_cptac@clinical.data$Tumor_Sample_Barcode,rownames(him_cptac)),]$clust
maf_cptac@clinical.data <- maf_cptac@clinical.data[order(maf_cptac@clinical.data$clust,decreasing = F),]

getSampleSummary(maf_meta)

plotmafSummary(maf = maf_meta, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

mafcolor = list(clust = c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6'))

# 12*16 horizontal
oncoplot(maf = maf_meta, 
         top = 50, 
         clinicalFeatures = 'clust',
         sortByAnnotation = T,
         annotationColor = mafcolor,
         groupAnnotationBySize = F,
         draw_titv = TRUE)

oncoplot(maf = maf_cptac,
         top = 50,
         draw_titv = TRUE)



anno_cols = list(clust = c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6'))

oncoplot(maf = maf_meta, 
         clinicalFeatures = 'clust', 
         sortByAnnotation = TRUE, 
         annotationColor = anno_cols)




dev.off()

write.table()

# This part finds genes with significant mutations across the three subtypes
# Initialize binary mutation matrix

mut_meta_samples <- names(table(maf_meta@data$Tumor_Sample_Barcode))

mut_bi = data.frame(all=unique(maf_meta@data$Hugo_Symbol))#162
for (i in 1:length(mut_meta_samples)) {
  mut_bi[,i+1] = ifelse(mut_bi$all %in% maf_meta@data$Hugo_Symbol[which(maf_meta@data$Tumor_Sample_Barcode==mut_meta_samples[i])],
                      1,
                      0)
}
rownames(mut_bi) = mut_bi$all
mut_bi = mut_bi[,-1]
colnames(mut_bi) <- mut_meta_samples

# Count number of samples mutated per gene per subtype using binary matrix

mut_subtype = data.frame(gene = rownames(mut_bi))

mut_bi_temp = mut_bi[,colnames(mut_bi) %in% rownames(filter(him_meta, clust == 'CS1'))] # 31
mut_subtype$CS1_count = apply(mut_bi_temp, 1, sum)
mut_subtype$CS1_ratio = apply(mut_bi_temp, 1, sum) / ncol(mut_bi_temp)

mut_bi_temp = mut_bi[,colnames(mut_bi) %in% rownames(filter(him_meta, clust == 'CS2'))] # 59
mut_subtype$CS2_count = apply(mut_bi_temp, 1, sum)
mut_subtype$CS2_ratio = apply(mut_bi_temp, 1, sum) / ncol(mut_bi_temp)

mut_bi_temp = mut_bi[,colnames(mut_bi) %in% rownames(filter(him_meta, clust == 'CS3'))] # 67
mut_subtype$CS3_count = apply(mut_bi_temp, 1, sum)
mut_subtype$CS3_ratio = apply(mut_bi_temp, 1, sum) / ncol(mut_bi_temp)

# Using fisher's exact test to determine significant p-value

for (i in 1:nrow(mut_subtype)){
  mut_subtype$p_overall[i] = fisher.test(matrix(c(mut_subtype$CS1_count[i],
                                                  33 - mut_subtype$CS1_count[i],
                                                  mut_subtype$CS2_count[i],
                                                  59 - mut_subtype$CS2_count[i],
                                                  mut_subtype$CS3_count[i],
                                                  67 - mut_subtype$CS3_count[i]), 
                                                nrow=2))$p.value
}

for (i in 1:nrow(mut_subtype)){
  mut_subtype$p_1v2[i] = fisher.test(matrix(c(mut_subtype$CS1_count[i],
                                              33 - mut_subtype$CS1_count[i],
                                              mut_subtype$CS2_count[i],
                                              59 - mut_subtype$CS2_count[i]), 
                                            nrow=2))$p.value
}

for (i in 1:nrow(mut_subtype)){
  mut_subtype$p_1v3[i] = fisher.test(matrix(c(mut_subtype$CS1_count[i],
                                              33 - mut_subtype$CS1_count[i],
                                              mut_subtype$CS3_count[i],
                                              67 - mut_subtype$CS3_count[i]), 
                                            nrow=2))$p.value
}

for (i in 1:nrow(mut_subtype)){
  mut_subtype$p_2v3[i] = fisher.test(matrix(c(mut_subtype$CS2_count[i],
                                              59 - mut_subtype$CS2_count[i],
                                              mut_subtype$CS3_count[i],
                                              67 - mut_subtype$CS3_count[i]), 
                                            nrow=2))$p.value
}

write.csv(mut_subtype, file='data/mutation_char.csv', row.names=F)


# Picking new genes and looking at pathway enrichment

mut_gene_list = c('COL12A1', 'AHNAK', 'DNAH5', 'DNAH2', 'MLL2', 'TP53')

de_list = list()

for (i in mut_gene_list) {
  samples = (maf_meta@data %>% filter(Hugo_Symbol == i))$Tumor_Sample_Barcode
  de_list[[i]] = limma_de(exp_meta, samples)
}

library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)


gsea_list <- list()

for (i in names(de_list)) {
  
  gene <- bitr(rownames(de_list[[i]]),
               fromType = "SYMBOL",
               toType =  "ENTREZID",
               OrgDb = org.Hs.eg.db)
  
  gene$lfc <- de_list[[i]][gene$SYMBOL,]$logFC
  gene <- gene[complete.cases(gene),]
  gene <- gene[order(gene$lfc, decreasing = T),]
  
  geneList <- gene$lfc
  names(geneList) <- gene$ENTREZID
  
  ego <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               ont          = "ALL",
               nPerm        = 1000,
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 1,
               verbose      = FALSE)
  
  gsea_list[[i]] <- ego
  
}

#View(gsea_list[[n]]@result)

# 5*6 horizontal
n = 'COL12A1'
gseaplot2(gsea_list[[n]], 
          geneSetID = c(match('extracellular matrix structural constituent', gsea_list[[n]]@result$Description),
                        match('extracellular structure organization', gsea_list[[n]]@result$Description),
                        match('extracellular matrix organization', gsea_list[[n]]@result$Description),
                        match('collagen-containing extracellular matrix', gsea_list[[n]]@result$Description)))

n = 'DNAH2'
gseaplot2(gsea_list[[n]], 
          geneSetID = c(match('extracellular matrix organization', gsea_list[[n]]@result$Description),
                        match('adaptive immune response', gsea_list[[n]]@result$Description),
                        match('collagen-containing extracellular matrix', gsea_list[[n]]@result$Description),
                        match('humoral immune response', gsea_list[[n]]@result$Description)))

n = 'DNAH5'
gseaplot2(gsea_list[[n]], 
          geneSetID = c(match('mitotic nuclear division', gsea_list[[n]]@result$Description),
                        match('nuclear chromosome segregation', gsea_list[[n]]@result$Description),
                        match('mitotic cell cycle phase transition', gsea_list[[n]]@result$Description),
                        match('positive regulation of cell cycle process', gsea_list[[n]]@result$Description)))

n = 'TP53'
gseaplot2(gsea_list[[n]], 
          geneSetID = c(match('mitotic nuclear division', gsea_list[[n]]@result$Description),
                        match('nuclear chromosome segregation', gsea_list[[n]]@result$Description),
                        match('mitotic cell cycle phase transition', gsea_list[[n]]@result$Description),
                        match('positive regulation of cell cycle process', gsea_list[[n]]@result$Description)))


# Plotting

plot_df_1 = mut_subtype[mut_subtype$gene %in% c('COL12A1', 'DNAH2'), c(1,3,5,7)] 
plot_df_1 = reshape2::melt(plot_df_1)

p1 = ggplot(plot_df_1, aes(x=gene, y=value, fill=variable)) +
  geom_bar(position='dodge', stat='identity') +
  scale_fill_manual(values = c(CS1_ratio='#E71D36', CS2_ratio='#FF9F1C', CS3_ratio='#2EC4B6')) +
  theme_classic()

plot_df_2 = mut_subtype[mut_subtype$gene %in% c('TP53', 'DNAH5'), c(1,3,5,7)] 
plot_df_2 = reshape2::melt(plot_df_2)

p2 = ggplot(plot_df_2, aes(x=gene, y=value, fill=variable)) +
  geom_bar(position='dodge', stat='identity') +
  scale_fill_manual(values = c(CS1_ratio='#E71D36', CS2_ratio='#FF9F1C', CS3_ratio='#2EC4B6')) +
  theme_classic()

plot_df_3 = mut_subtype[mut_subtype$gene %in% c('BRCA1', 'BRCA2'), c(1,3,5,7)] 
plot_df_3 = reshape2::melt(plot_df_3)

p3 = ggplot(plot_df_3, aes(x=gene, y=value, fill=variable)) +
  geom_bar(position='dodge', stat='identity') +
  scale_fill_manual(values = c(CS1_ratio='#E71D36', CS2_ratio='#FF9F1C', CS3_ratio='#2EC4B6')) +
  theme_classic()

# 4*12 horizontal
ggarrange(p1, p2, p3, nrow=1)


#####



