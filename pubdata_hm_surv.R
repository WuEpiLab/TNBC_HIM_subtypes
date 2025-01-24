library(maftools)
library(data.table)
library(org.Hs.eg.db)
library(dplyr)

source('functions.R')

him_orig <- readRDS('HIM_orig.rds')

# TCGA Prediction Part ####

# TCGA isolating TNBC
exp_tcga = read.table('D:/Data/TCGA/BRCA/HiSeqV2', sep='\t', header=T)
colnames(exp_tcga) = gsub('.', '', colnames(exp_tcga), fixed=T)
pam_tcga = tcga_pam50(exp_tcga)

write.csv(pam_tcga, 'data/BRCA_PAM50_TCGA.csv')

# Predicting TCGA into HIM subtypes

rownames(exp_tcga) = exp_tcga[,1]
exp_tcga = exp_tcga[,-1]
exp_tcga = exp_tcga[,colnames(exp_tcga) %in% rownames(pam_tcga %>% filter(subtype == 'Basal'))]

him_tcga <- HIM_predict_new(na.omit(exp_tcga), him_orig)

nrow(him_tcga %>% filter(clust == 'CS1')) # 38
nrow(him_tcga %>% filter(clust == 'CS2')) # 94
nrow(him_tcga %>% filter(clust == 'CS3')) # 75

write.csv(him_tcga, 'data/TNBC_HIM_TCGA.csv')

#####

# Metabric Prediction Part ####

# Predicting HIM subtypes for metabric data

exp_meta <- read.table("/public/home/frank/data/METABRIC/data_expression_median.txt", sep = '\t', header=T)
exp_meta <- exp_meta[!duplicated(exp_meta$Hugo_Symbol),]
exp_meta <- exp_meta[,-2]

pam_meta = tcga_pam50(exp_meta)
rownames(pam_meta) = gsub('.', '-', rownames(pam_meta), fixed=T)

write.csv(pam_meta, 'data/BRCA_PAM50_META.csv')

rownames(exp_meta) <- exp_meta$Hugo_Symbol
exp_meta <- exp_meta[,-1,]
colnames(exp_meta) <- gsub('.', '-', colnames(exp_meta), fixed=T)

exp_meta <- exp_meta[,colnames(exp_meta) %in% rownames(pam_meta %>% filter(subtype == 'Basal'))]

him_meta <- HIM_predict_new(na.omit(exp_meta), him_orig)

nrow(him_meta %>% filter(clust == 'CS1')) # 44
nrow(him_meta %>% filter(clust == 'CS2')) # 99
nrow(him_meta %>% filter(clust == 'CS3')) # 101

write.csv(him_meta, 'data/TNBC_HIM_META.csv')

# Adding a 4 subtype object for 4 subtype survival analysis

result_4c <- results[[4]][[3]]
result_4c <- result_4c[order(result_4c, decreasing=F)]

him_orig_4 <- list(expression = pt,
                 subtype = paste0('CS', result_4c),
                 samplelabels = names(result_4c))

him_meta_4 <- HIM_predict_new(na.omit(exp_meta), him_orig_4)


#####

# GSVA for TCGA and Metabric ####

library(GSVA)
library(GSEABase)

#pathways <- getGmt('files/h.all.v7.4.symbols.gmt')
pathways <- getGmt("files/c5.go.v7.5.1.symbols.gmt")

# Poisson is too slow
gsva_meta <- gsva(expr=as.matrix(exp_meta), gset.idx.list=pathways, 
                  method="gsva", kcdf="Gaussian", min.sz=10, max.sz=1000)

gsva_tcga <- gsva(expr=as.matrix(exp_tcga), gset.idx.list=pathways, 
                  method="gsva", kcdf="Gaussian", min.sz=10, max.sz=1000)

gsva_nc <- gsva(expr=as.matrix(exp_nc), gset.idx.list=pathways, 
                  method="gsva", kcdf="Gaussian", min.sz=10, max.sz=1000)

gsva_meta <- as.data.frame(gsva_meta)
gsva_meta <- gsva_meta[,rownames(him_meta[order(him_meta$clust),])]
rownames(gsva_meta) <- gsub('HALLMARK_','',rownames(gsva_meta))

gsva_tcga <- as.data.frame(gsva_tcga)
gsva_tcga <- gsva_tcga[,rownames(him_tcga[order(him_tcga$clust),])]
rownames(gsva_tcga) <- gsub('HALLMARK_','',rownames(gsva_tcga))

gsva_nc <- as.data.frame(gsva_nc)
gsva_nc <- gsva_nc[,rownames(him_nc[order(him_nc$clust),])]

# Drawing heatmap
library(pheatmap)

ann_colors <- list(
  clust = c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6'),
  subtype = c(UNS='#cc4e4f', BL1='#de8cdd', BL2='#4150b0', IM='#c98d47', M='#bdbdbd', MSL='#9131bd', LAR='#4cde2f')
)

hallmark_label <- read.csv('data/hallmark_labelling.csv', row.name = 1)

gsva_meta <- gsva_meta[rownames(hallmark_label),]
gsva_tcga <- gsva_tcga[rownames(hallmark_label),]

gsva_meta <- gsva_meta[rownames(gsva_meta) %in% gsea_diff_list,]
gsva_meta <- gsva_meta[gsea_diff_list[gsea_diff_list %in% rownames(gsva_meta)],]

gsva_tcga <- gsva_tcga[rownames(gsva_tcga) %in% gsea_diff_list,]
gsva_tcga <- gsva_tcga[gsea_diff_list[gsea_diff_list %in% rownames(gsva_tcga)],]

gsva_nc <- gsva_nc[rownames(gsva_nc) %in% gsea_diff_list,]
gsva_nc <- gsva_nc[gsea_diff_list,]


# Add estimate prediction for data

esti_tcga = as.data.frame(immunedeconv::deconvolute(exp_tcga, method='estimate'))
rownames(esti_tcga) = esti_tcga$cell_type
esti_tcga = esti_tcga[,-1]
esti_tcga = esti_tcga[,substr(colnames(esti_tcga), 11, 12) == '01']
colnames(esti_tcga) = substr(colnames(esti_tcga), 1, 10)

esti_tcga = esti_tcga[c('immune score', 'stroma score', 'tumor purity'),]


esti_meta = as.data.frame(immunedeconv::deconvolute(exp_meta, method='estimate'))
rownames(esti_meta) = esti_meta$cell_type
esti_meta = esti_meta[,-1]
colnames(esti_meta) = gsub('-', '', colnames(esti_meta))

esti_meta = esti_meta[c('immune score', 'stroma score', 'tumor purity'),]


{
him_tcga = read.csv('data/TNBC_HIM_TCGA.csv', row.names=1)
lehmann_tcga = read.csv('data/lehmann/lehmann_tcga.csv')

him_tcga = him_tcga[substr(rownames(him_tcga), 11, 12) == '01',]
rownames(him_tcga) = substr(rownames(him_tcga), 1, 10)
rownames(lehmann_tcga) = gsub('-', '', lehmann_tcga$patient)

annotations = merge(him_tcga, lehmann_tcga, by='row.names')
rownames(annotations) = annotations$Row.names
annotations = annotations[,c('clust', 'subtype')]
annotations = annotations[order(annotations$clust),]

gsva_tcga = gsva_tcga[,substr(colnames(gsva_tcga), 11, 12) == '01']
colnames(gsva_tcga) = substr(colnames(gsva_tcga), 1, 10)

} # Reading and processing TCGA subtype labels

gsva_tcga_2 = rbind(esti_tcga, gsva_tcga)

# w12 h9
pheatmap(as.matrix(gsva_tcga_2[,rownames(annotations)]), scale = 'row',
         gaps_col = c(sum(annotations$clust %in% 'CS1'), 
                      sum(annotations$clust %in% c('CS1', 'CS2'))),
         gaps_row = c(3, 23, 50),
         annotation_col = annotations,
         annotation_colors = ann_colors,
         cluster_rows = F, cluster_cols = F,
         show_colnames = F,
         main = 'TCGA HIM',
         breaks = seq(-3, 3, length.out = 100),
         color=colorRampPalette(c('#b6e8ff','white','#da3446'))(100))


{
  him_meta = read.csv('data/TNBC_HIM_META.csv', row.names=1)
  lehmann_meta = read.csv('data/lehmann/lehmann_meta.csv')
  
  lehmann_meta = lehmann_meta[!is.na(lehmann_meta$ESR1_mRNA),]
  
  rownames(him_meta) = gsub('-', '', rownames(him_meta))
  rownames(lehmann_meta) = gsub('_', '', lehmann_meta$sample)
  
  annotations = merge(him_meta, lehmann_meta, by='row.names')
  rownames(annotations) = annotations$Row.names
  annotations = annotations[,c('clust', 'TNBCtype')]
  annotations = annotations[order(annotations$clust),]
  colnames(annotations) = c('clust', 'subtype')
  
  colnames(gsva_meta) = gsub('-', '', colnames(gsva_meta))
  
} # Reading and processing METABRIC subtype labels

gsva_meta_2 = rbind(esti_meta, gsva_meta)

pheatmap(as.matrix(gsva_meta_2[,rownames(annotations)]), scale = 'row',
         gaps_col = c(sum(annotations$clust %in% 'CS1'), 
                      sum(annotations$clust %in% c('CS1', 'CS2'))),
         gaps_row = c(3, 23, 50),
         annotation_col = annotations,
         annotation_colors = ann_colors,
         cluster_rows = F, cluster_cols = F,
         show_colnames = F,
         main = 'META HIM',
         breaks = seq(-3, 3, length.out = 100),
         color=colorRampPalette(c('#b6e8ff','white','#da3446'))(100))





annotations <- data.frame(row.names = rownames(him_nc),
                          subtype = him_nc$clust)

pheatmap(as.matrix(gsva_nc), scale = 'row',
         gaps_col = c(sum(him_nc$clust %in% 'CS1'), 
                      sum(him_nc$clust %in% c('CS1', 'CS2'))),
         gaps_row = c(21, 51),
         annotation_col = annotations,
         annotation_colors = ann_colors,
         cluster_rows = F, cluster_cols = F,
         show_colnames = F,
         main = 'NC HIM',
         breaks = seq(-3, 3, length.out = 100),
         color=colorRampPalette(c('#b6e8ff','white','#da3446'))(100))





# Plotting overall distributions
library(ggplot2)

dist <- data.frame(subtype = c('CS1', 'CS2', 'CS3'),
                   tcga = c(nrow(him_tcga %>% filter(clust == 'CS1')),
                            nrow(him_tcga %>% filter(clust == 'CS2')),
                            nrow(him_tcga %>% filter(clust == 'CS3'))),
                   meta = c(nrow(him_meta %>% filter(clust == 'CS1')),
                            nrow(him_meta %>% filter(clust == 'CS2')),
                            nrow(him_meta %>% filter(clust == 'CS3'))))

ggplot(dist, aes(x=subtype, y=tcga, fill=subtype)) +
  geom_col() + 
  scale_fill_manual(values=c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6')) +
  theme_bw()

ggplot(dist, aes(x=subtype, y=meta, fill=subtype)) +
  geom_col() + 
  scale_fill_manual(values=c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6')) +
  theme_bw()


# Sankey plot with lehmann distribution ####

library(ggsankey)
library(ggplot2)

lehmann_colors = c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6',
                   UNS='#cc4e4f', BL1='#de8cdd', BL2='#4150b0', IM='#c98d47', M='#bdbdbd', MSL='#9131bd', LAR='#4cde2f')

lehmann_meta <- read.csv('data/lehmann_meta.csv', row.names = 1)
rownames(lehmann_meta) <- gsub('.', '-', rownames(lehmann_meta), fixed=T)


sankey_meta <- merge(lehmann_meta, him_meta, by=0)

sankey_meta_df <- sankey_meta %>% make_long(clust, subtype)

ggplot(sankey_meta_df, aes(x = x, 
                           next_x = next_x, 
                           node = node, 
                           next_node = next_node,
                           fill = factor(node))) +
  geom_sankey(flow.alpha = 0.3,
              node.color = 'black') +
  theme_sankey(base_size = 16) +
  scale_fill_manual(values = lehmann_colors) +
  ggtitle("Sankey METABRIC")


lehmann_tcga <- read.csv('data/lehmann_tcga.csv', row.names = 1)

lehmann_tcga <- lehmann_tcga[!is.na(lehmann_tcga$correlation),]

sankey_tcga <- merge(lehmann_tcga, him_tcga, by=0)

sankey_tcga_df <- sankey_tcga %>% make_long(clust, subtype)

ggplot(sankey_tcga_df, aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = factor(node))) +
  geom_sankey(flow.alpha = 0.3,
              node.color = 'black') +
  theme_sankey(base_size = 16) +
  scale_fill_manual(values = lehmann_colors) +
  ggtitle("Sankey TCGA")

#####

# Survival and prognostics ####

library(survminer)
library(survival)

# TCGA

surv_tcga <- read.table('data/survival_BRCA_survival.txt', sep='\t', header=T)
him_tcga <- read.csv('data/TNBC_HIM_TCGA.csv', row.names = 1)

surv_tcga$sample <- gsub('-', '', surv_tcga$sample, fixed=T)

surv_tcga <- surv_tcga %>% filter(sample %in% rownames(him_tcga))
rownames(surv_tcga) <- surv_tcga$sample

surv_tcga <- merge(surv_tcga, him_tcga, by=0)

# Change date from days to months

surv_tcga$OS.time <- surv_tcga$OS.time/30
surv_tcga$DSS.time <- surv_tcga$DSS.time/30
surv_tcga$DFI.time <- surv_tcga$DFI.time/30
surv_tcga$PFI.time <- surv_tcga$PFI.time/30

for (i in c(1:nrow(surv_tcga))) {
  if (surv_tcga$OS.time[i] > 120 & !is.na(surv_tcga$OS.time[i])) {
    surv_tcga$OS.time[i] = 120
    surv_tcga$OS[i] = 0
  }
  if (surv_tcga$DSS.time[i] > 120 & !is.na(surv_tcga$DSS.time[i])) {
    surv_tcga$DSS.time[i] = 120
    surv_tcga$DSS[i] = 0
  }
  if (surv_tcga$DFI.time[i] > 120 & !is.na(surv_tcga$DFI.time[i])) {
    surv_tcga$DFI.time[i] = 120
    surv_tcga$DFI[i] = 0
  }
  if (surv_tcga$PFI.time[i] > 120 & !is.na(surv_tcga$PFI.time[i])) {
    surv_tcga$PFI.time[i] = 120
    surv_tcga$PFI[i] = 0
  }
}

ggsurvplot(
  survfit(Surv(OS.time, OS) ~ clust, data = surv_tcga), # survfit object with calculated statistics.
  data = surv_tcga,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  palette = c('#E71D36', '#FF9F1C', '#2EC4B6'),
  # point estimates of survival curves.
  xlim = c(0,120),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  ylab = "TCGA Overall Survival",
  break.time.by = 24,     # break X axis in time intervals by 500.
  ggtheme = theme_classic(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

ggsurvplot(
  survfit(Surv(DSS.time, DSS) ~ clust, data = surv_tcga),  # survfit object with calculated statistics.
  data = surv_tcga,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  palette = c('#E71D36', '#FF9F1C', '#2EC4B6'),
  # point estimates of survival curves.
  xlim = c(0,3650),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  ylab = "TCGA Disease Specific Survival",
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_classic(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

ggsurvplot(
  survfit(Surv(DFI.time, DFI) ~ clust, data = surv_tcga), # survfit object with calculated statistics.
  data = surv_tcga,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  palette = c('#E71D36', '#FF9F1C', '#2EC4B6'),
  # point estimates of survival curves.
  xlim = c(0,3650),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  ylab = "TCGA Disease Free Interval",
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_classic(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

ggsurvplot(
  survfit(Surv(PFI.time, PFI) ~ clust, data = surv_tcga), # survfit object with calculated statistics.
  data = surv_tcga,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  palette = c('#E71D36', '#FF9F1C', '#2EC4B6'),
  # point estimates of survival curves.
  xlim = c(0,3650),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  ylab = "TCGA Progression Free Interval",
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_classic(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

pdf('survival_tcga.pdf', height=6, width=7)
dev.off()

# METABRIC

surv_meta <- read.table("/public/home/frank/data/METABRIC/data_clinical_patient.txt", sep='\t', header=T)

surv_meta <- surv_meta[,c('PATIENT_ID', 'OS_MONTHS', 'OS_STATUS', 'RFS_MONTHS', 'RFS_STATUS')]
surv_meta <- surv_meta %>% filter(PATIENT_ID %in% rownames(him_meta))
rownames(surv_meta) <- surv_meta$PATIENT_ID

surv_meta <- merge(surv_meta, him_meta, by=0)

surv_meta$OS_STATUS <- ifelse(surv_meta$OS_STATUS == '0:LIVING', 0, 1)
surv_meta$RFS_STATUS <- ifelse(surv_meta$RFS_STATUS == '0:Not Recurred', 0, 1)

for (i in c(1:nrow(surv_meta))) {
  if (surv_meta$OS_MONTHS[i] > 120) {
    surv_meta$OS_MONTHS[i] = 120
    surv_meta$OS_STATUS[i] = 0
  }
  if (surv_meta$RFS_MONTHS[i] > 120) {
    surv_meta$RFS_MONTHS[i] = 120
    surv_meta$RFS_STATUS[i] = 0
  }
}

ggsurvplot(
  survfit(Surv(OS_MONTHS, OS_STATUS) ~ clust, data = surv_meta),  # survfit object with calculated statistics.
  data = surv_meta,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  palette = c('#E71D36', '#FF9F1C', '#2EC4B6'),
  # point estimates of survival curves.
  xlim = c(0,120),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  ylab = "META Overall Survival",
  break.time.by = 24,     # break X axis in time intervals by 500.
  ggtheme = theme_classic(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

ggsurvplot(
  survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ clust, data = surv_meta),  # survfit object with calculated statistics.
  data = surv_meta,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  palette = c('#E71D36', '#FF9F1C', '#2EC4B6'),
  # point estimates of survival curves.
  xlim = c(0,120),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  ylab = "META RFS",
  break.time.by = 24,     # break X axis in time intervals by 500.
  ggtheme = theme_classic(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

# Metabric 4 subtypes

surv_meta <- read.table("/public/home/frank/data/METABRIC/data_clinical_patient.txt", sep='\t', header=T)

surv_meta <- surv_meta[,c('PATIENT_ID', 'OS_MONTHS', 'OS_STATUS', 'RFS_MONTHS', 'RFS_STATUS')]
surv_meta <- surv_meta %>% filter(PATIENT_ID %in% rownames(him_meta_4))
rownames(surv_meta) <- surv_meta$PATIENT_ID

surv_meta <- merge(surv_meta, him_meta_4, by=0)

surv_meta$OS_STATUS <- ifelse(surv_meta$OS_STATUS == '0:LIVING', 0, 1)
surv_meta$RFS_STATUS <- ifelse(surv_meta$RFS_STATUS == '0:Not Recurred', 0, 1)

for (i in c(1:nrow(surv_meta))) {
  if (surv_meta$OS_MONTHS[i] > 120) {
    surv_meta$OS_MONTHS[i] = 120
    surv_meta$OS_STATUS[i] = 0
  }
  if (surv_meta$RFS_MONTHS[i] > 120) {
    surv_meta$RFS_MONTHS[i] = 120
    surv_meta$RFS_STATUS[i] = 0
  }
}

ggsurvplot(
  survfit(Surv(OS_MONTHS, OS_STATUS) ~ clust, data = surv_meta),  # survfit object with calculated statistics.
  data = surv_meta,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # palette = c('#E71D36', '#FF9F1C', '#2EC4B6'),
  # point estimates of survival curves.
  xlim = c(0,120),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  ylab = "META Overall Survival",
  break.time.by = 24,     # break X axis in time intervals by 500.
  ggtheme = theme_classic(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

ggsurvplot(
  survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ clust, data = surv_meta),  # survfit object with calculated statistics.
  data = surv_meta,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # palette = c('#E71D36', '#FF9F1C', '#2EC4B6'),
  # point estimates of survival curves.
  xlim = c(0,120),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  ylab = "META RFS",
  break.time.by = 24,     # break X axis in time intervals by 500.
  ggtheme = theme_classic(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

#####

