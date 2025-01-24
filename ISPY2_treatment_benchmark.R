library(ggpubr)

ispy_exp <- read.table('data/ISPY2/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt', sep = '\t' , header = T)
ispy_exp <- ispy_exp[!duplicated(ispy_exp$X), ]

colnames(ispy_exp) <- gsub('X', '', colnames(ispy_exp))
colnames(ispy_exp) <- substr(colnames(ispy_exp), 1, 6)
rownames(ispy_exp) <- ispy_exp[,1]
ispy_exp <- ispy_exp[,-1]

ispy_meta <- read.csv('data/ISPY2/newISPY2_TableS2_info.csv')

ispy_exp <- ispy_exp[,as.character(ispy_meta[ispy_meta$Receptor.Subtype=='TN',]$Patient.Identifier)]

write.csv(ispy_exp, 'data/ISPY2/exp_filtered.csv')

# Reading ISPY2 proteomic data

ispy_prot_1 <- read.csv('data/ISPY2/GSE196093_ISPY2_RPPAdat_139endpts_RPPA1_raw_DeID.csv')
ispy_prot_2 <- read.csv('data/ISPY2/GSE196093_ISPY2_RPPAdat_139endpts_RPPA2_raw_DeID.csv')
ispy_prot_3 <- read.csv('data/ISPY2/GSE196093_ISPY2_RPPAdat_139endpts_RPPA3_raw_DeID.csv')


# Loading in Lehmann subtype for ISPY2 data
ispy_lehmann <- read.csv('data/ISPY2/91671054-2482-42a3-85a4-8b44ead7fe7d_result.csv')

ispy_lehmann$X <-gsub('X', '', ispy_lehmann$X)
rownames(ispy_lehmann) <- ispy_lehmann$X

# Loading in HIM subtype for ISPY2 data
ispy_him <- read.csv('data/ISPY2/PamrRes_3clust_newISPY2_combat.csv')

ispy_him$X <-gsub('X', '', ispy_him$X)
rownames(ispy_him) <- ispy_him$X

# TN data grouping and organization

ispy_tn <- ispy_meta[ispy_meta$Receptor.Subtype=='TN', c('Patient.Identifier','Arm..short.name.','pCR','Immune.','DRD.','I.SPY2.Subtypes','RPS.5')]
colnames(ispy_tn) <- c('patient', 'arm', 'pCR', 'immune', 'DRD', 'ispy_subtype1', 'ispy_subtype2')
ispy_tn$patient <- as.character(ispy_tn$patient)

ispy_tn$lehmann_subtype <- ispy_lehmann[ispy_tn$patient, 'subtype']

ispy_tn$him_subtype <- ispy_him[ispy_tn$patient, 'clust']


# Looking at IM subtype (immune therapy prediction)
ispy_label_im <- filter(ispy_lehmann, subtype=='IM')$X

ispy_rate <- function(ispy_meta, name_list) {
  
  ispy_meta <- filter(ispy_meta, Patient.Identifier %in% name_list)
  arms <- unique(ispy_meta$Arm..short.name.)
  
  result <- data.frame(row.names = arms,
                       arm = arms,
                       pCR_rate = 0,
                       sample_count = 0,
                       responsive = '',
                       unresponsive = '')
  
  for (i in arms) {
    temp_meta <- filter(ispy_meta, Arm..short.name. == i)
    result[i,2] <- nrow(filter(temp_meta, pCR==1))/nrow(temp_meta)
    result[i,3] <- nrow(temp_meta)
    result[i,4] <- paste(filter(temp_meta, pCR==1)$Patient.Identifier, collapse = ' ')
    result[i,5] <- paste(filter(temp_meta, pCR==0)$Patient.Identifier, collapse = ' ')
  }
  
  return(result)
}

ispy_rate_lehmann_IM <- ispy_rate(ispy_meta, ispy_label_im)
ispy_rate_lehmann_IM$pCR_rate = round(ispy_rate_lehmann_IM$pCR_rate*100, digits=2)


ggbarplot(ispy_rate_lehmann_IM, x = "arm", y = "pCR_rate", 
          fill = "arm", palette = c('#000000','#D1D1D1','#FDDED0','#FBB9A0','#FB9072',
                                    '#CA161E','#A50D16','#66000B'), #palette为配色方案，可以设置不同杂志的配色方案，如aaas是science,其他的方案可以百度
          width=0.7,label=T,  lab.col="black")+ 
  ylim(0,110)+
  labs(title = 'Lehmann IM',x='',y='% pCR')+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))


# Comparing responsive to non-responsive in pembro group

ispy_exp_filtered <- ispy_exp[,colnames(ispy_exp) %in% filter(ispy_meta, Arm..short.name.=='Pembro')$Patient.Identifier]
ispy_exp_filtered <- ispy_exp_filtered[,colnames(ispy_exp_filtered) %in% ispy_label_im]

ispy_exp_filtered['COL12A1',]
ispy_exp_filtered['EPCAM',]
ispy_exp_filtered['PECAM1',]



View(filter(ispy_him, X %in% colnames(ispy_exp_filtered)))

View(filter(ispy_him, clust=="CS1" & (X %in% as.character(filter(ispy_meta, Arm..short.name.=='Pembro')$Patient.identifier))))


library(ggsankey)

ispy_tn_pembro <- filter(ispy_tn, arm=='Pembro')

ispy_sankey <- ispy_tn %>% make_long(ispy_subtype2, him_subtype, lehmann_subtype)

ispy_sankey_colors = c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6',
                       UNS='#cc4e4f', BL1='#de8cdd', BL2='#4150b0', IM='#c98d47', M='#bdbdbd', MSL='#9131bd', LAR='#4cde2f',
                       `HER2-/Immune+`='#7575B9' ,`HER2-/Immune-/DRD.v3-`='#A1A45E', `HER2-/Immune-/DRD.v3+`='#ADB16E')


ggplot(ispy_sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = .6) +
  geom_alluvial_text(size = 3, color = "white") +
  #scale_fill_viridis_d(drop = FALSE) +
  scale_fill_manual(values = ispy_sankey_colors) +
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("ISPY2 data")

# Different plot

ispy_sankey <- ispy_tn %>% make_long(ispy_subtype2, him_subtype, arm)

ispy_sankey_colors = c(`HER2-/Immune+`='#7575B9' ,`HER2-/Immune-/DRD.v3-`='#A1A45E', `HER2-/Immune-/DRD.v3+`='#ADB16E',
                       CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6',
                       Ctr='#000000',  N='#f9f2ed', VC='#ecdaca', AMG386='#dab096', MK2206='#ca866a', Ganitumab='#8b2825', Ganetespib='#6c1d1c', Pembro='#3d0e10')

# landscape 9*6
ggplot(ispy_sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = .6) +
  geom_alluvial_text(size = 3, color = "white") +
  #scale_fill_viridis_d(drop = FALSE) +
  scale_fill_manual(values = ispy_sankey_colors) +
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("ISPY2 data")

# Plot sankey for each 

ispy_pcr_sankey <- ispy_tn_pembro %>% make_long(pCR, ispy_subtype2)
ispy_pcr_sankey <- ispy_tn_pembro %>% make_long(pCR, him_subtype)
ispy_pcr_sankey <- ispy_tn_pembro %>% make_long(pCR, lehmann_subtype)

ggplot(ispy_pcr_sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = .6) +
  geom_alluvial_text(size = 3, color = "white") +
  scale_fill_viridis_d(drop = FALSE) +
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("ISPY2 data")

# Calculating estimate score ---------------------------------------------------

ispy_ideconv <- immunedeconv::deconvolute_estimate(na.omit(ispy_exp))

library(estimate)

estimate::outputGCT(na.omit(ispy_exp), 'ispy_exp.gct')
estimate::estimateScore('ispy_exp.gct', 'ispy_estimate.txt', platform='illumina')

ispy_estimate <- read.table('ispy_estimate.txt', sep='\t', header=T, row.names=1)
ispy_estimate <- ispy_estimate[,-1]
colnames(ispy_estimate) <- gsub('X', '', ispy_estimate[1,])
ispy_estimate <- ispy_estimate[-1,]


ispy_estimate[, as.character(filter(ispy_tn_pembro, immune==1 & him_subtype=='CS2')$patient)]


filter(ispy_tn_pembro, immune==1 & him_subtype=='CS2')


ispy_ideconv <- as.data.frame(t(ispy_ideconv))

ispy_tn <- ispy_tn[,c(1:9)]
ispy_tn_pembro <- ispy_tn_pembro[,c(1:9)]


ispy_tn <- cbind(ispy_tn, ispy_ideconv[ispy_tn$patient,])
ispy_tn_pembro <- cbind(ispy_tn_pembro, ispy_ideconv[ispy_tn_pembro$patient,])

ispy_tn_pembro$pCR <- as.character(ispy_tn_pembro$pCR)

# landscape 3*4
ggplot(filter(ispy_tn_pembro, him_subtype%in%c('CS2', 'CS1')), aes(x=ImmuneScore, y=StromalScore, color=pCR==1)) +
  geom_point(size = 2) +
  scale_color_manual(values=c('red', 'black')) +
  theme_classic()

# landscape 3*4.3
ggplot(filter(ispy_tn_pembro, him_subtype%in%c('CS2', 'CS1')), aes(x=ImmuneScore, y=StromalScore, color=him_subtype)) +
  geom_point(size = 3) +
  scale_color_manual(values=c('red', 'black')) +
  theme_classic()

# landscape 3*4
ggplot(filter(ispy_tn_pembro, him_subtype%in%c('CS3', 'CS1')), aes(x=ImmuneScore, y=TumorPurity, color=pCR==1)) +
  geom_point(size = 2) +
  scale_color_manual(values=c('red', 'black')) +
  theme_classic()

# landscape 3*4.3
ggplot(filter(ispy_tn_pembro, him_subtype%in%c('CS3', 'CS1')), aes(x=ImmuneScore, y=TumorPurity, color=him_subtype)) +
  geom_point(size = 3) +
  scale_color_manual(values=c('red', 'black')) +
  theme_classic()

# portrait 5*4.5
ggplot(filter(ispy_tn_pembro, immune==1), aes(y=ImmuneScore, x=him_subtype, fill=him_subtype)) +
  geom_boxplot() +
  geom_jitter() +
  scale_fill_manual(values = c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6')) +
  theme_classic()
   



ispy_tn$pCR <- as.character(ispy_tn$pCR)

# landscape 3*4
ggplot(filter(ispy_tn, immune==1 & him_subtype=='CS2'), aes(x=ImmuneScore, y=StromalScore, color=pCR==1)) +
  geom_point() +
  scale_color_manual(values=c('red', 'black')) +
  theme_classic()

# portrait 5*4.5
ggplot(filter(ispy_tn, immune==1), aes(y=ImmuneScore, x=him_subtype, fill=him_subtype)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6')) +
  theme_classic()

# portrait 5*4.5
ggplot(filter(ispy_tn, immune==1), aes(y=ImmuneScore, x=him_subtype, fill=him_subtype)) +
  geom_boxplot() +
  geom_jitter() +
  scale_fill_manual(values = c(CS1='#E71D36', CS2='#FF9F1C', CS3='#2EC4B6')) +
  theme_classic() +
  stat_compare_means(test = 'wilcox.test', comparisons=list(c('CS1','CS2'), c('CS1','CS3'), c('CS2','CS3')), label = 'p.signif', show.legend = F)

ggsave('ispy_him_immune_boxplot.pdf', width=4.5, height=5)

