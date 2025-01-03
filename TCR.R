install.packages("immunarch")
library(immunarch)
data(immdata)
colnames(immdata$data$`A2-i129`)
repExplore(immdata$data, "lens") %>% vis()  # Visualise the length distribution of CDR3
repClonality(immdata$data, "homeo") %>% vis()  # Visualise the relative abundance of clonotypes
repOverlap(immdata$data) %>% vis()  # Build the heatmap of public clonotypes shared between repertoires
geneUsage(immdata$data[[1]]) %>% vis()  # Visualise the V-gene distribution for the first repertoire

.
├── vdj_v1_mm_c57bl6_pbmc_t_filtered_contig_annotations.csv <-- This contains the count data we want!
  ├── vdj_v1_mm_c57bl6_pbmc_t_consensus_annotations.csv 
├── vdj_v1_mm_c57bl6_pbmc_t_clonotypes.csv
├── vdj_v1_mm_c57bl6_pbmc_t_all_contig_annotations.csv 
├── vdj_v1_mm_c57bl6_pbmc_t_matrix.h5
├── vdj_v1_mm_c57bl6_pbmc_t_bam.bam.bai
├── vdj_v1_mm_c57bl6_pbmc_t_molecule_info.h5
├── vdj_v1_mm_c57bl6_pbmc_t_raw_feature_bc_matrix.tar.gz
├── vdj_v1_mm_c57bl6_pbmc_t_analysis.tar.gz

#data input
path <- file.path("D:", "TCRLEN003", fsep="\\")
immdata_10x <- repLoad(path)
$meta
immdata_10x$data$LEN003F
colnames(immdata_10x$data$LEN002)
library(immunarch)  # Load the package into R
data(immdata)  # Load the test dataset

#repertoire overlap
imm_ov1 <- repOverlap(immdata_10x$data, .method = "public", .verbose = F)
imm_ov2 <- repOverlap(immdata_10x$data, .method = "morisita", .verbose = F)

p1 <- vis(imm_ov1)
p2 <- vis(imm_ov2, .text.size = 2)
vis(imm_ov1, "heatmap2")
p1 + p2

#Diversity
library(immunarch); data(immdata)       # Load the package and the test dataset
?repDiversity

div_div <- repDiversity(immdata_10x$data, "inv.simp")
div_div

#Compute statistics and visualise them
# Chao1 diversity measure
div_chao <- repDiversity(immdata_10x$data, "chao1")

# Hill numbers
div_hill <- repDiversity(immdata$data, "hill")

# D50
div_d50 <- repDiversity(immdata$data, "d50")

# Ecological diversity measure
div_div <- repDiversity(immdata$data, "div")


p1 <- vis(div_chao)
p2 <- vis(div_chao, .by = c("Status", "Sex"), .meta = immdata$meta)
p3 <- vis(div_hill, .by = c("Status", "Sex"), .meta = immdata$meta)

p4 <- vis(div_d50)
p5 <- vis(div_d50, .by = "Status", .meta = immdata$meta)


p1 + p2

div_div <- repDiversity(immdata_10x$data, "div")
p6 <- vis(div_div)


imm_raref <- repDiversity(immdata_10x$data, "raref", .verbose = F)

imm_raref[1:5,]

p1 <- vis(imm_raref)
p2 <- vis(imm_raref, .by = "Status", .meta = immdata$meta)
p1 + p2

#geneusge
# Next four function calls are equal. "hs" is from the "alias" column.
imm_gu <- geneUsage(immdata_10x$data, "hs.trbd")
imm_gu
p1 <- vis(imm_gu[c(1, 2)])
p2 <- vis(imm_gu[c(1, 2)]) + coord_polar()
p3<- vis(imm_gu[c(1, 2)]) + coord_flip() + theme_bw()

library(patchwork)
p1 + p2 +p3 

imm_gu <- geneUsage(immdata_10x$data[c(1, 2)], "hs.trbd", .norm = T, .ambig = "exc")
vis(imm_gu)
imm_gu <- geneUsage(immdata_10x$data[c(1, 2)], "musmus.trbv")


##heliu
library(ggforce)
immdata_10x$data$NC%>%
  gather_set_data(c(6,7,5),.ambig = "exc") %>%
  ggplot(aes(x, id = id, split = y, value = 1))  +
  geom_parallel_sets(aes(fill = J.name), show.legend = FALSE, alpha = 0.3) +
  geom_parallel_sets_axes(axis.width = 0.1, color = "lightgrey", fill = "white") +
  geom_parallel_sets_labels(angle = 0) +
  theme_no_axes()


#clonetype
library(immunarch) 
exp_vol <- repExplore(immdata_10x$data, .method = "volume")
p1 <- vis(exp_vol)
p2 <- vis(exp_vol, .by = c("Status", "Sex"), .meta = immdata$meta)
p1 + p2

exp_len <- repExplore(immdata_10x$data, .method = "len", .col = "aa")
p1 <- vis(exp_len)

#top
imm_top <- repClonality(immdata_10x$data, .method = "top", .head = c(10, 100, 1000))
imm_top
imm_top %>% vis()
#rare
imm_rare <- repClonality(immdata_10x$data, .method = "rare", .bound =c(1, 3, 10))
imm_rare %>% vis()

imm_hom <- repClonality(immdata_10x$data,
                        .method = "homeo",
                        .clone.types = c(Rare=1e-05,Small = 1e-04, Medium =0.001, Large =0.01,Hyperexpand=1)
)
vis(imm_hom)
