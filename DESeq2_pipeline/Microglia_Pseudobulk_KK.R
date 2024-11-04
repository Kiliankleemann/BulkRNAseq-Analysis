#### -------- INSTALLING AND LOADING PACKAGES ------- ####
list.of.packages <- c("BiocGenerics","tximport","S4Vectors", "DESeq2", "biomaRt",'GseaVis',
                      "ggplot2", "ggsignif", "ggpubr", "sva", "devtools", "org.Hs.eg.db", 
                      "org.Mm.eg.db", "limma","stringr","KEGGREST","ggrepel", "openxlsx",
                      "fgsea","clusterProfiler","pheatmap","ggpubr","cowplot",'dplyr',
                      "RColorBrewer",'AnnotationDbi', 'tidyverse','pheatmap', 'dendextend',
                      'sva','caTools', 'quantmod','MASS', 'corrplot','mltools','sjPlot', 
                      'sjmisc','jtools','sjlabelled', 'data.table', 'webshot','Seurat')



new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if(length(new.packages)) BiocManager::install(new.packages)

# Special packages
#devtools::install_github("junjunlab/GseaVis")

# Packages loading
invisible(lapply(list.of.packages, library, character.only = TRUE))




#### -------- FILE PREPARATION ------- ####
# Load Seurat Object
df_mg <- readRDS("~/Library/CloudStorage/GoogleDrive-kiliankleemann@gmail.com/Shared drives/APOE Seurat Objects/Updated objects with pmd in metadata/6.2 Microglia Round 3_metanew.rds")
meta <- df_mg@meta.data


# Generate Metadata file for DESEQ
meta.mod <- meta %>% 
  dplyr::select(c(sample,diagnosis,group,sex, apoe, Tau_score,AB_score, snRNAseq_batch, pmd_min, percent.mito, percent.ribo, age)) %>% 
  group_by(sample) %>%
  mutate(percent.mit.aver = mean(percent.mito),
         percent.ribo.aver = mean(percent.ribo)) %>%
  ungroup() %>%
  remove_rownames %>% 
  dplyr::select(-c(percent.mito, percent.ribo)) %>%
  distinct() %>%
  mutate(Grouped = paste0(sample,"_",group))


write.xlsx(meta.mod, paste0('Metadata/metadata_allsamples_mg.xlsx'))

#Adjusting metadata to combine APOE34 and APOE44
#meta.mod$apoe <- ifelse(meta.mod$apoe == 44, 34, meta.mod$apoe)
meta.mod$group <- paste0(meta.mod$diagnosis,'_', meta.mod$sex,'_', meta.mod$apoe)


#Set up dataframe
sub_meta <- meta.mod
sub_meta <- sub_meta %>% 
  mutate(sex = as.factor(sex),
         apoe = as.factor(apoe),
         group_comb = as.factor(group),
         snRNAseq_batch = as.factor(snRNAseq_batch),
         age = scale(age, center = TRUE), 
         percent.mit.aver = scale(percent.mit.aver, center = TRUE),
         percent.ribo.aver = scale(percent.ribo.aver, center = TRUE), 
         diagnosis = as.factor(diagnosis),
         pmd_min = scale(pmd_min, center = TRUE))


##### ---- Checking for co-linearity all metadata snRNAseq_batch encoding ---- ####
#Correlation Plot
regression_data <- one_hot(as.data.table(sub_meta))

regression_data_select <- regression_data[,c(-1,-4,-11, -24)]
regression_data_select = as.data.frame(sapply(regression_data_select, as.numeric))

model_all <- lm(formula = diagnosis_AD ~ sex_m + apoe_33 + Tau_score  + 
                  snRNAseq_batch_AD10 + snRNAseq_batch_AD11 + snRNAseq_batch_AD4 + snRNAseq_batch_AD5 + snRNAseq_batch_AD6 + 
                  snRNAseq_batch_AD7 + snRNAseq_batch_AD8 + snRNAseq_batch_AD9 + pmd_min + percent.mit.aver + percent.ribo.aver + 
                  age, data = regression_data_select)  


summary(model_all)


#Corrlation plot
M <-cor(regression_data_select)
corrplot(M, type="lower",order="hclust", tl.col="black", tl.srt=45,sig.level = 0.01, insig = "blank")



#barplot for potential co-linearity
#Settings for barplots EDIT!!
percent.ribo.aver

lineWidth = 1
pointSize = 20
dir.create(paste0('plots/statistics_colinearity/metadata/'))
dir.create(paste0('plots/statistics_colinearity/metadata/barplots/'))

#Barplot metadata
goi_data3 <- meta.mod %>% dplyr::select(diagnosis,percent.ribo.aver) 
goi_data3$condition <- factor(goi_data3$diagnosis, levels=unique(goi_data3$diagnosis))
comparisons <- compare_means(
  data = goi_data3,
  formula = percent.ribo.aver ~ condition,
  method = "t.test") 
comparison_list_sign <- comparisons %>% mutate(comparison_list = map2(group1, group2,c)) %>% pull(comparison_list)
barplot <- ggplot(goi_data3, aes(x = diagnosis, y = percent.ribo.aver, fill = diagnosis)) + 
  geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0.4) +
  geom_point(aes(y = percent.ribo.aver)) +
  expand_limits(x = 0, y = 0) +
  theme(panel.spacing = unit(1, "lines")) +
  labs(x = NULL, y = c(paste('percent.mit.aver', "units ± SE"))) +
  ggtitle(paste0('percent.mit.aver'))+
  theme(
    text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(linewidth = lineWidth, colour = "black"),
    plot.title  = element_text(color="black", size=40, face="bold.italic"),
    axis.title  = element_text(size = pointSize, colour = "black"),
    axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
    axis.text.y  = element_text(size = pointSize , colour = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = pointSize , colour = "black"),
    # legend.key.height = unit(0.1, "cm"),
    # legend.key.width = unit(0.2, "cm"),
    axis.line = element_line(linewidth = lineWidth, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  stat_compare_means(comparisons = comparison_list_sign, method = 't.test', label = "p" ) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(expand = c(0, 0, .05, 0))

pdf(file = paste0('plots/statistics_colinearity/metadata/barplots/','percent.mit.aver', '_vs_diagnosis.pdf'), pointsize = 10, width = 5, height= 10)
print(barplot)
dev.off()

#####



#Aggregate by sample
psudobulk_mg <- AggregateExpression(df_mg, group.by = c("sample"),
                                 assays = "RNA",
                                 slot = "counts",
                                 return.seurat = FALSE)

# Extract the counts
psudobulk_mg <- as.matrix(psudobulk_mg$RNA)

#BATCH Correaction
count_batch <- sub_meta$snRNAseq_batch
count_batch
psudobulk_mg
data_adjusted_mg <- ComBat_seq(psudobulk_mg, batch = count_batch, group = NULL)



#### -------- ANALYSIS ALL SAMPLES - SIMPLE DESIGN - WALD - ALL COMPARISONS------- ####
#Setting experiment conditions 
file_prefix <- 'All_samples_simple_WALD_contrast_KK'
data <- data_adjusted_mg
experiment <- sub_meta
experiment

#Designs 
#dds <- DESeqDataSetFromMatrix(data, colData = experiment, design = ~ sex + age + pmd_min + apoe + diagnosis) 
dds <- DESeqDataSetFromMatrix(data, colData = experiment, design = ~ 0 + group + age + pmd_min)

## Run DESeq analysis to gather differential expression results
#Run DESeq (WALD - decault to extract the relevant comparisons
dds_run <- DESeq(dds)

####
#### -------- PCA ----- #####
#Transform counts for data visualization
vst <- vst(dds_run, blind=TRUE)

#Add nametags groups for pca 
groups_pca <- c('sex','apoe', 'snRNAseq_batch', 'diagnosis','group_comb')

z <- plotPCA(vst, intgroup=groups_pca, ntop = 500)
theme_PCA <- theme(aspect.ratio = 1, 
                   panel.background = element_blank(),
                   panel.border=element_rect(fill=NA),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   strip.background=element_blank(),
                   axis.text.x=element_text(colour="black"),
                   axis.text.y=element_text(colour="black"),
                   axis.ticks=element_line(colour="black"),
                   legend.key=element_blank(),
                   plot.margin=unit(c(1,1,1,1),"line"))

#Plot no labels
dir.create(paste0('plots/pca/Microglia/'))
dir.create(paste0('plots/pca/Microglia/',file_prefix))

pdf(file = paste0('plots/pca/Microglia/', file_prefix, '/PCA_top200_', 'SEX_labelled.pdf'), pointsize = 10)
ggplot(z$data,aes(x=z$data$PC1, y=z$data$PC2, )) + geom_point(aes(color = group_comb)) + theme_PCA + 
  labs(title = 'PCA (Top 200 variable genes)',x=paste(z$labels$x), y=paste(z$labels$y))
dev.off()

pdf(file = paste0('plots/pca/Microglia/', file_prefix, '/PCA_top200_', 'APOE_labelled.pdf'), pointsize = 10)
ggplot(z$data,aes(x=z$data$PC1, y=z$data$PC2, )) + geom_point(aes(color = apoe)) + theme_PCA + 
  labs(title = 'PCA (Top 200 variable genes)',x=paste(z$labels$x), y=paste(z$labels$y))
dev.off()

pdf(file = paste0('plots/pca/Microglia/', file_prefix, '/PCA_top200_', 'DIAGNOSIS_labelled.pdf'), pointsize = 10)
ggplot(z$data,aes(x=z$data$PC1, y=z$data$PC2, )) + geom_point(aes(color = diagnosis)) + theme_PCA + 
  labs(title = 'PCA (Top 200 variable genes)',x=paste(z$labels$x), y=paste(z$labels$y))
dev.off()

pdf(file = paste0('plots/pca/Microglia/', file_prefix, '/PCA_top200_', 'BATCH_labelled.pdf'), pointsize = 10)
ggplot(z$data,aes(x=z$data$PC1, y=z$data$PC2, )) + geom_point(aes(color = snRNAseq_batch)) + theme_PCA + 
  labs(title = 'PCA (Top 200 variable genes)',x=paste(z$labels$x), y=paste(z$labels$y))
dev.off()

####


# View names of estimated effects and set up comparisons
resultsNames(dds_run)

# meta data: to be used for complex contrasts
metadt <- data.frame(group = resultsNames(dds_run))

# define the complex contrasts
cont.list <- list( 
  #### case of 'even weighting', no double-delta
  # by default, the factors on each side are: 1 / uniqueN(condition)
  # eqv. to: 0.5*c(A + B) - 0.5*c(C + D)
 # AD vs CTRL
  "AD_vs_CTR" = data.frame(
    left      = c("groupAD_f_33", "groupAD_f_34","groupAD_f_44", "groupAD_m_33", "groupAD_m_34","groupAD_m_44"),
    right     = c("groupCTRL_f_33", "groupCTRL_f_34", "groupCTRL_m_33", "groupCTRL_m_34", NA, NA)
  ),
  "FEMALE:AD_vs_CTR" = data.frame(
    left      = c("groupAD_f_33", "groupAD_f_34", "groupAD_f_44"),
    right     = c("groupCTRL_f_33", "groupCTRL_f_34", NA)
  ),
  "MALE:AD_vs_CTR" = data.frame(
    left      = c("groupAD_m_33", "groupAD_m_34", "groupAD_m_44"),
    right     = c("groupCTRL_m_33", "groupCTRL_m_34",NA)
  ),
  "APOE34:AD_vs_CTR" = data.frame(
    left      = c("groupAD_f_34", "groupAD_m_34"),
    right     = c("groupCTRL_f_34", "groupCTRL_m_34")
  ),
  "FEMALE:APOE34:AD_vs_CTR" = data.frame(
    left      = c("groupAD_f_34"),
    right     = c("groupCTRL_f_34")
  ),
  "MALE:APOE34:AD_vs_CTR" = data.frame(
    left      = c("groupAD_m_34"),
    right     = c("groupCTRL_m_34")
  ),
  "APOE33:AD_vs_CTR" = data.frame(
    left      = c("groupAD_f_33", "groupAD_m_33"),
    right     = c("groupCTRL_f_33", "groupCTRL_m_33")
  ),
  "FEMALE:APOE33:AD_vs_CTR" = data.frame(
    left      = c("groupAD_f_33"),
    right     = c("groupCTRL_f_33")
  ),
  "MALE:APOE33:AD_vs_CTR" = data.frame(
    left      = c("groupAD_m_33"),
    right     = c("groupCTRL_m_33")
  ),
  "APOE34-APOE44:AD_vs_CTR" = data.frame(
    left      = c("groupAD_f_34", "groupAD_m_34","groupAD_f_44", "groupAD_m_44"),
    right     = c("groupCTRL_f_34", "groupCTRL_m_34",NA,NA)
  ),
  "FEMALE:APOE34-APOE44:AD_vs_CTR" = data.frame(
    left      = c("groupAD_f_34", "groupAD_m_34","groupAD_f_44", "groupAD_m_44"),
    right     = c("groupCTRL_f_34", "groupCTRL_m_34",NA,NA)
  ),  
  "MALE:APOE34-APOE44:AD_vs_CTR" = data.frame(
    left      = c("groupAD_f_34", "groupAD_m_44"),
    right     = c("groupCTRL_m_34",NA)
  ),

  #APOE4  vs APOE3
  "APOE34_vs_APOE33" = data.frame(
    left      = c("groupAD_f_34", "groupCTRL_f_34", "groupAD_m_34", "groupCTRL_m_34"),
    right     = c("groupAD_f_33", "groupCTRL_f_33", "groupAD_m_33", "groupCTRL_m_33")
  ),
  "FEMALE:APOE34_vs_APOE33" = data.frame(
    left      = c("groupCTRL_f_34", "groupAD_f_34"),
    right     = c("groupCTRL_f_33", "groupAD_f_33")
  ),
  "MALE:APOE34_vs_APOE33" = data.frame(
    left      = c("groupCTRL_m_34","groupAD_m_34"),
    right     = c("groupCTRL_m_33","groupAD_m_33")
  ),
  "AD:APOE34_vs_APOE33" = data.frame(
    left      = c("groupAD_f_34","groupAD_m_34"),
    right     = c("groupAD_f_33","groupAD_m_33")
  ),
  "AD:APOE44_vs_APOE33" = data.frame(
    left      = c("groupAD_f_44","groupAD_m_44"),
    right     = c("groupAD_f_33","groupAD_m_33")
  ),
  "AD:APOE34-APOE44_vs_APOE33" = data.frame(
    left      = c("groupAD_f_34","groupAD_m_34","groupAD_f_44","groupAD_m_44"),
    right     = c("groupAD_f_33","groupAD_m_33",NA,NA)
  ),
  "AD:APOE33-APOE34_vs_APOE44" = data.frame(
    left      = c("groupAD_f_34","groupAD_m_34","groupAD_f_33","groupAD_m_33"),
    right     = c("groupAD_f_44","groupAD_m_44",NA,NA)
  ),
  "CTRL:APOE34_vs_APOE33" = data.frame(
    left      = c("groupCTRL_f_34","groupCTRL_m_34"),
    right     = c("groupCTRL_f_33","groupCTRL_m_33")
  ),
  "FEMALE:CTRL:APOE34_vs_APOE33" = data.frame(
    left      = c("groupCTRL_f_34"),
    right     = c("groupCTRL_f_33")
  ),
  "MALE:CTRL:APOE34_vs_APOE33" = data.frame(
    left      = c("groupCTRL_m_34"),
    right     = c("groupCTRL_m_33")
  ),
  "FEMALE:AD:APOE34_vs_APOE33" = data.frame(
    left      = c("groupAD_f_34"),
    right     = c("groupAD_f_33")
  ),
  "MALE:AD:APOE34_vs_APOE33" = data.frame(
    left      = c("groupAD_m_34"),
    right     = c("groupAD_m_33")
  ),
  "FEMALE:AD:APOE44_vs_APOE33" = data.frame(
    left      = c("groupAD_f_44"),
    right     = c("groupAD_f_33")
  ),
  "MALE:AD:APOE44_vs_APOE33" = data.frame(
    left      = c("groupAD_m_44"),
    right     = c("groupAD_m_33")
  ),
  "FEMALE:AD:APOE34-APOE44_vs_APOE33" = data.frame(
    left      = c("groupAD_f_34","groupAD_f_44"),
    right     = c("groupAD_f_33", NA)
  ),
  "MALE:AD:APOE34-APOE44_vs_APOE33" = data.frame(
    left      = c("groupAD_m_34","groupAD_m_44"),
    right     = c("groupAD_m_33",NA)
  ),
  ### case of double-delta with 'even weighting'
  #eqv. to: c(A - C) - c(B - D) !!
  #AD vs CTRL comparing Female and Male
  "(FEMALE:AD_vs_CTRL)_vs_(MALE:AD_vs_CTRL)"  = data.frame(
    left      = c("groupAD_f_34", "groupAD_f_33","groupAD_f_44","groupCTRL_f_34", "groupCTRL_f_33"),
    right     = c("groupAD_m_34", "groupAD_m_33","groupAD_m_44","groupCTRL_m_34", "groupCTRL_m_33"),
    leftFC    = c(0.3333,0.3333,0.3333,-0.5,-0.5),
    rightFC   = c(0.3333,0.3333,0.3333,-0.5,-0.5)
  ),
  #AD vs CTRL comparing APOE34 and APOE3
  "(APOE34:AD_vs_CTRL)_vs_(APOE33:AD_vs_CTRL)"  = data.frame(
    left      = c("groupAD_f_34", "groupAD_m_34","groupCTRL_m_34", "groupCTRL_f_34"),
    right     = c("groupAD_f_33", "groupAD_m_33","groupCTRL_m_33", "groupCTRL_f_33"),
    leftFC    = c(0.5,0.5,-0.5,-0.5),
    rightFC   = c(0.5,0.5,-0.5,-0.5)
  ),
  #AD vs CTRL comparing APOE34/44 and APOE3
  "(APOE34-44:AD_vs_CTRL)_vs_(APOE33:AD_vs_CTRL)"  = data.frame(
    left      = c("groupAD_f_34", "groupAD_f_44","groupAD_m_34", "groupAD_m_44","groupCTRL_f_34", "groupCTRL_m_34"),
    right     = c("groupAD_f_33", "groupAD_m_33","groupCTRL_f_33", "groupCTRL_m_33",NA, NA),
    leftFC    = c(0.25,0.25,0.25,0.25,-0.5,-0.5),
    rightFC   = c(0.5,0.5,-0.5,-0.5,0,0)
  ),
  # ## SAME DIRECTION!
  # #AD vs CTRL comparing APOE34 and APOE3
  # "(APOE34:AD_vs_CTRL)_AND_(APOE33:AD_vs_CTRL)"  = data.frame(
  #   left      = c("groupAD_f_34", "groupAD_m_34","groupCTRL_m_34", "groupCTRL_f_34"),
  #   right     = c("groupAD_f_33", "groupAD_m_33","groupCTRL_m_33", "groupCTRL_f_33"),
  #   leftFC    = c(0.5,0.5,-0.5,-0.5),
  #   rightFC   = c(-0.5,-0.5,0.5,0.5)
  # ),
  # #AD vs CTRL comparing APOE34/44 and APOE3
  # "(APOE34-44:AD_vs_CTRL)_AND_(APOE33:AD_vs_CTRL)"  = data.frame(
  #   left      = c("groupAD_f_34", "groupAD_f_44","groupAD_m_34", "groupAD_m_44","groupCTRL_f_34", "groupCTRL_m_34"),
  #   right     = c("groupAD_f_33", "groupAD_m_33","groupCTRL_f_33", "groupCTRL_m_33",NA, NA),
  #   leftFC    = c(0.25,0.25,0.25,0.25,-0.5,-0.5),
  #   rightFC   = c(-0.5,-0.5,0.5,0.5,0,0)
  # )
  #Female APOE 34-44 AD vs CTRL vs female apoe33 AD vs CTRL
  "(FEMALE:APOE34-44:AD_vs_CTRL)_vs_(FEMALE:APOE33:AD_vs_CTRL)"  = data.frame(
    left      = c("groupAD_f_34", "groupAD_f_44","groupCTRL_f_34"),
    right     = c("groupAD_f_33", "groupCTRL_f_33",NA),
    leftFC    = c(0.5,0.5, -1),
    rightFC   = c(1, -1,0)
  ),
  #Male APOE34-44 AD vs CTRL vs Male APOE33 AD vs CTRL
  "(MALE:APOE34-44:AD_vs_CTRL)_vs_(MALE:APOE33:AD_vs_CTRL)"  = data.frame(
    left      = c("groupAD_m_34", "groupAD_m_44","groupCTRL_m_34"),
    right     = c("groupAD_m_33", "groupCTRL_m_33",NA),
    leftFC    = c(0.5,0.5, -1),
    rightFC   = c(1, -1,0)
  ),
  #Female APOE34 AD vs CTRL vs Female APOE33 AD vs CTRL
  "(FEMALE:APOE34:AD_vs_CTRL)_vs_(FEMALE:APOE33:AD_vs_CTRL)"  = data.frame(
    left      = c("groupAD_f_34","groupCTRL_f_34"),
    right     = c("groupAD_f_33", "groupCTRL_f_33"),
    leftFC    = c(1, -1),
    rightFC   = c(1, -1)
  ),
  #Male APOE34 AD vs CTRL vs Male APOE34 AD vs CTRL
  "(MALE:APOE34:AD_vs_CTRL)_vs_(MALE:APOE33:AD_vs_CTRL)"  = data.frame(
    left      = c("groupAD_m_34","groupCTRL_m_34"),
    right     = c("groupAD_m_33", "groupCTRL_m_33"),
    leftFC    = c(1, -1),
    rightFC   = c(1, -1)
  ),

  #Apoe34 vs APOE33 comparing CTRL females and males
  "(FEMALE:CTRL:APOE34_vs_APOE33)_vs_(MALE:CTRL:APOE34_vs_APOE33)"  = data.frame(
    left      = c("groupCTRL_f_34", "groupCTRL_f_33"),
    right     = c("groupCTRL_m_34", "groupCTRL_m_33"),
    leftFC    = c(1, -1),
    rightFC   = c(1, -1)
  ),
  #Apoe34 vs APOE33 comparing AD females and males
  "(FEMALE:AD:APOE34_vs_APOE33)_vs_(MALE:AD:APOE34_vs_APOE33)"  = data.frame(
    left      = c("groupAD_f_34", "groupAD_f_33"),
    right     = c("groupAD_m_34", "groupAD_m_33"),
    leftFC    = c(1, -1),
    rightFC   = c(1, -1)
  ),
  #Apoe34/44 vs APOE33 comparing AD females and males
  "(FEMALE:AD:APOE34-44_vs_APOE33)_vs_(MALE:AD:APOE34-44_vs_APOE33)"  = data.frame(
    left      = c("groupAD_f_34","groupAD_f_44", "groupAD_f_33"),
    right     = c("groupAD_m_34","groupAD_m_44", "groupAD_m_33"),
    leftFC    = c(0.5,0.5, -1),
    rightFC   = c(0.5,0.5, -1)
  ),
  #Apoe34 vs APOE33 comparing AD and CTRL only females
  "(FEMALE:AD:APOE34_vs_APOE33)_vs_(FEMALE:CTRL:APOE34_vs_APOE33)"  = data.frame(
    left      = c("groupAD_f_34", "groupAD_f_33"),
    right     = c("groupCTRL_f_34", "groupCTRL_f_33"),
    leftFC    = c(1, -1),
    rightFC   = c(1, -1)
  ),
  #Apoe34/44 vs APOE33 comparing AD and CTRL only females
  "(FEMALE:AD:APOE34-44_vs_APOE33)_vs_(FEMALE:CTRL:APOE34_vs_APOE33)"  = data.frame(
    left      = c("groupAD_f_34", "groupAD_f_44", "groupAD_f_33"),
    right     = c("groupCTRL_f_34", "groupCTRL_f_33",NA),
    leftFC    = c(0.5,0.5, -1),
    rightFC   = c(1, -1,0)
  ),
  #Apoe34 vs APOE33 comparing AD and CTRL only males
  "(MALE:AD:APOE34_vs_APOE33)_vs_(MALE:CTRL:APOE34_vs_APOE33)"  = data.frame(
    left      = c("groupAD_m_34", "groupAD_m_33"),
    right     = c("groupCTRL_m_34", "groupCTRL_m_33"),
    leftFC    = c(1, -1),
    rightFC   = c(1, -1)
  ),
  #Apoe34/44 vs APOE33 comparing AD and CTRL only males
  "(MALE:AD:APOE34-44_vs_APOE33)_vs_(MALE:CTRL:APOE34_vs_APOE33)"  = data.frame(
    left      = c("groupAD_m_34", "groupAD_m_44", "groupAD_m_33"),
    right     = c("groupCTRL_m_34","groupCTRL_m_33",NA),
    leftFC    = c(0.5,0.5, -1),
    rightFC   = c(1, -1,0)
  )
)

cont.list
#!!! EDIT EXPORT FILE AND DOUBLE CHECK COMPARISONS
cmplx.cont <- makeComplexContrasts(cmplx = cont.list,
                                   contrastDT = metadt,
                                   useComplexNames = FALSE,
                                   save2file = TRUE,
                                   filename = 'results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/ALL_COMPARISONS.xlsx',  #requires .xlsx at the end
                                   export.as = 'df')

cmplx.cont <- cmplx.cont[resultsNames(dds_run),]
cmplx.cont
dds_result <- data.frame()

for (i in 1:ncol(cmplx.cont)) {
  dds_result_total <- results(dds_run, pAdjustMethod = 'BH', independentFiltering = T, cooksCutoff = T, 
                            contrast = cmplx.cont[,i])
  dds_result_total$Contrast <- colnames(cmplx.cont)[i]
  dds_result_total <- as.data.frame(dds_result_total) %>% filter(pvalue < 0.05) %>% filter(baseMean > 1) %>% na.omit()
  dds_result_total <- dds_result_total %>% rownames_to_column('gene')
  dds_result <- rbind(dds_result,dds_result_total)
}

#Results
results_pval <- as.data.frame(table(dds_result$Contrast))
results_pval$group <- paste('pvalue_005')

dds_result_pval001 <- dds_result %>% filter(pvalue < 0.01)
results_pval001 <- as.data.frame(table(dds_result_pval001$Contrast))
results_pval001$group <- paste('pvalue_001')

dds_result_pval0001 <- dds_result %>% filter(pvalue < 0.001)
results_pval0001 <- as.data.frame(table(dds_result_pval0001$Contrast))
results_pval0001$group <- paste('pvalue_0001')

dds_result_padj01 <- dds_result %>% filter(padj < 0.1)
results_padj01 <- as.data.frame(table(dds_result_padj01$Contrast))
results_padj01$group <- paste('padj_01')

dds_result_padj005 <- dds_result %>% filter(padj < 0.05)
results_padj005 <- as.data.frame(table(dds_result_padj005$Contrast))
results_padj005$group <- paste('padj_005')

total_Degs_subset <- rbind(results_pval,results_pval001,results_pval0001,results_padj01,results_padj005)
write.xlsx(total_Degs_subset,"total_Degs_subset.xlsx")

dir.create(paste0('plots/barplots/Microglia/','DEGs_across_comparisons','/'))
barplot_dir <- paste0('plots/barplots/Microglia/','DEGs_across_comparisons','/',file_prefix)
barplot_dir
lineWidth=0.5
pointSize=20
barplot <- ggplot(total_Degs_subset, aes(x = Var1, y = Freq, fill = group)) + 
    geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
     expand_limits(x = 0, y = 0) +
    theme(panel.spacing = unit(1, "lines")) +
    labs(x = NULL, y = c("Number of DEGs")) +
    ggtitle(paste0('No_DEGs_cutoffs'))+
    theme(
      text = element_text(size = pointSize, colour = "black"),
      rect = element_blank(),
      line = element_line(linewidth = lineWidth, colour = "black"),
      plot.title  = element_text(color="black", size=40, face="bold.italic"),
      axis.title  = element_text(size = pointSize, colour = "black"),
      axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
      axis.text.y  = element_text(size = pointSize , colour = "black"),
      legend.position = "left",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = pointSize , colour = "black"),
      axis.line = element_line(linewidth = lineWidth, colour = "black"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    scale_fill_brewer(palette = "Paired") + scale_y_log10(expand = c(0, 0, .05, 0))
pdf(file = paste0(barplot_dir,'.pdf'), pointsize = 10, width = 30, height= 20)
print(barplot)
dev.off()

  
# View Total number of normalized counts per sample
dds_run <- estimateSizeFactors(dds_run)
DS_norm_counts <- counts(dds_run, normalized = TRUE)
colSums(DS_norm_counts)

#### -------- OUTPUT TABLES & RESULTS 
# Build significant gene table P - VALUE 0.05 and extract list of sorted DE genes.
sig_res <- dds_result %>%
  data.frame() %>%
  #rownames_to_column(var='Gene_Contrast') %>%  
  as_tibble %>%
  filter(pvalue < 0.05) %>%
  filter(baseMean > 1) %>%
  arrange(pvalue) 

# Build significant gene table P - ADJUSTED 0.05 and extract list of sorted DE genes.
sig_res_adj <- sig_res %>%
  filter(padj < 0.05) %>%
  arrange(padj) 

sorted_DEGenes <- sig_res$gene %>% unique()
sorted_DEGenes_padj <- sig_res_adj$gene %>% unique()

# Export significant gene count data
name_list <- c('gene', (paste0(dds_run$group, "_", dds_run$sample)))
name_list

sig_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes, gene)) 

sig_export_adj <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes_padj, gene))

# Setup Order of Counts
reordered_index <- c(grep("CTRL_m_33", names(sig_export), ignore.case = T),
                     grep("CTRL_m_34", names(sig_export), ignore.case = T),
                     grep("CTRL_f_33", names(sig_export), ignore.case = T),
                     grep("CTRL_f_34", names(sig_export), ignore.case = T),
                     grep("AD_m_33", names(sig_export), ignore.case = T),
                     grep("AD_m_34", names(sig_export), ignore.case = T),
                     grep("AD_m_44", names(sig_export), ignore.case = T),
                     grep("AD_f_33", names(sig_export), ignore.case = T),
                     grep("AD_f_34", names(sig_export), ignore.case = T),
                     grep("AD_f_44", names(sig_export), ignore.case = T))

sig_export <- sig_export %>%
  dplyr::select(c(1, reordered_index))

sig_export_adj <- sig_export_adj %>%
  dplyr::select(c(1, reordered_index))

# Column Names
check_cols <- c(colnames(sig_export))
check_cols

# Normalized counts
check_counts <- DS_norm_counts %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  dplyr::select(c(1, reordered_index)) %>%
  as_tibble %>% 
  `colnames<-`(check_cols) 


# Results
check_res <- dds_result 
file_prefix
# Output files
dir.create(paste0('results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/', file_prefix))
dir.create(paste0('results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/', file_prefix, '/# All_DEGs'))

results_dir <- paste0('results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/', file_prefix, '/# All_DEGs')
results_dir

#Statistics
write.xlsx(sig_res, file = paste0(results_dir, '/DEGene_statistics_pval05.xlsx'), overwrite = T)
write.xlsx(sig_res_adj, file = paste0(results_dir, '/DEGene_statistics_padj05.xlsx'), overwrite = T)
write.xlsx(check_res, file = paste0(results_dir, '/Statistics.xlsx'), overwrite = T)
#Counts
write.xlsx(sig_export, file = paste0(results_dir, '/DEGene_counts_pval05.xlsx'), overwrite = T)
write.xlsx(sig_export_adj, file = paste0(results_dir, '/DEGene_counts_padj05.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0(results_dir, '/DS_counts.xlsx'), overwrite = T)



##### HEATMAP 
# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_final <- sig_export %>% column_to_rownames('gene')
data_final_adj <- sig_export_adj %>% column_to_rownames('gene')

data_subset_50 <- data_final %>% head(n = 50)
data_subset_200 <- data_final %>% head(n = 200)
data_subset_400 <- data_final %>% head(n = 400)


#Grouped data final
data_grouped <- data.frame(
  CTRL_m_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("CTRL_m_33")))),
  CTRL_m_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("CTRL_m_34")))),
  CTRL_f_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("CTRL_f_33")))),
  CTRL_f_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("CTRL_f_34")))),
  AD_m_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_m_33")))),
  AD_m_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_m_34")))),
  AD_m_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_m_44")))),
  AD_f_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_34")))),
  AD_f_44 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_44")))))

data_subset_200_grouped <- data.frame(
  CTRL_m_33 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("CTRL_m_33")))),
  CTRL_m_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("CTRL_m_34")))),
  CTRL_f_33 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("CTRL_f_33")))),
  CTRL_f_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("CTRL_f_34")))),
  AD_m_33 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_m_33")))),
  AD_m_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_m_34")))),
  AD_m_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_m_44")))),
  AD_f_33 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_f_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_f_34")))),
  AD_f_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_f_44")))))


data_grouped_padj <- data.frame(
  CTRL_m_33 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("CTRL_m_33")))),
  CTRL_m_34 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("CTRL_m_34")))),
  CTRL_f_33 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("CTRL_f_33")))),
  CTRL_f_34 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("CTRL_f_34")))),
  AD_m_33 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_m_33")))),
  AD_m_34 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_m_34")))),
  AD_m_44 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_m_44")))),
  AD_f_33 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_f_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_f_34")))),
  AD_f_44 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_f_44")))))


# Normalize your data according to Z score using cal_z_score function.
data_norm <- t(apply(data_final, 1, cal_z_score))
data_final_adj_norm <- t(apply(data_final_adj, 1, cal_z_score))
data_subset_50 <- t(apply(data_subset_50, 1, cal_z_score))
data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))
data_grouped_norm_padj <- t(apply(data_grouped_padj, 1, cal_z_score))

#Export directory
dir.create(paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix))
dir.create(paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix, '/# All_DEGs'))
heatmap_dir <- paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix,'/# All_DEGs')

#Setting up Heatmap parameters
gaps = c(6,10,15,21,26,31,36,41,44)
clusters = 20
color_set <- c("navy", "white", "firebrick3")
all_degs_export_height <- nrow(data_final)
padj_degs_export_height <- nrow(data_final_adj)
## ------------- DENDROGRAM AND CLUSTER MANIPULATION (ADVANCED) -------------##
# Create heatmap
pdf(file = paste0(heatmap_dir, '/All_DEGs_padj05.pdf'), pointsize = 10, height = padj_degs_export_height/4)
phm_full <- pheatmap(data_final_adj_norm,
                     color = colorRampPalette(color_set)(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     cellheight = 9,
                     #cellwidth = 25,
                     scale = 'row'
)
dev.off()



pdf(file = paste0(heatmap_dir,'/All_DEGs_pval05.pdf'), pointsize = 10)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(color_set)(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row'
                     #cellheight = 3,
)
dev.off()


pdf(file = paste0(heatmap_dir, '/All_DEGs_labelled_pval05.pdf'), pointsize = 10, height = all_degs_export_height/4)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 8
)
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_50.pdf'), pointsize = 10, height = 10)
phm_50 <- pheatmap(data_subset_50,
                   color = colorRampPalette(color_set)(41),
                   breaks = seq(-2, 2, by = 0.1),
                   kmeans_k = NA,
                   cluster_rows = T,
                   cutree_row = clusters,
                   cluster_cols = F,
                   #cutree_cols = 4,
                   gaps_col = gaps,              
                   legend = TRUE,
                   show_rownames = T,
                   #cellwidth = 25,
                   cellheight = 9,
                   treeheight_col = 0,
                   treeheight_row = 50,
                   border_color = 'NA',
                   fontsize = 8,
                   scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_200.pdf'), pointsize = 10, height = 35)
phm_200 <- pheatmap(data_subset_200,
                    color = colorRampPalette(color_set)(41),
                    breaks = seq(-2, 2, by = 0.1),
                    kmeans_k = NA,
                    cluster_rows = T,
                    cutree_row = clusters,
                    cluster_cols = F,
                    #cutree_cols = 4,
                    gaps_col = gaps,              
                    legend = TRUE,
                    show_rownames = T,
                    #cellwidth = 25,
                    cellheight = 9,
                    treeheight_col = 0,
                    treeheight_row = 50,
                    border_color = 'NA',
                    fontsize = 8,
                    scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_400.pdf'), pointsize = 10, height = 60)
phm_full <- pheatmap(data_subset_400,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 9)
dev.off()


### GOURPED
pdf(file = paste0(heatmap_dir, '/All_DEGs_grouped.pdf'), pointsize = 10)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(40),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 20,
                     #cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_grouped.pdf'), pointsize = 10, height = nrow(data_grouped_norm)/7)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(40),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 20,
                     cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_padj05_grouped.pdf'), pointsize = 10, height= nrow(data_grouped_norm_padj)/6)
phm_full <- pheatmap(data_grouped_norm_padj,
                     color = colorRampPalette(c(color_set))(40),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()

#Top 200 grouped labelled
pdf(file = paste0(heatmap_dir,'/Top_200_grouped.pdf'), pointsize = 10, height= nrow(data_subset_200_grouped)/6)
phm_full <- pheatmap(data_subset_200_grouped,
                     color = colorRampPalette(c(color_set))(40),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 15,
                     cellheight = 9,
                     scale = 'row')
dev.off()


#####--------




##### ----- EXPORT RESULTS PER COMPARISON ----- ##### 
fgseaResTidy_curated_full <- data.frame()
KEGG_data_full <- data.frame()
GO_data_full <- data.frame()


for (comp_name in (colnames(cmplx.cont))) {
    dds_result <- results(dds_run, pAdjustMethod = 'BH', independentFiltering = T, cooksCutoff = T, 
                                contrast = cmplx.cont[, comp_name])
    dds_result$Contrast <- paste(comp_name)
    
  # Build significant gene table P - VALUE 0.05 and extract list of sorted DE genes.
  sig_res <- dds_result %>%
    data.frame() %>%
    rownames_to_column(var='gene') %>%  
    as_tibble %>%
    filter(pvalue < 0.05) %>%
    filter(baseMean > 1) %>%
    arrange(pvalue) 
  
  # Build significant gene table P - ADJUSTED 0.05 and extract list of sorted DE genes.
  sig_res_adj <- sig_res %>%
    filter(padj < 0.05) %>%
    arrange(padj) 
  
  sorted_DEGenes <- sig_res$gene
  sorted_DEGenes_padj <- sig_res_adj$gene
  
  # Export significant gene count data
  name_list <- c('gene', (paste0(dds_run$group, "_", dds_run$sample)))
  name_list
  
  sig_export <- DS_norm_counts %>%
    data.frame() %>%
    rownames_to_column(var = 'gene') %>%
    as_tibble() %>%
    `colnames<-`(name_list) %>%
    semi_join(sig_res) %>%
    dplyr::slice(match(sorted_DEGenes, gene)) 
  
  sig_export_adj <- DS_norm_counts %>%
    data.frame() %>%
    rownames_to_column(var = 'gene') %>%
    as_tibble() %>%
    `colnames<-`(name_list) %>%
    semi_join(sig_res) %>%
    dplyr::slice(match(sorted_DEGenes_padj, gene))
  
  
  #Selecting only columns specific for comparison
  selected_samples <- cmplx.cont %>% rownames_to_column('group_names') %>% 
    dplyr::select(c(group_names,comp_name)) 
  
  colnames(selected_samples) <- c('group_names','comparison_name')
  
  selected_samples_name <- selected_samples %>%
    dplyr::filter(comparison_name != 0) %>% pull('group_names') 
  
  selected_samples_name <-  gsub("group","",selected_samples_name)
 
  # Filter the data frame based on the index
  sig_export <- sig_export %>% dplyr::select(c(1, contains(selected_samples_name)))
  
  # Setup Order of Counts
  reordered_index <- c(grep("CTRL_m_33", names(sig_export), ignore.case = T),
                       grep("CTRL_m_34", names(sig_export), ignore.case = T),
                       grep("CTRL_f_33", names(sig_export), ignore.case = T),
                       grep("CTRL_f_34", names(sig_export), ignore.case = T),
                       grep("AD_m_33", names(sig_export), ignore.case = T),
                       grep("AD_m_34", names(sig_export), ignore.case = T),
                       grep("AD_m_44", names(sig_export), ignore.case = T),
                       grep("AD_f_33", names(sig_export), ignore.case = T),
                       grep("AD_f_34", names(sig_export), ignore.case = T),
                       grep("AD_f_44", names(sig_export), ignore.case = T))
  
  sig_export <- sig_export %>%
    dplyr::select(c(1, reordered_index))
  
  sig_export_adj <- sig_export_adj %>%
    dplyr::select(c(1, reordered_index))
  
  # Column Names
  check_cols <- c(colnames(sig_export))
  check_cols
  
  # Normalized counts
  check_counts <- DS_norm_counts %>%
    as.data.frame() %>%
    rownames_to_column('gene') %>%
    dplyr::select(c(1, reordered_index)) %>%
    as_tibble %>% 
    `colnames<-`(check_cols)
  
  # Results
  check_res <- dds_result %>%
    as.data.frame() %>%
    rownames_to_column('gene') 
  
  # Output files
  dir.create(paste0('results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/', file_prefix,'/', comp_name))
  
  results_dir <- paste0('results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/', file_prefix, '/', comp_name ,'/')

  #Statistics
  # write.xlsx(sig_res, file = paste0(results_dir,comp_name,'_DEGene_statistics_pval05.xlsx'), overwrite = T)
  # write.xlsx(sig_res_adj, file = paste0(results_dir,comp_name, '_DEGene_statistics_padj05.xlsx'), overwrite = T)
  # write.xlsx(check_res, file = paste0(results_dir,comp_name, '_Statistics.xlsx'), overwrite = T)
  #Counts
  # write.xlsx(sig_export, file = paste0(results_dir,comp_name,'_DEGene_counts_pval05.xlsx'), overwrite = T)
  # write.xlsx(sig_export_adj, file = paste0(results_dir,comp_name,'_DEGene_counts_padj05.xlsx'), overwrite = T)
  # write.xlsx(check_counts, file = paste0(results_dir,comp_name,'_DS_counts.xlsx'), overwrite = T)
  
  #Setting up organism data base for enrichment analysis
  organism_name <-  'hsa' #human: hsa
  OrgDb_name    <-  "org.Hs.eg.db" #mouse: 'org.Mm.eg.db'

  ##### --------- FGSEA analysis ####
  DEgenes <- sig_res %>% filter(pvalue < 0.05) %>%  pull('gene')
  DEgenes <- mapIds(org.Hs.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')
  
  DEG_Log2FC <-  sig_res %>% filter(pvalue < 0.05) %>%  pull('log2FoldChange')
  
  names(DEG_Log2FC) = DEgenes
  DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
  DEG_Log2FC
  
  #####CURATED pathways
  pathways.all.curated <- gmtPathways("/Users/kiliankleemann/Dropbox/Neutrophils Projects/Human_SK-5DC3/GSEA_database/c2.all.v2022.1.Hs.entrez.gmt") 
  fgseaRes <- fgsea(pathways=pathways.all.curated, stats=DEG_Log2FC)
  fgseaResTidy_curated <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    filter(pval<0.1)
  
  dir.create(paste0(results_dir,'pathway_analysis/'))
  dir.create(paste0(results_dir,'pathway_analysis/FGSEA/'))
  write.xlsx(fgseaResTidy_curated, paste0(results_dir,'pathway_analysis/FGSEA/',comp_name,'_curated_pathways.xlsx'))
  
  #Masterlist of KEGG results
  # fgseaResTidy_curated$comparison <- paste(comp_name)
  # fgseaResTidy_curated_full <- rbind(fgseaResTidy_curated_full,fgseaResTidy_curated)
  # 
  
  #Run KEGG 
  kegg_result <- gseKEGG(geneList = DEG_Log2FC,
                         organism     = organism_name,
                         pAdjustMethod = "BH",
                         nPermSimple = 1000,
                         pvalueCutoff = 0.5)
  
  kegg_result <- setReadable(kegg_result,OrgDb= 'org.Hs.eg.db', keyType = 'ENTREZID')
  KEGG_data <- as.data.frame(kegg_result) 
  KEGG_data$log_pval <- log(KEGG_data$pvalue,10)
  
  dir.create(paste0(results_dir,'pathway_analysis/KEGG/'))
  write.xlsx(KEGG_data, paste0(results_dir,'pathway_analysis/KEGG/',comp_name,'_pathways.xlsx'))
  
  #Masterlist of KEGG results
  # KEGG_data$comparison <- paste(comp_name)
  # KEGG_data_full <- rbind(KEGG_data_full,KEGG_data)
  # 
  #Run GO 
  go_result <- gseGO(geneList = DEG_Log2FC,
                     OrgDb     = OrgDb_name,
                     pAdjustMethod = "BH",
                     nPermSimple = 1000,
                     pvalueCutoff = 0.1)
  
  go_result <- setReadable(go_result, OrgDb= OrgDb_name, keyType = 'ENTREZID')
  GO_data <- as.data.frame(go_result) 
  GO_data$log_pval <- log(GO_data$pvalue,10)
  
  dir.create(paste0(results_dir,'pathway_analysis/GO/'))
  write.xlsx(GO_data, paste0(results_dir,'pathway_analysis/GO/',comp_name,'_pathways.xlsx'))
  
  #Masterlist of GO results
  # GO_data$comparison <- paste(comp_name)
  # GO_data_full <- rbind(GO_data_full,GO_data)
  # 
  ##### ------ BARPLOTS 
  #Settings for barplots
  # lineWidth = 1
  # pointSize = 20
  # dir.create(paste0('plots/barplots/Microglia/Aggregate_Sample/'))
  # dir.create(paste0('plots/barplots/Microglia/Aggregate_Sample/',file_prefix))
  # dir.create(paste0('plots/barplots/Microglia/Aggregate_Sample/',file_prefix,'/',comp_name))
  # 
  # barplot_dir <- paste0('plots/barplots/Microglia/Aggregate_Sample/',file_prefix,'/', comp_name, '/')
  # 
  # for (i in goi) { 
  #   goi_data3 <- goi_data2 %>% filter(gene == i) %>% pivot_longer(
  #     cols = -(1:1),
  #     values_to = c("DS_Counts"),
  #     names_to = c("condition", "replicate"),
  #     names_sep  = "-")
  #   
  #   goi_data3$condition <- factor(goi_data3$condition, levels=unique(goi_data3$condition))
  #   
  #   comparisons <- compare_means(
  #     data = goi_data3,
  #     formula = DS_Counts ~ condition,
  #     method = "t.test",
  #     p.adjust.method = "BH") %>% filter(p<0.05)
  #   
  #   comparison_list_sign <- comparisons %>% mutate(comparison_list = map2(group1, group2,c)) %>% pull(comparison_list)
  #   
  #   barplot <- ggplot(goi_data3, aes(x = condition, y = DS_Counts, fill = condition)) + 
  #     geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
  #     geom_errorbar(stat = 'summary', position = 'dodge', width = 0.4) +
  #     geom_point(aes(x = condition), position = position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width=0.9)) +
  #     expand_limits(x = 0, y = 0) +
  #     theme(panel.spacing = unit(1, "lines")) +
  #     labs(x = NULL, y = c("Normalized Counts ±SE")) +
  #     ggtitle(paste0(i))+
  #     theme(
  #       text = element_text(size = pointSize, colour = "black"),
  #       rect = element_blank(),
  #       line = element_line(linewidth = lineWidth, colour = "black"),
  #       plot.title  = element_text(color="black", size=40, face="bold.italic"),
  #       axis.title  = element_text(size = pointSize, colour = "black"),
  #       axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
  #       axis.text.y  = element_text(size = pointSize , colour = "black"),
  #       legend.position = "none",
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       legend.title = element_blank(),
  #       legend.text = element_text(size = pointSize , colour = "black"),
  #       axis.line = element_line(linewidth = lineWidth, colour = "black"),
  #       plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
  #     stat_compare_means(comparisons = comparison_list_sign, method = 't.test', label = "p.signif" ) +
  #     scale_fill_brewer(palette = "Paired") + scale_y_continuous(expand = c(0, 0, .05, 0))
  #   
  #   pdf(file = paste0(barplot_dir,  i,'.pdf'), pointsize = 10, width = 5, height= 10)
  #   print(barplot)
  #   dev.off()
  # }
  # 
  # ##### HEATMAP 
  # # Define a normalization function to calculate Z scores.
  # cal_z_score <- function(x){
  #   (x - mean(x)) / sd(x)
  # }
  # 
  # data_final <- sig_export %>% column_to_rownames('gene')
  # #data_final_adj <- sig_export_adj %>% column_to_rownames('gene')
  # data_subset_50 <- data_final %>% head(n = 50)
  # data_subset_200 <- data_final %>% head(n = 200)
  # data_subset_400 <- data_final %>% head(n = 400)
  # 
  # 
  # #Grouped data final
  # # data_grouped <- data.frame(
  # #   CTRL_m_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("CTRL_m_33")))),
  # #   CTRL_m_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("CTRL_m_34")))),
  # #   CTRL_f_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("CTRL_f_33")))),
  # #   CTRL_f_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("CTRL_f_34")))),
  # #   AD_m_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_m_33")))),
  # #   AD_m_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_m_34")))),
  # #   AD_f_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_33")))),
  # #   AD_f_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_34")))))
  # 
  # 
  # # Normalize your data according to Z score using cal_z_score function.
  # data_norm <- t(apply(data_final, 1, cal_z_score))
  # #data_final_adj_norm <- t(apply(data_final_adj, 1, cal_z_score))
  # data_subset_50 <- t(apply(data_subset_50, 1, cal_z_score))
  # data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
  # data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
  # #data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))
  # 
  # #Export directory
  # dir.create(paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix,'/', comp_name))
  # 
  # heatmap_dir <- paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix,'/',  comp_name)
  # 
  # #Setting up Heatmap parameters
  # gaps = FALSE
  # clusters = NA
  # color_set <- c("navy", "white", "firebrick3")
  # all_degs_export_height <- nrow(data_final)
  # padj_degs_export_height <- nrow(data_final_adj)
  # ## ------------- DENDROGRAM AND CLUSTER MANIPULATION (ADVANCED) -------------##
  # # Create heatmap
  # # pdf(file = paste0(heatmap_dir, '/All_DEGs_padj05.pdf'), pointsize = 10, height = padj_degs_export_height/4)
  # # phm_full <- pheatmap(data_final_adj_norm,
  # #                      color = colorRampPalette(color_set)(41),
  # #                      breaks = seq(-2, 2, by = 0.1),
  # #                      kmeans_k = NA,
  # #                      cluster_rows = T,
  # #                      cutree_row = clusters,
  # #                      cluster_cols = F,
  # #                      #cutree_cols = 4,
  # #                      gaps_col = gaps,             
  # #                      legend = TRUE,
  # #                      show_rownames = T,
  # #                      border_color = 'NA',
  # #                      fontsize = 8,
  # #                      treeheight_col = 0,
  # #                      treeheight_row = 10,
  # #                      cellheight = 9,
  # #                      #cellwidth = 25,
  # #                      scale = 'row'
  # # )
  # # dev.off()
  # 
  # pdf(file = paste0(heatmap_dir,'/All_DEGs_pval05.pdf'), pointsize = 10)
  # phm_full <- pheatmap(data_norm,
  #                      color = colorRampPalette(color_set)(41),
  #                      breaks = seq(-2, 2, by = 0.1),
  #                      kmeans_k = NA,
  #                      cluster_rows = T,
  #                      cutree_row = clusters,
  #                      cluster_cols = F,
  #                      #cutree_cols = 4,
  #                      gaps_col = gaps,             
  #                      legend = TRUE,
  #                      show_rownames = F,
  #                      border_color = 'NA',
  #                      fontsize = 8,
  #                      treeheight_col = 0,
  #                      treeheight_row = 10,
  #                      #cellwidth = 25,
  #                      scale = 'row'
  #                      #cellheight = 3,
  # )
  # dev.off()
  # 
  # 
  # pdf(file = paste0(heatmap_dir, '/All_DEGs_labelled_pvall05.pdf'), pointsize = 10, height = all_degs_export_height/4)
  # phm_full <- pheatmap(data_norm,
  #                      color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
  #                      breaks = seq(-2, 2, by = 0.1),
  #                      kmeans_k = NA,
  #                      cluster_rows = T,
  #                      cutree_row = clusters,
  #                      cluster_cols = F,
  #                      #cutree_cols = 4,
  #                      gaps_col = gaps,             
  #                      legend = TRUE,
  #                      show_rownames = T,
  #                      border_color = 'NA',
  #                      fontsize = 8,
  #                      treeheight_col = 0,
  #                      treeheight_row = 10,
  #                      #cellwidth = 25,
  #                      scale = 'row',
  #                      cellheight = 8
  # )
  # dev.off()
  # 
  # pdf(file = paste0(heatmap_dir, '/Top_50.pdf'), pointsize = 10, height = 10)
  # phm_50 <- pheatmap(data_subset_50,
  #                    color = colorRampPalette(color_set)(41),
  #                    breaks = seq(-2, 2, by = 0.1),
  #                    kmeans_k = NA,
  #                    cluster_rows = T,
  #                    cutree_row = clusters,
  #                    cluster_cols = F,
  #                    #cutree_cols = 4,
  #                    gaps_col = gaps,              
  #                    legend = TRUE,
  #                    show_rownames = T,
  #                    #cellwidth = 25,
  #                    cellheight = 9,
  #                    treeheight_col = 0,
  #                    treeheight_row = 50,
  #                    border_color = 'NA',
  #                    fontsize = 8,
  #                    scale = 'row')
  # dev.off()
  # 
  # pdf(file = paste0(heatmap_dir, '/Top_200.pdf'), pointsize = 10, height = 35)
  # phm_200 <- pheatmap(data_subset_200,
  #                     color = colorRampPalette(color_set)(41),
  #                     breaks = seq(-2, 2, by = 0.1),
  #                     kmeans_k = NA,
  #                     cluster_rows = T,
  #                     cutree_row = clusters,
  #                     cluster_cols = F,
  #                     #cutree_cols = 4,
  #                     gaps_col = gaps,              
  #                     legend = TRUE,
  #                     show_rownames = T,
  #                     #cellwidth = 25,
  #                     cellheight = 9,
  #                     treeheight_col = 0,
  #                     treeheight_row = 50,
  #                     border_color = 'NA',
  #                     fontsize = 8,
  #                     scale = 'row')
  # dev.off()
  # 
  # pdf(file = paste0(heatmap_dir, '/Top_400.pdf'), pointsize = 10, height = 60)
  # phm_full <- pheatmap(data_subset_400,
  #                      color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
  #                      breaks = seq(-2, 2, by = 0.1),
  #                      kmeans_k = NA,
  #                      cluster_rows = T,
  #                      cutree_row = clusters,
  #                      cluster_cols = F,
  #                      #cutree_cols = 3,
  #                      gaps_col = gaps,             
  #                      legend = TRUE,
  #                      show_rownames = T,
  #                      border_color = 'NA',
  #                      fontsize = 8,
  #                      treeheight_col = 0,
  #                      treeheight_row = 0,
  #                      #cellwidth = 25,
  #                      scale = 'row',
  #                      cellheight = 9)
  # dev.off()
  # 
  # 
  # # ### GOURPED
  # # pdf(file = paste0(heatmap_dir, '/All_DEGs_grouped.pdf'), pointsize = 10)
  # # phm_full <- pheatmap(data_grouped_norm,
  # #                      color = colorRampPalette(c(color_set))(41),
  # #                      breaks = seq(-2, 2, by = 0.1),
  # #                      kmeans_k = NA,
  # #                      cluster_rows = T,
  # #                      cutree_row = clusters,
  # #                      cluster_cols = F,
  # #                      #cutree_cols = 3,
  # #                      gaps_col = c(),             
  # #                      legend = TRUE,
  # #                      show_rownames = F,
  # #                      border_color = 'NA',
  # #                      fontsize = 8,
  # #                      treeheight_col = 0,
  # #                      treeheight_row = 0,
  # #                      #cellwidth = 20,
  # #                      #cellheight = 9,
  # #                      scale = 'row')
  # # dev.off()
  # # 
  # # 
  # # #Grouped labelled
  # # pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_grouped.pdf'), pointsize = 10, height = all_degs_export_height/4)
  # # phm_full <- pheatmap(data_grouped_norm,
  # #                      color = colorRampPalette(c(color_set))(41),
  # #                      breaks = seq(-2, 2, by = 0.1),
  # #                      kmeans_k = NA,
  # #                      cluster_rows = T,
  # #                      cutree_row = clusters,
  # #                      cluster_cols = F,
  # #                      #cutree_cols = 3,
  # #                      gaps_col = c(),             
  # #                      legend = TRUE,
  # #                      show_rownames = T,
  # #                      border_color = 'NA',
  # #                      fontsize = 8,
  # #                      treeheight_col = 0,
  # #                      treeheight_row = 0,
  # #                      cellwidth = 20,
  # #                      cellheight = 9,
  # #                      scale = 'row')
  # # dev.off()
}
















#### -------- ANALYSIS FEMALE AD SAMPLES - SIMPLE DESIGN - LRT - APOE COMPARISONS------- ####
#Setting experiment conditions !!!
file_prefix <- 'Female_AD_samples_simple_LRT_KK'
data <- data_adjusted_mg
experiment <- sub_meta %>% filter(sex == 'f' & diagnosis == 'AD')
experiment


data <-  as.data.frame(data) %>% dplyr::select(experiment$sample)
experiment

#Designs 
dds <- DESeqDataSetFromMatrix(data, colData = experiment, design = ~ 0 + age + pmd_min + group)

## Run DESeq analysis to gather differential expression results
#Run DESeq (LRT - for complex comparisons) 
dds_run <- DESeq(dds, test = 'LRT', reduced = ~ 1)
resultsNames(dds_run)

dds_result <- results(dds_run, pAdjustMethod = 'BH', independentFiltering = T, cooksCutoff = T)

# View Total number of normalized counts per sample
dds_run <- estimateSizeFactors(dds_run)
DS_norm_counts <- counts(dds_run, normalized = TRUE)
colSums(DS_norm_counts)

#### -------- OUTPUT TABLES & RESULTS 
# Build significant gene table P - VALUE 0.05 and extract list of sorted DE genes.
sig_res <- dds_result %>%
  data.frame() %>%
  rownames_to_column(var='gene') %>%  
  as_tibble %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) 

# Build significant gene table P - ADJUSTED 0.05 and extract list of sorted DE genes.
sig_res_adj <- sig_res %>%
  filter(padj < 0.05) %>%
  arrange(padj) 

sorted_DEGenes <- sig_res$gene %>% unique()
sorted_DEGenes_padj <- sig_res_adj$gene %>% unique()

# Export significant gene count data
name_list <- c('gene', (paste0(dds_run$group, "_", dds_run$sample)))
name_list

sig_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes, gene)) 

sig_export_adj <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes_padj, gene))

# Setup Order of Counts
reordered_index <- c(grep("CTRL_m_33", names(sig_export), ignore.case = T),
                     grep("CTRL_m_34", names(sig_export), ignore.case = T),
                     grep("CTRL_f_33", names(sig_export), ignore.case = T),
                     grep("CTRL_f_34", names(sig_export), ignore.case = T),
                     grep("AD_m_33", names(sig_export), ignore.case = T),
                     grep("AD_m_34", names(sig_export), ignore.case = T),
                     grep("AD_m_44", names(sig_export), ignore.case = T),
                     grep("AD_f_33", names(sig_export), ignore.case = T),
                     grep("AD_f_34", names(sig_export), ignore.case = T),
                     grep("AD_f_44", names(sig_export), ignore.case = T))

sig_export <- sig_export %>%
  dplyr::select(c(1, reordered_index))

sig_export_adj <- sig_export_adj %>%
  dplyr::select(c(1, reordered_index))

# Column Names
check_cols <- c(colnames(sig_export))
check_cols

# Normalized counts
check_counts <- DS_norm_counts %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  dplyr::select(c(1, reordered_index)) %>%
  as_tibble %>% 
  `colnames<-`(check_cols) #%>%
#filter(!str_detect(gene, unwanted_genes))


# Results
check_res <- dds_result 
file_prefix
# Output files
dir.create(paste0('results/Microglia/Aggregate_Sample/'))
dir.create(paste0('results/Microglia/Aggregate_Sample/', file_prefix))

#Statistics
write.xlsx(sig_res, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DEGene_statistics_pval05.xlsx'), overwrite = T)
write.xlsx(sig_res_adj, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DEGene_statistics_padj05.xlsx'), overwrite = T)
write.xlsx(check_res, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/Statistics.xlsx'), overwrite = T)
#Counts
write.xlsx(sig_export, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DEGene_counts_pval05.xlsx'), overwrite = T)
write.xlsx(sig_export_adj, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DEGene_counts_padj05.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DS_counts.xlsx'), overwrite = T)



##### HEATMAP 
# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_final <- sig_export %>% column_to_rownames('gene')
data_final_adj <- sig_export_adj %>% column_to_rownames('gene')

data_subset_50 <- data_final %>% head(n = 50)
data_subset_200 <- data_final %>% head(n = 200)
data_subset_400 <- data_final %>% head(n = 400)


#Grouped data final
data_grouped <- data.frame(
  AD_f_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_34")))),
  AD_m_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_44")))))

data_subset_200_grouped <- data.frame(
  AD_f_33 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_f_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_f_34")))),
  AD_f_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_f_44")))))


data_grouped_padj <- data.frame(
  AD_f_33 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_f_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_f_34")))),
  AD_f_44 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_f_44")))))


# Normalize your data according to Z score using cal_z_score function.
data_norm <- t(apply(data_final, 1, cal_z_score))
data_final_adj_norm <- t(apply(data_final_adj, 1, cal_z_score))
data_subset_50 <- t(apply(data_subset_50, 1, cal_z_score))
data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))
data_grouped_norm_padj <- t(apply(data_grouped_padj, 1, cal_z_score))

#Export directory
dir.create(paste0('plots/heatmaps/Microglia/'))
dir.create(paste0('plots/heatmaps/Microglia/',file_prefix))

heatmap_dir <- paste0('plots/heatmaps/Microglia/',file_prefix)

#Setting up Heatmap parameters
gaps = c(5,8)
clusters = 8
color_set <- c("navy", "white", "firebrick3")
all_degs_export_height <- nrow(data_final)
padj_degs_export_height <- nrow(data_final_adj)
## ------------- DENDROGRAM AND CLUSTER MANIPULATION (ADVANCED) -------------##
# Create heatmap
pdf(file = paste0(heatmap_dir, '/All_DEGs_padj05.pdf'), pointsize = 10, height = padj_degs_export_height/4)
phm_full <- pheatmap(data_final_adj_norm,
                     color = colorRampPalette(color_set)(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     cellheight = 9,
                     #cellwidth = 25,
                     scale = 'row'
)
dev.off()



pdf(file = paste0(heatmap_dir,'/All_DEGs_pval05.pdf'), pointsize = 10)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(color_set)(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row'
                     #cellheight = 3,
)
dev.off()


pdf(file = paste0(heatmap_dir, '/All_DEGs_labelled_pval05.pdf'), pointsize = 10, height = all_degs_export_height/4)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 8
)
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_50.pdf'), pointsize = 10, height = 10)
phm_50 <- pheatmap(data_subset_50,
                   color = colorRampPalette(color_set)(41),
                   breaks = seq(-2, 2, by = 0.1),
                   kmeans_k = NA,
                   cluster_rows = T,
                   cutree_row = clusters,
                   cluster_cols = F,
                   #cutree_cols = 4,
                   gaps_col = gaps,              
                   legend = TRUE,
                   show_rownames = T,
                   #cellwidth = 25,
                   cellheight = 9,
                   treeheight_col = 0,
                   treeheight_row = 50,
                   border_color = 'NA',
                   fontsize = 8,
                   scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_200.pdf'), pointsize = 10, height = 35)
phm_200 <- pheatmap(data_subset_200,
                    color = colorRampPalette(color_set)(41),
                    breaks = seq(-2, 2, by = 0.1),
                    kmeans_k = NA,
                    cluster_rows = T,
                    cutree_row = clusters,
                    cluster_cols = F,
                    #cutree_cols = 4,
                    gaps_col = gaps,              
                    legend = TRUE,
                    show_rownames = T,
                    #cellwidth = 25,
                    cellheight = 9,
                    treeheight_col = 0,
                    treeheight_row = 50,
                    border_color = 'NA',
                    fontsize = 8,
                    scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_400.pdf'), pointsize = 10, height = 60)
phm_full <- pheatmap(data_subset_400,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 9)
dev.off()


### GOURPED
pdf(file = paste0(heatmap_dir, '/All_DEGs_grouped.pdf'), pointsize = 10)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 20,
                     #cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_grouped.pdf'), pointsize = 10, height = nrow(data_grouped_norm)/7)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 20,
                     cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_padj05_grouped.pdf'), pointsize = 10, height= nrow(data_grouped_norm_padj)/6)
phm_full <- pheatmap(data_grouped_norm_padj,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir,'/Top_200_grouped.pdf'), pointsize = 10, height= nrow(data_subset_200_grouped)/6)
phm_full <- pheatmap(data_subset_200_grouped,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()















#### -------- ANALYSIS MALE AD SAMPLES - SIMPLE DESIGN - LRT - APOE COMPARISONS------- ####
#Setting experiment conditions 
file_prefix <- 'Male_AD_samples_simple_LRT_KK'
data <- data_adjusted
experiment <- sub_meta %>% filter(sex == 'm' & diagnosis == 'AD')
experiment

data <- data_adjusted 

data <-  as.data.frame(data) %>% dplyr::select(experiment$sample)
experiment

#Designs 
dds <- DESeqDataSetFromMatrix(data, colData = experiment, design = ~ 0 + age + pmd_min + group)

## Run DESeq analysis to gather differential expression results
#Run DESeq (LRT - for complex comparisons) 
dds_run <- DESeq(dds, test = 'LRT', reduced = ~ 1)
resultsNames(dds_run)

dds_result <- results(dds_run, pAdjustMethod = 'BH', independentFiltering = T, cooksCutoff = T)



# View Total number of normalized counts per sample
dds_run <- estimateSizeFactors(dds_run)
DS_norm_counts <- counts(dds_run, normalized = TRUE)
colSums(DS_norm_counts)

#### -------- OUTPUT TABLES & RESULTS 
# Build significant gene table P - VALUE 0.05 and extract list of sorted DE genes.
sig_res <- dds_result %>%
  data.frame() %>%
  rownames_to_column(var='gene') %>%  
  as_tibble %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) 

# Build significant gene table P - ADJUSTED 0.05 and extract list of sorted DE genes.
sig_res_adj <- sig_res %>%
  filter(padj < 0.05) %>%
  arrange(padj) 

sorted_DEGenes <- sig_res$gene %>% unique()
sorted_DEGenes_padj <- sig_res_adj$gene %>% unique()

# Export significant gene count data
name_list <- c('gene', (paste0(dds_run$group, "_", dds_run$sample)))
name_list

sig_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes, gene)) 

sig_export_adj <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes_padj, gene))

# Setup Order of Counts
reordered_index <- c(grep("CTRL_m_33", names(sig_export), ignore.case = T),
                     grep("CTRL_m_34", names(sig_export), ignore.case = T),
                     grep("CTRL_f_33", names(sig_export), ignore.case = T),
                     grep("CTRL_f_34", names(sig_export), ignore.case = T),
                     grep("AD_m_33", names(sig_export), ignore.case = T),
                     grep("AD_m_34", names(sig_export), ignore.case = T),
                     grep("AD_m_44", names(sig_export), ignore.case = T),
                     grep("AD_f_33", names(sig_export), ignore.case = T),
                     grep("AD_f_34", names(sig_export), ignore.case = T),
                     grep("AD_f_44", names(sig_export), ignore.case = T))

sig_export <- sig_export %>%
  dplyr::select(c(1, reordered_index))

sig_export_adj <- sig_export_adj %>%
  dplyr::select(c(1, reordered_index))

# Column Names
check_cols <- c(colnames(sig_export))
check_cols

# Normalized counts
check_counts <- DS_norm_counts %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  dplyr::select(c(1, reordered_index)) %>%
  as_tibble %>% 
  `colnames<-`(check_cols) #%>%
#filter(!str_detect(gene, unwanted_genes))


# Results
check_res <- dds_result 
file_prefix
# Output files
dir.create(paste0('results/Microglia/Aggregate_Sample/'))
dir.create(paste0('results/Microglia/Aggregate_Sample/', file_prefix))

#Statistics
write.xlsx(sig_res, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DEGene_statistics_pval05.xlsx'), overwrite = T)
write.xlsx(sig_res_adj, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DEGene_statistics_padj05.xlsx'), overwrite = T)
write.xlsx(check_res, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/Statistics.xlsx'), overwrite = T)
#Counts
write.xlsx(sig_export, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DEGene_counts_pval05.xlsx'), overwrite = T)
write.xlsx(sig_export_adj, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DEGene_counts_padj05.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DS_counts.xlsx'), overwrite = T)



##### HEATMAP 
# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_final <- sig_export %>% column_to_rownames('gene')
data_final_adj <- sig_export_adj %>% column_to_rownames('gene')

data_subset_50 <- data_final %>% head(n = 50)
data_subset_200 <- data_final %>% head(n = 200)
data_subset_400 <- data_final %>% head(n = 400)


#Grouped data final
data_grouped <- data.frame(
  AD_f_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_m_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_m_34")))),
  AD_m_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_m_44")))))

data_subset_200_grouped <- data.frame(
  AD_f_33 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_m_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_m_34")))),
  AD_f_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_m_44")))))


data_grouped_padj <- data.frame(
  AD_f_33 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_m_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_m_34")))),
  AD_f_44 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_m_44")))))


# Normalize your data according to Z score using cal_z_score function.
data_norm <- t(apply(data_final, 1, cal_z_score))
data_final_adj_norm <- t(apply(data_final_adj, 1, cal_z_score))
data_subset_50 <- t(apply(data_subset_50, 1, cal_z_score))
data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))
data_grouped_norm_padj <- t(apply(data_grouped_padj, 1, cal_z_score))

#Export directory
dir.create(paste0('plots/heatmaps/Microglia/'))
dir.create(paste0('plots/heatmaps/Microglia/',file_prefix))

heatmap_dir <- paste0('plots/heatmaps/Microglia/',file_prefix)

#Setting up Heatmap parameters
gaps = c(5,10)
clusters = 8
color_set <- c("navy", "white", "firebrick3")
all_degs_export_height <- nrow(data_final)
padj_degs_export_height <- nrow(data_final_adj)
## ------------- DENDROGRAM AND CLUSTER MANIPULATION (ADVANCED) -------------##
# Create heatmap
pdf(file = paste0(heatmap_dir, '/All_DEGs_padj05.pdf'), pointsize = 10, height = padj_degs_export_height/4)
phm_full <- pheatmap(data_final_adj_norm,
                     color = colorRampPalette(color_set)(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     cellheight = 9,
                     #cellwidth = 25,
                     scale = 'row'
)
dev.off()



pdf(file = paste0(heatmap_dir,'/All_DEGs_pval05.pdf'), pointsize = 10)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(color_set)(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row'
                     #cellheight = 3,
)
dev.off()


pdf(file = paste0(heatmap_dir, '/All_DEGs_labelled_pval05.pdf'), pointsize = 10, height = all_degs_export_height/4)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 8
)
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_50.pdf'), pointsize = 10, height = 10)
phm_50 <- pheatmap(data_subset_50,
                   color = colorRampPalette(color_set)(41),
                   breaks = seq(-2, 2, by = 0.1),
                   kmeans_k = NA,
                   cluster_rows = T,
                   cutree_row = clusters,
                   cluster_cols = F,
                   #cutree_cols = 4,
                   gaps_col = gaps,              
                   legend = TRUE,
                   show_rownames = T,
                   #cellwidth = 25,
                   cellheight = 9,
                   treeheight_col = 0,
                   treeheight_row = 50,
                   border_color = 'NA',
                   fontsize = 8,
                   scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_200.pdf'), pointsize = 10, height = 35)
phm_200 <- pheatmap(data_subset_200,
                    color = colorRampPalette(color_set)(41),
                    breaks = seq(-2, 2, by = 0.1),
                    kmeans_k = NA,
                    cluster_rows = T,
                    cutree_row = clusters,
                    cluster_cols = F,
                    #cutree_cols = 4,
                    gaps_col = gaps,              
                    legend = TRUE,
                    show_rownames = T,
                    #cellwidth = 25,
                    cellheight = 9,
                    treeheight_col = 0,
                    treeheight_row = 50,
                    border_color = 'NA',
                    fontsize = 8,
                    scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_400.pdf'), pointsize = 10, height = 60)
phm_full <- pheatmap(data_subset_400,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 9)
dev.off()


### GOURPED
pdf(file = paste0(heatmap_dir, '/All_DEGs_grouped.pdf'), pointsize = 10)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 20,
                     #cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_grouped.pdf'), pointsize = 10, height = nrow(data_grouped_norm)/7)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 20,
                     cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_padj05_grouped.pdf'), pointsize = 10, height= nrow(data_grouped_norm_padj)/6)
phm_full <- pheatmap(data_grouped_norm_padj,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir,'/Top_200_grouped.pdf'), pointsize = 10, height= nrow(data_subset_200_grouped)/6)
phm_full <- pheatmap(data_subset_200_grouped,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()




















#### -------- ANALYSIS FEMALE SAMPLES - SIMPLE DESIGN - LRT - APOE / DIAGNOSIS COMPARISONS------- ####
#Setting experiment conditions 
file_prefix <- 'Female_APOE33_APOE34_samples_simple_LRT_KK'
data <- data_adjusted
experiment <- sub_meta %>% filter(sex == 'f' & apoe != '44')
experiment

data <- data_adjusted 

data <-  as.data.frame(data) %>% dplyr::select(experiment$sample)
experiment

#Designs 
dds <- DESeqDataSetFromMatrix(data, colData = experiment, design = ~ 0 + age + pmd_min + apoe * diagnosis)

## Run DESeq analysis to gather differential expression results
#Run DESeq (LRT - for complex comparisons) 
dds_run <- DESeq(dds, test = 'LRT', reduced = ~ 1)
resultsNames(dds_run)

dds_result <- results(dds_run, pAdjustMethod = 'BH', independentFiltering = T, cooksCutoff = T)



# View Total number of normalized counts per sample
dds_run <- estimateSizeFactors(dds_run)
DS_norm_counts <- counts(dds_run, normalized = TRUE)
colSums(DS_norm_counts)

#### -------- OUTPUT TABLES & RESULTS 
# Build significant gene table P - VALUE 0.05 and extract list of sorted DE genes.
sig_res <- dds_result %>%
  data.frame() %>%
  rownames_to_column(var='gene') %>%  
  as_tibble %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) 

# Build significant gene table P - ADJUSTED 0.05 and extract list of sorted DE genes.
sig_res_adj <- sig_res %>%
  filter(padj < 0.05) %>%
  arrange(padj) 

sorted_DEGenes <- sig_res$gene %>% unique()
sorted_DEGenes_padj <- sig_res_adj$gene %>% unique()

# Export significant gene count data
name_list <- c('gene', (paste0(dds_run$group, "_", dds_run$sample)))
name_list

sig_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes, gene)) 

sig_export_adj <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes_padj, gene))

# Setup Order of Counts
reordered_index <- c(grep("CTRL_m_33", names(sig_export), ignore.case = T),
                     grep("CTRL_m_34", names(sig_export), ignore.case = T),
                     grep("CTRL_f_33", names(sig_export), ignore.case = T),
                     grep("CTRL_f_34", names(sig_export), ignore.case = T),
                     grep("AD_m_33", names(sig_export), ignore.case = T),
                     grep("AD_m_34", names(sig_export), ignore.case = T),
                     grep("AD_m_44", names(sig_export), ignore.case = T),
                     grep("AD_f_33", names(sig_export), ignore.case = T),
                     grep("AD_f_34", names(sig_export), ignore.case = T),
                     grep("AD_f_44", names(sig_export), ignore.case = T))

sig_export <- sig_export %>%
  dplyr::select(c(1, reordered_index))

sig_export_adj <- sig_export_adj %>%
  dplyr::select(c(1, reordered_index))

# Column Names
check_cols <- c(colnames(sig_export))
check_cols

# Normalized counts
check_counts <- DS_norm_counts %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  dplyr::select(c(1, reordered_index)) %>%
  as_tibble %>% 
  `colnames<-`(check_cols) #%>%
#filter(!str_detect(gene, unwanted_genes))


# Results
check_res <- dds_result 
file_prefix
# Output files
dir.create(paste0('results/Microglia/Aggregate_Sample/'))
dir.create(paste0('results/Microglia/Aggregate_Sample/', file_prefix))

#Statistics
write.xlsx(sig_res, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DEGene_statistics_pval05.xlsx'), overwrite = T)
write.xlsx(sig_res_adj, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DEGene_statistics_padj05.xlsx'), overwrite = T)
write.xlsx(check_res, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/Statistics.xlsx'), overwrite = T)
#Counts
write.xlsx(sig_export, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DEGene_counts_pval05.xlsx'), overwrite = T)
write.xlsx(sig_export_adj, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DEGene_counts_padj05.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0('results/Microglia/Aggregate_Sample/', file_prefix, '/DS_counts.xlsx'), overwrite = T)



##### HEATMAP 
# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_final <- sig_export %>% column_to_rownames('gene')
data_final_adj <- sig_export_adj %>% column_to_rownames('gene')

data_subset_50 <- data_final %>% head(n = 50)
data_subset_200 <- data_final %>% head(n = 200)
data_subset_400 <- data_final %>% head(n = 400)


#Grouped data final
data_grouped <- data.frame(
  AD_f_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_34")))),
  AD_m_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_44")))))

data_subset_200_grouped <- data.frame(
  AD_f_33 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_f_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_f_34")))),
  AD_f_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_f_44")))))


data_grouped_padj <- data.frame(
  AD_f_33 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_f_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_f_34")))),
  AD_f_44 = rowMeans(dplyr::select(data_final_adj, dplyr::contains(("AD_f_44")))))


# Normalize your data according to Z score using cal_z_score function.
data_norm <- t(apply(data_final, 1, cal_z_score))
data_final_adj_norm <- t(apply(data_final_adj, 1, cal_z_score))
data_subset_50 <- t(apply(data_subset_50, 1, cal_z_score))
data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))
data_grouped_norm_padj <- t(apply(data_grouped_padj, 1, cal_z_score))

#Export directory
dir.create(paste0('plots/heatmaps/Microglia/'))
dir.create(paste0('plots/heatmaps/Microglia/',file_prefix))

heatmap_dir <- paste0('plots/heatmaps/Microglia/',file_prefix)

#Setting up Heatmap parameters
gaps = c(5,8)
clusters = 8
color_set <- c("navy", "white", "firebrick3")
all_degs_export_height <- nrow(data_final)
padj_degs_export_height <- nrow(data_final_adj)
## ------------- DENDROGRAM AND CLUSTER MANIPULATION (ADVANCED) -------------##
# Create heatmap
pdf(file = paste0(heatmap_dir, '/All_DEGs_padj05.pdf'), pointsize = 10, height = padj_degs_export_height/4)
phm_full <- pheatmap(data_final_adj_norm,
                     color = colorRampPalette(color_set)(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     cellheight = 9,
                     #cellwidth = 25,
                     scale = 'row'
)
dev.off()



pdf(file = paste0(heatmap_dir,'/All_DEGs_pval05.pdf'), pointsize = 10)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(color_set)(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row'
                     #cellheight = 3,
)
dev.off()


pdf(file = paste0(heatmap_dir, '/All_DEGs_labelled_pval05.pdf'), pointsize = 10, height = all_degs_export_height/4)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 8
)
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_50.pdf'), pointsize = 10, height = 10)
phm_50 <- pheatmap(data_subset_50,
                   color = colorRampPalette(color_set)(41),
                   breaks = seq(-2, 2, by = 0.1),
                   kmeans_k = NA,
                   cluster_rows = T,
                   cutree_row = clusters,
                   cluster_cols = F,
                   #cutree_cols = 4,
                   gaps_col = gaps,              
                   legend = TRUE,
                   show_rownames = T,
                   #cellwidth = 25,
                   cellheight = 9,
                   treeheight_col = 0,
                   treeheight_row = 50,
                   border_color = 'NA',
                   fontsize = 8,
                   scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_200.pdf'), pointsize = 10, height = 35)
phm_200 <- pheatmap(data_subset_200,
                    color = colorRampPalette(color_set)(41),
                    breaks = seq(-2, 2, by = 0.1),
                    kmeans_k = NA,
                    cluster_rows = T,
                    cutree_row = clusters,
                    cluster_cols = F,
                    #cutree_cols = 4,
                    gaps_col = gaps,              
                    legend = TRUE,
                    show_rownames = T,
                    #cellwidth = 25,
                    cellheight = 9,
                    treeheight_col = 0,
                    treeheight_row = 50,
                    border_color = 'NA',
                    fontsize = 8,
                    scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_400.pdf'), pointsize = 10, height = 60)
phm_full <- pheatmap(data_subset_400,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 9)
dev.off()


### GOURPED
pdf(file = paste0(heatmap_dir, '/All_DEGs_grouped.pdf'), pointsize = 10)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 20,
                     #cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_grouped.pdf'), pointsize = 10, height = nrow(data_grouped_norm)/7)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 20,
                     cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_padj05_grouped.pdf'), pointsize = 10, height= nrow(data_grouped_norm_padj)/6)
phm_full <- pheatmap(data_grouped_norm_padj,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir,'/Top_200_grouped.pdf'), pointsize = 10, height= nrow(data_subset_200_grouped)/6)
phm_full <- pheatmap(data_subset_200_grouped,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()







#### Specific Heatmap
##### HEATMAP 
# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

spec_heatmap <- read_excel("results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/All_samples_simple_WALD_contrast_KK/(APOE34:AD_vs_CTRL)_vs_(APOE33:AD_vs_CTRL)/(APOE34:AD_vs_CTRL)_vs_(APOE33:AD_vs_CTRL)_DEGene_counts_pval05.xlsx")

data_final <- spec_heatmap %>% column_to_rownames('gene')
data_final_adj <- sig_export_adj %>% column_to_rownames('gene')

data_subset_50 <- data_final %>% head(n = 50)
data_subset_200 <- data_final %>% head(n = 200)
data_subset_400 <- data_final %>% head(n = 400)


#Grouped data final
data_grouped <- data.frame(
  CTRL_33 = rowMeans(dplyr::select(data_final, matches('CTRL.*33'))),
  AD_33 = rowMeans(dplyr::select(data_final, matches('AD.*33'))),
  CTRL_34 = rowMeans(dplyr::select(data_final, matches('CTRL.*34'))),
  AD_34 = rowMeans(dplyr::select(data_final, matches('AD.*34'))))
  
  
  
data_subset_200_grouped <- data.frame(
  CTRL_33 = rowMeans(dplyr::select(data_subset_200, matches('CTRL.*33'))),
  AD_33 = rowMeans(dplyr::select(data_subset_200, matches('AD.*33'))),
  CTRL_34 = rowMeans(dplyr::select(data_subset_200, matches('CTRL.*34'))),
  AD_34 = rowMeans(dplyr::select(data_subset_200, matches('AD.*34'))))


data_subset_50_grouped <- data.frame(
  CTRL_33 = rowMeans(dplyr::select(data_subset_50, matches('CTRL.*33'))),
  AD_33 = rowMeans(dplyr::select(data_subset_50, matches('AD.*33'))),
  CTRL_34 = rowMeans(dplyr::select(data_subset_50, matches('CTRL.*34'))),
  AD_34 = rowMeans(dplyr::select(data_subset_50, matches('AD.*34'))))


# Normalize your data according to Z score using cal_z_score function.
data_norm <- t(apply(data_final, 1, cal_z_score))
data_final_adj_norm <- t(apply(data_final_adj, 1, cal_z_score))
data_subset_50 <- t(apply(data_subset_50, 1, cal_z_score))
data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))
data_subset_50_grouped <- t(apply(data_subset_50_grouped, 1, cal_z_score))
data_subset_200_grouped <- t(apply(data_subset_200_grouped, 1, cal_z_score))

#Export directory
dir.create(paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix))
dir.create(paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix, '/(APOE34:AD_vs_CTRL)_vs_(APOE33:AD_vs_CTRL)_specific'))
heatmap_dir <- paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix,'/(APOE34:AD_vs_CTRL)_vs_(APOE33:AD_vs_CTRL)_specific')

#Setting up Heatmap parameters
gaps = c(6,10,15,21,26,31)
clusters = 10
color_set <- c("navy", "white", "firebrick3")
all_degs_export_height <- nrow(data_final)
padj_degs_export_height <- nrow(data_final_adj)
## ------------- DENDROGRAM AND CLUSTER MANIPULATION (ADVANCED) -------------##
### GOURPED
pdf(file = paste0(heatmap_dir, '/All_DEGs_grouped.pdf'), pointsize = 10)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(33),
                     #breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 20,
                     #cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_grouped.pdf'), pointsize = 10, height = nrow(data_grouped_norm)/7)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(33),
                     #breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 20,
                     cellheight = 9,
                     scale = 'row')
dev.off()


#Top 200 grouped labelled
pdf(file = paste0(heatmap_dir,'/Top_200_grouped.pdf'), pointsize = 10, height= nrow(data_subset_200_grouped)/6)
phm_full <- pheatmap(data_subset_200_grouped,
                     color = colorRampPalette(c(color_set))(33),
                     #breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()


#Top 50 grouped labelled
pdf(file = paste0(heatmap_dir,'/Top_50_grouped.pdf'), pointsize = 10, height= nrow(data_subset_200_grouped)/6)
phm_full <- pheatmap(data_subset_50_grouped,
                     color = colorRampPalette(c(color_set))(33),
                     #breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()


##### HEATMAP 
# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

spec_heatmap <- read_excel("results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/All_samples_simple_WALD_contrast_KK/(APOE34:AD_vs_CTRL)_vs_(APOE33:AD_vs_CTRL)/(APOE34:AD_vs_CTRL)_vs_(APOE33:AD_vs_CTRL)_DEGene_counts_pval05.xlsx")

data_final <- spec_heatmap %>% column_to_rownames('gene')
data_final_adj <- sig_export_adj %>% column_to_rownames('gene')

data_subset_50 <- data_final %>% head(n = 50)
data_subset_200 <- data_final %>% head(n = 200)
data_subset_400 <- data_final %>% head(n = 400)


#Grouped data final
data_grouped <- data.frame(
  CTRL_33 = rowMeans(dplyr::select(data_final, matches('CTRL.*33'))),
  AD_33 = rowMeans(dplyr::select(data_final, matches('AD.*33'))),
  CTRL_34 = rowMeans(dplyr::select(data_final, matches('CTRL.*34'))),
  AD_34 = rowMeans(dplyr::select(data_final, matches('AD.*34'))))



data_subset_200_grouped <- data.frame(
  CTRL_33 = rowMeans(dplyr::select(data_subset_200, matches('CTRL.*33'))),
  AD_33 = rowMeans(dplyr::select(data_subset_200, matches('AD.*33'))),
  CTRL_34 = rowMeans(dplyr::select(data_subset_200, matches('CTRL.*34'))),
  AD_34 = rowMeans(dplyr::select(data_subset_200, matches('AD.*34'))))


data_subset_50_grouped <- data.frame(
  CTRL_33 = rowMeans(dplyr::select(data_subset_50, matches('CTRL.*33'))),
  AD_33 = rowMeans(dplyr::select(data_subset_50, matches('AD.*33'))),
  CTRL_34 = rowMeans(dplyr::select(data_subset_50, matches('CTRL.*34'))),
  AD_34 = rowMeans(dplyr::select(data_subset_50, matches('AD.*34'))))


# Normalize your data according to Z score using cal_z_score function.
data_norm <- t(apply(data_final, 1, cal_z_score))
data_final_adj_norm <- t(apply(data_final_adj, 1, cal_z_score))
data_subset_50 <- t(apply(data_subset_50, 1, cal_z_score))
data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))
data_subset_50_grouped <- t(apply(data_subset_50_grouped, 1, cal_z_score))
data_subset_200_grouped <- t(apply(data_subset_200_grouped, 1, cal_z_score))

#Export directory
dir.create(paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix))
dir.create(paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix, '/(APOE34:AD_vs_CTRL)_vs_(APOE33:AD_vs_CTRL)_specific'))
heatmap_dir <- paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix,'/(APOE34:AD_vs_CTRL)_vs_(APOE33:AD_vs_CTRL)_specific')

#Setting up Heatmap parameters
gaps = c(6,10,15,21,26,31)
clusters = 10
color_set <- c("navy", "white", "firebrick3")
all_degs_export_height <- nrow(data_final)
padj_degs_export_height <- nrow(data_final_adj)
## ------------- DENDROGRAM AND CLUSTER MANIPULATION (ADVANCED) -------------##
### GOURPED
pdf(file = paste0(heatmap_dir, '/All_DEGs_grouped.pdf'), pointsize = 10)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(33),
                     #breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 20,
                     #cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_grouped.pdf'), pointsize = 10, height = nrow(data_grouped_norm)/7)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(33),
                     #breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 20,
                     cellheight = 9,
                     scale = 'row')
dev.off()


#Top 200 grouped labelled
pdf(file = paste0(heatmap_dir,'/Top_200_grouped.pdf'), pointsize = 10, height= nrow(data_subset_200_grouped)/6)
phm_full <- pheatmap(data_subset_200_grouped,
                     color = colorRampPalette(c(color_set))(33),
                     #breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()


#Top 50 grouped labelled
pdf(file = paste0(heatmap_dir,'/Top_50_grouped.pdf'), pointsize = 10, height= nrow(data_subset_200_grouped)/6)
phm_full <- pheatmap(data_subset_50_grouped,
                     color = colorRampPalette(c(color_set))(33),
                     #breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()











#### -------- ANALYSIS AD SAMPLES APOE DOSAGE - SIMPLE DESIGN - WALD - ALL COMPARISONS------- ####
#Setting experiment conditions 
file_prefix <- 'All_samples_simple_WALD_contrast_APOE_DOSSAGE_AD_KK'
data <- data_adjusted_mg
experiment <- sub_meta
experiment

#Designs 
#dds <- DESeqDataSetFromMatrix(data, colData = experiment, design = ~ sex + age + pmd_min + apoe + diagnosis) 
dds <- DESeqDataSetFromMatrix(data, colData = experiment, design = ~ 0 + group + age + pmd_min)

## Run DESeq analysis to gather differential expression results
#Run DESeq (WALD - decault to extract the relevant comparisons
dds_run <- DESeq(dds)

# View names of estimated effects and set up comparisons
resultsNames(dds_run)

# meta data: to be used for complex contrasts
metadt <- data.frame(group = resultsNames(dds_run))

# define the complex contrasts
cont.list <- list( 
  "AD:APOE34_vs_APOE33" = data.frame(
    left      = c("groupAD_f_34","groupAD_m_34"),
    right     = c("groupAD_f_33","groupAD_m_33")
  ),
  "AD:APOE44_vs_APOE33" = data.frame(
    left      = c("groupAD_f_44","groupAD_m_44"),
    right     = c("groupAD_f_33","groupAD_m_33")
  ),
  "AD:APOE34-APOE44_vs_APOE33" = data.frame(
    left      = c("groupAD_f_34","groupAD_m_34","groupAD_f_44","groupAD_m_44"),
    right     = c("groupAD_f_33","groupAD_m_33",NA,NA)
  ),
  "AD:APOE33-APOE34_vs_APOE44" = data.frame(
    left      = c("groupAD_f_34","groupAD_m_34","groupAD_f_33","groupAD_m_33"),
    right     = c("groupAD_f_44","groupAD_m_44",NA,NA)
  )
)

cont.list
cmplx.cont <- makeComplexContrasts(cmplx = cont.list,
                                   contrastDT = metadt,
                                   useComplexNames = FALSE,
                                   save2file = TRUE,
                                   filename = 'results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/AD_APOE_dosage_WALD.xlsx',  #requires .xlsx at the end
                                   export.as = 'df')

cmplx.cont <- cmplx.cont[resultsNames(dds_run),]
cmplx.cont
dds_result <- data.frame()

for (i in 1:ncol(cmplx.cont)) {
  dds_result_total <- results(dds_run, pAdjustMethod = 'BH', independentFiltering = T, cooksCutoff = T, 
                              contrast = cmplx.cont[,i])
  dds_result_total$Contrast <- colnames(cmplx.cont)[i]
  dds_result_total <- as.data.frame(dds_result_total) %>% filter(pvalue < 0.05) %>% filter(baseMean > 1) %>% na.omit()
  dds_result_total <- dds_result_total %>% rownames_to_column('gene')
  dds_result <- rbind(dds_result,dds_result_total)
}

#Results
results_pval <- as.data.frame(table(dds_result$Contrast))
results_pval$group <- paste('pvalue_005')

dds_result_pval001 <- dds_result %>% filter(pvalue < 0.01)
results_pval001 <- as.data.frame(table(dds_result_pval001$Contrast))
results_pval001$group <- paste('pvalue_001')

dds_result_pval0001 <- dds_result %>% filter(pvalue < 0.001)
results_pval0001 <- as.data.frame(table(dds_result_pval0001$Contrast))
results_pval0001$group <- paste('pvalue_0001')

dds_result_padj01 <- dds_result %>% filter(padj < 0.1)
results_padj01 <- as.data.frame(table(dds_result_padj01$Contrast))
results_padj01$group <- paste('padj_01')

dds_result_padj005 <- dds_result %>% filter(padj < 0.05)
results_padj005 <- as.data.frame(table(dds_result_padj005$Contrast))
results_padj005$group <- paste('padj_005')

total_Degs_subset <- rbind(results_pval,results_pval001,results_pval0001,results_padj01,results_padj005)


dir.create(paste0('plots/barplots/Microglia/','DEGs_across_comparisons','/'))
barplot_dir <- paste0('plots/barplots/Microglia/','DEGs_across_comparisons','/',file_prefix)

lineWidth=0.5
pointSize=20
barplot <- ggplot(total_Degs_subset, aes(x = Var1, y = Freq, fill = group)) + 
  geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
  expand_limits(x = 0, y = 0) +
  theme(panel.spacing = unit(1, "lines")) +
  labs(x = NULL, y = c("Number of DEGs")) +
  ggtitle(paste0('No_DEGs_cutoffs'))+
  theme(
    text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(linewidth = lineWidth, colour = "black"),
    plot.title  = element_text(color="black", size=40, face="bold.italic"),
    axis.title  = element_text(size = pointSize, colour = "black"),
    axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
    axis.text.y  = element_text(size = pointSize , colour = "black"),
    legend.position = "left",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = pointSize , colour = "black"),
    axis.line = element_line(linewidth = lineWidth, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
  scale_fill_brewer(palette = "Paired") + scale_y_log10(expand = c(0, 0, .05, 0))
pdf(file = paste0(barplot_dir,'.pdf'), pointsize = 10, width = 30, height= 20)
print(barplot)
dev.off()


# View Total number of normalized counts per sample
dds_run <- estimateSizeFactors(dds_run)
DS_norm_counts <- counts(dds_run, normalized = TRUE)
colSums(DS_norm_counts)

#### -------- OUTPUT TABLES & RESULTS 
# Build significant gene table P - VALUE 0.05 and extract list of sorted DE genes.
sig_res <- dds_result %>%
  data.frame() %>%
  #rownames_to_column(var='Gene_Contrast') %>%  
  as_tibble %>%
  filter(pvalue < 0.05) %>%
  filter(baseMean > 1) %>%
  arrange(pvalue) 

# Build significant gene table P - ADJUSTED 0.05 and extract list of sorted DE genes.
sig_res_adj <- sig_res %>%
  filter(padj < 0.05) %>%
  arrange(padj) 

sorted_DEGenes <- sig_res$gene %>% unique()
sorted_DEGenes_padj <- sig_res_adj$gene %>% unique()

# Export significant gene count data
name_list <- c('gene', (paste0(dds_run$group, "_", dds_run$sample)))
name_list

sig_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes, gene)) 

sig_export_adj <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes_padj, gene))

# Setup Order of Counts
reordered_index <- c(grep("AD_m_33", names(sig_export), ignore.case = T),
                     grep("AD_f_33", names(sig_export), ignore.case = T),
                     grep("AD_m_34", names(sig_export), ignore.case = T),
                     grep("AD_f_34", names(sig_export), ignore.case = T),
                     grep("AD_m_44", names(sig_export), ignore.case = T),
                     grep("AD_f_44", names(sig_export), ignore.case = T))

sig_export <- sig_export %>%
  dplyr::select(c(1, reordered_index))

sig_export_adj <- sig_export_adj %>%
  dplyr::select(c(1, reordered_index))

# Column Names
check_cols <- c(colnames(sig_export))
check_cols

# Normalized counts
check_counts <- DS_norm_counts %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  dplyr::select(c(1, reordered_index)) %>%
  as_tibble %>% 
  `colnames<-`(check_cols) #%>%
#filter(!str_detect(gene, unwanted_genes))


# Results
check_res <- dds_result 
file_prefix
# Output files
dir.create(paste0('results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/', file_prefix))
dir.create(paste0('results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/', file_prefix, '/# All_DEGs'))

results_dir <- paste0('results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/', file_prefix, '/# All_DEGs')
results_dir

#Statistics
write.xlsx(sig_res, file = paste0(results_dir, '/DEGene_statistics_pval05.xlsx'), overwrite = T)
write.xlsx(sig_res_adj, file = paste0(results_dir, '/DEGene_statistics_padj05.xlsx'), overwrite = T)
write.xlsx(check_res, file = paste0(results_dir, '/Statistics.xlsx'), overwrite = T)
#Counts
write.xlsx(sig_export, file = paste0(results_dir, '/DEGene_counts_pval05.xlsx'), overwrite = T)
write.xlsx(sig_export_adj, file = paste0(results_dir, '/DEGene_counts_padj05.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0(results_dir, '/DS_counts.xlsx'), overwrite = T)



##### HEATMAP 
# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_final <- sig_export %>% column_to_rownames('gene')
data_final_adj <- sig_export_adj %>% column_to_rownames('gene')

data_subset_50 <- data_final %>% head(n = 50)
data_subset_200 <- data_final %>% head(n = 200)
data_subset_400 <- data_final %>% head(n = 400)


#Grouped data final
data_grouped <- data.frame(
  AD_33 = rowMeans(dplyr::select(data_final, matches('AD.*33'))),
  AD_34 = rowMeans(dplyr::select(data_final, matches("AD.*34"))),
  AD_44 = rowMeans(dplyr::select(data_final, matches("AD.*44"))))

data_subset_50_grouped <- data.frame(
  AD_33 = rowMeans(dplyr::select(data_subset_50, matches('AD.*33'))),
  AD_34 = rowMeans(dplyr::select(data_subset_50, matches("AD.*34"))),
  AD_44 = rowMeans(dplyr::select(data_subset_50, matches("AD.*44"))))

data_subset_200_grouped <- data.frame(
  AD_33 = rowMeans(dplyr::select(data_subset_200, matches('AD.*33'))),
  AD_34 = rowMeans(dplyr::select(data_subset_200, matches("AD.*34"))),
  AD_44 = rowMeans(dplyr::select(data_subset_200, matches("AD.*44"))))




# Normalize your data according to Z score using cal_z_score function.
data_norm <- t(apply(data_final, 1, cal_z_score))
data_final_adj_norm <- t(apply(data_final_adj, 1, cal_z_score))
data_subset_50 <- t(apply(data_subset_50, 1, cal_z_score))
data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))
data_grouped_50_norm <- t(apply(data_subset_50_grouped, 1, cal_z_score))
data_grouped_200_norm <- t(apply(data_subset_200_grouped, 1, cal_z_score))


#Export directory
dir.create(paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix))
dir.create(paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix, '/# All_DEGs'))
heatmap_dir <- paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix,'/# All_DEGs')
heatmap_dir
#Setting up Heatmap parameters
gaps = c(6,10,15,21)
clusters = 5
color_set <- c("navy", "white", "firebrick3")
all_degs_export_height <- nrow(data_final)
padj_degs_export_height <- nrow(data_final_adj)
## ------------- DENDROGRAM AND CLUSTER MANIPULATION (ADVANCED) -------------##
# Create heatmap
pdf(file = paste0(heatmap_dir, '/All_DEGs_padj05.pdf'), pointsize = 10, height = padj_degs_export_height/4)
phm_full <- pheatmap(data_final_adj_norm,
                     color = colorRampPalette(color_set)(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     cellheight = 9,
                     #cellwidth = 25,
                     scale = 'row'
)
dev.off()



pdf(file = paste0(heatmap_dir,'/All_DEGs_pval05.pdf'), pointsize = 10)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(color_set)(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row'
                     #cellheight = 3,
)
dev.off()


pdf(file = paste0(heatmap_dir, '/All_DEGs_labelled_pval05.pdf'), pointsize = 10, height = all_degs_export_height/4)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 8
)
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_50.pdf'), pointsize = 10, height = 10)
phm_50 <- pheatmap(data_subset_50,
                   color = colorRampPalette(color_set)(41),
                   breaks = seq(-2, 2, by = 0.1),
                   kmeans_k = NA,
                   cluster_rows = T,
                   cutree_row = clusters,
                   cluster_cols = F,
                   #cutree_cols = 4,
                   gaps_col = gaps,              
                   legend = TRUE,
                   show_rownames = T,
                   #cellwidth = 25,
                   cellheight = 9,
                   treeheight_col = 0,
                   treeheight_row = 50,
                   border_color = 'NA',
                   fontsize = 8,
                   scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_200.pdf'), pointsize = 10, height = 35)
phm_200 <- pheatmap(data_subset_200,
                    color = colorRampPalette(color_set)(41),
                    breaks = seq(-2, 2, by = 0.1),
                    kmeans_k = NA,
                    cluster_rows = T,
                    cutree_row = clusters,
                    cluster_cols = F,
                    #cutree_cols = 4,
                    gaps_col = gaps,              
                    legend = TRUE,
                    show_rownames = T,
                    #cellwidth = 25,
                    cellheight = 9,
                    treeheight_col = 0,
                    treeheight_row = 50,
                    border_color = 'NA',
                    fontsize = 8,
                    scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_400.pdf'), pointsize = 10, height = 60)
phm_full <- pheatmap(data_subset_400,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 9)
dev.off()


### GOURPED
pdf(file = paste0(heatmap_dir, '/All_DEGs_grouped.pdf'), pointsize = 10)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 20,
                     #cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_grouped.pdf'), pointsize = 10, height = nrow(data_grouped_norm)/7)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 20,
                     cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_padj05_grouped.pdf'), pointsize = 10, height= nrow(data_grouped_norm_padj)/6)
phm_full <- pheatmap(data_grouped_norm_padj,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()

#Top 50 grouped labelled
pdf(file = paste0(heatmap_dir,'/Top_50_grouped.pdf'), pointsize = 10, height= nrow(data_subset_200_grouped)/6)
phm_full <- pheatmap(data_subset_50_grouped,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()

#Top 200 grouped labelled
pdf(file = paste0(heatmap_dir,'/Top_200_grouped.pdf'), pointsize = 10, height= nrow(data_subset_200_grouped)/6)
phm_full <- pheatmap(data_subset_200_grouped,
                     color = colorRampPalette(c(color_set))(33),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()










#### -------- ANALYSIS AD SAMPLES APOE DOSAGE - SIMPLE DESIGN - WALD - ALL COMPARISONS - MALE vs FEMALE ------- ####
#Setting experiment conditions 
file_prefix <- 'All_samples_simple_WALD_contrast_APOE_DOSSAGE_AD_Female_vs_Male_KK'

#Metadata
experiment<- sub_meta
experiment

#Comparison dataframe
data <- data_adjusted_mg
data

#Designs 
dds <- DESeqDataSetFromMatrix(data, colData = experiment, design = ~ 0 + group + age + pmd_min)

## Run DESeq analysis to gather differential expression results
#Run DESeq (WALD - decault to extract the relevant comparisons
dds_run <- DESeq(dds)

# View names of estimated effects and set up comparisons
resultsNames(dds_run)

# meta data: to be used for complex contrasts
metadt <- data.frame(group = resultsNames(dds_run))

# define the complex contrasts
cont.list <- list( 
  #Simple comparisons
  # "FEMALE:AD:APOE34_vs_APOE33" = data.frame(
  #   left      = c("groupAD_f_34"),
  #   right     = c("groupAD_f_33")
  # ),
  # "MALE:AD:APOE34_vs_APOE33" = data.frame(
  #   left      = c("groupAD_m_34"),
  #   right     = c("groupAD_m_33")
  # ),
  # "FEMALE:AD:APOE44_vs_APOE33" = data.frame(
  #   left      = c("groupAD_f_44"),
  #   right     = c("groupAD_f_33")
  # ),
  # "MALE:AD:APOE44_vs_APOE33" = data.frame(
  #   left      = c("groupAD_m_44"),
  #   right     = c("groupAD_m_33")
  # ),
  # "FEMALE:AD:APOE34-APOE44_vs_APOE33" = data.frame(
  #   left      = c("groupAD_f_34","groupAD_f_44"),
  #   right     = c("groupAD_f_33",NA)
  # ),
  # "MALE:AD:APOE34-APOE44_vs_APOE33" = data.frame(
  #   left      = c("groupAD_m_34","groupAD_m_44"),
  #   right     = c("groupAD_m_33",NA)
  # ),
  #Apoe34 vs APOE33 comparing AD females and males
  "(FEMALE:AD:APOE34_vs_APOE33)_vs_(MALE:AD:APOE34_vs_APOE33)"  = data.frame(
    left      = c("groupAD_f_34", "groupAD_f_33"),
    right     = c("groupAD_m_34", "groupAD_m_33"),
    leftFC    = c(1, -1),
    rightFC   = c(1, -1)
  ),
  #Apoe34 vs APOE33 comparing AD females and males
  "(FEMALE:AD:APOE44_vs_APOE33)_vs_(MALE:AD:APOE44_vs_APOE33)"  = data.frame(
    left      = c("groupAD_f_44", "groupAD_f_33"),
    right     = c("groupAD_m_44", "groupAD_m_33"),
    leftFC    = c(1, -1),
    rightFC   = c(1, -1)
  ),
  #Apoe34/44 vs APOE33 comparing AD females and males
  "(FEMALE:AD:APOE34-44_vs_APOE33)_vs_(MALE:AD:APOE34-44_vs_APOE33)"  = data.frame(
    left      = c("groupAD_f_34","groupAD_f_44", "groupAD_f_33"),
    right     = c("groupAD_m_34","groupAD_m_44", "groupAD_m_33"),
    leftFC    = c(0.5,0.5, -1),
    rightFC   = c(0.5,0.5, -1)
  ),
  #Apoe34/44 vs APOE33 comparing AD females and males
  "(FEMALE:AD:APOE33-34_vs_APOE44)_vs_(MALE:AD:APOE33-34_vs_APOE44)"  = data.frame(
    left      = c("groupAD_f_34","groupAD_f_33", "groupAD_f_44"),
    right     = c("groupAD_m_34","groupAD_m_33", "groupAD_m_44"),
    leftFC    = c(0.5,0.5, -1),
    rightFC   = c(0.5,0.5, -1)
  )
)

cont.list
cmplx.cont <- makeComplexContrasts(cmplx = cont.list,
                                   contrastDT = metadt,
                                   useComplexNames = FALSE,
                                   save2file = TRUE,
                                   filename = 'results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/AD_APOE_dosage_Sex_differences_WALD.xlsx',  #requires .xlsx at the end
                                   export.as = 'df')

cmplx.cont <- cmplx.cont[resultsNames(dds_run),]
cmplx.cont
dds_result <- data.frame()

for (i in 1:ncol(cmplx.cont)) {
  dds_result_total <- results(dds_run, pAdjustMethod = 'BH', independentFiltering = T, cooksCutoff = T, 
                              contrast = cmplx.cont[,i])
  dds_result_total$Contrast <- colnames(cmplx.cont)[i]
  dds_result_total <- as.data.frame(dds_result_total) %>% filter(pvalue < 0.05) %>% filter(baseMean > 1) %>% na.omit()
  dds_result_total <- dds_result_total %>% rownames_to_column('gene')
  dds_result <- rbind(dds_result,dds_result_total)
}

#Results
results_pval <- as.data.frame(table(dds_result$Contrast))
results_pval$group <- paste('pvalue_005')

dds_result_pval001 <- dds_result %>% filter(pvalue < 0.01)
results_pval001 <- as.data.frame(table(dds_result_pval001$Contrast))
results_pval001$group <- paste('pvalue_001')

dds_result_pval0001 <- dds_result %>% filter(pvalue < 0.001)
results_pval0001 <- as.data.frame(table(dds_result_pval0001$Contrast))
results_pval0001$group <- paste('pvalue_0001')

dds_result_padj01 <- dds_result %>% filter(padj < 0.1)
results_padj01 <- as.data.frame(table(dds_result_padj01$Contrast))
results_padj01$group <- paste('padj_01')

dds_result_padj005 <- dds_result %>% filter(padj < 0.05)
results_padj005 <- as.data.frame(table(dds_result_padj005$Contrast))
results_padj005$group <- paste('padj_005')

total_Degs_subset <- rbind(results_pval,results_pval001,results_pval0001,results_padj01,results_padj005)

#Write results
write.xlsx(total_Degs_subset,paste0('results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/DEGs_per_comparison/AD_APOE_dosage_Sex_differences_WALD.xlsx'))

dir.create(paste0('plots/barplots/Microglia/','DEGs_across_comparisons','/'))
barplot_dir <- paste0('plots/barplots/Microglia/','DEGs_across_comparisons','/',file_prefix)
barplot_dir
lineWidth=0.5
pointSize=20
barplot <- ggplot(total_Degs_subset, aes(x = Var1, y = Freq, fill = group)) + 
  geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
  expand_limits(x = 0, y = 0) +
  theme(panel.spacing = unit(1, "lines")) +
  labs(x = NULL, y = c("Number of DEGs")) +
  ggtitle(paste0('No_DEGs_cutoffs'))+
  theme(
    text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(linewidth = lineWidth, colour = "black"),
    plot.title  = element_text(color="black", size=40, face="bold.italic"),
    axis.title  = element_text(size = pointSize, colour = "black"),
    axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
    axis.text.y  = element_text(size = pointSize , colour = "black"),
    legend.position = "left",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = pointSize , colour = "black"),
    axis.line = element_line(linewidth = lineWidth, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
  scale_fill_brewer(palette = "Paired") + scale_y_log10(expand = c(0, 0, .05, 0))
pdf(file = paste0(barplot_dir,'.pdf'), pointsize = 10, width = 30, height= 20)
print(barplot)
dev.off()


# View Total number of normalized counts per sample
dds_run <- estimateSizeFactors(dds_run)
DS_norm_counts <- counts(dds_run, normalized = TRUE)
colSums(DS_norm_counts)

#### -------- OUTPUT TABLES & RESULTS 
# Build significant gene table P - VALUE 0.05 and extract list of sorted DE genes.
sig_res <- dds_result %>%
  data.frame() %>%
  #rownames_to_column(var='Gene_Contrast') %>%  
  as_tibble %>%
  filter(pvalue < 0.05) %>%
  filter(baseMean > 1) %>%
  arrange(pvalue) 

# Build significant gene table P - ADJUSTED 0.05 and extract list of sorted DE genes.
sig_res_adj <- sig_res %>%
  filter(padj < 0.05) %>%
  arrange(padj) 

sorted_DEGenes <- sig_res$gene %>% unique()
sorted_DEGenes_padj <- sig_res_adj$gene %>% unique()

# Export significant gene count data
name_list <- c('gene', (paste0(dds_run$group, "_", dds_run$sample)))
name_list

sig_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes, gene)) 

sig_export_adj <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes_padj, gene))

# Setup Order of Counts
reordered_index <- c(grep("AD_m_33", names(sig_export), ignore.case = T),
                     grep("AD_f_33", names(sig_export), ignore.case = T),
                     grep("AD_m_34", names(sig_export), ignore.case = T),
                     grep("AD_f_34", names(sig_export), ignore.case = T),
                     grep("AD_m_44", names(sig_export), ignore.case = T),
                     grep("AD_f_44", names(sig_export), ignore.case = T))

sig_export <- sig_export %>%
  dplyr::select(c(1, reordered_index))

sig_export_adj <- sig_export_adj %>%
  dplyr::select(c(1, reordered_index))

# Column Names
check_cols <- c(colnames(sig_export))
check_cols

# Normalized counts
check_counts <- DS_norm_counts %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  dplyr::select(c(1, reordered_index)) %>%
  as_tibble %>% 
  `colnames<-`(check_cols) #%>%
#filter(!str_detect(gene, unwanted_genes))


# Results
check_res <- dds_result 
file_prefix
# Output files
dir.create(paste0('results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/', file_prefix))
dir.create(paste0('results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/', file_prefix, '/# All_DEGs'))

results_dir <- paste0('results/Microglia/CCA_Round3/Pseudobulk/Aggregate_Sample/', file_prefix, '/# All_DEGs')
results_dir

#Statistics
write.xlsx(sig_res, file = paste0(results_dir, '/DEGene_statistics_pval05.xlsx'), overwrite = T)
write.xlsx(sig_res_adj, file = paste0(results_dir, '/DEGene_statistics_padj05.xlsx'), overwrite = T)
write.xlsx(check_res, file = paste0(results_dir, '/Statistics.xlsx'), overwrite = T)
#Counts
write.xlsx(sig_export, file = paste0(results_dir, '/DEGene_counts_pval05.xlsx'), overwrite = T)
write.xlsx(sig_export_adj, file = paste0(results_dir, '/DEGene_counts_padj05.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0(results_dir, '/DS_counts.xlsx'), overwrite = T)



##### HEATMAP 
# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_final <- sig_export %>% column_to_rownames('gene')
data_final_adj <- sig_export_adj %>% column_to_rownames('gene')

data_subset_50 <- data_final %>% head(n = 50)
data_subset_200 <- data_final %>% head(n = 200)
data_subset_400 <- data_final %>% head(n = 400)


#Grouped data final
data_grouped <- data.frame(
    AD_m_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_m_33")))),
    AD_m_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_m_34")))),
    AD_m_44 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_m_44")))),
    AD_f_33 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_33")))),
    AD_f_34 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_34")))),
    AD_f_44 = rowMeans(dplyr::select(data_final, dplyr::contains(("AD_f_44")))))

data_subset_50_grouped <- data.frame(
  AD_m_33 = rowMeans(dplyr::select(data_subset_50, dplyr::contains(("AD_m_33")))),
  AD_m_34 = rowMeans(dplyr::select(data_subset_50, dplyr::contains(("AD_m_34")))),
  AD_m_44 = rowMeans(dplyr::select(data_subset_50, dplyr::contains(("AD_m_44")))),
  AD_f_33 = rowMeans(dplyr::select(data_subset_50, dplyr::contains(("AD_f_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_subset_50, dplyr::contains(("AD_f_34")))),
  AD_f_44 = rowMeans(dplyr::select(data_subset_50, dplyr::contains(("AD_f_44")))))

data_subset_200_grouped <- data.frame(
  AD_m_33 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_m_33")))),
  AD_m_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_m_34")))),
  AD_m_44 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_m_44")))),
  AD_f_33 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_f_33")))),
  AD_f_34 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_f_34")))),
  AD_f_44 = rowMeans(dplyr::select(data_subset_200, dplyr::contains(("AD_f_44")))))



# Normalize your data according to Z score using cal_z_score function.
data_norm <- t(apply(data_final, 1, cal_z_score))
data_final_adj_norm <- t(apply(data_final_adj, 1, cal_z_score))
data_subset_50 <- t(apply(data_subset_50, 1, cal_z_score))
data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))
data_grouped_50_norm <- t(apply(data_subset_50_grouped, 1, cal_z_score))
data_grouped_200_norm <- t(apply(data_subset_200_grouped, 1, cal_z_score))


#Export directory
dir.create(paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix))
dir.create(paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix, '/# All_DEGs'))
heatmap_dir <- paste0('plots/heatmaps/Microglia/Pseudobulk/Aggregate_Sample/',file_prefix,'/# All_DEGs')
heatmap_dir
#Setting up Heatmap parameters
gaps = c(5,10,15,19,24)
clusters = 10
color_set <- c("navy", "white", "firebrick3")
all_degs_export_height <- nrow(data_final)
padj_degs_export_height <- nrow(data_final_adj)
## ------------- DENDROGRAM AND CLUSTER MANIPULATION (ADVANCED) -------------##
# Create heatmap
# pdf(file = paste0(heatmap_dir, '/All_DEGs_padj05.pdf'), pointsize = 10, height = padj_degs_export_height/4)
# phm_full <- pheatmap(data_final_adj_norm,
#                      color = colorRampPalette(color_set)(41),
#                      breaks = seq(-2, 2, by = 0.1),
#                      kmeans_k = NA,
#                      cluster_rows = T,
#                      cutree_row = clusters,
#                      cluster_cols = F,
#                      #cutree_cols = 4,
#                      gaps_col = gaps,             
#                      legend = TRUE,
#                      show_rownames = T,
#                      border_color = 'NA',
#                      fontsize = 8,
#                      treeheight_col = 0,
#                      treeheight_row = 10,
#                      cellheight = 9,
#                      #cellwidth = 25,
#                      scale = 'row'
# )
# dev.off()



pdf(file = paste0(heatmap_dir,'/All_DEGs_pval05.pdf'), pointsize = 10)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(color_set)(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row'
                     #cellheight = 3,
)
dev.off()


pdf(file = paste0(heatmap_dir, '/All_DEGs_labelled_pval05.pdf'), pointsize = 10, height = all_degs_export_height/4)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 8
)
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_50.pdf'), pointsize = 10, height = 10)
phm_50 <- pheatmap(data_subset_50,
                   color = colorRampPalette(color_set)(41),
                   breaks = seq(-2, 2, by = 0.1),
                   kmeans_k = NA,
                   cluster_rows = T,
                   cutree_row = clusters,
                   cluster_cols = F,
                   #cutree_cols = 4,
                   gaps_col = gaps,              
                   legend = TRUE,
                   show_rownames = T,
                   #cellwidth = 25,
                   cellheight = 9,
                   treeheight_col = 0,
                   treeheight_row = 50,
                   border_color = 'NA',
                   fontsize = 8,
                   scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_200.pdf'), pointsize = 10, height = 35)
phm_200 <- pheatmap(data_subset_200,
                    color = colorRampPalette(color_set)(41),
                    breaks = seq(-2, 2, by = 0.1),
                    kmeans_k = NA,
                    cluster_rows = T,
                    cutree_row = clusters,
                    cluster_cols = F,
                    #cutree_cols = 4,
                    gaps_col = gaps,              
                    legend = TRUE,
                    show_rownames = T,
                    #cellwidth = 25,
                    cellheight = 9,
                    treeheight_col = 0,
                    treeheight_row = 50,
                    border_color = 'NA',
                    fontsize = 8,
                    scale = 'row')
dev.off()

pdf(file = paste0(heatmap_dir, '/Top_400.pdf'), pointsize = 10, height = 60)
phm_full <- pheatmap(data_subset_400,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 9)
dev.off()


### GOURPED
pdf(file = paste0(heatmap_dir, '/All_DEGs_grouped.pdf'), pointsize = 10)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(40),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 20,
                     #cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_grouped.pdf'), pointsize = 10, height = nrow(data_grouped_norm)/7)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(40),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 20,
                     cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
# pdf(file = paste0(heatmap_dir,'/All_DEGs_labelled_padj05_grouped.pdf'), pointsize = 10, height= nrow(data_grouped_norm_padj)/6)
# phm_full <- pheatmap(data_grouped_norm_padj,
#                      color = colorRampPalette(c(color_set))(33),
#                      breaks = seq(-2, 2, by = 0.1),
#                      #annotation_row = sig_res$Contrast,
#                      kmeans_k = NA,
#                      cluster_rows = T,
#                      cutree_row = clusters,
#                      cluster_cols = F,
#                      #cutree_cols = 3,
#                      gaps_col = c(),             
#                      legend = TRUE,
#                      show_rownames = T,
#                      border_color = 'NA',
#                      fontsize = 8,
#                      treeheight_col = 0,
#                      treeheight_row = 0,
#                      cellwidth = 30,
#                      cellheight = 9,
#                      scale = 'row')
# dev.off()

#Top 50 grouped labelled
pdf(file = paste0(heatmap_dir,'/Top_50_grouped.pdf'), pointsize = 10, height= nrow(data_subset_200_grouped)/6)
phm_full <- pheatmap(data_subset_50_grouped,
                     color = colorRampPalette(c(color_set))(40),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()

#Top 200 grouped labelled
pdf(file = paste0(heatmap_dir,'/Top_200_grouped.pdf'), pointsize = 10, height= nrow(data_subset_200_grouped)/6)
phm_full <- pheatmap(data_subset_200_grouped,
                     color = colorRampPalette(c(color_set))(40),
                     breaks = seq(-2, 2, by = 0.1),
                     #annotation_row = sig_res$Contrast,
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 30,
                     cellheight = 9,
                     scale = 'row')
dev.off()













#### ---------CLUSTER EXTRACTION ####################################
# Build cluster tables from dendrogram based on k
gene_clusters <- data.frame(sort(cutree(phm_full$tree_row, k=4))) %>%
  rownames_to_column(var="gene") 

temp = gene_clusters$sort.cutree.phm_full.tree_row..k...4..

gene_clusters <- gene_clusters %>% 
  mutate(cluster = temp) %>%
  select(gene, cluster)
gene_clusters

write.xlsx(gene_clusters, file=paste0('results/', file_prefix, 'clustergenes.xlsx'))




