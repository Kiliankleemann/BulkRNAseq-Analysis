#Immune Deconvolution of human bulk RNAseq data
install.packages("remotes")
remotes::install_github("grst/immunedeconv")
library(immunedeconv)
#load expression matrix
imm_deconv_data <- read.xlsx(paste0('results/', 'All_data_', 'DS_counts.xlsx'), rowNames = T)
imm_deconv_data <- read.xlsx(paste0('results/', 'All_data_', 'DEGene_counts_pval05.xlsx'), rowNames = T)

#Pseudobulk RNA-seq analysis
immune_dev_res <- immunedeconv::deconvolute(imm_deconv_data, "quantiseq")

#Immunodeconvolution 
immune_dev_res %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rev(levels(immune_dev_res)))
