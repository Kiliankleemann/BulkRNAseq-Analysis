#Crosstalk Analysis based on NichenetR database
#Database from NichenetR
#Cell Crosstalk 
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
write.xlsx(lr_network,paste0('lr_network.xlsx'))


#Load datasets for analysis
crosstalk_prefix <- 'PREFIX'

data <- read.xlsx(paste0('results/',crosstalk_prefix, 'DEGene_statistics_pval05.xlsx'))

sig_genes <- data %>% pull(gene)
sig_genes_network <- lr_network %>% filter(from %in% sig_genes)

#Import Interaction cell genes
interaction_celltype <- 'MICROGLIA'
#Receiver Cell Expression of receptors
geneExpression <- read.xlsx("-DIRECTORY-")
geneExpression_genes <- geneExpression_Microglia %>% filter(AVERAGE > 10) %>% pull(gene)


#Final list 
sig_genes_network_final <- sig_genes_network %>% filter(to %in% geneExpression_genes)
sig_genes_network_final$gene <- sig_genes_network_final$from

sig_genes_network_final <- merge(sig_genes_network_final,data,by=c("gene")) 
sig_genes_network_final$Receptor_Ligand <- paste0(sig_genes_network_final$from, "_",sig_genes_network_final$to)
sig_genes_network_final <- sig_genes_network_final[!duplicated(sig_genes_network_final$Receptor_Ligand),]


sig_genes_network_final <- sig_genes_network_final %>% group_by(gene) %>% mutate(Count = n())
sig_genes_network_final$Combined <- sig_genes_network_final$Count * sig_genes_network_final$log2FoldChange
write.xlsx(sig_genes_network_final, paste0('results/crosstalk/',crosstalk_prefix,'Neutorphil_to_',interaction_celltype,'.xlsx'))

#Edit overlapping receptors 
sig_genes_network_final <- sig_genes_network_final %>% 
  mutate(to = ifelse(to == 'PIK3CB','PIK3BR',to)) 

sig_genes_network_final <- sig_genes_network_final %>% 
  mutate(to = ifelse(to == 'ICAM1','ICAM1_R',to))

sig_genes_network_final <- sig_genes_network_final %>% 
  mutate(to = ifelse(to == 'ITGB1','ITGB1_R',to))

sig_genes_network_final <- sig_genes_network_final %>% 
  mutate(to = ifelse(to == 'ITGB2','ITGB2_R',to))

#########---------Making Circlos Plots---------########
#Select Data
df <- as.data.frame(sig_genes_network_final)
df <- df   %>% 
  select(c(`from`, `to`, `Combined`))
df <- df %>% arrange(desc(`Combined`))
# 
# 
# #use dyplr to remove small values 
# df_up <- filter(df, `Log2FC` > 0) %>% arrange(desc(`Log2FC`))
# df_down <- filter(df,`Log2FC` < 0)


#Setting Colour scheme
combined_max <- max(df$Combined)
combined_min <- min(df$Combined)
col_fun = colorRamp2(c(combined_min,0,combined_max), c("blue","whitesmoke","red"))
col_fun(seq(combined_min,combined_max, by = 0.01))

# grid.col = c(Lgals3 = "darkblue", Cxcl2 = "blue", Tnfsf13 = "lightblue", Cxcl16='cyan', Spp1 = 'blue2', Apoe = 'deepskyblue', Chad = 'dodgerblue',Slit2 = 'cadetblue',
#              Nptn = "forestgreen", Itgb1 = "seagreen2", Cadm1 = "palegreen", Colec12 = "olivedrab", Flt4 = "limegreen", Ackr3 = "green1", Tgfbr2 = 'khaki1', Fas = 'lightgreen', Tfrc = 'darkgreen', S1pr1 = 'darkolivegreen4', Lrp1 ='darkseagreen2', Robo4 = 'darkseagreen4', Sdc3 = 'darkslategrey', Scarb1 = 'chartreuse3')

#Customize graphic parameters before initializing the circlos: 
circos.clear()
circos.par(canvas.xlim = c(-1.5, 1.5), 
           canvas.ylim = c(-1.5, 1.5),
           track.margin= c(0.01, 0.01),
           start.degree = -90,     #rotating circlos plot # of °
           #gap.degree = 0.8,
           "track.height" = 0.1)

#circos.initialize(): allocates sectors on the circle.
#circos.track(): creates plotting regions for cells in one single track.

#Assigning grid and annotation regions / size
chordDiagram(df, big.gap = 10, small.gap = 2,
             annotationTrack = c('grid','names'), 
             col=col_fun, 
             annotationTrackHeight = mm_h(2), 
             preAllocateTracks = list(track.height = mm_h(4)),h.ratio=0.4,
             transparency = 0.2,  
             directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             link.sort = TRUE, link.decreasing = TRUE,
             link.zindex = rank(df[[3]]))
#abline(v = 0, lty = 2, col = "#00000080")
#circos.par

#Assign Annotations 
circos.track(track.index = 1, panel.fun = function(x,y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1]+ mm_y(5), CELL_META$sector.index,
              cex = 0.7,
              facing = 'clockwise',
              niceFacing = T, 
              adj = c(0,0.9))
}, bg.border = NA)
title(paste0('Crosstalk_Neutrophils_', crosstalk_prefix, 'OLAH_MG_Baseline_expression'))

#adding legend 
lgd_links = Legend(at = c(combined_min,0,combined_max), col_fun = col_fun, 
                   title_position = "topleft", title = "Regulatory Potential")
draw(lgd_links, x = unit(6, "mm"), y = unit(6, "mm"), just = c("left", "bottom"))

#highlighting sectors 
highlight.sector(df$from, track.index = 1, col = 'navyblue', text = 'NEUTROPHIL LIGANDS', cex = 0.7, text.col = 'white', facing = 'bending.inside', niceFacing = T)
highlight.sector(df$to, track.index = 1, col = 'cyan3', text = paste0(interaction_celltype,'RECEPTORS'), cex = 0.7, text.col = 'white', facing = 'bending.inside', niceFacing = T)
