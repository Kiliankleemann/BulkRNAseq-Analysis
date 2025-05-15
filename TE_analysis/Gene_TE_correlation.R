#### -------- INSTALLING AND LOADING PACKAGES ------- ####
list.of.packages <- c("BiocGenerics","tximport","S4Vectors", "GenomicRanges","DESeq2", "biomaRt","data.table","readr",
                      "ggplot2", "ggsignif", "ggpubr", "sva", "devtools", "rtracklayer", "org.Hs.eg.db", "ggnewscale","VennDiagram",
                      "org.Mm.eg.db", "limma","stringr","KEGGREST","ggrepel", "openxlsx", 'splitstackshape',
                      "fgsea","clusterProfiler","pheatmap","ggpubr","cowplot",'dplyr','janitor',"purrr",
                      "RColorBrewer","Gviz", 'AnnotationDbi', 'tidyverse','pheatmap', 'dendextend',"factoextra")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if(length(new.packages)) BiocManager::install(new.packages)


# Packages loading
invisible(lapply(list.of.packages, library, character.only = TRUE))


###### ----- SETTING WORK DIRECTORY -----#####
setwd("/home/kilian/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Aref_DNA_damage_Oct_2023/Kilian_results/")

#Set prefix
file_prefix <- "TELocal_18h_vs_ctrl"
Pval_cutoff <-  "Padj05"
  
#Create Directories
dir.create('results/TE_gene_correlation')
dir.create('plots/TE_gene_correlation/')


dir.create(paste0('results/TE_gene_correlation/',file_prefix))
dir.create(paste0('plots/TE_gene_correlation/',file_prefix))

###### ----- FILE PREPARATION -----#####
# Load differential gene information
DEG_data <- read.xlsx(paste0("results/",file_prefix,"/DEGene_statistics_pval05.xlsx"), sheet = 1)

# Filter for significantly up/downregulated genes
de_genes_filtered <- DEG_data %>%
  filter(!is.na(padj)) %>%
  filter(abs(log2FoldChange) >= 1, padj < 0.05) %>%
  filter(!grepl(':',gene))


# Load GTF annotation file (e.g., Ensembl) 
gtf <- import("/Users/kiliankleemann/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Immune_priming_TE/Database_GTF_files_TE/mm10.refGene.gtf")  # adjust file name as needed
# Filter for gene entries
gene_ranges <- gtf[gtf$type == "transcript"]

# Keep relevant metadata
gene_ranges$gene_name <- mcols(gene_ranges)$gene_name


# Remove duplicates (some genes may have multiple transcripts; take the longest or representative)
gene_df <- as.data.frame(gene_ranges) %>%
  group_by(gene_name) %>%
  summarise(
    seqnames = unique(seqnames)[1],
    start = min(start),
    end = max(end),
    strand = unique(strand)[1]
  )

# Convert to GRanges
genes_gr <- GRanges(
  seqnames = gene_df$seqnames,
  ranges = IRanges(start = gene_df$start, end = gene_df$end),
  strand = gene_df$strand,
  gene_name = gene_df$gene_name
)

#Filter relevant genes
genes_gr_filtered <- genes_gr[mcols(genes_gr)$gene_name %in% de_genes_filtered$gene]


#TE gene ranges file
# Load differential TE information.
tes <- read.xlsx(paste0("results/",file_prefix,"/DETE_locus_statistics_pval05.xlsx"), sheet = 1) %>% filter(padj < 0.05)

te_gr <- GRanges(
  seqnames = tes$Chromosome,
  ranges = IRanges(start = tes$Start, end = tes$Stop),
  te_id = tes$TE
)

# Set maximum gap
maxgap <- 100000

# Find overlaps within maxgap
hits <- findOverlaps(genes_gr_filtered, te_gr, maxgap = maxgap, ignore.strand = TRUE)

# Extract gene–TE proximity pairs
result <- data.frame(
  gene = mcols(genes_gr_filtered[queryHits(hits)])$gene_name,
  gene_chr  = as.character(seqnames(genes_gr_filtered[queryHits(hits)])),
  gene_start = start(genes_gr_filtered[queryHits(hits)]),
  gene_end   = end(genes_gr_filtered[queryHits(hits)]),
  TE     = mcols(te_gr[subjectHits(hits)])$te_id,
  te_chr    = as.character(seqnames(te_gr[subjectHits(hits)])),
  te_start  = start(te_gr[subjectHits(hits)]),
  te_end    = end(te_gr[subjectHits(hits)])
)
#Joining with dif. expression data
result <- result %>%
  left_join(DEG_data, by = "gene") %>%
  rename(
    gene_log2FC = log2FoldChange,
    gene_pvalue = pvalue
  )

# Step 3: Join TE DE information
result <- result %>%
  left_join(tes, by = "TE") %>%
  rename(
    te_log2FC = log2FoldChange,
    te_pvalue = pvalue
  )

#Export Results
write.xlsx(result, file = paste0('results/TE_gene_correlation/', file_prefix, '/Overlap_padj05.xlsx'), overwrite = T)



#Scatterplot
pointSize = 15
lineWidth = 0.5
scatter <- ggplot(result, aes(x= gene_log2FC, y= te_log2FC)) +
  geom_smooth(data = result, method = "lm", se=TRUE, color="red", formula = y ~ x) +
  stat_cor(data = result, color = 'red', label.x.npc = "left", label.y.npc = "top") +
  geom_point(size=0.5, color = 'grey90') +
  # geom_point(data = cmb_dt_filter_padj_up_up, size=1, color = 'firebrick') +
  # annotate("text", x = 1.2, y = 1.5, label= paste(pct_up_up,"%"), size = 5, color = 'firebrick')  +
  # geom_point(data = cmb_dt_filter_padj_up_dw, size=1, color = 'purple3') +
  # annotate("text", x = 1.2, y = -2.2, label= paste(pct_up_dw,"%"), size = 5, color = 'purple3')  +
  # geom_point(data = cmb_dt_filter_padj_dw_dw, size=1, color = 'cyan3') +
  # annotate("text", x = -1, y = -2.2, label= paste(pct_dw_dw,"%"), size = 5, color = 'cyan3')  +
  # geom_point(data = cmb_dt_filter_padj_dw_up, size=1, color = 'navy') +
  # annotate("text", x = -1, y = 1.5, label= paste("0 %"), size = 5, color = 'navy')  +
  theme(aspect.ratio = 1, 
        panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title = element_text(size = pointSize, colour = "black"),
        axis.text.x=element_text(colour="black",size =pointSize),
        axis.text.y=element_text(colour="black",size =pointSize),
        axis.ticks=element_line(colour="black", size = lineWidth),
        plot.margin=unit(c(1,1,1,1),"line"),
        legend.position = "right",
        plot.title = element_text(colour="black",size =pointSize)) +
  # xlim(-5,5)+
  # ylim(0,20)+
  geom_vline(xintercept = 0,colour="grey", linetype = "longdash") +
  geom_hline(yintercept = 0, colour="grey", linetype = "longdash") +
  labs(title = str_wrap("LPS gene & TE correlation",60),
       x = "Log2FC Genes",
       y = "Log2FC TEs (loci)") #+
  # geom_text_repel(data = result,
  #                 color = 'black', size = 5,
  #                 min.segment.length = 0.5,fontface = 'italic',
  #                 #label.padding = unit(0.2,"lines"),
  #                 # nudge_y = 1,
  #                 # nudge_x = 2.5,
  #                 segment.angle = 20,
  #                 max.overlaps = 15) 
scatter 


ggsave(
  paste0('plots/TE_gene_correlation/',file_prefix, '/','Scatter_gene_TE_', Pval_cutoff, '.pdf'),
  plot = scatter,
  width = 100,
  height = 100,
  units = c("mm"),
  dpi = 300)


#Plot number of overlapping features:

# Get unique TE IDs that are in proximity to DE genes
tes_in_overlap <- unique(result$TE)

# Label TEs as overlapping or not
tes_overlap_df <- tes %>%
  mutate(overlap_status = ifelse(TE %in% tes_in_overlap, "Near DE Gene", "Not Near DE Gene"))

# Count categories
te_counts <- tes_overlap_df %>%
  count(overlap_status)

pieplot = ggplot(te_counts, aes(x = "", y = n, fill = overlap_status)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "Proximity of TEs to DE Genes") +
  theme_void() +
  geom_text(aes(label = paste0(round(n / sum(n) * 100), "%")),
            position = position_stack(vjust = 0.5)) +
  scale_fill_brewer(palette = "Set2")

ggsave(
  paste0('plots/TE_gene_correlation/',file_prefix, '/','Pieplot_TE_close_far_genes_', Pval_cutoff, '.pdf'),
  plot = pieplot,
  width = 100,
  height = 100,
  units = c("mm"),
  dpi = 300)


#Pieplot of TE families close and far away
# Define overlap group
tes_overlap_annot <- tes %>%
  mutate(overlap_status = ifelse(TE %in% unique(result$TE), "Near DEG", "Far DEG"))

# Check column name for family: adjust if needed (e.g., "Repeat_family", "Family", etc.)
family_col <- "Class"  # change this to your actual column name

te_family_pct <- tes_overlap_annot %>%
  group_by(overlap_status, !!sym(family_col)) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(overlap_status) %>%
  mutate(percent = 100 * n / sum(n))

stackbar = ggplot(te_family_pct, aes(x = overlap_status, y = percent, fill = .data[[family_col]])) +
  geom_bar(stat = "identity") +
  labs(title = "TE Family Composition by Proximity to DE Genes",
       x = "", y = "Percentage of TEs") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

ggsave(
  paste0('plots/TE_gene_correlation/',file_prefix, '/','Barplot_Percentage_TE_class_close_far_genes_', Pval_cutoff, '.pdf'),
  plot = stackbar,
  width = 100,
  height = 100,
  units = c("mm"),
  dpi = 300)


#Vizualize location for gene of interest
#Find out top genes according to number of DETEs deregulated:
#Count number of nearby DE TEs per gene
top_genes <- result %>%
  count(gene, sort = TRUE) %>%
  top_n(20, n)  # adjust number as needed
top_gene_list <- top_genes$gene

#Check downregulated genes
downregulated_genes <- DEG_data %>%
  filter(log2FoldChange < 0) %>%
  pull(gene)
downreg_result <- result %>%
  filter(gene %in% downregulated_genes)
top_down_genes <- downreg_result %>%
  count(gene, sort = TRUE) %>%
  top_n(10, n)  # top 10 with most nearby TEs
top_down_gene_list <- top_down_genes$gene

#Set output path for
output_path <- paste0("plots/TE_gene_correlation/", file_prefix, "/Gviz_100000/")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

#Set function
plot_gene_te_region <- function(gene_name, genes_gr, DEG_data, tes, file_prefix, flank = 100000) {
  # Get gene GRanges
  gene_gr <- genes_gr[genes_gr$gene_name == gene_name]
  
  if (length(gene_gr) == 0) {
    warning(paste("Gene not found:", gene_name))
    return(NULL)
  }
  
  chr <- as.character(seqnames(gene_gr))
  window_start <- start(gene_gr) - flank
  window_end   <- end(gene_gr) + flank
  
  # Filter nearby TEs
  tes_window <- tes %>%
    filter(Chromosome == chr,
           Start >= window_start,
           Stop <= window_end)
  
  te_window_gr <- GRanges(
    seqnames = tes_window$Chromosome,
    ranges = IRanges(start = tes_window$Start, end = tes_window$Stop),
    te_id = tes_window$TE,
    log2FC = tes_window$log2FoldChange
  )
  
  # Assign color to TEs
  te_colors <- ifelse(te_window_gr$log2FC > 0, "firebrick", "navy")
  
  # Gene color based on log2FC
  gene_log2FC <- DEG_data %>%
    filter(gene == gene_name) %>%
    pull(log2FoldChange)
  
  gene_color <- ifelse(gene_log2FC > 0, "firebrick", "navy")
  
  # Create Gviz tracks
  axisTrack <- GenomeAxisTrack()
  geneTrack <- AnnotationTrack(
    range = gene_gr,
    name = "Gene",
    genome = "mm10",
    fill = gene_color,
    col = NA,
    stacking = "squish",
    feature = gene_name,
    showFeatureId = TRUE,
    fontcolor.feature = "black"
  )
  teTrack <- AnnotationTrack(
    range = te_window_gr,
    genome = "mm10",
    name = "DE TEs",
    fill = te_colors,
    col = NA,
    stacking = "full",
    feature = te_window_gr$te_id,
    showFeatureId = TRUE,
    cex.feature = 0.7,
    fontcolor.feature = "black"
  )
  
  # Plot and save
  output_path2 <- paste0(output_path, gene_name, ".pdf")
  pdf(output_path2, width = 8, height = 4)
  plotTracks(
    list(axisTrack, geneTrack, teTrack),
    from = window_start,
    to = window_end,
    chromosome = chr,
    main = paste("Gene & TE map:", gene_name),
    cex.title = 1.2,
    cex.axis = 0.8
  )
  dev.off()
}

# Set list of genes to plot
genes_to_plot <- c(top_gene_list)  # <-- replace with your own

# Loop through and plot each gene
for (gene in genes_to_plot) {
  plot_gene_te_region(
    gene_name = gene,
    genes_gr = genes_gr,
    DEG_data = DEG_data,
    tes = tes,
    file_prefix = file_prefix
  )
}



##### ----- PERMUTATION TEST ----- #####
maxgap <- 100000  # window around gene

# Get observed overlaps
observed_hits <- findOverlaps(genes_gr_filtered, te_gr, maxgap = maxgap, ignore.strand = TRUE)
observed_te_ids <- unique(mcols(te_gr[subjectHits(observed_hits)])$te_id)
observed_count <- length(observed_te_ids)

set.seed(42)
n_perm <- 1000
perm_counts <- numeric(n_perm)

# List of all non-DE gene names
all_gene_names <- unique(mcols(genes_gr)$gene_name)
non_de_genes <- setdiff(all_gene_names, mcols(genes_gr_filtered)$gene_name)

for (i in seq_len(n_perm)) {
  # Randomly sample same number of genes as in DE set
  sampled_genes <- sample(non_de_genes, length(genes_gr_filtered))
  
  sampled_gr <- genes_gr[mcols(genes_gr)$gene_name %in% sampled_genes]
  
  # Find TEs near random genes
  random_hits <- findOverlaps(sampled_gr, te_gr, maxgap = maxgap, ignore.strand = TRUE)
  unique_tes <- unique(mcols(te_gr[subjectHits(random_hits)])$te_id)
  perm_counts[i] <- length(unique_tes)
}

# Calculate p-value: how often do permuted counts >= observed?
pval <- mean(perm_counts >= observed_count)
x_max <- max(observed_count, max(perm_counts)) + 20  # or use 600 for round number

# Plot
# Set output file path
out_file <- paste0('plots/TE_gene_correlation/', file_prefix, '/Histogram_permutation_DETE_DEG_', Pval_cutoff, '.pdf')
pdf(out_file, width = 5, height = 4)
hist(perm_counts,
     breaks = 50,
     main = "Permutation Test: DE TE Enrichment Near DE Genes",
     xlab = "Number of TEs near random genes",
     col = "lightgrey",
     xlim = c(min(perm_counts), x_max)) 
abline(v = observed_count, col = "red", lwd = 2)
legend("topright",legend = paste("Observed =", observed_count,"\nP =", signif(pval, 3)),text.col = "red")
dev.off()





##### ----- DENSITY PLOT ----- #####
# For DE TEs
dist_to_deg <- distanceToNearest(te_gr, genes_gr_filtered, ignore.strand = TRUE)
de_distances <- mcols(dist_to_deg)$distance

# For DE TEs vs all genes (null)
dist_to_all <- distanceToNearest(te_gr, genes_gr, ignore.strand = TRUE)
null_distances <- mcols(dist_to_all)$distance

# Plot
df <- data.frame(
  distance = c(de_distances, null_distances),
  group = rep(c("To DE Genes", "To All Genes"), c(length(de_distances), length(null_distances)))
)

densityplot = ggplot(df, aes(x = distance, fill = group)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +
  labs(title = "Distance from DE TEs to Nearest Gene", x = "Distance (log10 bp)", y = "Density") +
  theme_minimal()


ggsave(
  paste0('plots/TE_gene_correlation/',file_prefix,'/','Densityplot_TE_class_close_far_genes_', Pval_cutoff, '.pdf'),
  plot = densityplot,
  width = 100,
  height = 100,
  units = c("mm"),
  dpi = 300)


#####
# Parameters
flank <- 100000     # 100kb upstream/downstream of gene
bin_size <- 10000   # 10kb bins

# Step 1: Create windows around DE genes
gene_windows <- promoters(genes_gr_filtered, upstream = flank, downstream = flank)

# Step 2: Bin windows into fixed intervals relative to TSS
bins_list <- lapply(seq_along(gene_windows), function(i) {
  gene <- gene_windows[i]
  tss <- ifelse(strand(gene) == "+", start(gene), end(gene))
  rel_bins <- seq(-flank, flank - bin_size, by = bin_size)
  
  bin_starts <- if (strand(gene) == "+") tss + rel_bins else tss - rev(rel_bins)
  bin_ends   <- bin_starts + bin_size - 1
  
  GRanges(seqnames = seqnames(gene),
          ranges = IRanges(start = bin_starts, end = bin_ends),
          bin = rel_bins + bin_size / 2)
})

all_bins <- do.call(c, bins_list)

# Step 3: Count overlaps with DE TEs
overlaps <- countOverlaps(all_bins, te_gr)
all_bins$overlap_count <- overlaps

# Step 4: Aggregate across all genes by bin position
df <- as.data.frame(mcols(all_bins))
df_summary <- df %>%
  group_by(bin) %>%
  summarise(mean_count = mean(overlap_count))

# Step 5: Plot
ggplot(df_summary, aes(x = bin / 1000, y = mean_count)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Distance from TSS (kb)", y = "Mean TE Overlap Count",
       title = "TE Density Around DE Genes (±100kb)") +
  theme_minimal()

