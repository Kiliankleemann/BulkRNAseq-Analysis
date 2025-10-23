#### -------- INSTALLING AND LOADING PACKAGES ------- ####
list.of.packages <- c("BiocGenerics","tximport","S4Vectors", "DESeq2", "biomaRt","data.table",
                      "ggplot2", "ggsignif", "ggpubr", "sva", "devtools", "org.Hs.eg.db", "ggnewscale",
                      "org.Mm.eg.db", "limma","stringr","KEGGREST","ggrepel", "openxlsx", 'splitstackshape',
                      "fgsea","clusterProfiler","pheatmap","ggpubr","cowplot",'dplyr','janitor',
                      "RColorBrewer",'AnnotationDbi', "igraph","ggraph","enrichplot",
                      'tidyverse','pheatmap', 'dendextend',"factoextra")


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if(length(new.packages)) BiocManager::install(new.packages)

# Packages loading
invisible(lapply(list.of.packages, library, character.only = TRUE))

###### ----- FILE PREPARATION -----#####
file_prefix <- "TETranscript_multi_all_limma"

###### ----- SETTING WORK DIRECTORY -----#####
setwd("/Users/kiliankleemann/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Human_RTI_study/Bulk_RNA_seq")

#Setting up directories for analysis
dir.create(paste('results/Pathway_analysis/'))
dir.create(paste('results/Pathway_analysis/FGSEA_enrichment/'))
dir.create(paste('results/Pathway_analysis/KEGG_enrichment/'))
dir.create(paste('results/Pathway_analysis/GO_enrichment/'))

dir.create(paste('Pathway_analysis/'))
dir.create(paste('Pathway_analysis/FGSEA_enrichment/'))
dir.create(paste('Pathway_analysis/KEGG_enrichment/'))
dir.create(paste('Pathway_analysis/GO_enrichment/'))

#Setting up organism data base for enrichment analysis
organism_name <-  'hsa' #human: hsa
OrgDb_name    <-  "org.Hs.eg.db" #human: 'org.Hs.eg.db'

#Setting up filter for comparison
comparison_prefix <- "V3_vs_V2"

###### ----- LOADING DATA FOR GSEA ####
pathway_data <- read.xlsx(paste0('results/',file_prefix, '/limma_AllResults.xlsx')) %>% 
  filter(!grepl(':',gene)) %>% filter(comparison == comparison_prefix)

DEgenes <- pathway_data %>% pull('gene') #%>% filter(P.Value < 0.05)
#DEgenes <- mapIds(org.Mm.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')

DEG_Log2FC <-  pathway_data  %>%  pull('t') #%>% filter(P.Value < 0.05)

names(DEG_Log2FC) = DEgenes
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
DEG_Log2FC

#Reversing for opposite plotting of GSEA plot:
reverse_direction <- TRUE   # set FALSE for normal orientation
if (reverse_direction) {
  DEG_Log2FC <- -DEG_Log2FC
}
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)

###### --------- CURATED pathways ------- #######
pathways.all.curated <- gmtPathways("GSEA_database/c2.all.v2025.1.Hs.symbols.gmt") 
fgseaRes <- fgsea(pathways=pathways.all.curated, 
                  stats=DEG_Log2FC,
                  minSize   = 15,
                  maxSize   = 500)

fgseaResTidy_curated <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(padj<0.05)

dir.create(paste0('results/Pathway_analysis/FGSEA_enrichment/',file_prefix))
write.xlsx(fgseaResTidy_curated, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_prefix,'/Curated_pathways_',comparion_prefix,'.xlsx'))

###### --------- CANONICAL pathways ------- #######
pathways.all.canonical <- gmtPathways("GSEA_database/c2.cp.v2025.1.Hs.symbols.gmt") 
fgseaRes <- fgsea(pathways=pathways.all.canonical, 
                  stats=DEG_Log2FC,
                  minSize   = 15,
                  maxSize   = 500)

fgseaResTidy_canonical <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(padj<0.05)

dir.create(paste0('results/Pathway_analysis/FGSEA_enrichment/',file_prefix))
write.xlsx(fgseaResTidy_curated, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_prefix,'/Canonical_pathways_',comparion_prefix,'.xlsx'))


###### --------- REACTOME pathways ------- #######
pathways.all.reactome <- gmtPathways("GSEA_database/c2.cp.reactome.v2025.1.Hs.symbols.gmt") 
fgseaRes <- fgsea(pathways=pathways.all.reactome, 
                  stats=DEG_Log2FC,
                  minSize   = 15,
                  maxSize   = 500)

fgseaResTidy_reactome <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(padj<0.05)

dir.create(paste0('results/Pathway_analysis/FGSEA_enrichment/',file_prefix))
write.xlsx(fgseaResTidy_reactome, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_prefix,'/Reactome_pathways_',comparison_prefix,'.xlsx'))


####### --------- BARPLOT  ----- #####
Canonical_pathway_plot_data <- read.xlsx(paste0('results/Pathway_analysis/FGSEA_enrichment/',file_prefix,'/Canonical_pathways_',comparion_prefix,'.xlsx'), sheet =2)

Canonical_pathway_plot_data$log_padj <- -log10(Canonical_pathway_plot_data$padj)
FGSEA_data <- as.data.frame(Canonical_pathway_plot_data)

#Up
col_max <- 'firebrick'
col_min <- 'navyblue'

#Fgsea plot theme
lineWidth = 0.5
pointSize = 10
theme_gsea_bar <- theme(
  text = element_text(size = pointSize, colour = "black"),
  rect = element_blank(),
  line = element_line(size = lineWidth, colour = "black"),
  plot.title  = element_text(color="black", size = pointSize),
  axis.title  = element_text(size = pointSize, colour = "black"),
  axis.text.x  = element_text(size = pointSize, colour = "black"),
  axis.text.y  = element_text(size = pointSize, colour = "black"),
  axis.ticks.x = element_line(size = lineWidth, colour = "black"),
  axis.ticks.y = element_line(size = lineWidth, colour = "black"),
  axis.line = element_line(size = lineWidth, colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title = element_text(size = pointSize , colour = "black"),
  legend.text = element_text(size = pointSize, colour = "black"),
  legend.key.height = unit(0.5, "cm"),
  legend.key.width = unit(0.5, "cm"))

#Barplot
barplot_curated <- ggplot(FGSEA_data, aes(NES,reorder(pathway, NES))) +
  geom_col(aes(fill= NES), width = 0.5,color = 'black') +
  geom_point(aes(size = log_padj))+
  theme_gsea_bar +
  scale_fill_gradient2(low = col_min, high = col_max) +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= str_wrap("Curated pathways NES from GSEA",width =30)) + 
  geom_vline(xintercept = 0, color ='black') +
  scale_x_continuous(expand = c(0, 0.1, 0.1, 0))

barplot_curated

dev.off()

dir.create(paste0('plots/pathways/'))
dir.create(paste0('plots/pathways/FGSEA_pathways/'))
dir.create(paste0('plots/pathways/FGSEA_pathways/',file_prefix))
ggsave(paste0('plots/pathways/FGSEA_pathways/', file_prefix,'/Barplot_canonical_pathways_',comparion_prefix,'.pdf'),
       barplot_curated,
       height = 10,
       width = 15,
       dpi = 300)


####### --------- CONNECTOME ------ ######
# Clean and extract gene sets ! Make sure to simplify names *before* computing similarity
gsea <- fgseaResTidy_reactome
gsea <- gsea %>%
  rename(name = pathway)

# Simplify names
gsea$name <- gsea$name %>%
  gsub("^REACTOME_", "", .) %>%
  gsub("^HALLMARK_", "", .) %>%
  gsub("^KEGG_", "", .) %>%
  gsub("^GO_", "", .) %>%
  gsub("^WP_", "", .) %>%
  gsub("_PATHWAY$", "", .) %>%
  gsub("_", " ", .) %>%
  stringr::str_to_title()

# Extract leading edge genes
gene_sets <- strsplit(as.character(gsea$leadingEdge), "/|,")
names(gene_sets) <- gsea$name  # <— now uses cleaned names

# Compute Jaccard similarity
pathways <- names(gene_sets)
n <- length(pathways)
similarity_matrix <- matrix(0, n, n, dimnames = list(pathways, pathways))

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    overlap <- length(intersect(gene_sets[[i]], gene_sets[[j]]))
    union <- length(unique(c(gene_sets[[i]], gene_sets[[j]])))
    similarity_matrix[i, j] <- similarity_matrix[j, i] <- ifelse(union > 0, overlap / union, 0)
  }
}

# Build edge list based on similarity threshold
threshold <- 0.1  # adjust for desired connection density
edges <- as.data.frame(as.table(similarity_matrix)) %>%
  filter(Freq > threshold & Var1 != Var2) %>%  # exclude self-edges
  rename(from = Var1, to = Var2, weight = Freq)

# Create graph object
g <- graph_from_data_frame(edges, vertices = gsea, directed = FALSE)

#Adjust padj value
gsea$logP <- -log10(pmax(gsea$padj, 1e-300))

# Match vertex data to graph
V(g)$logP <- gsea$logP[match(V(g)$name, gsea$name)]
V(g)$NES  <- gsea$NES[match(V(g)$name, gsea$name)]

# Plot
set.seed(42)
plot_g <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = weight, alpha = weight), color = "grey80") +
  geom_node_point(aes(size = logP, color = NES)) +
  scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
  scale_size_continuous(name = expression(-log[10](padj)),
                        range = c(3, 12)) +
  scale_edge_width(range = c(0.2, 2)) +
  theme_void() +
  theme(legend.position = "right")

# Save
ggsave(
  paste0("plots/pathways/FGSEA_pathways/", file_prefix, "/GGraph_canonical_pathways_no_label_", comparison_prefix, ".pdf"),
  plot_g,
  height = 7,
  width = 7,
  dpi = 300
)


plot_g <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = weight, alpha = weight), color = "grey80") +
  geom_node_point(aes(size = logP, color = NES)) +
  geom_node_text(aes(label = name),repel = TRUE, size = 1, max.overlaps = Inf) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  scale_size_continuous(name = expression(-log[10](padj)),
                        range = c(3, 12)) +
  scale_edge_width(range = c(0.2, 2)) +
  theme_void() +
  theme(legend.position = "right")

# Save
ggsave(
  paste0("plots/pathways/FGSEA_pathways/", file_prefix, "/GGraph_canonical_pathways_all_label_", comparison_prefix, ".pdf"),
  plot_g,
  height = 10,
  width = 10,
  dpi = 300
)

# Load manually selected top labels
top_labels <- read.xlsx(
  paste0("results/Pathway_analysis/FGSEA_enrichment/", file_prefix, "/Reactome_pathways_", comparison_prefix, ".xlsx"),
  sheet = 2) %>% pull("Labels")

# Simplify top labels
top_labels <- top_labels %>%
  gsub("^REACTOME_", "", .) %>%
  gsub("^HALLMARK_", "", .) %>%
  gsub("^KEGG_", "", .) %>%
  gsub("^GO_", "", .) %>%
  gsub("^WP_", "", .) %>%
  gsub("_PATHWAY$", "", .) %>%
  gsub("_", " ", .) %>%
  stringr::str_to_title()

# Plot
set.seed(42)
plot_g <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = weight, alpha = weight), color = "grey80") +
  geom_node_point(aes(size = logP, color = NES)) +
  geom_node_text(aes(label = ifelse(name %in% top_labels, name, "")),
                 repel = TRUE, size = 3, max.overlaps = Inf) +
  scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
  scale_size_continuous(name = expression(-log[10](padj)),
                        range = c(3, 12)) +
  scale_edge_width(range = c(0.2, 2)) +
  theme_void() +
  theme(legend.position = "right")

# Save
ggsave(
  paste0("plots/pathways/FGSEA_pathways/", file_prefix, "/GGraph_canonical_pathways_selected_label_", comparison_prefix, ".pdf"),
  plot_g,
  height = 7,
  width = 8,
  dpi = 300
)


####### --------- HEATMAP - FGSEA EXTRACTING SPECIFIC PATHWAYS ---- ####
toMatch <- c("REACTOME_DNA_REPAIR", "REACTOME_HIV_INFECTION","REACTOME_FCGAMMA_RECEPTOR_FCGR_DEPENDENT_PHAGOCYTOSIS")

pathway_data <- fgseaResTidy_reactome #%>%
  filter(grepl(paste(toMatch, collapse = "|"), pathway, ignore.case = TRUE))

pathway_ID <- fgseaResTidy_reactome$pathway

#Extract genes and make heatmap for pathways matching term
cluster_no = 2
gaps = c(11,22)
data <- read.xlsx(paste0('results/',file_prefix, '/limma_DEG_counts_pval05.xlsx'))

dir.create(paste0("plots/pathways/FGSEA_pathways/", file_prefix, "/Heatmap_pval05/"))
outdir <- paste0("plots/pathways/FGSEA_pathways/", file_prefix, "/Heatmap_pval05/")

dev.off()
for (id in pathway_ID) {
  pathway_data_select <- pathway_data %>% filter(pathway == id)
  pathway_data_genes <- unique(unlist(pathway_data_select$leadingEdge))
  
  # subset expression data
  pathway_data_selected <- data %>%
    filter(gene %in% pathway_data_genes)

  if (nrow(pathway_data_selected) == 0) {
    message(paste("No matching genes found for", id))
    next
  }
  
  pathway_data_selected <- pathway_data_selected %>%
    tibble::column_to_rownames("gene")
  
  phm_kegg <- pheatmap(pathway_data_selected,
                       color = colorRampPalette(color_set)(40),
                       breaks = seq(-2, 2, by = 0.1),
                       kmeans_k = NA,
                       cluster_rows = F,
                       cutree_row = cluster_no,
                       cluster_cols = F,
                       gaps_col = gaps,             
                       legend = TRUE,
                       show_rownames = T,
                       border_color = 'NA',
                       fontsize = 10,
                       treeheight_col = 50,
                       treeheight_row = 50,
                       cellwidth = 10,
                       scale = 'row',
                       cellheight = 10,
                       main = id
  )
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    paste0(outdir, "Heatmap_", gsub("[^A-Za-z0-9]", "_", id), ".pdf"),
    phm_kegg,
    height = 10,
    width = 7,
    dpi = 300
  )
  
}


####### --------- DOTPLOT for overlapping Gene-sets ------ #####
# Import pathway results
Reactome_V3_vs_V2 <- read.xlsx(
  paste0("results/Pathway_analysis/FGSEA_enrichment/", file_prefix, "/Reactome_pathways_V3_vs_V2.xlsx"))

Reactome_V4_vs_V3 <- read.xlsx(
  paste0("results/Pathway_analysis/FGSEA_enrichment/", file_prefix, "/Reactome_pathways_V4_vs_V3.xlsx"))

#Clean names
clean_pathway_names <- function(df) {
  df %>%
    rename(name = pathway) %>%
    mutate(name = name %>%
             gsub("^REACTOME_|^HALLMARK_|^KEGG_|^GO_", "", .) %>%
             gsub("_PATHWAY$", "", .) %>%
             gsub("_", " ", .) %>%
             stringr::str_to_title())
}
fgsea_A <- clean_pathway_names(Reactome_V3_vs_V2)
fgsea_B <- clean_pathway_names(Reactome_V4_vs_V3)

#Merge data and select
fgsea_A$group <- "V3_vs_V2"
fgsea_B$group <- "V4_vs_V3"
merged_data_pathway <- rbind(fgsea_A, fgsea_B)

merged_data_pathway_select <- merge(fgsea_A,fgsea_B, by = "name") %>% pull("name")

merged_data_pathway_overlap <- merged_data_pathway %>% filter(name %in% merged_data_pathway_select)
merged_data_pathway_overlap$logpadj <- -log(merged_data_pathway_overlap$padj)
merged_data_pathway_overlap$NES <- as.numeric(merged_data_pathway_overlap$NES)

selected_pathways <- c(
  "Mitochondrial Translation",
  "Neutrophil Degranulation",
  "Sars Cov Infections",
  "Hiv Infection",
  "Beta Defensins",
  "Autophagy",
  "Eukaryotic Translation Initiation",
  "Signaling By The B Cell Receptor Bcr", 
  "Interferon Signaling",
  "Fc Epsilon Receptor Fceri Signaling",
  "Interleukin 1 Family Signaling",
  "Interleukin 17 Signaling",
  "Chromatin Organization")

merged_data_pathway_overlap_select <- merged_data_pathway_overlap %>% filter(name %in% selected_pathways)

lineWidth = 0.5
pointSize = 10

col_max <- max(merged_data_pathway_overlap_select$NES)
#Dotplot
dotplot <- ggplot(merged_data_pathway_overlap_select, aes(group,reorder(name, NES), fill = NES)) +
  geom_point(aes(size = logpadj),  shape = 21, color = 'black') +
  scale_fill_gradient2(limits=c(-col_max, col_max), low="navyblue", mid="whitesmoke", high = "firebrick", na.value = 'navyblue') +
  theme(
    text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(size = lineWidth, colour = "black"),
    plot.title  = element_text(color="black", size = pointSize),
    axis.title  = element_text(size = pointSize, colour = "black"),
    axis.text.x  = element_text(size = pointSize, colour = "black"),
    axis.text.y  = element_text(size = pointSize, colour = "black"),
    axis.ticks.x = element_line(size = lineWidth, colour = "black"),
    axis.ticks.y = element_line(size = lineWidth, colour = "black"),
    axis.line = element_line(size = lineWidth, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = pointSize , colour = "black"),
    legend.text = element_text(size = pointSize, colour = "black"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm")) +
  labs(x="Group comparison", y="Pathways overlapping",
       title= str_wrap("Reactome pathways NES from GSEA",width =30)) + 
  geom_vline(xintercept = 0, color ='black') 

dotplot

# Save
ggsave(
  paste0("plots/pathways/FGSEA_pathways/", file_prefix, "/Dotplot_reactome_pathways_overlap_selected", ".pdf"),
  dotplot,
  height = 5,
  width = 5,
  dpi = 300
)



####### --------- GSEA plot for specific pathways ------- #######
# ---------------- USER INPUT ----------------
comparison_prefix <- "V2_vs_V3"
# Pathway and highlights
geneset_name      <- "REACTOME_INTERFERON_SIGNALING"
highlight_genes   <- c("ISG15", "IRF7", "IFI27")

# Style parameters
heatmap_height <- 0.1   # height of heatmap/ticks
lineWidth      <- 0.5   # line thickness for axes, ticks, hline
textSize       <- 20    # base text size
ylim_min <- 0   # adjust to your preferred scale
ylim_max <- 0.4

# ---- Extract adjusted p-value (and optionally NES) ----
padj_value <- fgseaResTidy_reactome %>%
  filter(pathway == geneset_name) %>%
  pull(padj)

nes_value <- fgseaResTidy_reactome %>%
  filter(pathway == geneset_name) %>%
  pull(NES) 

# ---- Data preparation ----------------
ranks <- sort(DEG_Log2FC, decreasing = TRUE)
geneset <- intersect(pathways.all.reactome[[geneset_name]], names(ranks))

# Enrichment curve data from fgsea
enrichment_data <- plotEnrichment(geneset, ranks)
df_es <- ggplot_build(enrichment_data)$data[[1]] |> 
  dplyr::select(x, y)

# Rank metric dataframe
df_rank <- data.frame(rank = 1:length(ranks), score = ranks)
pathway_positions <- match(geneset, names(ranks))

# Highlight genes and positions
highlight_genes <- intersect(highlight_genes, geneset)
highlight_pos   <- match(highlight_genes, names(ranks))

# Label dataframe for gene text
label_df <- data.frame(
  gene = highlight_genes,
  x = highlight_pos,
  y = rep(heatmap_height / 2 + 0.03, length(highlight_pos))
)
# ---- Make 11-bin, contrast-enhanced band colors ----
lim <- quantile(df_rank$score, probs = c(0.01, 0.99), na.rm = TRUE)
df_rank$score_clipped <- pmax(pmin(df_rank$score, lim[2]), lim[1])
df_rank$score_scaled  <- sign(df_rank$score_clipped) * abs(df_rank$score_clipped)^(1/2)

range_limit  <- max(abs(df_rank$score_scaled), na.rm = TRUE)
break_points <- seq(-range_limit, range_limit, length.out = 12)  # 11 bins
df_rank$score_bin <- cut(df_rank$score_scaled, breaks = break_points, include.lowest = TRUE)
bin_levels <- levels(df_rank$score_bin)

color_bins <- colorRampPalette(c("navy", "white", "firebrick"))(length(bin_levels))
names(color_bins) <- bin_levels

# ---- Build a left-to-right raster for the band ----
band_cols <- color_bins[as.character(df_rank$score_bin)]            # vector of colors per rank
band_img  <- as.raster(matrix(band_cols, nrow = 1))                  # 1 row, many columns (left→right)

# A tiny dataset to carry the legend (invisible tile)
legend_df <- data.frame(
  rank = seq_along(bin_levels),
  score_bin = factor(bin_levels, levels = bin_levels),
  y = 0
)

# ---- Enhanced scaling (as before) ----
lim <- quantile(df_rank$score, probs = c(0.01, 0.99), na.rm = TRUE)
df_rank$score_clipped <- pmax(pmin(df_rank$score, lim[2]), lim[1])
df_rank$score_scaled  <- sign(df_rank$score_clipped) * abs(df_rank$score_clipped)^(1/2)

range_limit  <- max(abs(df_rank$score_scaled), na.rm = TRUE)
break_points <- seq(-range_limit, range_limit, length.out = 12)
df_rank$score_bin <- cut(df_rank$score_scaled, breaks = break_points, include.lowest = TRUE)
bin_levels <- levels(df_rank$score_bin)

color_bins <- colorRampPalette(c("navy", "white", "firebrick"))(length(bin_levels))
names(color_bins) <- bin_levels

# ---- Build raster for narrow band around y=0 ----
band_cols <- color_bins[as.character(df_rank$score_bin)]
band_img  <- as.raster(matrix(band_cols, nrow = 1))  # 1 row, many columns (left→right)

# Legend data — tiny dummy tile just for legend, placed outside plot range
legend_df <- data.frame(
  rank = seq_along(bin_levels),
  score_bin = factor(bin_levels, levels = bin_levels),
  y = max(df_es$y) + 10   # safely outside view
)

# ---- Y range for ES curve (tight scaling) ----
y_range <- range(df_es$y, na.rm = TRUE)
y_margin <- diff(y_range) * 0.05   # small top/bottom padding

# ---- Plot ####
p <- ggplot(df_es, aes(x = x, y = y)) +
  # (A) Raster band around y=0
  annotation_raster(
    band_img,
    xmin = min(df_rank$rank),
    xmax = max(df_rank$rank),
    ymin = -heatmap_height/2,
    ymax =  heatmap_height/2,
    interpolate = FALSE
  ) +
  
  # (B) Legend generator (opaque but off-screen)
  geom_tile(
    data = legend_df,
    aes(x = rank, y = y, fill = score_bin),
    show.legend = TRUE, inherit.aes = FALSE
  ) +
  
  scale_x_continuous(expand = c(0, 0, 0, 0)) +
  coord_cartesian(ylim = c(ylim_min, ylim_max)) +
  
  # (C) Pathway ticks
  geom_segment(
    data = data.frame(x = pathway_positions),
    aes(x = x, xend = x,
        y = -heatmap_height / 4,
        yend =  heatmap_height / 4),
    color = "black", alpha = 0.5, linewidth = lineWidth * 0.3
  ) +
  
  # (D) Highlighted ticks
  geom_segment(
    data = data.frame(x = highlight_pos),
    aes(x = x, xend = x,
        y = -heatmap_height / 4,
        yend =  heatmap_height / 4),
    color = "goldenrod", linewidth = lineWidth * 0.8, alpha = 0.9
  ) +
  
  # (E) Enrichment curve and baseline
  geom_line(color = "black", linewidth = lineWidth) +
  geom_hline(yintercept = 0, color = "grey80", linetype = "dashed", linewidth = lineWidth) +
  
  # (F) Labels for highlighted genes
  geom_text_repel(
    data = label_df,
    aes(x = x, y = y, label = gene),
    color = "goldenrod",
    size = textSize * 0.25,
    angle = 0,
    vjust = 0.1,
    nudge_y = 0.02,
    segment.size = 0,
    direction = "y",
    box.padding = 0.05,
    max.overlaps = Inf
  ) +
  
  # (G) Add padj and NES text annotation
  annotate(
    "text",
    x = max(df_rank$rank) * 0.98,
    y = ylim_max - (ylim_max - ylim_min) * 0.05,
    label = paste0("NES = ", nes_value, "\nFDR = ", padj_value),
    hjust = 1,
    vjust = 1,
    size = textSize * 0.25,
    color = "black"
  ) +
  
  # (H) Discrete color legend
  scale_fill_manual(
    values = color_bins,
    drop = FALSE,
    name = "Rank metric"
  ) +
  
  # (I) Theme and axis styling
  theme_minimal(base_size = textSize) +
  labs(
    title = geneset_name,
    x = "Gene Rank (sorted by statistic)",
    y = "Running Enrichment Score"
  ) +
  theme(
    plot.title = element_text(size = textSize, color = "black", face = "plain", hjust = 0.5),
    axis.title = element_text(color = "black", face = "plain"),
    axis.text  = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.line  = element_line(color = "black", linewidth = lineWidth),
    axis.ticks = element_line(color = "black", linewidth = lineWidth),
    legend.position = "right"
  )

print(p)

# ---- Export (flat RGB, crisp colors) ----
ggsave(
  filename = file.path(outdir, paste0("GSEAplot_", geneset_name, "_", comparison_prefix, ".pdf")),
  plot = p,
  height = 5, width = 9, dpi = 300,
  device = cairo_pdf,
  bg = "white")





















###### ---- KEGG ENRICHMENT ---- #####
pathway_data <- read.xlsx(paste0('results/',file_prefix, '/DEGene_statistics_pval05.xlsx'))

DEgenes <- pathway_data %>% filter(pvalue < 0.01)  %>% pull('gene')
DEgenes <- mapIds(org.Mm.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')

DEG_Log2FC <-  pathway_data %>% filter(pvalue < 0.01) %>%  pull('log2FoldChange')

names(DEG_Log2FC) = DEgenes
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
DEG_Log2FC

#Run KEGG 
kegg_result <- gseKEGG(geneList = DEG_Log2FC,
                       organism     = organism_name,
                       pAdjustMethod = "BH",
                       minGSSize = 15,
                       maxGSSize = 500,
                       nPermSimple = 1000,
                       pvalueCutoff = 0.05)

kegg_result <- setReadable(kegg_result,OrgDb= 'org.Mm.eg.db', keyType = 'ENTREZID')
#head(kegg_result@geneSets)

#GSEA plot
# gseaNb(object = kk2,
#        geneSetID = 'hsa05169',
#        subPlot = 2,
#        addGene = T,
#        kegg = T,
#        markTopgene = T)
# 
# ggsave(paste0('plots/pathways/KEGG_pathways/', file_export,'.pdf'),
#        gsea_KEGG,
#        dpi = 300)


#Barplot
KEGG_data <- as.data.frame(kegg_result) %>% filter(NES > 0.5 | NES < -1.5)
KEGG_data$log_pval <- log(KEGG_data$pvalue,10)

#Filter names 
KEGG_data$Description <- gsub('Mus','', KEGG_data$Description)
KEGG_data$Description <- gsub('musculus','', KEGG_data$Description)
KEGG_data$Description <- gsub('house','', KEGG_data$Description)
KEGG_data$Description <- gsub('mouse','', KEGG_data$Description)
KEGG_data$Description <- gsub('\\( )','', KEGG_data$Description)
KEGG_data$Description <- gsub('\\-','', KEGG_data$Description)
KEGG_data$Description <- gsub('    ','', KEGG_data$Description)

barplot_Kegg <- ggplot(KEGG_data, aes(NES,reorder(Description, NES))) +
  geom_col(aes(fill= NES), width = 0.5,color = 'black') +
  geom_point(aes(size = log_pval))+
  theme_gsea_bar +
  scale_fill_gradient2(low = col_min, high = col_max) +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= str_wrap("Curated pathways NES from GSEA",width =30)) + 
  geom_vline(xintercept = 0, color ='black') +
  scale_x_continuous(expand = c(0, 0, 0.1, 0))

barplot_Kegg

dev.off()

dir.create(paste0('plots/pathways/'))
dir.create(paste0('plots/pathways/KEGG_pathways/'))
dir.create(paste0('plots/pathways/KEGG_pathways/',file_prefix))
ggsave(paste0('plots/pathways/KEGG_pathways/', file_prefix,'/Barplot_pathways.pdf'),
       barplot_Kegg,
       height = 20,
       width = 30,
       dpi = 300)
dir.create(paste0('Pathway_analysis/KEGG_enrichment/',file_prefix))
write.xlsx(KEGG_data, paste0('Pathway_analysis/KEGG_enrichment/',file_prefix,'/KEGG_pathways.xlsx'))


###### ---- EXTRACTING SPECIFIC PATHWAYS - HEATMAP ---- ####
toMatch <- c("Transposable")
pathway_data <- dplyr::filter(kegg_result,Description %in% unique(grep(paste(toMatch,collapse="|"), kegg_result$Description, value=TRUE))) 
pathway_ID <- pathway_data@result$ID


#Extract genes and make heatmap for pathways matching term
cluster_no = 2
gaps = c(4)
data <- read.xlsx(paste0('results/',file_prefix, '/DEGene_counts_pval05.xlsx'))

remove_species_info <- function(input_strings) {
  cleaned_strings <- gsub('Mus','', input_strings)
  cleaned_strings <- gsub('musculus','', cleaned_strings)
  cleaned_strings <- gsub('house','', cleaned_strings)
  cleaned_strings <- gsub('mouse','', cleaned_strings)
  cleaned_strings <- gsub('Mus','', cleaned_strings)
  cleaned_strings <- gsub('\\( )','', cleaned_strings)
  cleaned_strings <- gsub('\\-','', cleaned_strings)
  cleaned_strings <- gsub('    ','', cleaned_strings)
  return(cleaned_strings)
}

for (id in pathway_ID) {
  pathway_data_genes <- pathway_data[[id]]
  pathway_name <- pathway_data[id] %>% pull(Description)
  pathway_name <- remove_species_info(pathway_name)
  pathway_data_selected <- data %>% filter(gene %in% pathway_data_genes) %>% column_to_rownames('gene')
  phm_kegg <- pheatmap(pathway_data_selected,
                       color = colorRampPalette(color_set)(40),
                       breaks = seq(-2, 2, by = 0.1),
                       kmeans_k = NA,
                       cluster_rows = F,
                       cutree_row = cluster_no,
                       cluster_cols = F,
                       gaps_col = gaps,             
                       legend = TRUE,
                       show_rownames = T,
                       border_color = 'NA',
                       fontsize = 10,
                       treeheight_col = 50,
                       treeheight_row = 50,
                       cellwidth = 10,
                       scale = 'row',
                       cellheight = 10,
                       main = pathway_name
  )
  
  ggsave(paste0('plots/pathways/KEGG_pathways/',file_prefix, '/Heatmap_',pathway_name , '_', 'pval05.pdf'),
         phm_kegg,
         height = 10,
         width = 5,
         dpi = 300)
  
}





###### ---- GO ENRICHMENT ---- #####
pathway_data <- read.xlsx(paste0('results/',file_prefix, '/DEGene_statistics_pval05.xlsx'))

DEgenes <- pathway_data %>% filter(pvalue < 0.05)  %>% pull('gene')
DEgenes <- mapIds(org.Mm.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')

DEG_Log2FC <-  pathway_data %>% filter(pvalue < 0.05) %>%  pull('log2FoldChange')

names(DEG_Log2FC) = DEgenes
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
DEG_Log2FC

#Run GO 
go_result <- gseGO(geneList = DEG_Log2FC,
                       OrgDb     = OrgDb_name,
                       pAdjustMethod = "BH",
                       minGSSize = 15,
                       maxGSSize = 500,
                       nPermSimple = 1000,
                       pvalueCutoff = 0.05)

go_result <- setReadable(go_result, OrgDb= OrgDb_name, keyType = 'ENTREZID')


#Barplot
GO_data <- as.data.frame(go_result) %>% filter(NES > 1.8 | NES < -2.4)
GO_data$log_pval <- log(GO_data$pvalue,10)

barplot_Go <- ggplot(GO_data, aes(NES,reorder(Description, NES))) +
  geom_col(aes(fill= NES), width = 0.5,color = 'black') +
  geom_point(aes(size = log_pval))+
  theme_gsea_bar +
  scale_fill_gradient2(low = col_min, high = col_max) +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= str_wrap("Curated pathways NES from GSEA",width =30)) + 
  geom_vline(xintercept = 0, color ='black') +
  scale_x_continuous(expand = c(0, 0, 0.1, 0))

barplot_Go

dev.off()
dir.create(paste0('plots/pathways/GO_pathways/'))
dir.create(paste0('plots/pathways/GO_pathways/',file_prefix))
ggsave(paste0('plots/pathways/GO_pathways/', file_prefix,'/Barplot_pathways.pdf'),
       barplot_Go,
       height = 20,
       width = 30,
       dpi = 300)

dir.create(paste0('Pathway_analysis/GO_enrichment/',file_prefix))
write.xlsx(GO_data, paste0('Pathway_analysis/GO_enrichment/',file_prefix,'/GO_pathways.xlsx'))

###### ---- EXTRACTING SPECIFIC PATHWAYS  ---- ####
toMatch <- c("immune")
go_result
pathway_data <- dplyr::filter(go_result,Description %in% unique(grep(paste(toMatch,collapse="|"), go_result$Description, value=TRUE))) 
pathway_ID <- pathway_data@result$ID
pathway_ID

cluster_no = 2
gaps = c(4)
data <- read.xlsx(paste0('results/',file_prefix, '/DEGene_counts_pval05.xlsx'))


for (id in pathway_ID) {
  pathway_data_genes <- pathway_data[[id]]
  pathway_name <- pathway_data[id] %>% pull(Description)
  pathway_data_selected <- data %>% filter(gene %in% pathway_data_genes) %>% column_to_rownames('gene')
  phm_kegg <- pheatmap(pathway_data_selected,
                       color = colorRampPalette(color_set)(40),
                       breaks = seq(-2, 2, by = 0.1),
                       kmeans_k = NA,
                       cluster_rows = F,
                       cutree_row = cluster_no,
                       cluster_cols = F,
                       gaps_col = gaps,             
                       legend = TRUE,
                       show_rownames = T,
                       border_color = 'NA',
                       fontsize = 10,
                       treeheight_col = 50,
                       treeheight_row = 50,
                       cellwidth = 10,
                       scale = 'row',
                       cellheight = 10,
                       main = pathway_name
  )
  
  ggsave(paste0('plots/pathways/GO_pathways/',file_prefix, '/Heatmap_',pathway_name , '_', 'pval05.pdf'),
         phm_kegg,
         height = 30,
         units = 'cm',
         width = 12,
         dpi = 300)
  
}







###### ---- HEATMAP CLUSTER ENRICHMENT ANALYSIS ---- ####
lineWidth = 0.1
pointSize = 12
data <- read.xlsx(paste0('results/', file_prefix, 'clustergenes.xlsx'))
clusters <- unique(data$cluster)

for (i in clusters) {
  # Get genes for this cluster and convert to ENTREZID
  cluster_genes <- data %>% filter(cluster == i) %>% pull(gene)
  cluster_genes <- mapIds(org.Hs.eg.db, cluster_genes, 'ENTREZID', 'SYMBOL')
  cluster_genes <- na.omit(cluster_genes)
  
  # Perform enrichment
  cluster_go <- enrichGO(gene = cluster_genes,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         readable = TRUE)
  
  # Skip empty results
  if (is.null(cluster_go) || nrow(cluster_go@result) == 0) {
    cat("No enrichment found for cluster", i, "\n")
    next
  }
  
  # Process results
  p1 <- cluster_go@result %>%
    filter(pvalue < 0.05, Count > 1)
  p1$log_p <- -log10(p1$pvalue)
  color_threshold <- max(p1$log_p) + 1
  
  # Show top N terms (adjustable)
  top_n <- 15
  plot_data <- p1 %>% slice_min(order_by = pvalue, n = min(nrow(p1), top_n))
  
  # Make plot
  dotplot_GO <- ggplot(plot_data, aes(x = Count, y = reorder(Description, log_p))) + 
    geom_point(aes(size = Count, color = log_p), stroke = 1.5, shape = 21, fill = "white") +
    scale_color_gradient2(limits = c(0, color_threshold), low = "white", mid = "white", high = "purple") +
    theme_minimal(base_size = 14) +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title  = element_text(color="black", size = pointSize),
          axis.title  = element_text(size = pointSize, colour = "black"),
          axis.text.x  = element_text(size = pointSize, colour = "black"),
          axis.text.y  = element_text(size = pointSize, colour = "black"),
          axis.ticks.x = element_line(size = lineWidth, colour = "black"),
          axis.ticks.y = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize, colour = "black"),
          legend.key.height = unit(0.3, "cm"),
          legend.key.width = unit(0.3, "cm")) +
    xlab("Gene Count") + ylab(NULL) +
    labs(title = str_wrap(paste0("GO Biological Process – Cluster ", i), 40))
  
  # Save plot with height scaled to number of terms
  ggsave(paste0("Pathway_analysis/GO_enrichment/", file_prefix, "_cluster", i, "_GO_pathway.pdf"),
         dotplot_GO,
         width = 10,
         height = max(4, 0.5 * nrow(plot_data)),  # Adjust height dynamically
         dpi = 300)
}



