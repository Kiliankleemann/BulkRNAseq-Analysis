#Script to investigate repeatmasker before and after filtering
# Load required library
library(data.table)
library(ggplot2)
library(dplyr)

# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)

# Function to read GTF file
read_gtf <- function(file_path) {
  gtf <- fread(file_path, header = FALSE)
  colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  return(gtf)
}

# Function to count records by family ID
count_by_family <- function(gtf) {
  gtf[, family_id := gsub(".*Family:", "", gsub(";.*", "", attribute))]
  family_count <- gtf[, .(count = .N), by = family_id]
  return(family_count)
}

# Function to calculate percentage of each family
calculate_percentage <- function(family_count) {
  total_records <- sum(family_count$count)
  family_count <- family_count %>%
    mutate(percentage = count / total_records * 100)
  return(family_count)
}

# Function to plot pie chart
plot_pie_chart <- function(family_percentage) {
  ggplot(family_percentage, aes(x = "", y = percentage, fill = family_id)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    labs(title = "Percentage of Each Family") +
    theme_void() +
    theme(legend.position = "bottom")
}

# Main function
main <- function(file_path) {
  gtf <- read_gtf(file_path)
  cat("Total number of records in GTF file:", nrow(gtf), "\n\n")
  
  cat("Summary by Family ID:\n")
  family_count <- count_by_family(gtf)
  print(family_count)
  
  cat("\nPercentage of Each Family:\n")
  family_percentage <- calculate_percentage(family_count)
  print(family_percentage)
  
  cat("\n")
  plot_pie_chart(family_percentage)
}

# Replace 'file_path' with the path to your GTF file
file_path <- "/media/kilian/OS/GTF_files_TEtranscript/GRCm38_GENCODE_rmsk_TE.gtf"
main(file_path)

