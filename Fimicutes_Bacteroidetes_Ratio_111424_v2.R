# First, analyze at phylum level
ps_phylum <- tax_glom(ps_rel, taxrank="Phylum")
# Convert to relative abundance
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x/sum(x) * 100)

# Extract abundance matrix and taxonomy
otu_mat <- as(otu_table(ps_phylum_rel), "matrix")
if(taxa_are_rows(ps_phylum_rel)) {
  otu_mat <- t(otu_mat)
}
tax_mat <- as(phyloseq::tax_table(ps_phylum_rel), "matrix")
metadata <- data.frame(sample_data(ps_phylum_rel))

# Get Firmicutes and Bacteroidetes abundances
firmicutes_idx <- which(tax_mat[, "Phylum"] == "Firmicutes")
bacteroidetes_idx <- which(tax_mat[, "Phylum"] == "Bacteroidetes")

# Calculate F/B ratio
fb_ratio <- data.frame(
  Sample = rownames(otu_mat),
  Firmicutes = otu_mat[, firmicutes_idx],
  Bacteroidetes = otu_mat[, bacteroidetes_idx],
  stringsAsFactors = FALSE
)
fb_ratio$FB_Ratio <- fb_ratio$Firmicutes / fb_ratio$Bacteroidetes

# Add metadata by matching sample names
fb_ratio$Disease.status <- metadata[match(fb_ratio$Sample, rownames(metadata)), "Disease.status"]
fb_ratio$CellType <- metadata[match(fb_ratio$Sample, rownames(metadata)), "CellType"]

# Print first few rows to check
print(head(fb_ratio))

# Statistical testing
cell_types <- unique(fb_ratio$CellType)
stats_results <- data.frame()

for(cell_type in cell_types) {
  subset_data <- fb_ratio[fb_ratio$CellType == cell_type,]
  test <- wilcox.test(FB_Ratio ~ Disease.status, data=subset_data)
  
  stats_results <- rbind(stats_results, data.frame(
    CellType = cell_type,
    p_value = test$p.value,
    Control_median = median(subset_data$FB_Ratio[subset_data$Disease.status == "Control"]),
    Celiac_median = median(subset_data$FB_Ratio[subset_data$Disease.status == "Celiac"])
  ))
}
stats_results$adj_p_value <- p.adjust(stats_results$p_value, method="BH")

# Visualiz