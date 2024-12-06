# Load required libraries
library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)

# Define color schemes
control_color <- "#00AED7"  # Light blue
celiac_color <- "#FD9347"   # Orange

# Function to perform ordination analysis for a specific comparison
perform_ordination_analysis <- function(ps, group1, group2, type = "PCA", output_prefix) {
  # Filter samples for the specific comparison
  subset_samples <- sample_data(ps)$Group %in% c(group1, group2)
  ps_subset <- prune_samples(subset_samples, ps)
  
  # Create a simplified group designation
  sample_data(ps_subset)$SimpleGroup <- factor(
    ifelse(grepl("Control", sample_data(ps_subset)$Group), "Control", "Celiac"),
    levels = c("Control", "Celiac")
  )
  
  # Transform counts to relative abundance
  ps_rel <- transform_sample_counts(ps_subset, function(x) x/sum(x))
  
  if(type == "PCA") {
    # Perform PCA
    ord <- ordinate(ps_rel, method = "RDA", distance = "euclidean")
    
    # Create PCA plot
    p <- plot_ordination(ps_rel, ord, color = "SimpleGroup") +
      geom_point(size = 3) +
      stat_ellipse(type = "norm") +
      scale_color_manual(values = c(Control = control_color, Celiac = celiac_color)) +
      ggtitle(paste("PCA -", sub("_.*", "", group1), "Comparison")) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )
    
  } else {  # PCoA
    # Calculate Bray-Curtis distances and perform PCoA
    ord <- ordinate(ps_rel, method = "PCoA", distance = "bray")
    
    # Create PCoA plot
    p <- plot_ordination(ps_rel, ord, color = "SimpleGroup") +
      geom_point(size = 3) +
      stat_ellipse(type = "norm") +
      scale_color_manual(values = c(Control = control_color, Celiac = celiac_color)) +
      ggtitle(paste("PCoA -", sub("_.*", "", group1), "Comparison")) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )
  }
  
  # Save plot
  filename <- paste0(output_prefix, "_", sub("_.*", "", group1), "_comparison.pdf")
  ggsave(filename, p, width = 8, height = 6)
  
  # Perform PERMANOVA
  otu_mat <- t(otu_table(ps_rel)@.Data)
  metadata <- data.frame(sample_data(ps_rel))
  dist_mat <- vegdist(otu_mat, method = "bray")
  
  perm_result <- adonis2(dist_mat ~ SimpleGroup, data = metadata)
  
  # Save PERMANOVA results
  write.csv(perm_result, 
            paste0(output_prefix, "_PERMANOVA_", sub("_.*", "", group1), "_comparison.csv"))
  
  return(list(plot = p, permanova = perm_result))
}

# Perform analyses for each comparison
# 1. Presort comparison
presort_groups <- c("Presort_Control_1", "Presort_Celiac_1")
pca_presort <- perform_ordination_analysis(ps, presort_groups[1], presort_groups[2], 
                                           "PCA", "ordination_pca")
pcoa_presort <- perform_ordination_analysis(ps, presort_groups[1], presort_groups[2], 
                                            "PCoA", "ordination_pcoa")

# 2. IGpos comparison
igpos_groups <- c("IGpos_Control_1", "IGpos_Celiac_1")
pca_igpos <- perform_ordination_analysis(ps, igpos_groups[1], igpos_groups[2], 
                                         "PCA", "ordination_pca")
pcoa_igpos <- perform_ordination_analysis(ps, igpos_groups[1], igpos_groups[2], 
                                          "PCoA", "ordination_pcoa")

# 3. IGneg comparison
igneg_groups <- c("IGneg_Control_1", "IGneg_Celiac_1")
pca_igneg <- perform_ordination_analysis(ps, igneg_groups[1], igneg_groups[2], 
                                         "PCA", "ordination_pca")
pcoa_igneg <- perform_ordination_analysis(ps, igneg_groups[1], igneg_groups[2], 
                                          "PCoA", "ordination_pcoa")

# Create combined plots for publication
library(gridExtra)

# Combine PCA plots
combined_pca <- grid.arrange(
  pca_presort$plot + ggtitle("Presort"),
  pca_igpos$plot + ggtitle("IG+"),
  pca_igneg$plot + ggtitle("IG-"),
  ncol = 3
)
ggsave("combined_pca_plots.pdf", combined_pca, width = 18, height = 6)

# Combine PCoA plots
combined_pcoa <- grid.arrange(
  pcoa_presort$plot + ggtitle("Presort"),
  pcoa_igpos$plot + ggtitle("IG+"),
  pcoa_igneg$plot + ggtitle("IG-"),
  ncol = 3
)
ggsave("combined_pcoa_plots.pdf", combined_pcoa, width = 18, height = 6)

# Save all PERMANOVA results in a single file
all_permanova <- list(
  Presort = pca_presort$permanova,
  IGpos = pca_igpos$permanova,
  IGneg = pca_igneg$permanova
)

write.csv(do.call(rbind, lapply(names(all_permanova), function(name) {
  data.frame(
    Comparison = name,
    all_permanova[[name]]
  )
})), "all_permanova_results.csv")

