# Load required libraries
library(ggplot2)
library(phyloseq)
library(dplyr)
library(vegan)
library(patchwork)

# Function to create alpha diversity plot
create_alpha_plot <- function(ps_obj) {
  # Calculate alpha diversity metrics
  alpha_div <- estimate_richness(ps_obj, measures = "Observed")
  
  # Add metadata information
  alpha_div$Disease.status <- sample_data(ps_obj)$Disease.status
  alpha_div$Disease.status <- factor(alpha_div$Disease.status, levels = c("Control", "Celiac"))
  
  # Create the box plot
  p <- ggplot(alpha_div, aes(x = Disease.status, y = Observed, fill = Disease.status)) +
    geom_boxplot(outlier.shape = 19, outlier.size = 1.5) +
    scale_fill_manual(values = c("Control" = "#00AED7", "Celiac" = "#FD9347")) +
    labs(
      title = "Observed ASV",
      y = "Alpha Diversity\n(ASV)",
      x = ""
    ) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
  
  # Perform statistical test
  stat_test <- wilcox.test(Observed ~ Disease.status, data = alpha_div, exact = FALSE)
  
  # Add p-value annotation if significant
  if (stat_test$p.value < 0.05) {
    sig_symbol <- if(stat_test$p.value < 0.001) "***"
    else if(stat_test$p.value < 0.01) "**"
    else "*"
    
    y_pos <- max(alpha_div$Observed) * 1.1
    p <- p + 
      geom_segment(aes(x = 1, xend = 2, y = y_pos, yend = y_pos)) +
      geom_text(aes(x = 1.5, y = y_pos * 1.05), label = sig_symbol)
  }
  
  return(list(plot = p, statistics = data.frame(
    test = "Wilcoxon",
    statistic = stat_test$statistic,
    p_value = stat_test$p.value
  )))
}

# Function to create Bray-Curtis analysis with PERMANOVA
create_beta_analysis <- function(ps_obj) {
  # Calculate Bray-Curtis distances
  bc_dist <- phyloseq::distance(ps_obj, method = "bray")
  
  # Get metadata
  meta <- data.frame(sample_data(ps_obj))
  
  # Perform PERMANOVA
  permanova_result <- adonis2(bc_dist ~ Disease.status, data = meta, permutations = 999)
  
  # Create distance matrix
  dist_matrix <- as.matrix(bc_dist)
  
  # Function to extract between-group distances
  get_between_group_distances <- function(dist_matrix, meta) {
    distances <- data.frame(
      Distance = numeric(),
      Disease.status = character()
    )
    
    # Get distances within Control group
    control_samples <- rownames(meta)[meta$Disease.status == "Control"]
    if(length(control_samples) >= 2) {
      control_dist <- dist_matrix[control_samples, control_samples]
      control_dist <- control_dist[upper.tri(control_dist)]
      distances <- rbind(distances, 
                         data.frame(Distance = control_dist,
                                    Disease.status = "Control"))
    }
    
    # Get distances within Celiac group
    celiac_samples <- rownames(meta)[meta$Disease.status == "Celiac"]
    if(length(celiac_samples) >= 2) {
      celiac_dist <- dist_matrix[celiac_samples, celiac_samples]
      celiac_dist <- celiac_dist[upper.tri(celiac_dist)]
      distances <- rbind(distances, 
                         data.frame(Distance = celiac_dist,
                                    Disease.status = "Celiac"))
    }
    
    distances$Disease.status <- factor(distances$Disease.status, 
                                       levels = c("Control", "Celiac"))
    return(distances)
  }
  
  # Get distances
  bc_distances <- get_between_group_distances(dist_matrix, meta)
  
  # Create box plot
  p <- ggplot(bc_distances, aes(x = Disease.status, y = Distance, fill = Disease.status)) +
    geom_boxplot(outlier.shape = 19, outlier.size = 1.5) +
    scale_fill_manual(values = c("Control" = "#00AED7", "Celiac" = "#FD9347")) +
    labs(
      title = "Bray-Curtis Dissimilarity",
      subtitle = sprintf("PERMANOVA R² = %.3f, p = %.3f", 
                         permanova_result$R2[1], 
                         permanova_result$`Pr(>F)`[1]),
      y = "Distance",
      x = ""
    ) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 10, hjust = 0),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
  
  # Perform within-group statistical test
  stat_test <- wilcox.test(Distance ~ Disease.status, data = bc_distances, exact = FALSE)
  
  # Add p-value annotation if significant
  if (stat_test$p.value < 0.05) {
    sig_symbol <- if(stat_test$p.value < 0.001) "***"
    else if(stat_test$p.value < 0.01) "**"
    else "*"
    
    y_pos <- max(bc_distances$Distance) * 1.1
    p <- p + 
      geom_segment(aes(x = 1, xend = 2, y = y_pos, yend = y_pos)) +
      geom_text(aes(x = 1.5, y = y_pos * 1.05), label = sig_symbol)
  }
  
  # Test for homogeneity of dispersion
  disp <- betadisper(bc_dist, meta$Disease.status)
  disp_test <- permutest(disp, permutations = 999)
  
  return(list(
    plot = p,
    statistics = list(
      permanova = permanova_result,
      within_group_test = data.frame(
        test = "Wilcoxon",
        statistic = stat_test$statistic,
        p_value = stat_test$p.value
      ),
      dispersion_test = disp_test
    )
  ))
}

# Create both analyses
alpha_results <- create_alpha_plot(ps_filtered_prev)
beta_results <- create_beta_analysis(ps_filtered_prev)

# Combine plots vertically
combined_plot <- alpha_results$plot / beta_results$plot +
  plot_layout(heights = c(1, 1))

# Save the combined plot
ggsave("alpha_beta_diversity_boxplots.pdf", combined_plot, width = 6, height = 8, dpi = 300)

# Save individual plots
ggsave("alpha_diversity_boxplot.pdf", alpha_results$plot, width = 6, height = 4, dpi = 300)
ggsave("beta_diversity_boxplot.pdf", beta_results$plot, width = 6, height = 4, dpi = 300)

# Save statistical results
statistical_summary <- list(
  alpha_diversity = alpha_results$statistics,
  beta_diversity = beta_results$statistics
)

# Write statistical results to a file
sink("diversity_analysis_statistics.txt")
cat("Alpha Diversity Analysis:\n")
print(alpha_results$statistics)
cat("\nBeta Diversity Analysis:\n")
cat("\nPERMANOVA Results:\n")
print(beta_results$statistics$permanova)
cat("\nWithin-group Distance Test:\n")
print(beta_results$statistics$within_group_test)
cat("\nDispersion Homogeneity Test:\n")
print(beta_results$statistics$dispersion_test)
sink()
