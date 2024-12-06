# Load required libraries
library(ggplot2)
library(phyloseq)
library(dplyr)
library(vegan)
library(patchwork)

# Function to create alpha diversity plot for a specific sort type
create_alpha_plot <- function(ps_obj, sort_type) {
  # Filter for specific sort type
  ps_subset <- subset_samples(ps_obj, Sort %in% sort_type)
  
  # Calculate alpha diversity metrics
  alpha_div <- estimate_richness(ps_subset, measures = c("Observed", "Shannon", "Simpson"))
  
  # Add metadata information
  alpha_div$Disease.status <- sample_data(ps_subset)$Disease.status
  alpha_div$Sort <- sample_data(ps_subset)$Sort
  alpha_div$Disease.status <- factor(alpha_div$Disease.status, levels = c("Control", "Celiac"))
  
  # Create the box plot
  p <- ggplot(alpha_div, aes(x = Disease.status, y = Observed, fill = Disease.status)) +
    geom_boxplot(outlier.shape = 19, outlier.size = 1.5) +
    geom_jitter(width = 0.2, alpha = 0.5) +  # Add individual points
    scale_fill_manual(values = c("Control" = "#00AED7", "Celiac" = "#FD9347")) +
    labs(
      title = paste("Observed ASV -", sort_type),
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
  
  # Perform statistical tests for each diversity measure
  measures <- c("Observed", "Shannon", "Simpson")
  stats_results <- data.frame()
  
  for(measure in measures) {
    # Wilcoxon test
    test <- wilcox.test(alpha_div[[measure]] ~ alpha_div$Disease.status, exact = FALSE)
    
    # Calculate effect size (Cliff's Delta)
    control_vals <- alpha_div[[measure]][alpha_div$Disease.status == "Control"]
    celiac_vals <- alpha_div[[measure]][alpha_div$Disease.status == "Celiac"]
    
    n1 <- length(control_vals)
    n2 <- length(celiac_vals)
    dominance <- sum(outer(control_vals, celiac_vals, ">")) - 
      sum(outer(control_vals, celiac_vals, "<"))
    cliff_delta <- dominance / (n1 * n2)
    
    stats_results <- rbind(stats_results, data.frame(
      sort_type = sort_type,
      measure = measure,
      statistic = test$statistic,
      p_value = test$p.value,
      cliff_delta = cliff_delta,
      control_mean = mean(control_vals),
      control_sd = sd(control_vals),
      celiac_mean = mean(celiac_vals),
      celiac_sd = sd(celiac_vals)
    ))
  }
  
  # Save statistical test results
  write.csv(stats_results, 
            paste0("alpha_diversity_statistics_", sort_type, ".csv"), 
            row.names = FALSE)
  
  # Add p-value annotation if significant
  if (stats_results$p_value[1] < 0.05) {  # For Observed diversity
    sig_symbol <- if(stats_results$p_value[1] < 0.001) "***"
    else if(stats_results$p_value[1] < 0.01) "**"
    else "*"
    
    y_pos <- max(alpha_div$Observed) * 1.1
    p <- p + 
      geom_segment(aes(x = 1, xend = 2, y = y_pos, yend = y_pos)) +
      geom_text(aes(x = 1.5, y = y_pos * 1.05), 
                label = sig_symbol)
  }
  
  return(list(plot = p, stats = stats_results))
}

# Function to create Bray-Curtis analysis for a specific sort type
create_bray_curtis_analysis <- function(ps_obj, sort_type) {
  # Filter for specific sort type
  ps_subset <- subset_samples(ps_obj, Sort %in% sort_type)
  
  # Calculate Bray-Curtis distances
  bc_dist <- phyloseq::distance(ps_subset, method = "bray")
  
  # Get metadata
  meta <- data.frame(sample_data(ps_subset))
  
  # Perform PERMANOVA
  permanova_result <- adonis2(bc_dist ~ Disease.status, data = meta, permutations = 999)
  
  # Perform betadisper test for homogeneity of dispersion
  beta_disper <- betadisper(bc_dist, meta$Disease.status)
  disper_test <- permutest(beta_disper, permutations = 999)
  
  # Extract within-group and between-group distances
  dist_matrix <- as.matrix(bc_dist)
  
  # Calculate within-group distances
  control_samples <- meta$Disease.status == "Control"
  celiac_samples <- meta$Disease.status == "Celiac"
  
  within_control <- dist_matrix[control_samples, control_samples][upper.tri(dist_matrix[control_samples, control_samples])]
  within_celiac <- dist_matrix[celiac_samples, celiac_samples][upper.tri(dist_matrix[celiac_samples, celiac_samples])]
  between_groups <- dist_matrix[control_samples, celiac_samples]
  
  # Prepare data for plotting
  plot_data <- data.frame(
    Distance = c(within_control, within_celiac),
    Group = factor(c(rep("Control", length(within_control)),
                     rep("Celiac", length(within_celiac))),
                   levels = c("Control", "Celiac"))
  )
  
  # Create box plot
  p <- ggplot(plot_data, aes(x = Group, y = Distance, fill = Group)) +
    geom_boxplot(outlier.shape = 19, outlier.size = 1.5) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values = c("Control" = "#00AED7", "Celiac" = "#FD9347")) +
    labs(
      title = paste("Bray-Curtis Dissimilarity -", sort_type),
      subtitle = paste("PERMANOVA R² =", round(permanova_result$R2[1], 3),
                       ", p =", format(permanova_result$`Pr(>F)`[1], digits = 3)),
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
  
  # Save statistical results
  stats_results <- data.frame(
    sort_type = sort_type,
    test = c("PERMANOVA", "Dispersion"),
    statistic = c(permanova_result$F[1], disper_test$F[1]),
    R2 = c(permanova_result$R2[1], NA),
    p_value = c(permanova_result$`Pr(>F)`[1], disper_test$tab$`Pr(>F)`[1])
  )
  
  write.csv(stats_results,
            paste0("beta_diversity_statistics_", sort_type, ".csv"),
            row.names = FALSE)
  
  return(list(plot = p, stats = stats_results))
}

# Process each sort type
sort_types <- c("IGneg", "IGpos")
all_results <- list()

for (sort_type in sort_types) {
  # Run analyses
  alpha_results <- create_alpha_plot(ps_filtered_prev, sort_type)
  beta_results <- create_bray_curtis_analysis(ps_filtered_prev, sort_type)
  
  # Store results
  all_results[[sort_type]] <- list(
    alpha = alpha_results,
    beta = beta_results
  )
  
  # Create combined plot for this sort type
  combined_plot <- alpha_results$plot / beta_results$plot +
    plot_layout(heights = c(1, 1))
  
  # Save plots
  ggsave(
    paste0("diversity_analysis_", sort_type, ".pdf"),
    combined_plot,
    width = 6,
    height = 8,
    dpi = 300
  )
}

# Create final combined figure
combined_figure <- (all_results$IGneg$alpha$plot + all_results$IGpos$alpha$plot) /
  (all_results$IGneg$beta$plot + all_results$IGpos$beta$plot) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Diversity Analyses by Sort Type",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

# Save combined figure
ggsave(
  "all_diversity_analyses.pdf",
  combined_figure,
  width = 12,
  height = 10,
  dpi = 300
)

# Create comprehensive statistical summary
alpha_stats <- do.call(rbind, lapply(all_results, function(x) x$alpha$stats))
beta_stats <- do.call(rbind, lapply(all_results, function(x) x$beta$stats))

write.csv(alpha_stats, "alpha_diversity_statistics_summary.csv", row.names = FALSE)
write.csv(beta_stats, "beta_diversity_statistics_summary.csv", row.names = FALSE)
