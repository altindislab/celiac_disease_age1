# Load required libraries
suppressPackageStartupMessages({
  library(phyloseq)
  library(MicrobiotaProcess)
  library(tidyverse)
  library(ggplot2)
  library(gridExtra)
})

# Set up color scheme
group_colors <- c(
  "Control" = "#00AED7",
  "Celiac" = "#FD9347"
)

# Create output directories
dir.create("lefse_results", showWarnings = FALSE)
dir.create("lefse_plots", showWarnings = FALSE)

# Process each cell type
cell_types <- c("Presort", "IGpos", "IGneg")

# Initialize an empty list to store results
all_cell_results <- list()

# Function to process differential analysis results
process_diff_analysis <- function(deres, cell_type) {
  if (!is.null(deres) && nrow(data.frame(deres)) > 0) {
    # Convert deres to data frame and add cell type
    diff_table <- data.frame(deres)
    diff_table$CellType <- cell_type
    
    # Create differential abundance boxplot
    tryCatch({
      p1 <- ggdiffbox(deres, 
                      box_notch = FALSE) +
        scale_fill_manual(values = group_colors) +
        scale_color_manual(values = group_colors) +
        theme_bw() +
        labs(title = paste(cell_type, "- Differential Features")) +
        theme(axis.text.y = element_text(size = 8))
      
      ggsave(paste0("lefse_plots/", tolower(cell_type), "_diffbox.pdf"),
             p1, width = 10, height = 8)
    }, error = function(e) {
      warning(sprintf("Error creating boxplot for %s: %s", cell_type, e$message))
    })
    
    # Create cladogram
    tryCatch({
      p2 <- ggdiffclade(deres, 
                        alpha = 0.3,
                        linewd = 0.15,
                        skpointsize = 0.6,
                        layout = "radial",
                        taxlevel = 3,
                        removeUnknown = TRUE,
                        reduce = TRUE) +
        scale_fill_manual(values = group_colors) +
        scale_color_manual(values = group_colors)
      
      ggsave(paste0("lefse_plots/", tolower(cell_type), "_cladogram.pdf"),
             p2, width = 12, height = 12)
    }, error = function(e) {
      warning(sprintf("Error creating cladogram for %s: %s", cell_type, e$message))
    })
    
    return(diff_table)
  }
  return(NULL)
}

# Main analysis loop
for (cell_type in cell_types) {
  cat(sprintf("\nProcessing differential analysis for %s...\n", cell_type))
  
  # Subset data for cell type
  ps_subset <- tryCatch({
    subset_samples(ps_rel, CellType == cell_type)
  }, error = function(e) {
    cat(sprintf("Error subsetting data for %s: %s\n", cell_type, e$message))
    return(NULL)
  })
  
  if (is.null(ps_subset)) next
  
  # Run differential analysis
  deres <- tryCatch({
    diff_analysis(ps_subset, 
                  classgroup = "Disease.status",
                  mlfun = "lda",
                  filtermod = "pvalue",
                  firstcomfun = "kruskal_test",
                  firstalpha = 0.05,
                  strictmod = FALSE,
                  secondcomfun = "wilcox_test",
                  subclmin = 3,
                  subclwilc = FALSE,
                  secondalpha = 0.05,
                  lda = 2)
  }, error = function(e) {
    cat(sprintf("Error in differential analysis for %s: %s\n", cell_type, e$message))
    return(NULL)
  })
  
  # Process results
  if (!is.null(deres)) {
    results <- process_diff_analysis(deres, cell_type)
    if (!is.null(results)) {
      all_cell_results[[cell_type]] <- results
      
      # Save individual results
      saveRDS(deres, file = paste0("lefse_results/", tolower(cell_type), "_lefse_results.rds"))
      write.csv(results,
                paste0("lefse_results/", tolower(cell_type), "_results.csv"),
                row.names = FALSE)
    }
  } else {
    cat(sprintf("No significant results found for %s\n", cell_type))
  }
}

# Process combined results
if (length(all_cell_results) > 0) {
  # Combine all results into one dataframe
  all_results <- bind_rows(all_cell_results)
  
  # Save combined summary
  write.csv(all_results,
            "lefse_results/combined_lefse_summary.csv",
            row.names = FALSE)
  
  # Print available columns
  cat("\nAvailable columns in the dataset:\n")
  print(names(all_results))
  
  # Create summary statistics
  summary_table <- all_results %>%
    group_by(CellType) %>%
    summarise(
      Total_Features = n(),
      Features_in_Control = sum(Disease.status == "Control", na.rm = TRUE),
      Features_in_Celiac = sum(Disease.status == "Celiac", na.rm = TRUE),
      Mean_LDA_Score = mean(LDAmean, na.rm = TRUE),
      Mean_FDR = mean(fdr, na.rm = TRUE)
    )
  
  # Save summary statistics
  write.csv(summary_table,
            "lefse_results/summary_statistics.csv",
            row.names = FALSE)
  
  # Print top features for each cell type
  cat("\nTop features by LDA score for each cell type:\n")
  for (ct in unique(all_results$CellType)) {
    cat(sprintf("\n%s:\n", ct))
    
    tryCatch({
      # Filter for this cell type
      ct_data <- all_results %>%
        filter(CellType == ct) %>%
        arrange(desc(abs(LDAmean)))
      
      # Select columns that exist
      cols_to_select <- c("f", "Disease.status", "LDAmean", "pvalue", "fdr")
      cols_present <- intersect(cols_to_select, names(ct_data))
      
      # Select top 10 features
      top_features <- ct_data %>%
        select(!!!syms(cols_present)) %>%
        head(10)
      
      if(nrow(top_features) > 0) {
        # Format the output for better readability
        formatted_features <- top_features %>%
          mutate(
            LDAmean = round(LDAmean, 3),
            pvalue = format.pval(pvalue, digits = 3),
            fdr = format.pval(fdr, digits = 3)
          )
        
        print(formatted_features)
      } else {
        cat("No features found.\n")
      }
    }, error = function(e) {
      cat(sprintf("Error processing top features for %s: %s\n", ct, e$message))
    })
  }
  
  # Create visualization plots
  tryCatch({
    # Summary plot
    summary_plot <- ggplot(all_results, 
                           aes(x = reorder(f, LDAmean),
                               y = LDAmean,
                               fill = Disease.status)) +
      geom_bar(stat = "identity") +
      facet_wrap(~CellType, scales = "free_y") +
      coord_flip() +
      theme_bw() +
      theme(axis.text.y = element_text(size = 8)) +
      labs(x = "Features",
           y = "LDA Score",
           title = "Differential Analysis Summary",
           fill = "Disease Status") +
      scale_fill_manual(values = group_colors)
    
    ggsave("lefse_plots/combined_summary.pdf", 
           summary_plot, width = 12, height = 8)
    
    # Volcano plots for each cell type
    for (ct in unique(all_results$CellType)) {
      p_volcano <- all_results %>%
        filter(CellType == ct) %>%
        ggplot(aes(x = LDAmean, y = -log10(fdr))) +
        geom_point(aes(color = Disease.status), alpha = 0.6) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        scale_color_manual(values = group_colors) +
        theme_bw() +
        labs(title = paste(ct, "- Volcano Plot"),
             x = "LDA Score",
             y = "-log10(FDR)")
      
      ggsave(paste0("lefse_plots/", tolower(ct), "_volcano.pdf"),
             p_volcano, width = 8, height = 6)
    }
    
    # FDR distribution plot
    p_fdr <- ggplot(all_results, aes(x = fdr, fill = Disease.status)) +
      geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
      facet_wrap(~CellType) +
      scale_fill_manual(values = group_colors) +
      theme_bw() +
      labs(title = "FDR Distribution",
           x = "FDR",
           y = "Count")
    
    ggsave("lefse_plots/fdr_distribution.pdf", 
           p_fdr, width = 10, height = 6)
    
  }, error = function(e) {
    cat(sprintf("Error creating visualization plots: %s\n", e$message))
  })
  
  # Print statistical insights
  cat("\nStatistical Insights:\n")
  for (ct in unique(all_results$CellType)) {
    cat(sprintf("\n%s:\n", ct))
    
    # Get subset for this cell type
    ct_data <- all_results %>%
      filter(CellType == ct)
    
    # Number of significant features (with more detailed info)
    sig_features <- ct_data %>%
      filter(fdr < 0.05) %>%
      nrow()
    
    # Mean LDA score
    mean_lda <- mean(abs(ct_data$LDAmean), na.rm = TRUE)
    
    # Distribution between groups for significant features
    sig_dist <- ct_data %>%
      filter(fdr < 0.05) %>%
      count(Disease.status) %>%
      spread(Disease.status, n, fill = 0)
    
    cat(sprintf("- Significant features (FDR < 0.05): %d\n", sig_features))
    cat(sprintf("- Mean absolute LDA score: %.2f\n", mean_lda))
    cat("- Distribution of significant features:\n")
    if(nrow(sig_dist) > 0) {
      print(sig_dist)
    } else {
      cat("  No significant features found\n")
    }
    
    # Add additional statistics
    cat("\nTop 5 features by LDA score:\n")
    ct_data %>%
      arrange(desc(abs(LDAmean))) %>%
      select(f, LDAmean, pvalue, fdr) %>%
      head(5) %>%
      mutate(across(where(is.numeric), round, 3)) %>%
      print()
  }
} else {
  cat("\nNo results were generated for any cell type.\n")
}

cat("\nDifferential analysis complete. Results are saved in 'lefse_results' and 'lefse_plots' directories.\n")
