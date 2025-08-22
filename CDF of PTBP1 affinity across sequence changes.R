# Clear environment and load required libraries
rm(list = ls())
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(reshape2)
library(dplyr)

# Read PTBP1 SELEX data (col1: sequence, col2: read count)
selex_data <- read.table("./mean_combined_count.txt", header = FALSE)
colnames(selex_data) <- c("raw_sequence", "read_count")

# Data preprocessing: extract sequence features
selex_data <- selex_data %>%
  mutate(
    # Extract and combine first three and last three characters
    sequence = paste(substr(raw_sequence, 1, 3), substr(raw_sequence, 5, 7), sep = ""),
    # Extract the middle character as label
    label = substr(raw_sequence, 4, 4)
  )

# Sort data by sequence and get unique sequence names
sorted_data <- selex_data[order(selex_data$sequence), ]
unique_sequences <- unique(sorted_data$sequence)

# Create sub_file directory if it doesn't exist
if (!dir.exists("./sub_file/")) {
  dir.create("./sub_file/")
}

# Split data by unique sequences and save as separate files
for (seq_name in unique_sequences) {
  subset_data <- sorted_data[sorted_data$sequence == seq_name, ]
  file_path <- paste0("./sub_file/", seq_name, ".txt")
  write.table(
    subset_data, 
    file = file_path, 
    sep = "\t", 
    row.names = FALSE, 
    col.names = TRUE,
    quote = FALSE
  )
}

# Process each sub-file to calculate fold changes
# Create output directory if it doesn't exist
if (!dir.exists("./sub_file_FC/")) {
  dir.create("./sub_file_FC/")
}

# Calculate fold changes for each substitution type
for (filename in dir("./sub_file/")) {
  # Read data from current file
  data <- read.table(paste0("./sub_file/", filename), header = TRUE)
  
  # Extract counts for each nucleotide
  G_count <- data[data$label == "G", "read_count"]
  A_count <- data[data$label == "A", "read_count"]
  C_count <- data[data$label == "C", "read_count"]
  T_count <- data[data$label == "T", "read_count"]
  
  # Calculate fold changes for each possible substitution
  # A as original nucleotide
  A_to_G_FC <- G_count / A_count
  A_to_C_FC <- C_count / A_count
  A_to_T_FC <- T_count / A_count
  
  # G as original nucleotide
  G_to_A_FC <- A_count / G_count
  G_to_C_FC <- C_count / G_count
  G_to_T_FC <- T_count / G_count
  
  # C as original nucleotide
  C_to_A_FC <- A_count / C_count
  C_to_G_FC <- G_count / C_count
  C_to_T_FC <- T_count / C_count
  
  # T as original nucleotide
  T_to_A_FC <- A_count / T_count
  T_to_C_FC <- C_count / T_count
  T_to_G_FC <- G_count / T_count
  
  # Extract sequence identifier from filename
  seq_identifier <- substr(filename, 1, 6)
  
  # Create summary data frame
  fc_summary <- data.frame(
    sequence = seq_identifier,
    A_to_G_FC, A_to_C_FC, A_to_T_FC,
    G_to_A_FC, G_to_C_FC, G_to_T_FC,
    C_to_A_FC, C_to_G_FC, C_to_T_FC,
    T_to_A_FC, T_to_C_FC, T_to_G_FC
  )
  
  # Save fold change results
  output_file <- paste0("./sub_file_FC/", seq_identifier, "_FC.csv")
  write.csv(fc_summary, output_file, row.names = FALSE)
}

# The "FC_result.csv" file was obtained by merging all files under "./sub_file_FC/"
fc_combined <- read.csv("./FC_result.csv", header = TRUE)

# Prepare data for plotting
fc_long <- melt(fc_combined)
fc_type_data <- fc_long[, 2:3]
colnames(fc_type_data) <- c("type", "FC")
fc_type_data$FC <- log2(fc_type_data$FC)  # Convert to log2 scale

# Create background data for comparison
background_data <- fc_type_data %>%
  mutate(type = "Background")

# Combine data for plotting
plot_data <- bind_rows(fc_type_data, background_data)

# Generate CDF plots
ggplot(plot_data, aes(x = FC)) +
  stat_ecdf(aes(colour = type), linewidth = 0.7) +   
  scale_y_continuous() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~ type, ncol = 4) 

# Prepare groups for statistical comparison (exclude Background)
all_groups <- unique(plot_data$type)
all_groups <- all_groups[all_groups != "Background"]

# Filter data for statistical tests
filtered_data <- plot_data %>%
  filter(type %in% c(all_groups, "Background")) %>%
  droplevels()

# Function to calculate Cohen's d effect size
cohens_d <- function(x, y) {
  mean_diff <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
  sd_pooled <- sqrt(
    ((length(x) - 1) * var(x, na.rm = TRUE) + (length(y) - 1) * var(y, na.rm = TRUE)) / 
      (length(x) + length(y) - 2)
  )
  return(mean_diff / sd_pooled)
}

# Perform statistical comparisons with strict settings
test_results <- lapply(all_groups, function(group) {
  # Extract data for current group and background
  group_data <- filtered_data %>% filter(type == group) %>% pull(FC)
  background_data <- filtered_data %>% filter(type == "Background") %>% pull(FC)
  
  # Test for variance equality
  var_test <- var.test(group_data, background_data)
  var_equal <- var_test$p.value > 0.01  # Stringent variance equality criterion
  
  # Perform Welch's t-test (more conservative)
  test_result <- t.test(
    group_data, 
    background_data, 
    var.equal = FALSE,  # Force Welch correction
    conf.level = 0.99   # Higher confidence level
  )
  
  # Calculate effect size
  effect_size <- cohens_d(group_data, background_data)
  
  # Return formatted results
  return(data.frame(
    group1 = group,
    group2 = "Background",
    statistic = as.numeric(test_result$statistic),
    p_value = test_result$p.value,
    effect_size = effect_size,
    var_equal = var_equal,
    stringsAsFactors = FALSE
  ))
})

# Combine results into data frame
strict_comparisons <- do.call(rbind, test_results)

# Apply multiple testing correction and significance thresholds
strict_comparisons <- strict_comparisons %>%
  mutate(
    # Apply Bonferroni correction
    p_adj = p.adjust(p_value, method = "bonferroni"),
    
    # Define strict significance levels
    significance = case_when(
      p_adj < 1e-6 & abs(effect_size) > 1.5 ~ "******",
      p_adj < 1e-5 & abs(effect_size) > 1.2 ~ "*****",
      p_adj < 1e-4 & abs(effect_size) > 1.0 ~ "****",
      p_adj < 1e-3 & abs(effect_size) > 0.8 ~ "***",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)  # Sort by adjusted p-value

# Display comparison results
cat("Comparison results between all groups and Background (strict criteria):\n")
print(strict_comparisons %>% 
        select(group1, p_adj, effect_size, significance))

# Handle zero p-values and calculate absolute effect sizes
min_non_zero_p <- min(strict_comparisons$p_adj[strict_comparisons$p_adj > 0])
strict_comparisons <- strict_comparisons %>%
  mutate(
    p_adj_corrected = ifelse(p_adj == 0, min_non_zero_p / 10, p_adj),
    neg_log10_padj = -log10(p_adj_corrected),
    abs_effect_size = abs(effect_size)  
  )

# Generate bubble plot visualization
ggplot(strict_comparisons, aes(x = group1, y = significance)) +
  geom_point(aes(color = abs_effect_size, size = neg_log10_padj), alpha = 0.7) +
  scale_size_continuous(
    range = c(2, 8),
    name = "-log10(Adjusted p-value)",
    breaks = c(10, 50, 100, max(strict_comparisons$neg_log10_padj)),
    labels = c("10", "50", "100", expression(">="~100))
  ) +
  scale_color_gradient2(
    low = "#5d6bb2", 
    mid = "white", 
    high = "red", 
    midpoint = mean(strict_comparisons$abs_effect_size),
    name = "Absolute Effect Size"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Group",
    y = "Significance",
    title = "Absolute Effect Size (Color) vs. Adjusted p-value (Size)"
  ) +
  annotate("text", x = 1, y = "ns", 
           label = "Note: p=0 shown as largest bubbles", 
           hjust = 0, vjust = 0, size = 3, color = "gray50")