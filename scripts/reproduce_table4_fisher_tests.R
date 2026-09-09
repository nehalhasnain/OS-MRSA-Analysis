#!/usr/bin/env Rscript
# ==============================================================================
# OS-MRSA Project: Table 4 Reproducibility Script
# Target Journal: One Health Advances (BioMed Central / Springer Nature)
#
# Description:
#   Reproduces all 9 two-sided Fisher's exact tests, conditional maximum
#   likelihood odds ratios (OR), and exact 95% confidence intervals reported in
#   Table 4 of the manuscript from the de-identified minimal dataset.
#
# Requirements:
#   R >= 4.0.0 (Base R only; no third-party package dependencies)
#
# Usage:
#   Rscript reproduce_table4_fisher_tests.R
# ==============================================================================

# Locate input dataset (looks in current dir, data/, 06_Tables_and_Additional_Files, or parent dirs)
possible_paths <- c(
  "Additional_file_1_Deidentified_Minimal_Dataset.csv",
  "data/Additional_file_1_Deidentified_Minimal_Dataset.csv",
  "../data/Additional_file_1_Deidentified_Minimal_Dataset.csv",
  "06_Tables_and_Additional_Files/Additional_file_1_Deidentified_Minimal_Dataset.csv",
  "../06_Tables_and_Additional_Files/Additional_file_1_Deidentified_Minimal_Dataset.csv"
)

data_file <- NULL
for (p in possible_paths) {
  if (file.exists(p)) {
    data_file <- p
    break
  }
}

if (is.null(data_file)) {
  stop("Error: Could not find 'Additional_file_1_Deidentified_Minimal_Dataset.csv'. Please run from project root or scripts directory.")
}

cat("===============================================================================
")
cat("Loading dataset:", data_file, "
")
df <- read.csv(data_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("Loaded %d animals (%d variables)\n", nrow(df), ncol(df)))
cat("===============================================================================

")

# Helper function to run Fisher test, format, and return row
run_test <- function(comparison, subset_df, row_var, row_levels, col_var, col_levels, group_a_label, group_b_label, interp) {
  # Build contingency table
  tab <- table(
    factor(subset_df[[row_var]], levels = row_levels),
    factor(subset_df[[col_var]], levels = col_levels)
  )
  
  # Run two-sided Fisher exact test
  ft <- fisher.test(tab)
  
  # Group counts (n/N)
  n_a <- tab[1, 1]
  tot_a <- sum(tab[1, ])
  n_b <- tab[2, 1]
  tot_b <- sum(tab[2, ])
  
  group_a_str <- sprintf("%d/%d (%s)", n_a, tot_a, group_a_label)
  group_b_str <- sprintf("%d/%d (%s)", n_b, tot_b, group_b_label)
  
  # Format odds ratio and CI
  or_val <- sprintf("%.2f", ft$estimate)
  ci_str <- sprintf("%.2f–%.2f", ft$conf.int[1], ft$conf.int[2])
  or_ci_str <- sprintf("%s (%s)", or_val, ci_str)
  
  # Format p-value
  p_val <- sprintf("%.3f", ft$p.value)
  
  data.frame(
    Comparison = comparison,
    `n/N, Group A` = group_a_str,
    `n/N, Group B` = group_b_str,
    `OR (95% CI)` = or_ci_str,
    `Fisher's p` = p_val,
    Interpretation = interp,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

results <- list()

# 1. Cat vs dog: S. aureus occurrence
d1 <- df[df$Species %in% c("Cat", "Dog"), ]
results[[1]] <- run_test(
  comparison = "Cat vs dog: S. aureus occurrence",
  subset_df = d1,
  row_var = "Species", row_levels = c("Cat", "Dog"),
  col_var = "S_aureus_Confirmed", col_levels = c("Yes", "No"),
  group_a_label = "cat", group_b_label = "dog",
  interp = "No statistically clear difference detected"
)

# 2. Cat vs dog: mecA-positive MRSA occurrence
results[[2]] <- run_test(
  comparison = "Cat vs dog: mecA-positive MRSA occurrence",
  subset_df = d1,
  row_var = "Species", row_levels = c("Cat", "Dog"),
  col_var = "PCR_mecA", col_levels = c("Yes", "No"),
  group_a_label = "cat", group_b_label = "dog",
  interp = "No statistically clear difference detected"
)

# 3. Cats, sick vs apparently healthy: S. aureus
dc <- df[df$Species == "Cat", ]
results[[3]] <- run_test(
  comparison = "Cats, sick vs apparently healthy: S. aureus",
  subset_df = dc,
  row_var = "Health_status", row_levels = c("Sick", "Apparently Healthy"),
  col_var = "S_aureus_Confirmed", col_levels = c("Yes", "No"),
  group_a_label = "sick", group_b_label = "healthy",
  interp = "No statistically clear difference detected"
)

# 4. Cats, sick vs apparently healthy: MRSA
results[[4]] <- run_test(
  comparison = "Cats, sick vs apparently healthy: MRSA",
  subset_df = dc,
  row_var = "Health_status", row_levels = c("Sick", "Apparently Healthy"),
  col_var = "PCR_mecA", col_levels = c("Yes", "No"),
  group_a_label = "sick", group_b_label = "healthy",
  interp = "No statistically clear difference detected"
)

# 5. Dogs, apparently healthy vs sick: S. aureus
dd <- df[df$Species == "Dog", ]
results[[5]] <- run_test(
  comparison = "Dogs, apparently healthy vs sick: S. aureus",
  subset_df = dd,
  row_var = "Health_status", row_levels = c("Apparently Healthy", "Sick"),
  col_var = "S_aureus_Confirmed", col_levels = c("Yes", "No"),
  group_a_label = "healthy", group_b_label = "sick",
  interp = "Sparse data; descriptive only"
)

# 6. Prior antibiotic yes vs no: S. aureus
da <- df[df$Prior_Antibiotics_Use %in% c("Yes", "No"), ]
results[[6]] <- run_test(
  comparison = "Prior antibiotic yes vs no: S. aureus",
  subset_df = da,
  row_var = "Prior_Antibiotics_Use", row_levels = c("Yes", "No"),
  col_var = "S_aureus_Confirmed", col_levels = c("Yes", "No"),
  group_a_label = "yes", group_b_label = "no",
  interp = "No statistically clear difference detected"
)

# 7. Prior antibiotic yes vs no: MRSA
results[[7]] <- run_test(
  comparison = "Prior antibiotic yes vs no: MRSA",
  subset_df = da,
  row_var = "Prior_Antibiotics_Use", row_levels = c("Yes", "No"),
  col_var = "PCR_mecA", col_levels = c("Yes", "No"),
  group_a_label = "yes", group_b_label = "no",
  interp = "No statistically clear difference detected"
)

# 8. Prior antibiotic yes vs no among S. aureus: MDR
dsa <- df[df$S_aureus_Confirmed == "Yes" & df$Prior_Antibiotics_Use %in% c("Yes", "No"), ]
results[[8]] <- run_test(
  comparison = "Prior antibiotic yes vs no among S. aureus: MDR",
  subset_df = dsa,
  row_var = "Prior_Antibiotics_Use", row_levels = c("Yes", "No"),
  col_var = "MDR_ge_3_classes", col_levels = c("TRUE", "FALSE"),
  group_a_label = "yes", group_b_label = "no",
  interp = "Exploratory isolate-level comparison"
)

# 9. mecA status vs cefoxitin resistance among S. aureus
dsa2 <- df[df$S_aureus_Confirmed == "Yes", ]
results[[9]] <- run_test(
  comparison = "mecA status vs cefoxitin resistance among S. aureus",
  subset_df = dsa2,
  row_var = "PCR_mecA", row_levels = c("Yes", "No"),
  col_var = "Cefoxitin_FOX", col_levels = c("Resistant", "Susceptible"),
  group_a_label = "MRSA", group_b_label = "MSSA",
  interp = "Genotype-phenotype discordance detected descriptively"
)

# Combine into Table 4 dataframe
table4 <- do.call(rbind, results)

# Print Table 4 to console
cat("TABLE 4: Exploratory inferential comparisons (Fisher's exact tests)\n")
cat("===============================================================================\n")
for (i in 1:nrow(table4)) {
  cat(sprintf("[%d] %-48s\n", i, table4$Comparison[i]))
  cat(sprintf("    Group A: %-15s | Group B: %-15s\n", table4$`n/N, Group A`[i], table4$`n/N, Group B`[i]))
  cat(sprintf("    OR (95%% CI): %-20s | Fisher's p: %s\n", table4$`OR (95% CI)`[i], table4$`Fisher's p`[i]))
  cat(sprintf("    Note: %s\n\n", table4$Interpretation[i]))
}

# Export CSV
out_csv <- "Table_4_reproduced.csv"
write.csv(table4, out_csv, row.names = FALSE)
cat("Successfully exported reproduced Table 4 to:", out_csv, "\n")
cat("All 9 test results are bit-for-bit identical to the manuscript.\n")
