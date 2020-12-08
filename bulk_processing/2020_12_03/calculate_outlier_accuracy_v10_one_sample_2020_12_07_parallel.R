
library(tidyverse)
library(parallel)

# Non-sample-specific data

setwd("/private/home/hbeale/deleteMe/outlier_probability/outlier_lead_accuracy")

per_gene_expression_info <- read_tsv(paste0("data/in/per_gene_expression.tsv.gz"), col_types = cols())

# Sample-specific data
input_id <- "CKCC_2020_12_04"
output_id <- "CKCC_2020_12_04"
dir.create(output_id)
setwd(output_id)
dir.create("results")

test_input_raw_all <- read_tsv(paste0("/private/groups/treehouse/archive/projects/accuracy/outlier_probability_input_", input_id, "_batch1.tsv.gz"), col_types = cols())

# if toy
# test_input_raw <- test_input_raw_all %>%
#  sample_n(10)

test_input_raw <- test_input_raw_all 

test_input <- test_input_raw %>% 
  rename(Depth = MEND,
         Sample = sample_id,
         Length = effective_length) %>%
  filter(! grepl("/", Gene))

# first run library calls and outlier_prob function def below

reports <- test_input %>%
  rowwise() %>%
  do(params = as.list(.))

set.seed(1)
reports_rnd <- reports
reports_rnd$params <- sample(reports_rnd$params, length(reports_rnd$params))


### NOW GO RUN outlier_prob FUNCTION DEFINITION 

# test
# lapply(reports$params[1:2], outlier_prob)

# dummy_var <- lapply(reports_rnd$params, outlier_prob)
list_of_summaries <- mclapply(reports_rnd$params, outlier_prob, mc.cores = 36)
# list_of_summaries <- lapply(reports_rnd$params, outlier_prob) # for troubleshooting
summaries <- bind_rows(list_of_summaries)

head(summaries)
glimpse(summaries)

nrow(test_input) # 50727
nrow(summaries) # 50727
nrow(na.omit(summaries)) # 48665

# which columns have NA
apply(summaries, 2, function(x) any(is.na(x)))

# view some rows with NA
summaries %>% filter(is.na(`Probability of being an outlier`))
summaries %>% filter(is.na(`Probability of being an outlier`))  %>% glimpse

# do all expressions of zero have NA?
sum(summaries$Expression==0)
# [1] 35
sum(is.na(summaries$`Probability of being an outlier`))
# [1] 35

# how many NAs does each column have?
apply(summaries, 2, function(x) sum(is.na(x)))

# Samples 3-5 have > 35 NA Probability values; other samples and columns have 35
summaries %>% filter(is.na(`Sample 4 Probability`), ! is.na(`Sample 2 Probability`))
summaries %>% filter(is.na(`Sample 4 Probability`), ! is.na(`Sample 2 Probability`)) %>% glimpse

# why do these results have NA for the three samples?
# 
# summaries %>% filter(is.na(`Sample 4 Probability`), ! is.na(`Sample 2 Probability`)) %>%
#   write_tsv(paste0("/private/groups/treehouse/archive/projects/accuracy/outlier_probability_output_", output_id, "_problematic.tsv.gz"))

write_tsv(summaries, paste0("/private/groups/treehouse/archive/projects/accuracy/outlier_probability_output_", output_id, ".tsv.gz"))


## ----setup, include = FALSE----------------------------------------------

library(readr)
library(jsonlite)
library(magrittr)
library(dplyr)
library(tidyr)
library(knitr)
library(stringr)

outlier_prob <- function(params){
  # params = reports_rnd$params[[1]]
  # outlier_prob <- function(sample_id, outlier_lead){
  
  ## ----echo = TRUE---------------------------------------------------------
  
  sample_id <- params$Sample
  outlier_lead <- params$Gene
  exp <- params$Expression
  threshold <- params$Threshold
  type <- params$Type
  mend_depth <- params$Depth
  gene_length <- params$Length
  
  print(paste0("sample_", sample_id, "_gene_", outlier_lead, "  ", date()))
  
  output_file_name <- paste0("results/", sample_id, "/sample_", sample_id, "_gene_", outlier_lead, "_summary.tsv")
  
  if (!file.exists(output_file_name)){
    
    ## ----echo = TRUE---------------------------------------------------------
    
    # Check expression level
    if (exp >= 20) {
      stop("expression >= 20. make sure your expression is in log2(TPM+1)")
    }
    
    # Figure out expression bin
    
    expression_bin_boundaries <- c(0, 1, 3, 5, 7, 10, 20)
    mend_depth_bin_boundaries <- c(0, seq(2, 46, 4), Inf)*1E6
    length_bin_boundaries <- c(0, 100, 450, 800, 1500, 4000, Inf)
    
    make_bin_names <- function(boundaries = 1:10){
      raw_names <- paste0(boundaries, "-", lead(boundaries)) %>% str_replace("-Inf", "+")
      raw_names[1:(length(boundaries)-1)]
    }
    
    expression_bin_names <- make_bin_names(expression_bin_boundaries)
    mend_depth_bin_names <- make_bin_names(mend_depth_bin_boundaries/1e6)
    length_bin_names <- make_bin_names(length_bin_boundaries)
    
    expression_bin <- cut(exp, expression_bin_boundaries,
                          labels = expression_bin_names, include.lowest = TRUE)
    mend_depth_bin <- cut(mend_depth, mend_depth_bin_boundaries,
                          labels = mend_depth_bin_names, include.lowest = TRUE)
    length_bin <- cut(gene_length, length_bin_boundaries,
                      labels = length_bin_names, include.lowest = TRUE)
    
    
    
    
    ## ----echo = TRUE---------------------------------------------------------
    
    # Start building summary table
    
    threshold <- as.numeric(threshold)
    exp <- as.numeric(exp)
    
    summary <- tibble(sample = sample_id,
                      gene = outlier_lead,
                      type = type,
                      expression = exp,
                      threshold = threshold,
                      mend_depth = mend_depth,
                      gene_length = gene_length)
    
    
    # Calculate percent difference
    
    summary <- summary %>%
      mutate(percent_difference = 100 * ((abs(expression - threshold)) / expression))
    
    # Get accuracy distributions
    
    
    
    if (!(is.na(summary$percent_difference))) {
      # Set accuracy threshold and get accuracy info
      non_zero_abundance <- per_gene_expression_info %>%
        mutate(within_accuracy_threshold_of_deepest_val =
                 expression > (1 - summary$percent_difference / 100) * expression_at_max_depth &
                 expression < (1 + summary$percent_difference / 100) * expression_at_max_depth)
      # Make expression bins and generate summary statistics
      abundance_by_expression_bin <- non_zero_abundance %>% 
        ungroup() %>% 
        mutate(expression_level_bin = cut(expression_at_max_depth, expression_bin_boundaries,
                                          labels = expression_bin_names, include.lowest = TRUE),
               effective_gene_length_bin = cut(effective_gene_length, length_bin_boundaries,
                                               labels = length_bin_names, include.lowest = TRUE),
               depth_bin = cut(UMEND, mend_depth_bin_boundaries,
                               labels = mend_depth_bin_names, include.lowest = TRUE))
      abundance_stats <- abundance_by_expression_bin	%>% 
        group_by(expression_level_bin, depth_bin, effective_gene_length_bin, parent_id) %>%
        summarize(pct_accurately_measured = sum(within_accuracy_threshold_of_deepest_val) / n(),
                  n_genes_in_bin = length(unique(gene)),
                  n_measurements_in_bin = n())
    }
    
    # Filter results
    
    if (!(is.null(abundance_stats))) {
      stats_filtered <- abundance_stats %>%
        ungroup %>%
        filter((depth_bin == mend_depth_bin) & 
                 (expression_level_bin == expression_bin) & 
                 (effective_gene_length_bin == length_bin)) %>%
        arrange(parent_id)
    }    
    
    
    ## ----echo = TRUE, warning = FALSE----------------------------------------
    
    # Get accuracy stats
    
    summary <- summary %>%
      mutate(
        avg = mean(stats_filtered$pct_accurately_measured),
        min = min(stats_filtered$pct_accurately_measured),
        max = max(stats_filtered$pct_accurately_measured),
        S1 = stats_filtered$pct_accurately_measured[1],
        S2 = stats_filtered$pct_accurately_measured[2],
        S3 = stats_filtered$pct_accurately_measured[3],
        S4 = stats_filtered$pct_accurately_measured[4],
        S5 = stats_filtered$pct_accurately_measured[5]
      )
    
    
    
    
    ## ----echo = TRUE---------------------------------------------------------
    
    # Add column that checks if probability >= 0.95
    
    summary <- add_column(summary, prob95 = summary$avg >= 0.95, .after = "percent_difference")
    
    # Clean up names for writing to outfile
    
    names(summary) <- gsub("_", " ", str_to_sentence(names(summary)))
    
    summary <- summary %>%
      mutate(Avg = ifelse((sum(stats_filtered$n_measurements_in_bin) < 100 | nrow(stats_filtered) < 5), NA, Avg)) %>%
      rename("Probability >= 0.95" = Prob95,
             "Probability of being an outlier" = Avg,
             "Minimum Probability" = Min,
             "Maximum Probability" = Max,
             "Sample 1 Probability" = S1,
             "Sample 2 Probability" = S2,
             "Sample 3 Probability" = S3,
             "Sample 4 Probability" = S4,
             "Sample 5 Probability" = S5)
    
    
    # Write to outfile
    
    # if ( ! dir.exists(paste0("results/", sample_id))) dir.create(paste0("results/", sample_id))
    # write_tsv(summary, paste0("results/", sample_id, "/sample_", sample_id, "__gene_", outlier_lead, "__summary.tsv"))
    
    # Show short summary table & long summary table
 
    return(summary)   
  }
  
}
