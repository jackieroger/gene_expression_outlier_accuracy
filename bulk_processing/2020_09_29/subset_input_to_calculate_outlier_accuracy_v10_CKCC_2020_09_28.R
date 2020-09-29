library(tidyverse)
library(janitor) # just for data viewing, not processing
# Non-sample-specific data

setwd("/private/home/hbeale/deleteMe/outlier_probability/outlier_lead_accuracy")


# Sample-specific data
batch_input_id <- "CKCC_2020_07_13"
batch_output_id <- "CKCC_2020_09_28"

test_input_raw_all <- read_tsv(paste0("/private/groups/treehouse/archive/projects/accuracy/outlier_probability_input_", batch_input_id, ".tsv"), col_types = cols())

test_input_raw <- test_input_raw_all 

test_input <- test_input_raw %>% 
  rename(Depth = MEND,
         Sample = sample_id,
         Length = effective_length) %>%
  filter(! grepl("/", Gene))

# first run library calls and outlier_prob function def below

    
## ----echo = TRUE---------------------------------------------------------

expression_bin_boundaries <- c(0, 1, 3, 5, 7, 10, 20)
mend_depth_bin_boundaries <- c(0, seq(2, 46, 4), Inf)*1E6
length_bin_boundaries <- c(0, 100, 250, 400, 600, 900, 1200, 1500, 2000, 3000, 4000, Inf)

system.time(test_input_with_expr_bin <- test_input_raw_all %>%
  mutate(expr_bin = cut(Expression, expression_bin_boundaries, 
                        include.lowest = TRUE)))

system.time(test_input_with_mend_bin <- test_input_with_expr_bin %>%
              mutate(mend_bin = cut(MEND, mend_depth_bin_boundaries)))

system.time(test_input_with_length_bin <- test_input_with_mend_bin %>%
              mutate(length_bin = cut(effective_length, length_bin_boundaries, 
                                      include.lowest = TRUE)))

tab_by_bin <- test_input_with_length_bin %>% 
  group_by(expr_bin, mend_bin, length_bin) %>%
  mutate(n_entries = n(),
         n_to_sample = ifelse(n_entries < 100, n_entries, 100)) %>%
  sample_n(n_to_sample)

# how many bins could there be (max)
length(expression_bin_boundaries) * length(mend_depth_bin_boundaries) * length(length_bin_boundaries)


# how many bins have fewer than 100 values?
tab_by_bin %>%
  group_by(expr_bin, mend_bin, length_bin) %>%
  summarize(n_post_sample_entries = n(),
         has_100_values = n_post_sample_entries==100) %>%
  tabyl(has_100_values)

  
# test
tab_by_bin %>%
  group_by(expr_bin, mend_bin, length_bin) %>%
  mutate(n_post_sample_entries = n()) %>%
  pull(n_post_sample_entries) %>%
  summary

not_selected <- test_input_with_length_bin %>%
  anti_join(tab_by_bin, 
            by=c("sample_id", "Gene"))


output1_name <- paste0("/private/groups/treehouse/archive/projects/accuracy/outlier_probability_input_", batch_output_id, "_batch1.tsv.gz")
output2_name <- paste0("/private/groups/treehouse/archive/projects/accuracy/outlier_probability_input_", batch_output_id, "_batch2.tsv.gz")

tab_by_bin %>%
  ungroup %>%
  select(colnames(test_input_raw_all)) %>%
  write_tsv(output1_name)
  

not_selected %>%
  ungroup %>%
  select(colnames(test_input_raw_all)) %>%
  write_tsv(output2_name)

# use as input for the next step:
# "/private/groups/treehouse/archive/projects/accuracy/outlier_probability_input_CKCC_2020_09_28_batch1.tsv.gz"
