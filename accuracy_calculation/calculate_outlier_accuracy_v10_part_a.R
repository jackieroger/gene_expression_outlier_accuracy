library(tidyverse)


# Non-sample-specific data

per_gene_expression_info <- read_tsv(paste0("data/in/per_gene_expression.tsv.gz"), col_types = cols())


# Sample-specific data

test_input_raw <- read_tsv(paste0("test_input.tsv"), col_types = cols())

reports <- test_input_raw %>%
  rowwise() %>%
  do(params = as.list(.))

test_input <- unite(test_input_raw, file_name, remove = FALSE) %>%
  mutate(output_file = paste0(file_name, ".html"))  

reports$output_file <- test_input$output_file  

reports %>% 
  purrr::pwalk(rmarkdown::render, input = "calculate_outlier_accuracy_v10_part_b.Rmd")
