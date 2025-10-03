#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript create_samplesheet.R <path_to_fastq_folders> <samplesheet_name>", call.=FALSE)
}

folder  <- args[1] # "/mnt/N/NGS/3-Sekvenseringsbiblioteker/2022/Illumina_RunXXX/Run820_Virus/Run820/"
outfile <- args[2] # "2023.01.19-HCV_Run829.csv"

# Validate input folder exists
if (!dir.exists(folder)) {
  stop(paste("Error: Input folder does not exist:", folder), call.=FALSE)
}

# Validate output directory exists
output_dir <- dirname(outfile)
if (!dir.exists(output_dir)) {
  stop(paste("Error: Output directory does not exist:", output_dir), call.=FALSE)
}

# Get the fastq files
fastq <- list.files(folder,
           recursive = TRUE,
           full.names = TRUE,
           pattern = "\\.fastq\\.gz$|\\.fq\\.gz$")

# Filter for HCV files only
fastq <- fastq[str_detect(fastq, "HCV")]

# Check if any HCV files were found
if (length(fastq) == 0) {
  stop("Error: No HCV fastq files found in the specified folder", call.=FALSE)
}

cat(sprintf("Found %d HCV fastq files\n", length(fastq)))

R1 <- sort(fastq[grep("_R1_", fastq)])
R2 <- sort(fastq[grep("_R2_", fastq)])

# Check if we have both R1 and R2 files
if (length(R1) == 0) {
  stop("Error: No R1 files found", call.=FALSE)
}
if (length(R2) == 0) {
  stop("Error: No R2 files found", call.=FALSE)
}
if (length(R1) != length(R2)) {
  stop(sprintf("Error: Unequal number of R1 (%d) and R2 (%d) files", length(R1), length(R2)), call.=FALSE)
}

cat(sprintf("Found %d R1 and %d R2 files\n", length(R1), length(R2)))

df <- as_tibble(cbind(R1, R2))

# Check that the R1 and R2 files are correctly paired
tmp <- df %>%
  mutate(tmpR1 = gsub("_.*", "", basename(R1)),
         tmpR2 = gsub("_.*", "", basename(R2))) %>%
  select(tmpR1, tmpR2)

if (identical(tmp$tmpR1, tmp$tmpR2)) {
  df <- df %>%
    mutate(sample_id = gsub("_.*", "", basename(R1))) %>%
    mutate(
      original_sample_id = sample_id,
      # Remove invalid characters from sample names (keep alphanumeric characters and dashes)
      sample_id = str_replace_all(sample_id, "[^A-Za-z0-9-]", "")
    ) %>%
    select("sample" = sample_id,
           "fastq_1" = R1,
           "fastq_2" = R2,
           original_sample_id)
  
  # Check for empty sample names after cleaning
  empty_samples <- df %>% filter(sample == "")
  if (nrow(empty_samples) > 0) {
    cat("Error: Some samples have no valid characters left after cleaning:\n")
    for (i in seq_len(nrow(empty_samples))) {
      cat(sprintf("  '%s' resulted in empty sample name\n", empty_samples$original_sample_id[i]))
    }
    stop("Please check your file naming convention", call.=FALSE)
  }
  
  # Check for duplicate sample names
  duplicate_samples <- df %>% 
    group_by(sample) %>% 
    filter(n() > 1) %>%
    ungroup()
  
  if (nrow(duplicate_samples) > 0) {
    cat("Error: Duplicate sample names detected after cleaning:\n")
    duplicate_names <- unique(duplicate_samples$sample)
    for (name in duplicate_names) {
      originals <- df %>% filter(sample == name) %>% pull(original_sample_id)
      cat(sprintf("  Sample '%s' appears %d times (from: %s)\n", 
                  name, length(originals), paste(originals, collapse = ", ")))
    }
    stop("Please ensure unique sample names", call.=FALSE)
  }
  
  # Check if any sample names were modified and print a message
  modified_samples <- df %>%
    filter(sample != original_sample_id) %>%
    select(original_sample_id, sample)
  
  if (nrow(modified_samples) > 0) {
    cat("Warning: Sample names contained invalid characters and were modified:\n")
    for (i in seq_len(nrow(modified_samples))) {
      cat(sprintf("  '%s' -> '%s'\n", 
                  modified_samples$original_sample_id[i], 
                  modified_samples$sample[i]))
    }
    cat("Invalid characters (periods, spaces, and other symbols except dashes) have been removed.\n\n")
  }
  
  # Remove the helper column before writing output
  df <- df %>% select(-original_sample_id)
  
} else {
  stop("Error: R1 and R2 files not correctly paired", call.=FALSE)
}

# Attempt to write the output file
tryCatch({
  write_csv(df, outfile)
  cat(sprintf("Successfully created samplesheet with %d samples: %s\n", nrow(df), outfile))
}, error = function(e) {
  stop(paste("Error writing output file:", e$message), call.=FALSE)
})

