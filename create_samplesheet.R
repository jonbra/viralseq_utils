#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript create_samplesheet.R <path_to_fastq_folders> <samplesheet_name> [agent]", call.=FALSE)
}

folder  <- args[1]
outfile <- args[2]
# optional: if provided, filter fastq files by agent (case-insensitive); otherwise keep all
agens   <- if (length(args) >= 3) args[3] else "ALL"

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

# Optional agent filtering (case-insensitive). If agens == "ALL", skip filtering.
if (toupper(agens) != "ALL") {
  # Match the agent against the file name (case-insensitive). Matching the
  # basename rather than the full path avoids false matches from the parent
  # directory name (e.g. a temp dir called "fastq_hcv").
  fastq <- fastq[str_detect(basename(fastq), regex(agens, ignore_case = TRUE))]
  if (length(fastq) == 0) {
    stop(paste0("Error: No fastq files matching agent '", agens, "' found in the specified folder"), call.=FALSE)
  }
  cat(sprintf("Found %d fastq files matching agent '%s'\n", length(fastq), agens))
} else {
  cat(sprintf("Found %d fastq files\n", length(fastq)))
}

# Accept both Illumina-style _R1_/_R2_ and SRA-style _1/_2 pair labels.
# Files without these labels are treated as single-end files.
fastq_base <- basename(fastq)
is_R1 <- str_detect(fastq_base, regex("(?:^|_)R?1(?:_|\\.)", ignore_case = TRUE))
is_R2 <- str_detect(fastq_base, regex("(?:^|_)R?2(?:_|\\.)", ignore_case = TRUE))

fastq_info <- tibble(path = fastq, file = fastq_base) %>%
  mutate(
    read_type = case_when(
      is_R1 ~ "R1",
      is_R2 ~ "R2",
      TRUE ~ "SE"
    ),
    sample_id_raw = str_replace(
      file,
      regex("(?:_R?1(?:_[^_]+)?|_R?2(?:_[^_]+)?|_1(?:_[^_]+)?|_2(?:_[^_]+)?)?\\.f(?:ast)?q\\.gz$", ignore_case = TRUE),
      ""
    )
  )

# For paired-end labels, ensure we do not have duplicate R1 or R2 for the same sample.
duplicate_reads <- fastq_info %>%
  filter(read_type != "SE") %>%
  group_by(sample_id_raw, read_type) %>%
  filter(n() > 1) %>%
  ungroup()

if (nrow(duplicate_reads) > 0) {
  dup_names <- unique(paste0(duplicate_reads$sample_id_raw, " (", duplicate_reads$read_type, ")"))
  shown <- paste(head(dup_names, 10), collapse = ", ")
  extra <- if (length(dup_names) > 10) sprintf(" (and %d more)", length(dup_names) - 10) else ""
  stop(paste0("Error: Found multiple files for the same sample/read direction: ", shown, extra), call.=FALSE)
}

r1 <- fastq_info %>%
  filter(read_type == "R1") %>%
  transmute(sample_id_raw, fastq_1 = path)

r2 <- fastq_info %>%
  filter(read_type == "R2") %>%
  transmute(sample_id_raw, fastq_2 = path)

paired_df <- full_join(r1, r2, by = "sample_id_raw")
single_end_df <- fastq_info %>%
  filter(read_type == "SE") %>%
  transmute(sample_id_raw, fastq_1 = path, fastq_2 = "")

# If only one mate exists for a sample, keep it as single-end.
orphaned_pairs <- paired_df %>%
  filter(is.na(fastq_1) | is.na(fastq_2))

if (nrow(orphaned_pairs) > 0) {
  shown <- paste(head(orphaned_pairs$sample_id_raw, 10), collapse = ", ")
  extra <- if (nrow(orphaned_pairs) > 10) sprintf(" (and %d more)", nrow(orphaned_pairs) - 10) else ""
  cat(sprintf("Warning: %d sample(s) had only one mate file; keeping them as single-end: %s%s\n",
              nrow(orphaned_pairs), shown, extra))
  single_end_df <- bind_rows(
    single_end_df,
    orphaned_pairs %>%
      transmute(sample_id_raw, fastq_1 = if_else(is.na(fastq_1), fastq_2, fastq_1), fastq_2 = "")
  )
}

overlap_samples <- intersect(paired_df$sample_id_raw, single_end_df$sample_id_raw)
if (length(overlap_samples) > 0) {
  shown <- paste(head(overlap_samples, 10), collapse = ", ")
  extra <- if (length(overlap_samples) > 10) sprintf(" (and %d more)", length(overlap_samples) - 10) else ""
  cat(sprintf("Warning: %d sample(s) have both paired-end and single-end files; using paired-end and ignoring single-end duplicates: %s%s\n",
              length(overlap_samples), shown, extra))
  single_end_df <- single_end_df %>%
    filter(!sample_id_raw %in% overlap_samples)
}

paired_df <- paired_df %>%
  filter(!is.na(fastq_1) & !is.na(fastq_2))

cat(sprintf("Found %d paired and %d single-end fastq samples\n", nrow(paired_df), nrow(single_end_df)))

df <- bind_rows(paired_df, single_end_df) %>%
  mutate(
    original_sample_id = sample_id_raw,
    # Remove invalid characters from sample names (keep alphanumeric characters and dashes)
    sample = str_replace_all(sample_id_raw, "[^A-Za-z0-9-]", "")
  ) %>%
  select(sample, fastq_1, fastq_2, original_sample_id) %>%
  arrange(sample)

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

# Attempt to write the output file
tryCatch({
  write_csv(df, outfile)
  cat(sprintf("Successfully created samplesheet with %d samples: %s\n", nrow(df), outfile))
}, error = function(e) {
  stop(paste("Error writing output file:", e$message), call.=FALSE)
})

