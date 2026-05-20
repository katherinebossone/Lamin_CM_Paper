library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)


# Define files
file_paths <- list(
  WT_rep1 = "P0_5_HiC/indiv/P0_5_LB1Fplus_8450_2_5kb.tsv",
  WT_rep2 = "P0_5_HiC/indiv/P0_5_LB1Fplus_8450_4_5kb.tsv",
  KO_rep1 = "P0_5_HiC/indiv/P0_5_LB1Fd_8450_7_5kb.tsv",
  KO_rep2 = "P0_5_HiC/indiv/P0_5_LB1Fd_8450_8_5kb.tsv"
)

# Read in files and rename columns
column_names <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "count")

data_list <- lapply(names(file_paths), function(name) {
  read_tsv(file_paths[[name]], col_names = FALSE) %>%
    setNames(c(column_names, paste0("extra", 1:(ncol(.) - length(column_names))))) %>%
    select(all_of(column_names)) %>%
    mutate(sample = name,
           count = as.numeric(count))
})

data <- bind_rows(data_list)

data <- data %>%
  mutate(genotype = case_when(
    grepl("WT", sample) ~ "WT",
    grepl("KO", sample) ~ "KO",
    TRUE ~ "Other"
  ))

# Scale counts to equal total reads per sample
total_counts <- data %>%
  group_by(sample) %>%
  summarise(total = sum(count), .groups = "drop")

median_total <- median(total_counts$total)
scale_factors <- setNames(median_total / total_counts$total, total_counts$sample)

data <- data %>%
  mutate(scaled_count = count * scale_factors[sample])

# Filter for autosomal chromosomes
valid_chroms <- paste0("chr", c(1:19))

inter_data <- data %>%
  filter(chr1 %in% valid_chroms & chr2 %in% valid_chroms) %>%
  mutate(pair = ifelse(chr1 < chr2,
                       paste(chr1, chr2, sep = "_"),
                       paste(chr2, chr1, sep = "_")))

# === Sum scaled counts per chromosome pair and sample ===
pair_summary <- inter_data %>%
  group_by(genotype, sample, pair) %>%
  summarise(total = sum(scaled_count), .groups = "drop")

# === Average across replicates per genotype ===
genotype_avg <- pair_summary %>%
  group_by(genotype, pair) %>%
  summarise(avg_contact = mean(total), .groups = "drop")

genotype_avg <- genotype_avg %>%
  separate(pair, into = c("chrA", "chrB"), sep = "_")

# Duplicate to make symmetric matrix
genotype_avg_sym <- bind_rows(
  genotype_avg,
  genotype_avg %>% filter(chrA != chrB) %>%
    mutate(tmp = chrA, chrA = chrB, chrB = tmp) %>%
    select(-tmp)
)

heatmap_data <- genotype_avg_sym %>%
  pivot_wider(
    names_from = genotype,
    values_from = avg_contact
  ) %>%
  mutate(log2FC = log2(KO / WT) * 100)


# Ensure correct chromosome order
chrom_order <- paste0("chr", c(1:19))

heatmap_data <- heatmap_data %>%
  filter(chrA %in% chrom_order & chrB %in% chrom_order) %>%
  mutate(
    chrA = factor(chrA, levels = chrom_order, ordered = TRUE),
    chrB = factor(chrB, levels = chrom_order, ordered = TRUE)
  )

# Plot
ggplot(heatmap_data, aes(x = chrA, y = chrB, fill = log2FC)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "log2FC(KO / WT) * 100") +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Chromosome", y = "Chromosome")
