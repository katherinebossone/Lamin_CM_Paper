library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomeInfoDb)


# Import HiLands BED and RNA bigwig
hilands <- read.table("state_merged_sorted.bed", header = FALSE)
colnames(hilands) <- c("chr", "start", "end", "state")

hilands_gr <- GRanges(
  seqnames = hilands$chr,
  ranges   = IRanges(start = hilands$start, end = hilands$end),
  state    = hilands$state
)

clusters_gr <- import("LB1Fplus_clusters_fixed.bw", format = "bigwig")
seqlevelsStyle(clusters_gr) <- "UCSC"

# Import GTF
gtf_df <- read.delim("../genome/mm10-iGenomes.gtf", header = FALSE, stringsAsFactors = FALSE)
colnames(gtf_df) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Keep exons
exons_df <- gtf_df[gtf_df$type == "exon", ]

parse_attributes <- function(attr_string) {
  kvs <- strsplit(attr_string, ";")[[1]]
  kvs <- trimws(kvs)
  gene_id <- NA
  gene_name <- NA
  for (kv in kvs) {
    if (startsWith(kv, "gene_id ")) gene_id <- sub('gene_id "(.*)"', '\\1', kv)
    if (startsWith(kv, "gene_name ")) gene_name <- sub('gene_name "(.*)"', '\\1', kv)
  }
  return(c(gene_id, gene_name))
}

parsed <- t(sapply(exons_df$attributes, parse_attributes))
colnames(parsed) <- c("gene_id", "gene_name")
exons_df$gene_id <- parsed[, "gene_id"]
exons_df$gene_name <- parsed[, "gene_name"]

# Convert to GRanges
exons <- GRanges(
  seqnames = exons_df$chr,
  ranges   = IRanges(start = exons_df$start, end = exons_df$end),
  strand   = exons_df$strand,
  gene_id  = exons_df$gene_id,
  gene_name = exons_df$gene_name
)

# Merge exons per gene and compute gene length in kb
genes_list <- split(exons, mcols(exons)$gene_id)
genes_grl <- GRangesList(lapply(genes_list, GenomicRanges::reduce))
genes_gr <- unlist(genes_grl)
genes_gr$gene_id <- rep(names(genes_list), elementNROWS(genes_grl))
genes_gr$length_kb <- width(genes_gr) / 1000


# Keep autosomal chromosomes in both objects
keep_chr <- paste0("chr", c(1:19))

genes_gr <- keepSeqlevels(genes_gr, keep_chr, pruning.mode="coarse")
clusters_gr <- keepSeqlevels(clusters_gr, keep_chr, pruning.mode="coarse")
hilands_gr <- keepSeqlevels(hilands_gr, keep_chr, pruning.mode="coarse")
hits <- findOverlaps(genes_gr, clusters_gr)
df <- data.frame(
  gene_id = genes_gr$gene_id[queryHits(hits)],
  count   = mcols(clusters_gr)$score[subjectHits(hits)]
)

# Sum counts per gene
gene_counts <- df %>%
  group_by(gene_id) %>%
  summarise(raw_counts = sum(count, na.rm = TRUE)) %>%
  ungroup()

# Add gene lengths
gene_counts <- gene_counts %>%
  left_join(data.frame(gene_id = genes_gr$gene_id, length_kb = genes_gr$length_kb), by = "gene_id") %>%
  distinct(gene_id, .keep_all = TRUE)  # remove duplicates

# Compute TPM per gene
# TPM = (counts / length_kb) / sum(counts/length_kb) * 1e6
gene_counts <- gene_counts %>%
  mutate(rpk = raw_counts / length_kb) %>%
  mutate(TPM = rpk / sum(rpk) * 1e6)

# Map genes to HiLands bins
genes_gr$TPM <- gene_counts$TPM[match(genes_gr$gene_id, gene_counts$gene_id)]
hits_bins <- findOverlaps(hilands_gr, genes_gr)

df_bins <- data.frame(
  hiland_id = queryHits(hits_bins),
  TPM = genes_gr$TPM[subjectHits(hits_bins)]
)

# Sum TPM per bin
bin_TPM <- df_bins %>%
  group_by(hiland_id) %>%
  summarise(sum_TPM = sum(TPM, na.rm = TRUE)) %>%
  ungroup()

# Add HiLands state
bin_TPM$state <- hilands_gr$state[bin_TPM$hiland_id]

# Plot
bin_TPM$state <- factor(bin_TPM$state, levels = c("P", "B", "G", "Y", "O", "R"))

ggplot(bin_TPM, aes(x = state, y = log10(sum_TPM), fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = c(
    "P"="#a020f0","B"="#0000ff","G"="#006400",
    "Y"="#ffd700","O"="#ffa500","R"="#ff0000"
  )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  ylab("log10(TPM)") +
  xlab("HiLands state")
