library(data.table)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(rtracklayer)

# Load Hi-C data
load_hic_data <- function(file_path) {
  hic <- fread(file_path, header = FALSE)
  setnames(hic, c("chr1","start1","end1","chr2","start2","end2","Count_raw","Count"))
  hic <- hic[chr1 == chr2 & chr1 %in% keep_chrs]
  return(hic)
}

keep_chrs <- paste0("chr", 1:19)

# Load H3K4me3 promoters
genes <- fread("../homer/H3K4me3_all_promoters_with_genes.bed")
colnames(genes) <- c("chr","start","end","interval","score","strand",
                     "chr_prom","start_prom","end_prom","gene_name")

promoters <- GRanges(
  seqnames = genes$chr,
  ranges = IRanges(start = pmax(0, genes$start - 5000),
                   end = genes$end + 5000),
  gene_name = genes$gene_name
)

upstream_flank <- GRanges(
  seqnames = genes$chr,
  ranges = IRanges(start = pmax(0, genes$start - 5e6),
                   end = genes$start - 1),
  gene_name = genes$gene_name
)

downstream_flank <- GRanges(
  seqnames = genes$chr,
  ranges = IRanges(start = genes$end + 1,
                   end = genes$end + 5e6),
  gene_name = genes$gene_name
)

# Load HiLands states
hilands <- fread("state_location_fixed.bed")
hilands_gr <- GRanges(
  seqnames = hilands$chr,
  ranges = IRanges(start = hilands$start, end = hilands$end),
  state = hilands$state
)

# Load DE genes and define Non-DE
de_genes <- read.table("markers_combined_clusters_all_vent_LB1Fd_vs_LB1Fplus.tsv",
                       header = TRUE, stringsAsFactors = FALSE)

sig_up <- de_genes %>% filter(p_val_adj < 0.05, avg_log2FC > 0.26)
sig_down <- de_genes %>% filter(p_val_adj < 0.05, avg_log2FC < -0.26)

de_gene_symbols <- union(sig_up$symbol, sig_down$symbol)
non_de_genes <- de_genes %>% filter(p_val_adj > 0.05, avg_log2FC > -0.26, )

# Promoter GRanges subsets
prom_up <- promoters[promoters$gene_name %in% sig_up$symbol]
prom_down <- promoters[promoters$gene_name %in% sig_down$symbol]
prom_non <- promoters[promoters$gene_name %in% non_de_genes$symbol]

# Flank GRanges subsets
up_flank_up <- upstream_flank[upstream_flank$gene_name %in% sig_up$symbol]
down_flank_up <- downstream_flank[downstream_flank$gene_name %in% sig_up$symbol]

up_flank_down <- upstream_flank[upstream_flank$gene_name %in% sig_down$symbol]
down_flank_down <- downstream_flank[downstream_flank$gene_name %in% sig_down$symbol]

up_flank_non <- upstream_flank[upstream_flank$gene_name %in% non_de_genes$symbol]
down_flank_non <- downstream_flank[downstream_flank$gene_name %in% non_de_genes$symbol]

# Assign HiLands states to Hi-C bins
assign_hilands <- function(hic_df, hilands_gr) {
  bins1 <- GRanges(seqnames = hic_df$chr1,
                   ranges = IRanges(hic_df$start1, hic_df$end1))
  bins2 <- GRanges(seqnames = hic_df$chr2,
                   ranges = IRanges(hic_df$start2, hic_df$end2))
  
  overlap1 <- findOverlaps(bins1, hilands_gr)
  overlap2 <- findOverlaps(bins2, hilands_gr)
  
  hic_df$state1 <- NA
  hic_df$state2 <- NA
  hic_df$state1[queryHits(overlap1)] <- hilands_gr$state[subjectHits(overlap1)]
  hic_df$state2[queryHits(overlap2)] <- hilands_gr$state[subjectHits(overlap2)]
  
  return(hic_df)
}

# Scale Hi-C counts
scale_counts <- function(hic, target_total) {
  hic$Count_scaled <- hic$Count * (target_total / sum(hic$Count, na.rm=TRUE))
  return(hic)
}

# Filter promoter–flank contacts
filter_promoter_contacts <- function(hic_df, promoters, upstream_flank, downstream_flank) {
  bins1 <- GRanges(seqnames = hic_df$chr1,
                   ranges = IRanges(hic_df$start1, hic_df$end1))
  bins2 <- GRanges(seqnames = hic_df$chr2,
                   ranges = IRanges(hic_df$start2, hic_df$end2))
  
  hits1 <- findOverlaps(bins1, promoters)
  hits2 <- findOverlaps(bins2, c(upstream_flank, downstream_flank))
  keep1 <- intersect(queryHits(hits1), queryHits(hits2))
  
  hits1b <- findOverlaps(bins1, c(upstream_flank, downstream_flank))
  hits2b <- findOverlaps(bins2, promoters)
  keep2 <- intersect(queryHits(hits1b), queryHits(hits2b))
  
  keep <- unique(c(keep1, keep2))
  return(hic_df[keep, ])
}

# Aggregate contacts per gene by flank state
aggregate_contacts <- function(genes_gr, hic_df) {
  results <- list()
  
  for (i in seq_along(genes_gr)) {
    gene <- genes_gr[i]
    gene_name <- gene$gene_name
    chr <- as.character(seqnames(gene))
    
    overlaps <- hic_df %>%
      filter(chr1 == chr | chr2 == chr) %>%
      filter(
        (start1 <= end(gene) & end1 >= start(gene)) |
          (start2 <= end(gene) & end2 >= start(gene))
      )
    
    overlaps <- overlaps %>%
      mutate(
        flank_state = case_when(
          start1 <= end(gene) & end1 >= start(gene) ~ state2,
          start2 <= end(gene) & end2 >= start(gene) ~ state1,
          TRUE ~ NA_character_
        )
      )
    
    agg <- overlaps %>%
      group_by(flank_state) %>%
      summarise(contact_sum = sum(Count_scaled, na.rm=TRUE),
                .groups = "drop") %>%
      mutate(gene = gene_name)
    
    results[[i]] <- agg
  }
  
  return(rbindlist(results))
}

# Process replicates
process_replicate <- function(file_path, promoters, up_flank, down_flank, target_total) {
  hic <- load_hic_data(file_path)
  hic <- assign_hilands(hic, hilands_gr)
  hic <- scale_counts(hic, target_total)
  hic <- filter_promoter_contacts(hic, promoters, up_flank, down_flank)
  counts <- aggregate_contacts(promoters, hic)
  return(counts)
}

wt_files  <- c("P0_5_HiC/indiv/P0_5_LB1Fplus_8450_2_5kb.tsv",
               "P0_5_HiC/indiv/P0_5_LB1Fplus_8450_4_5kb.tsv")
mut_files <- c("P0_5_HiC/indiv/P0_5_LB1Fd_8450_7_5kb.tsv",
               "P0_5_HiC/indiv/P0_5_LB1Fd_8450_8_5kb.tsv")

all_files <- c(wt_files, mut_files)
target_total <- mean(sapply(all_files, function(f) sum(load_hic_data(f)$Count, na.rm=TRUE)))

# Calculate counts for promoters and flanking states
run_group <- function(files, promoters, up_flank, down_flank, regulation, condition) {
  lst <- lapply(files, process_replicate,
                promoters, up_flank, down_flank,
                target_total)
  for (i in seq_along(lst)) {
    lst[[i]]$rep <- paste0(condition, i)
    lst[[i]]$condition <- condition
    lst[[i]]$regulation <- regulation
  }
  return(lst)
}

# Usable if sufficient memory
# wt_up   <- run_group(wt_files, prom_up, up_flank_up, down_flank_up, "Up", "WT")
# mut_up  <- run_group(mut_files, prom_up, up_flank_up, down_flank_up, "Up", "Mut")
# wt_down <- run_group(wt_files, prom_down, up_flank_down, down_flank_down, "Down", "WT")
# mut_down<- run_group(mut_files, prom_down, up_flank_down, down_flank_down, "Down", "Mut")
# wt_non <- run_group(wt_files, prom_non, up_flank_non, down_flank_non, "NonDE", "WT")
# mut_non<- run_group(mut_files, prom_non, up_flank_non, down_flank_non, "NonDE", "Mut")

# Break into chunks if insufficient memory
chunk_size <- 4000
chunks <- split(prom_non, ceiling(seq_along(prom_non)/chunk_size))

#--------------------------non-------------------------------------
results_mut_non <- lapply(chunks, function(chunk) {
  genes_in_chunk <- chunk$gene_name
  run_group(
    mut_files,
    chunk,
    up_flank_non[up_flank_non$gene_name %in% genes_in_chunk],
    down_flank_non[down_flank_non$gene_name %in% genes_in_chunk],
    "NonDE",
    "Mut"
  )
})

results_wt_non <- lapply(chunks, function(chunk) {
  genes_in_chunk <- chunk$gene_name
  run_group(
    wt_files,
    chunk,
    up_flank_non[up_flank_non$gene_name %in% genes_in_chunk],
    down_flank_non[down_flank_non$gene_name %in% genes_in_chunk],
    "NonDE",
    "Wt"
  )
})

results_mut_non1 <- lapply(seq_along(results_mut_non[[1]]), function(r) {
  rbindlist(lapply(results_mut_non, `[[`, r))
})
results_wt_non1 <- lapply(seq_along(results_wt_non[[1]]), function(r) {
  rbindlist(lapply(results_wt_non, `[[`, r))
})

#--------------------------UP-------------------------------------
chunks <- split(prom_up, ceiling(seq_along(prom_up)/chunk_size))
results_mut_up <- lapply(chunks, function(chunk) {
  genes_in_chunk <- chunk$gene_name
  run_group(
    mut_files,
    chunk,
    up_flank_up[up_flank_up$gene_name %in% genes_in_chunk],
    down_flank_up[down_flank_up$gene_name %in% genes_in_chunk],
    "Up",
    "Mut"
  )
})

results_wt_up <- lapply(chunks, function(chunk) {
  genes_in_chunk <- chunk$gene_name
  run_group(
    wt_files,
    chunk,
    up_flank_up[up_flank_up$gene_name %in% genes_in_chunk],
    down_flank_up[down_flank_up$gene_name %in% genes_in_chunk],
    "Up",
    "Wt"
  )
})

results_mut_up1 <- lapply(seq_along(results_mut_up[[1]]), function(r) {
  rbindlist(lapply(results_mut_up, `[[`, r))
})
results_wt_up1 <- lapply(seq_along(results_wt_up[[1]]), function(r) {
  rbindlist(lapply(results_wt_up, `[[`, r))
})

#--------------------------DOWN-------------------------------------
chunks <- split(prom_down, ceiling(seq_along(prom_down)/chunk_size))
results_mut_down <- lapply(chunks, function(chunk) {
  genes_in_chunk <- chunk$gene_name
  run_group(
    mut_files,
    chunk,
    up_flank_down[up_flank_down$gene_name %in% genes_in_chunk],
    down_flank_down[down_flank_down$gene_name %in% genes_in_chunk],
    "Down",
    "Mut"
  )
})

results_wt_down <- lapply(chunks, function(chunk) {
  genes_in_chunk <- chunk$gene_name
  run_group(
    wt_files,
    chunk,
    up_flank_down[up_flank_down$gene_name %in% genes_in_chunk],
    down_flank_down[down_flank_down$gene_name %in% genes_in_chunk],
    "Down",
    "Wt"
  )
})

results_mut_down1 <- lapply(seq_along(results_mut_down[[1]]), function(r) {
  rbindlist(lapply(results_mut_down, `[[`, r))
})
results_wt_down1 <- lapply(seq_along(results_wt_down[[1]]), function(r) {
  rbindlist(lapply(results_wt_down, `[[`, r))
})

all_counts <- rbindlist(c(results_wt_up, results_mut_up, results_wt_down, results_mut_down, results_wt_non, results_mut_non))

# Calculate average and fold-change
avg_counts <- all_counts %>%
  group_by(gene, flank_state, regulation, condition) %>%
  summarise(mean_contact = mean(contact_sum, na.rm=TRUE), .groups="drop")

final_counts <- avg_counts %>%
  pivot_wider(names_from = condition, values_from = mean_contact) %>%
  mutate(fold_change = log2((Mut)/ (Wt))) %>%
  filter(is.finite(fold_change))

# Wilcoxon tests

wilcox_up_non <- final_counts %>%
  group_by(flank_state) %>%
  summarise(
    p_value = wilcox.test(fold_change[regulation=="Up"],
                          fold_change[regulation=="NonDE"])$p.value,
    .groups="drop"
  )

wilcox_down_non <- final_counts %>%
  group_by(flank_state) %>%
  summarise(
    p_value = wilcox.test(fold_change[regulation=="Down"],
                          fold_change[regulation=="NonDE"])$p.value,
    .groups="drop"
  )

## Plot fold change for HiLand flank states ##
final_counts <- final_counts %>% drop_na(flank_state)

custom_colors <- c("Up" = "aliceblue", "Down" = "darkgray", "NonDE" = "lightgrey")

final_counts$flank_state <- factor(final_counts$flank_state, 
                                   levels = c("P", "B", "G", "Y", "O", "R"))

# Plot
ggplot(final_counts, aes(x = flank_state, y = fold_change, fill = regulation)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA) +
  geom_jitter(aes(color = regulation),
              position = position_jitterdodge(jitter.width = 0.2),
              alpha = 0.1, size = 0.1, color = "black") +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Flanking HiLands State",
       y = "log2(Mut/WT)") +
  scale_fill_manual(values = custom_colors) + 
  scale_color_manual(values = custom_colors) +  
  theme_minimal()
