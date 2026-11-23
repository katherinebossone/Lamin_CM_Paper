library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(dplyr)  

# Promoters
keep_chrs <- paste0("chr", 1:19)

genes <- fread("../homer/H3K4me3_all_promoters_with_genes.bed")
colnames(genes) <- c("chr","start","end","interval","score","strand",
                      "chr_prom","start_prom","end_prom","gene_name")

promoters <- GRanges(
  seqnames = genes$chr,
  ranges   = IRanges(start = pmax(0, genes$start - 5000),
                      end   = genes$end + 5000),
  gene_name = genes$gene_name
)
promoters <- keepSeqlevels(promoters, keep_chrs, pruning.mode = "coarse")

# DE Genes
de_genes <- read.table("markers_combined_clusters_all_vent_LB1Fd_vs_LB1Fplus_stop.tsv",
                        header = TRUE, stringsAsFactors = FALSE)

sig_up   <- de_genes %>% filter(p_val_adj < 0.05, avg_log2FC >  0.26)
sig_down <- de_genes %>% filter(p_val_adj < 0.05, avg_log2FC < -0.26)
non_de_genes <- de_genes %>% filter(p_val_adj > 0.05, avg_log2FC > -0.26, avg_log2FC < 0.26)

# cScore tracks
atac_het  <- import("P0_5_LB1Fplus_5kb_cscore_modified.bedgraph")
atac_null <- import("P0_5_LB1Fd_5kb_cscore_modified.bedgraph")
score(atac_het)  <- -1 * score(atac_het) # if needed to flip so pos = A, neg = B
score(atac_null) <- -1 * score(atac_null) # if needed to flip so pos = A, neg = B

common_sl <- intersect(seqlevels(promoters), seqlevels(atac_het))
promoters <- keepSeqlevels(promoters, common_sl, pruning.mode = "coarse")
atac_het  <- keepSeqlevels(atac_het,  common_sl, pruning.mode = "coarse")
atac_null <- keepSeqlevels(atac_null, common_sl, pruning.mode = "coarse")

get_weighted_cscore <- function(regions, track) {
  hits <- findOverlaps(regions, track)
  if (length(hits) == 0) return(rep(NA_real_, length(regions)))
  ov <- pintersect(regions[queryHits(hits)], track[subjectHits(hits)])
  wdf <- tibble(
    region_id = queryHits(hits),
    score_val = score(track)[subjectHits(hits)],
    weight    = width(ov)
  )
  agg <- wdf %>%
    dplyr::group_by(region_id) %>%
    dplyr::summarise(mean_score = weighted.mean(score_val, weight), .groups = "drop")
  out <- rep(NA_real_, length(regions))
  out[agg$region_id] <- agg$mean_score
  out
}

# Promoter-level cScore + DE annotation
prom_df <- tibble(
  gene_name   = promoters$gene_name,
  chr         = as.character(seqnames(promoters)),
  start       = start(promoters),
  end         = end(promoters),
  cscore_het  = get_weighted_cscore(promoters, atac_het),
  cscore_null = get_weighted_cscore(promoters, atac_null)
) %>%
  filter(!is.na(cscore_het), !is.na(cscore_null)) %>%
  mutate(
    delta      = cscore_null - cscore_het,
    state_het  = ifelse(cscore_het  > 0, "A", "B"),
    state_null = ifelse(cscore_null > 0, "A", "B"),
    transition = paste0(state_het, "->", state_null)
  )

# Collapse to one row per gene (a gene may have multiple H3K4me3 promoter calls / isoforms) by averaging cscore across its promoter intervals
prom_df <- prom_df %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(
    cscore_het  = mean(cscore_het),
    cscore_null = mean(cscore_null),
    delta       = mean(delta),
    .groups = "drop"
  ) %>%
  mutate(
    state_het  = ifelse(cscore_het  > 0, "A", "B"),
    state_null = ifelse(cscore_null > 0, "A", "B"),
    transition = paste0(state_het, "->", state_null)
  )

# Attach DE category (Up / Down / Non-DE) + full DE stats
de_annot <- de_genes %>%
  transmute(
    gene_name = symbol,
    avg_log2FC,
    p_val_adj,
    de_category = case_when(
      symbol %in% sig_up$symbol       ~ "Up",
      symbol %in% sig_down$symbol     ~ "Down",
      symbol %in% non_de_genes$symbol ~ "Non-DE",
      TRUE ~ NA_character_
    )
  )

expr_compartment_df <- prom_df %>%
  inner_join(de_annot, by = "gene_name") %>%
  filter(!is.na(de_category)) %>%
  mutate(de_category = factor(de_category, levels = c("Non-DE","Down","Up")))


# Statistics
# Delta cScore differences at promoters by DE category
kw <- kruskal.test(delta ~ de_category, data = expr_compartment_df)
print(kw)
pairwise <- pairwise.wilcox.test(expr_compartment_df$delta, expr_compartment_df$de_category,
                                  p.adjust.method = "BH")
print(pairwise)

# log2FC vs delta Cscore (all DE-tested genes)
sp <- cor.test(expr_compartment_df$avg_log2FC, expr_compartment_df$delta, method = "spearman")
print(sp)

# Upregulated genes enrichment test for B->A switching relative to Non-DE/Down
switch_tab <- expr_compartment_df %>%
  dplyr::count(de_category, transition) %>%
  dplyr::group_by(de_category) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup()
print(switch_tab)

chisq_tab <- table(expr_compartment_df$de_category, expr_compartment_df$transition)
print(chisq.test(chisq_tab))

# Up vs Non-DE, B->A switch vs not
updf <- expr_compartment_df %>% filter(de_category %in% c("Up","Non-DE")) %>%
  mutate(switched_BtoA = transition == "B->A")
print(fisher.test(table(updf$de_category, updf$switched_BtoA)))

# Plot
p1 <- ggplot(expr_compartment_df, aes(x = de_category, y = delta, fill = de_category)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = expression(Delta*"Cscore at promoter (null - het)"),
       title = "Promoter compartment change vs transcriptional response") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")
ggsave("delta_cscore_by_DE_category.pdf", p1, width = 6, height = 5)

p2 <- ggplot(expr_compartment_df, aes(x = delta, y = avg_log2FC)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  labs(x = expression(Delta*"Cscore at promoter (null - het)"),
       y = "avg_log2FC (null vs het)",
       title = sprintf("Spearman rho = %.2f, p = %.2e", sp$estimate, sp$p.value)) +
  theme_bw(base_size = 13)
ggsave("logFC_vs_delta_cscore_scatter.pdf", p2, width = 6, height = 5)

p3 <- ggplot(switch_tab, aes(x = de_category, y = frac, fill = transition)) +
  geom_col(position = "stack") +
  labs(x = NULL, y = "Fraction of genes",
       title = "Promoter compartment transitions by DE category") +
  theme_bw(base_size = 13)
