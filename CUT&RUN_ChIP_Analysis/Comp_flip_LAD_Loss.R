library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(dplyr)


# Input

het_bed  <- "../LB1het_LAC_2state_S2.bed"
null_bed <- "../LB1null_LAC_2state_S2.bed"

het  <- import(het_bed, format = "BED")
null <- import(null_bed, format = "BED")
start(het)  <- start(het)  - 1
start(null) <- start(null) - 1

cscore_het  <- import("P0_5_LB1Fplus_5kb_cscore_modified.bedgraph")
cscore_null <- import("P0_5_LB1Fd_5kb_cscore_modified.bedgraph")

# Cscore sign convention needs flipping so that positive = A compartment, negative = B compartment (if needed)
score(cscore_het)  <- -1 * score(cscore_het)
score(cscore_null) <- -1 * score(cscore_null)

common_sl <- intersect(seqlevels(het), seqlevels(cscore_het))
het   <- keepSeqlevels(het,   common_sl, pruning.mode = "coarse")
null  <- keepSeqlevels(null,  common_sl, pruning.mode = "coarse")
cscore_het  <- keepSeqlevels(cscore_het,  common_sl, pruning.mode = "coarse")
cscore_null <- keepSeqlevels(cscore_null, common_sl, pruning.mode = "coarse")

# Define LAD categories
het_r  <- reduce(het)
null_r <- reduce(null)

lost_LAD   <- GenomicRanges::setdiff(het_r, null_r)   # LAD in het, lost in null
gained_LAD <- GenomicRanges::setdiff(null_r, het_r)   # LAD gained in null
stable_LAD <- GenomicRanges::intersect(het_r, null_r) # LAD in both

all_LAD  <- reduce(c(het_r, null_r))
# Background: everything covered by the cscore tracks that is never a LAD
genome_gr <- reduce(c(granges(cscore_het), granges(cscore_null)))
never_LAD <- GenomicRanges::setdiff(genome_gr, all_LAD)

cat(sprintf("Lost: %d regions (%.1f Mb)\n", length(lost_LAD),  sum(width(lost_LAD))/1e6))
cat(sprintf("Stable: %d regions (%.1f Mb)\n", length(stable_LAD), sum(width(stable_LAD))/1e6))
cat(sprintf("Gained: %d regions (%.1f Mb)\n", length(gained_LAD), sum(width(gained_LAD))/1e6))
cat(sprintf("Never-LAD: %d regions (%.1f Mb)\n", length(never_LAD), sum(width(never_LAD))/1e6))

# 2. Width-weighted mean cScore per region
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

build_df <- function(regions, category) {
  tibble(
    category   = category,
    chr        = as.character(seqnames(regions)),
    start      = start(regions),
    end        = end(regions),
    width      = width(regions),
    cscore_het  = get_weighted_cscore(regions, cscore_het),
    cscore_null = get_weighted_cscore(regions, cscore_null)
  )
}

df <- bind_rows(
  build_df(lost_LAD,   "Lost LAD"),
  build_df(stable_LAD, "Stable LAD"),
  build_df(gained_LAD, "Gained LAD"),
  build_df(never_LAD,  "Never LAD")
) %>%
  filter(!is.na(cscore_het), !is.na(cscore_null)) %>%
  mutate(
    delta      = cscore_null - cscore_het,
    state_het  = ifelse(cscore_het  > 0, "A", "B"),
    state_null = ifelse(cscore_null > 0, "A", "B"),
    transition = paste0(state_het, "->", state_null),
    category   = factor(category, levels = c("Never LAD","Stable LAD","Gained LAD","Lost LAD"))
  )

# Statistics
# Delta Cscore differences across categories
kw <- kruskal.test(delta ~ category, data = df)
print(kw)

pairwise <- pairwise.wilcox.test(df$delta, df$category, p.adjust.method = "BH")
print(pairwise)

# Within "Lost LAD" regions
lost_df <- df %>% filter(category == "Lost LAD")
wt <- wilcox.test(lost_df$cscore_het, lost_df$cscore_null, paired = TRUE)
print(wt)

# Proportion of B->A switching by category
switch_tab <- df %>%
  dplyr::count(category, transition) %>%
  dplyr::group_by(category) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup()
print(switch_tab)

chisq_tab <- table(df$category, df$transition)
print(chisq.test(chisq_tab))

# Plot
p1 <- ggplot(df, aes(x = category, y = delta, fill = category)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = expression(Delta*"Cscore (null - het)"),
       title = "Compartment score change by LAD category") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave("delta_cscore_by_LAD_category.pdf", p1, width = 6, height = 5)

p2 <- ggplot(switch_tab, aes(x = category, y = frac, fill = transition)) +
  geom_col(position = "stack") +
  labs(x = NULL, y = "Fraction of regions",
       title = "Compartment state transitions by LAD category") +
  theme_bw(base_size = 13)
