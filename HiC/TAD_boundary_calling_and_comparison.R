library(data.table)
library(zoo)
library(GenomicRanges)
library(VennDiagram)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Parameters
bin_size       <- 20000
smooth_k       <- 3
window_size    <- 8
tolerance_bins <- 2
tolerance      <- tolerance_bins * bin_size
flank_bins     <- 25

# Load FAN-C insulation files
load_insulation <- function(file) {
  ins <- fread(file, header = FALSE)
  setnames(ins, c("chr", "start", "end", "dot1", "insulation_score", "dot2"))
  ins[, insulation_score := as.numeric(insulation_score)]
  setorder(ins, chr, start)
  ins[, insulation_smooth := zoo::rollmean(insulation_score, k = smooth_k, fill = NA), by = chr]
  ins
}

# Call TAD boundaries
call_boundaries <- function(ins, min_score = -0.3) {
  boundary_list <- list()
  for (chr_name in unique(ins$chr)) {
    chr_data <- ins[chr == chr_name]
    chr_data[, boundary := FALSE]
    n <- nrow(chr_data)
    if (n <= (2 * window_size)) next
    for (i in (window_size + 1):(n - window_size)) {
      window_vals <- chr_data$insulation_smooth[(i - window_size):(i + window_size)]
      center_val  <- chr_data$insulation_smooth[i]
      if (!is.na(center_val) &&
          all(center_val < window_vals[-(window_size + 1)], na.rm = TRUE))
        chr_data$boundary[i] <- TRUE
    }
    boundary_list[[chr_name]] <- chr_data[boundary == TRUE &
                                            insulation_smooth < min_score]
  }
  boundaries <- rbindlist(boundary_list)
  GRanges(seqnames = boundaries$chr,
          ranges   = IRanges(start = boundaries$start, end = boundaries$end),
          score    = boundaries$insulation_smooth)
}

# Extract ±flank_bins insulation profile per boundary
extract_profiles <- function(gr, ins_dt, flank = flank_bins) {
  nc  <- 2 * flank + 1
  mat <- matrix(NA_real_, nrow = length(gr), ncol = nc)
  for (i in seq_along(gr)) {
    chr_name  <- as.character(seqnames(gr[i]))
    center_bp <- start(gr[i])
    chr_ins   <- ins_dt[chr == chr_name]
    setorder(chr_ins, start)
    idx <- which.min(abs(chr_ins$start - center_bp))
    s <- idx - flank; e <- idx + flank
    if (s >= 1 && e <= nrow(chr_ins)) {
      vals <- chr_ins$insulation_score[s:e]
      if (length(vals) == nc) mat[i, ] <- vals
    }
  }
  mat
}

## Run ##
wt_file  <- "Documents/LB1Fplus_fanc.insulation_1mb.bed"
mut_file <- "Documents/LB1Fd_fanc.insulation_1mb.bed"

wt_ins  <- load_insulation(wt_file)
mut_ins <- load_insulation(mut_file)

wt_gr  <- call_boundaries(wt_ins)
mut_gr <- call_boundaries(mut_ins)

# Categorize
ol <- findOverlaps(wt_gr, mut_gr, maxgap = tolerance)
wt_gr$category  <- "lost"
mut_gr$category <- "gained"
wt_gr$category[unique(queryHits(ol))]    <- "conserved"
mut_gr$category[unique(subjectHits(ol))] <- "conserved"

cat("Lmnb1+/- boundaries:\n");  print(table(wt_gr$category))
cat("\nLmnb1-/- boundaries:\n"); print(table(mut_gr$category))

# Group GRanges objects
gr_I   <- wt_gr[wt_gr$category == "conserved"]
gr_II  <- wt_gr[wt_gr$category == "lost"]
gr_III <- mut_gr[mut_gr$category == "gained"]

n_I   <- length(gr_I)
n_II  <- length(gr_II)
n_III <- length(gr_III)

# Venn diagram
wt_total  <- length(wt_gr)
mut_total <- length(mut_gr)
conserved <- sum(wt_gr$category == "conserved")

venn.plot <- draw.pairwise.venn(
  area1      = wt_total,
  area2      = mut_total,
  cross.area = conserved,
  category   = c("Lmnb1+/-", "Lmnb1-/-"),
  fill       = c("#87CEEB", "#E64B35"),
  lty        = "blank",
  cex        = 2,
  cat.cex    = 1.5,
  cat.pos    = c(-20, 20)
)

# Extract insulation profiles
clip2 <- function(m, v = 2) { m[m > v] <- v; m[m < -v] <- -v; m }

mI_wt    <- clip2(extract_profiles(gr_I,   wt_ins))
mI_mut   <- clip2(extract_profiles(gr_I,   mut_ins))
mII_wt   <- clip2(extract_profiles(gr_II,  wt_ins))
mII_mut  <- clip2(extract_profiles(gr_II,  mut_ins))
mIII_wt  <- clip2(extract_profiles(gr_III, wt_ins))
mIII_mut <- clip2(extract_profiles(gr_III, mut_ins))

# Order rows by central dip depth
center_dip <- function(mat, flank = flank_bins, half_win = 3) {
  ctr <- flank + 1
  idx <- (ctr - half_win):(ctr + half_win)
  rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
}

ord_I   <- order(center_dip(mI_wt))
ord_II  <- order(center_dip(mII_wt))
ord_III <- order(center_dip(mIII_mut))

mI_wt    <- mI_wt[ord_I,];    mI_mut   <- mI_mut[ord_I,]
mII_wt   <- mII_wt[ord_II,];  mII_mut  <- mII_mut[ord_II,]
mIII_wt  <- mIII_wt[ord_III,]; mIII_mut <- mIII_mut[ord_III,]

col_fn  <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
wt_ttl  <- expression(italic(Lmnb1)^"+/-")
mut_ttl <- expression(italic(Lmnb1)^"-/-")

# Make heatmap
mk_ht <- function(mat, nm, col_title = NULL, show_lgd = FALSE,
                  row_title = "", show_col_title = TRUE) {
  Heatmap(
    mat,
    name            = nm,
    col             = col_fn,
    column_title    = if (!is.null(col_title) && show_col_title) col_title else character(0),
    column_title_gp = gpar(fontsize = 11),
    row_title       = row_title,
    row_title_rot   = 0,
    row_title_gp    = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows    = FALSE, cluster_columns = FALSE,
    show_row_names  = FALSE, show_column_names = FALSE,
    show_heatmap_legend = show_lgd,
    use_raster      = TRUE, raster_quality = 2,
    border          = FALSE,
    width           = unit(3, "cm"),
    heatmap_legend_param = list(
      title     = "Insulation\nScore",
      at        = c(-2, -1, 0, 1, 2),
      direction = "vertical",
      title_gp  = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    )
  )
}

hI   <- mk_ht(mI_wt,   "I_wt",   wt_ttl, TRUE,  "I")   + mk_ht(mI_mut,   "I_mut",   mut_ttl, FALSE, " ")
hII  <- mk_ht(mII_wt,  "II_wt",  NULL,   FALSE, "II")  + mk_ht(mII_mut,  "II_mut",  NULL,    FALSE, " ")
hIII <- mk_ht(mIII_wt, "III_wt", NULL,   FALSE, "III") + mk_ht(mIII_mut, "III_mut", NULL,    FALSE, " ")

# Assemble heatmaps into one
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 1,
                                           heights = unit(c(n_I, n_II, n_III) / (n_I + n_II + n_III), "npc"))))

pushViewport(viewport(layout.pos.row = 1))
draw(hI, newpage = FALSE, padding = unit(c(0, 22, 0, 5), "mm"))
popViewport()

pushViewport(viewport(layout.pos.row = 2))
draw(hII, newpage = FALSE, show_heatmap_legend = FALSE,
     padding = unit(c(0, 22, 0, 5), "mm"))
popViewport()

pushViewport(viewport(layout.pos.row = 3))
draw(hIII, newpage = FALSE, show_heatmap_legend = FALSE,
     padding = unit(c(0, 22, 0, 5), "mm"))
popViewport()
