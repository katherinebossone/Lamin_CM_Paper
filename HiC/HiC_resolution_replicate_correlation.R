library(strawr)
library(reshape2)
library(Matrix)
library(tidyverse)

# Define samples and resolutions
samples <- c(
  "Rep1_Het" = "P0_5_LB1Fplus_8450_2.hic",
  "Rep2_Het" = "P0_5_LB1Fplus_8450_4.hic",
  "Rep1_Null" = "P0_5_LB1Fd_8450_7.hic",
  "Rep2_Null" = "P0_5_LB1Fd_8450_8.hic"
)
resolutions <- c(50000, 250000, 500000)

# Chromosomes to include
chromosomes <- paste0("chr", 1:19)  

# Functions to extract matrices
get_matrix_chr <- function(fname, chrom, res){
  df <- strawr::straw(
    norm = "KR",
    fname = fname,
    chr1loc = chrom,
    chr2loc = chrom,
    unit = "BP",
    binsize = res
  )
  
  if(nrow(df) == 0) return(NULL)
  
  # Columns in strawr v0.0.92 are x, y, counts
  mat <- acast(df, x ~ y, value.var = "counts", fill = NA)
  return(mat)
}

get_genome_matrix <- function(fname, chromosomes, res){
  mats <- lapply(chromosomes, function(chr) get_matrix_chr(fname, chr, res))
  mats <- mats[!sapply(mats, is.null)]
  if(length(mats) == 0) return(NULL)
  
  genome_mat <- as.matrix(Matrix::bdiag(mats))
  return(genome_mat)
}

# Define replicate pairings
rep_pairs <- list(
  c("Rep1_Het", "Rep2_Het"),
  c("Rep1_Null", "Rep2_Null")
)

# Calculate correlations
results <- data.frame()

for(res in resolutions){
  
  # Load genome matrices for all samples
  matrices <- lapply(samples, function(f) get_genome_matrix(f, chromosomes, res))
  
  # Compute correlation for each replicate pair
  for(pair in rep_pairs){
    mat1 <- matrices[[pair[1]]]
    mat2 <- matrices[[pair[2]]]
    
    flat1 <- as.vector(mat1)
    flat2 <- as.vector(mat2)
    mask <- !is.na(flat1) & !is.na(flat2)
    
    corr <- cor(flat1[mask], flat2[mask])
    
    results <- rbind(results, data.frame(
      Pair = paste(pair, collapse = " vs "),
      Resolution = res,
      Correlation = corr
    ))
  }
}

# -----------------------------
# 5. Line plot
# -----------------------------
ggplot(results, aes(x = Resolution/1000, y = Correlation, color = Pair)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0.99, 1)) +  # Set y-axis from 0 to 1
  labs(
    x = "Resolution (kb)",
    y = "Pearson Correlation",
    title = "Hi-C Replicate Correlation Across Resolutions"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
