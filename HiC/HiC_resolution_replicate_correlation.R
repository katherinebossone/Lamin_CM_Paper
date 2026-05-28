library(reticulate)
library(tidyverse)
library(ggplot2)

# Load environment
use_condaenv("r-reticulate2", required = TRUE)

cooler <- import("cooler")
np <- import("numpy")

# Load files and resolutions
samples <- c(
  "Rep1_Het"  = "P0_5_LB1Fplus_8450_2.mcool",
  "Rep2_Het"  = "P0_5_LB1Fplus_8450_4.mcool",
  "Rep1_Null" = "P0_5_LB1Fd_8450_7.mcool",
  "Rep2_Null" = "P0_5_LB1Fd_8450_8.mcool"
)

rep_pairs <- list(
  c("Rep1_Het", "Rep2_Het"),
  c("Rep1_Null", "Rep2_Null")
)

chromosomes <- paste0("chr", 1:19)
resolutions <- c(10000, 50000, 250000)

get_chr_matrix <- function(mcool_file, chr, res){
  path <- paste0(mcool_file, "::/resolutions/", res)
  c <- cooler$Cooler(path)
  mat <- c$matrix(balance = TRUE)$fetch(chr)
  matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat))
}

# Calculate correlations
results <- data.frame()

for(res in resolutions){
  cat("\nResolution:", res, "\n")
  
  for(pair in rep_pairs){
    cat(" Pair:", pair[1], "vs", pair[2], "\n")
    
    vec1_all <- c()
    vec2_all <- c()
    
    file1 <- samples[[pair[1]]]
    file2 <- samples[[pair[2]]]
    
    for(chr in chromosomes){
      cat("  ", chr, "\n")
      
      mat1 <- get_chr_matrix(file1, chr, res)
      mat2 <- get_chr_matrix(file2, chr, res)
      
      ut <- upper.tri(mat1, diag = TRUE)
      keep <- is.finite(mat1[ut]) & is.finite(mat2[ut])
      vec1_all <- c(vec1_all, mat1[ut][keep])
      vec2_all <- c(vec2_all, mat2[ut][keep])
    }
    
    corr <- cor(
      log1p(vec1_all),
      log1p(vec2_all),
      use = "complete.obs"
    )
    
    cat("  Correlation:", round(corr, 4), "\n")
    
    results <- rbind(
      results,
      data.frame(
        Pair = paste(pair, collapse = " vs "),
        Resolution = res,
        Correlation = corr
      )
    )
  }
}

# Plot
ggplot(results,
       aes(x = Resolution / 1000,
           y = Correlation,
           color = Pair,
           group = Pair)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = c(10, 50, 250),
    labels = c("10", "50", "250"),
    trans = "log10"
  ) +
  scale_y_continuous(limits = c(0.9, 1)) +
  labs(
    x = "Resolution (kb)",
    y = "Pearson Correlation",
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

