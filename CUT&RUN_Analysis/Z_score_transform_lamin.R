library(rtracklayer)
library(GenomicRanges)

# Import bigWig
data <- import("LB1het_LAC_merged.bw")
# data <- import("LB1null_LAC_merged.bw")

# Extract signal
signal <- data$score

# Compute z-score
avg   <- mean(signal, na.rm = TRUE)
stdev <- sd(signal, na.rm = TRUE)

data$score <- (signal - avg) / stdev

# Export as bigWig
export(data, "~/Desktop/Nkx2_5Cre_Lmnb1Fplus_lamin-A_zscore.bw", format = "bigWig")
#export(data, "~/Desktop/Nkx2_5Cre_Lmnb1Fminus_lamin-A_zscore.bw", format = "bigWig")
