library(edgeR)

# Import data
counts <- read.delim("Desktop/RNA/P7.5_RNA/counts_only.txt", header=T, row.names=1)
metadata <- read.delim("Desktop/RNA/P7.5_RNA/Names.txt", header=T)

# Define group and create DGE object
if (ncol(metadata) > 2) {
  groupinfo <- metadata[,2]
  for (i in 3:ncol(metadata)) {
    groupinfo <- paste(groupinfo, metadata[,i], sep=".")
  }
  Group <- factor(groupinfo)
} else {
  Group <- factor(metadata[,2])
}

dge <- DGEList(counts=counts, group=Group)
dge$samples$lib.size <- colSums(dge$counts)

# Filter lowly expressed genes
cpmVal     <- 0.26
gThreshold <- 4

cpm.count <- cpm(dge$counts)
keep      <- rowSums(cpm.count > cpmVal) >= gThreshold

dge <- dge[keep, keep.lib.sizes=FALSE]

# Normalize
dge <- calcNormFactors(dge, method="TMM")

# Make matrix
f      <- droplevels(dge$samples$group)
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)
rownames(design) <- rownames(dge$samples)

# Estimate dispersion and fit GLM model
dge <- estimateDisp(dge, design)

plotBCV(dge, main="Biological Coefficient of Variation")

fit <- glmFit(dge, design)

# Define comparison
group1 <- "Het"
group2 <- "Null"

stopifnot(group1 %in% colnames(design))
stopifnot(group2 %in% colnames(design))

comp <- makeContrasts(
  contrasts = paste(group2, group1, sep="-"),
  levels    = design
)

# DE test
lrt <- glmLRT(fit, contrast=comp)

# Filter and annotate
fcval <- 1.2
fdr   <- 0.05

tp     <- topTags(lrt, n=Inf, adjust.method="BH", sort.by="none")
filter <- decideTestsDGE(lrt, adjust.method="BH", p.value=fdr, lfc=log2(fcval))
filter_vec <- as.integer(filter)
tp$table$filter     <- as.factor(filter_vec)
tp$table$Expression <- ifelse(filter_vec ==  1, "Upregulated",
                              ifelse(filter_vec == -1, "Downregulated",
                                     "Not Significant"))

# Summary
cat("Upregulated:  ", sum(filter_vec ==  1), "\n")
cat("Downregulated:", sum(filter_vec == -1), "\n")
cat("Not Sig:      ", sum(filter_vec ==  0), "\n")

plotMDS(dge, col=as.numeric(f), main="MDS Plot")
legend("topright", legend=levels(f), col=1:nlevels(f), pch=16)

plotSmear(lrt, de.tags=rownames(tp$table)[filter_vec != 0],
          main=paste("Smear Plot:", group2, "vs", group1))
abline(h=c(-log2(fcval), log2(fcval)), col="blue", lty=2)

# Save file of genes
outfile <- paste0("Downloads/edgeR_", group2, "_vs_", group1, ".csv")
write.csv(tp$table, outfile, row.names=TRUE)

# Volcano plot

df <- tp$table

# Thresholds
fdr_thresh <- 0.05
lfc_thresh <- 0.26

# Color by significance
df$color <- "gray"
df$color[df$FDR < fdr_thresh & df$logFC >  lfc_thresh] <- "red"
df$color[df$FDR < fdr_thresh & df$logFC < -lfc_thresh] <- "blue"
df$neg_log10_FDR <- -log10(df$FDR)

ggplot(df, aes(x = logFC, y = neg_log10_FDR, color = color)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_identity() +
  
  geom_hline(yintercept = -log10(fdr_thresh),
             linetype = "22", linewidth = 0.4, color = "black") +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
             linetype = "22", linewidth = 0.4, color = "black") +

  scale_x_continuous(breaks = pretty(df$logFC, n = 4),
                     minor_breaks = NULL) +
  scale_y_continuous(breaks = pretty(df$neg_log10_FDR, n = 4),
                     minor_breaks = NULL) +
  
  labs(x = "log2 Fold Change",
       y = "-log10(FDR)",
       title = paste("Volcano Plot:", group2, "vs", group1)) +
  
  theme(
    panel.grid.major = element_line(color = "#d0d2d3", linewidth = 0.5),
    panel.grid.minor = element_line(color = "#d0d2d3", linewidth = 0.25),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )
