library("depmixS4")

# Load files
gr <- import("~/Desktop/Nkx2_5Cre_Lmnb1Fplus_lamin-A_zscore.bw", format = "BigWig")
# gr <- import("~/Desktop/Nkx2_5Cre_Lmnb1Fminus_lamin-A_zscore.bw", format = "BigWig")
dataset <- as.data.frame(gr)
log2 <- dataset$score

#Run HMM
mod <- depmix(log2 ~1, family = gaussian(), nstates = 2, data = dataset)
set.seed(123)
fm <- fit(mod, verbose = FALSE, emcontrol=em.control(crit ="absolute"))
print(fm)
summary(fm)
probs <- posterior(fm)
head(probs)

#Separate states
state <- probs$state
S1 <- probs$S1
S2 <- probs$S2
dataset$V5 <- state
dataset$V6 <- S1
dataset$V7 <- S2

write.table(dataset, file = "~/Desktop/LB1Fplus_LAC_2stateHMM.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Create bed (Switch files for mutant)
state <- dataset[dataset$V5==1,]
export <- state[1:3]
write.table(export, file = "~/Desktop/LB1Fplus_LAC_2state_S1plus.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
state <- dataset[dataset$V5==2,]
export <- state[1:3]
write.table(export, file = "~/Desktop/LB1Fplus_LAC_2state_S2plus.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Create bedgraphs (Switch files for mutant)
state <- dataset[dataset$V5==1,]
export <- state[1:4]
write.table(export, file = "~/Desktop/LB1Fplus_LAC_2state_S1plus.bedgraph", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
state <- dataset[dataset$V5==2,]
export <- state[1:4]
write.table(export, file = "~/Desktop/LB1Fplus_LAC_2state_S2plus.bedgraph", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

gr <- import("Desktop/bigwig_lamins/LB1null_LAC_merged.bw", format = "BigWig")
dataset <- as.data.frame(gr)
