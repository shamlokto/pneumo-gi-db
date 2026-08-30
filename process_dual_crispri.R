#!/usr/bin/env Rscript
# Simplified Dual-CRISPRi GI score calculation
# Based on veeninglab/Dual-CRISPRi Dual-CRISPRi.Rmd

library(DESeq2)
library(reshape2)
library(dplyr)
library(tidyr)
library(stringr)

# Set working directory
setwd("/workspace")

# Load sgRNA targets
D39V_targets <- read.table("targets_operon.txt", header=T)[c(1,26,27)]
colnames(D39V_targets) <- c("sgRNA", "Locus.Tag", "Gene.Name")
D39V_targets$sgRNA <- paste0("sgRNA", as.numeric(gsub("sgRNA", "", D39V_targets$sgRNA)))

cat("Loaded", nrow(D39V_targets), "sgRNA targets\n")

# ============================================================================
# Process 869x869 library
# ============================================================================
cat("\n=== Processing 869x869 library ===\n")

raw869 <- read.csv2("veeninglab-Dual-CRISPRi-b9f977f/Read Counts/raw869_NextSeq_NovaSeq.csv", sep=",", header=T, row.names=1)
cat("Raw data:", nrow(raw869), "rows,", ncol(raw869), "columns\n")

# Aggregating the lanes for each replicate
raw869_noLane <- melt(raw869) %>%
  separate('variable', into=c("Sequencer", "Treatment", "Replicate", "Lane"), sep="[_]")

raw869_noLane <- aggregate(.~SG1+SG2+Sequencer+Treatment+Replicate, raw869_noLane[-c(6)], sum)
raw869_noLane <- dcast(raw869_noLane, SG1+SG2~Sequencer+Treatment+Replicate)

raw869_ref <- raw869_noLane[raw869_noLane$SG1 == raw869_noLane$SG2,]
raw869_comb <- raw869_noLane[raw869_noLane$SG1 != raw869_noLane$SG2,]

cat("Reference (single sgRNA):", nrow(raw869_ref), "rows\n")
cat("Combinations (dual sgRNA):", nrow(raw869_comb), "rows\n")

# Create DDS object
raw869_colData <- data.frame(s=colnames(raw869_noLane[3:14]), row.names=colnames(raw869_noLane[3:14])) %>%
  separate(col="s", into=c("Sequencer", "Treatment", "Repl"))

row.names(raw869_comb) <- paste0(raw869_comb$SG1, "_", raw869_comb$SG2)
row.names(raw869_ref) <- raw869_ref$SG1

DDS_old <- DESeqDataSetFromMatrix(countData = raw869_comb[-c(1,2)],
                                  colData = raw869_colData,
                                  design = ~ Treatment + Sequencer) %>%
  estimateSizeFactors()

DDS_old_single <- DESeqDataSetFromMatrix(countData = raw869_ref[-c(1,2)],
                                  colData = raw869_colData,
                                  design = ~ Treatment + Sequencer) %>%
  estimateSizeFactors()

# Run DESeq2
cat("Running DESeq2...\n")
alpha <- 0.05
LFCt <- 1

DDS_old = DESeq(DDS_old)
DDS_old_single = DESeq(DDS_old_single)

DDS_old_res <- results(DDS_old, name = "Treatment_WithIPTG_vs_NoIPTG", lfcThreshold = LFCt, alpha = alpha)
DDS_old_single_res <- results(DDS_old_single, name = "Treatment_WithIPTG_vs_NoIPTG", lfcThreshold = LFCt, alpha = alpha)

df_old_res_shr <- as.data.frame(lfcShrink(DDS_old, res = DDS_old_res, type = "apeglm", coef = "Treatment_WithIPTG_vs_NoIPTG"))
df_old_single_res_shr <- as.data.frame(lfcShrink(DDS_old_single, res = DDS_old_single_res, type = "apeglm", coef = "Treatment_WithIPTG_vs_NoIPTG"))

df_old_single_res_shr$sgRNA <- row.names(df_old_single_res_shr)
df_old_res_shr$sgRNA <- row.names(df_old_res_shr)
df_old_res_shr <- separate(df_old_res_shr, 'sgRNA', into=c("SG1", "SG2"), sep="[_]", remove = T)

df_old_res_shr$avReadsNI <- rowMeans(raw869_comb[grepl("NoIPTG",  names(raw869_comb))])
df_old_res_shr$avReadsI <- rowMeans(raw869_comb[grepl("WithIPTG",  names(raw869_comb))])

df_old_res_shr$SG1.targets <- D39V_targets$Gene.Name[match(df_old_res_shr$SG1, D39V_targets$sgRNA)]
df_old_res_shr$SG2.targets <- D39V_targets$Gene.Name[match(df_old_res_shr$SG2, D39V_targets$sgRNA)]
df_old_res_shr$SG1.locus <- D39V_targets$Locus.Tag[match(df_old_res_shr$SG1, D39V_targets$sgRNA)]
df_old_res_shr$SG2.locus <- D39V_targets$Locus.Tag[match(df_old_res_shr$SG2, D39V_targets$sgRNA)]

df_old_res_shr$targets <- paste0(df_old_res_shr$SG1.targets, " | ", df_old_res_shr$SG2.targets)

df_old_res_shr$SG1.refLog2FC <- df_old_single_res_shr$log2FoldChange[match(df_old_res_shr$SG1, df_old_single_res_shr$sgRNA)]
df_old_res_shr$SG1.reflfcSE <- df_old_single_res_shr$lfcSE[match(df_old_res_shr$SG1, df_old_single_res_shr$sgRNA)]
df_old_res_shr$SG1.refpadj <- df_old_single_res_shr$padj[match(df_old_res_shr$SG1, df_old_single_res_shr$sgRNA)]
df_old_res_shr$SG2.refLog2FC <- df_old_single_res_shr$log2FoldChange[match(df_old_res_shr$SG2, df_old_single_res_shr$sgRNA)]
df_old_res_shr$SG2.reflfcSE <- df_old_single_res_shr$lfcSE[match(df_old_res_shr$SG2, df_old_single_res_shr$sgRNA)]
df_old_res_shr$SG2.refpadj <- df_old_single_res_shr$padj[match(df_old_res_shr$SG2, df_old_single_res_shr$sgRNA)]

df_old_res_shr <- df_old_res_shr %>%
  mutate(pairs = case_when(
    SG1.refpadj < alpha & SG2.refpadj < alpha ~ "E-E",
    SG1.refpadj > alpha & SG2.refpadj > alpha ~ "NE-NE",
    SG1.refpadj < alpha | SG2.refpadj < alpha ~ "NE-E"
  ))

df_old_res_shr <- transform(df_old_res_shr, expectedSum=SG1.refLog2FC+SG2.refLog2FC)
df_old_res_shr$expectedReads <- df_old_res_shr$avReadsNI*2^df_old_res_shr$expectedSum

# Score Epsilon: ε = Wxy − E(Wxy)
df_old_res_shr <- df_old_res_shr %>% mutate(
  epsilonSum = ifelse(avReadsI < 4 & expectedReads < 2, NA, log2FoldChange-expectedSum)
)

df_old_res_shr <- df_old_res_shr %>% mutate(
  interactionSum = 
    case_when(
      avReadsI < 4 & expectedReads < 2 ~ NA,
      pairs == "NE-NE" ~ ifelse(epsilonSum < -1 & padj < 0.05, "Negative", "Neutral"),
      pairs == "NE-E" & epsilonSum > 1.5 & baseMean > 10 ~"Positive",
      pairs == "NE-E" & epsilonSum < -1 & padj < 0.05 ~ "Negative",
      pairs == "E-E" & epsilonSum > 1.5 & baseMean > 10 ~"Positive",
      pairs == "E-E" & epsilonSum < -1 & padj < 0.05 ~ "Negative",
      TRUE ~ "Neutral"
    )
)

df_old_res_shr$Group <- factor(df_old_res_shr$interactionSum, levels=c("Positive", "Neutral", "Negative", "Synthetic Lethal"))
df_old_res_shr$Group[which(df_old_res_shr$Group == "Negative" & df_old_res_shr$log2FoldChange < -7)] <- "Synthetic Lethal"

# Summary
cat("\n869x869 Results:\n")
print(table(df_old_res_shr$interactionSum, useNA="ifany"))

# Export
write.table(df_old_res_shr, "/workspace/Result_869x869.txt", quote=F, row.names=F, sep="\t")
cat("\nSaved: /workspace/Result_869x869.txt\n")

# ============================================================================
# Process 19x1499 library
# ============================================================================
cat("\n=== Processing 19x1499 library ===\n")

sgRNAList <- c(
  "sgRNA3","sgRNA348","sgRNA493","sgRNA593","sgRNA768",
  "sgRNA785","sgRNA812","sgRNA932","sgRNA187","sgRNA294",
  "sgRNA440","sgRNA500","sgRNA628","sgRNA758","sgRNA788",
  "sgRNA870","sgRNA1029","sgRNA1357","sgRNA1500"
)

rawMini <- read.csv2("veeninglab-Dual-CRISPRi-b9f977f/Read Counts/raw19x1499_NovaSeq.csv", sep=",")
cat("Raw data:", nrow(rawMini), "rows,", ncol(rawMini), "columns\n")

# Aggregating the lanes for each replicate
rawMini <- separate(melt(rawMini), variable, into=c("Treatment", "Repl", "Lane"), sep="[_]")
rawMini <- aggregate(.~SG1+SG2+Treatment+Repl, rawMini[-c(5)], sum)
rawMini <- dcast(rawMini, SG1+SG2~Treatment+Repl)

# Formatting input
rawMini_format <- data.frame(row.names=paste0("sgRNA", 1:1499))

for(i in sgRNAList) {
  rawMini_format[paste(i,colnames(rawMini[-c(1,2)]), sep="_")] <- 0
}

library(pbapply)
pbapply(rawMini, 1, function(row) {
  v <- c()
  if(row[[1]] %in% sgRNAList) {
    v <- paste(row[[1]], colnames(rawMini[-c(1,2)]), sep="_")
    for(colName in v) {
      rawMini_format[row[[2]], colName] <<- as.numeric(row[[gsub(paste0(row[[1]], "_"), "", colName)]])
    }
  }
  if(row[[2]] %in% sgRNAList) {
    v <- paste(row[[2]], colnames(rawMini[-c(1,2)]), sep="_")
    for(colName in v) {
      rawMini_format[row[[1]], colName] <<- as.numeric(row[[gsub(paste0(row[[2]], "_"), "", colName)]])
    }
  }
})
rawMini_format[is.na(rawMini_format)] <- 0

# Normalising
lapply(sgRNAList, function(sg) {
  rawMini_format[sg, grepl(paste0(sg,"_"), colnames(rawMini_format))] <<- round(mean(as.numeric(rawMini_format[sg, !grepl(paste0(sg,"_"), colnames(rawMini_format))])))
})
rawMini_format[is.na(rawMini_format)] <- 0

# Removing the last line (sgRNA1500 => luciferase)
rawMini_format <- rawMini_format[1:1499,]

# DDS object
coldataMini <- data.frame(s=colnames(rawMini_format), row.names=colnames(rawMini_format))
coldataMini <- separate(coldataMini, col="s", into=c("Background", "Treatment", "Repl"))

DDS_mini <- DESeqDataSetFromMatrix(countData = rawMini_format, 
                                  colData = coldataMini, 
                                  design = ~ Background * Treatment + Repl)

DDS_mini <- estimateSizeFactors(DDS_mini)

# DESeq2
cat("Running DESeq2...\n")
alpha = 0.05
LFCt = 1
DDS_mini$Background = relevel(DDS_mini$Background, ref="sgRNA1500")

DDS_mini = DESeq(DDS_mini)

DDS_mini_res <- pblapply(paste0("Background", sgRNAList[1:18]), function(sg) {
  st <- as.character(paste0(sg, ".TreatmentWithIPTG"))
  res <- results(DDS_mini, lfcThreshold = LFCt, name=st)
  res_shr <- as.data.frame(lfcShrink(DDS_mini, res = res, type = "apeglm", coef = st))
  res_shr$sgRNA=row.names(res_shr)
  return(res_shr)
})
names(DDS_mini_res) <- sgRNAList[1:18]

df_mini <- bind_rows(DDS_mini_res, .id="Background")

# Reference (single sgRNA)
rawMini_refDF <- rawMini_format[grepl("sgRNA1500", colnames(rawMini_format))]

DDS_mini_ref_colData <- data.frame(s=colnames(rawMini_refDF), row.names=colnames(rawMini_refDF))
DDS_mini_ref_colData <- separate(DDS_mini_ref_colData, col="s", into=c("Background", "Treatment", "Repl"))[-c(1)]
DDS_mini_ref <- DESeqDataSetFromMatrix(countData = rawMini_refDF,
                                  colData = DDS_mini_ref_colData,
                                  design = ~ Treatment + Repl)
DDS_mini_ref <- estimateSizeFactors(DDS_mini_ref)
DDS_mini_ref <- DESeq(DDS_mini_ref)
DDS_mini_ref_res <- results(DDS_mini_ref, lfcThreshold = LFCt, name="Treatment_WithIPTG_vs_NoIPTG", alpha=alpha)
df_mini_ref <- as.data.frame(lfcShrink(DDS_mini_ref, res = DDS_mini_ref_res, type = "apeglm", coef = "Treatment_WithIPTG_vs_NoIPTG"))
df_mini_ref$sgRNA = row.names(df_mini_ref)
df_mini_ref$targets <- D39V_targets$Gene.Name[match(df_mini_ref$sgRNA, D39V_targets$sgRNA)]

df_mini$bg.targets <- D39V_targets$Gene.Name[match(df_mini$Background, D39V_targets$sgRNA)]
df_mini$targets <- D39V_targets$Gene.Name[match(df_mini$sgRNA, D39V_targets$sgRNA)]
df_mini$bg.locus <- D39V_targets$Locus.Tag[match(df_mini$Background, D39V_targets$sgRNA)]
df_mini$locus <- D39V_targets$Locus.Tag[match(df_mini$sgRNA, D39V_targets$sgRNA)]

# Export
write.table(df_mini, "/workspace/Result_19x1499.txt", quote=F, row.names=F, sep="\t")
cat("\nSaved: /workspace/Result_19x1499.txt\n")

cat("\n=== Done ===\n")
