# Load required libraries
library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(wateRmelon)
library(future)
library(gtools)
library(matrixStats)
library(data.table)
library(MASS)
library(sandwich)
library(lmtest)
library(parallel)
library(R.utils)
library(readxl)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(qqman)
library(ENmix)
library(DMRcate)

# Load input data
inputfile1 <- readRDS("cor_info.rds")
inputfile <- inputfile1 %>%
  select("CpGsite", "Estimate", "StdError", "pval") %>%
  rename(CpGID = CpGsite, coef = Estimate, se = StdError, pvalue = pval)

# Annotate CpG data
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotdata <- ann850k[c("chr", "pos", "UCSC_RefGene_Name", "Relation_to_Island")]

# Merge input data with annotation
inputfileBed <- merge(inputfile, annotdata, by = "row.names")
inputfileBed$pos <- as.integer(inputfileBed$pos)
inputfileBed$end <- inputfileBed$pos + 1
inputfileBed$stat <- inputfileBed$coef / inputfileBed$se
inputfileBed$CpGID <- as.character(inputfileBed$CpGID)
inputfileBed$CHR <- as.character(inputfileBed$chr)
inputfileBed$CHR1 <- gsub("chr", "", inputfileBed$CHR)

# Prepare Meta data
Meta <- data.frame(
  chr = inputfileBed$CHR1,
  start = inputfileBed$pos,
  end = inputfileBed$end,
  p = inputfileBed$pvalue,
  probe = inputfileBed$CpGID
) %>%
  na.omit()

saveRDS(Meta, "Meta.rds")

# Perform combp analysis
combp(Meta, dist.cutoff = 1000, bin.size = 310, seed = 0.01,
      region_plot = TRUE, mht_plot = TRUE, nCores = 10, verbose = TRUE)

# Load results and filter
ddt <- read.csv("resu_combp.csv", sep = ",", header = TRUE) %>%
  filter(sidak < 0.05, nprobe > 1)

# Save significant results
write.csv(ddt, "resu_combp.sig.csv", row.names = FALSE, quote = FALSE)

# Prepare data for BED file
ddt$chr <- paste0("chr", ddt$chr)
ddt$newname <- paste0(ddt$chr, "_", ddt$start, "_", ddt$end)
input <- ddt[, c(1:3, 9)]
write.table(input, "DMR.combp.bed", sep = "\t", row.names = FALSE, quote = FALSE)

# Load annotation data
krt <- read.csv("DMR.combp.anno.txt", sep = "\t", header = TRUE)
colnames(krt) <- c("Region", "Gene")

# Clean gene names
krt$Gene <- gsub("\\s*\\([^\\)]+\\)", "", as.character(krt$Gene))
krt$Gene <- gsub(" ", "", as.character(krt$Gene))
krt$Gene <- gsub(",", ";", as.character(krt$Gene))

# Annotate DMR results
for (i in seq_along(ddt$newname)) {
  reg <- ddt$newname[i]
  indx <- which(krt$Region %in% reg)
  
  if (length(indx) == 0) {
    ddt$newname[i] <- "NoGene"
  } else {
    ddt$newname[i] <- krt$Gene[indx]
  }
}

colnames(ddt)[8] <- "Genesymbol"
write.xlsx(ddt, "DMR.combp.Genesymbol.xlsx")

# DMRcate analysis
inputfile1 <- readRDS("cor_info.rds") %>%
  na.omit() %>%
  select("CpGsite", "Estimate", "StdError", "pval") %>%
  rename(CpGID = CpGsite, coef = Estimate, se = StdError, pvalue = pval)

# Merge with annotation
annotdata <- ann850k[c("chr", "pos", "UCSC_RefGene_Name", "Relation_to_Island")]
inputfileBed <- merge(inputfile, annotdata, by = "row.names")
inputfileBed$pos <- as.integer(inputfileBed$pos)
inputfileBed$end <- inputfileBed$pos + 1
inputfileBed$stat <- inputfileBed$coef / inputfileBed$se
inputfileBed$CpGID <- as.character(inputfileBed$CpGID)
inputfileBed$CHR <- as.character(inputfileBed$chr)
inputfileBed$CHR1 <- gsub("chr", "", inputfileBed$CHR)

# Create GRanges object for DMRcate
df2 <- data.frame(seqnames = inputfileBed$CHR, start = inputfileBed$pos, end = inputfileBed$end, strand = "*")
gr100 <- makeGRangesFromDataFrame(df2)
gr100$stat <- inputfileBed$stat
gr100$diff <- inputfileBed$coef
gr100$ind.fdr <- inputfileBed$pvalue
gr100$is.sig <- gr100$ind.fdr < 0.01

# DMRcate analysis
grl <- GRangesList(gr100)
class(grl) <- "CpGannotated"
grl@ranges <- grl@unlistData
grl <- changeFDR(grl, 0.01)

dmrcoutput <- dmrcate(grl, lambda = 1000, C = 2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")

saveRDS(results.ranges, "DMR.DMRcate.rds", compress = "xz")

# Filter DMR results
ddt <- readRDS("DMR.DMRcate.rds") %>%
  as.data.frame() %>%
  filter(no.cpgs > 1, min_smoothed_fdr < 0.05)

# Save filtered DMR results
saveRDS(ddt, "DMR.DMRcate.sign.rds", compress = "xz")

# Prepare DMR for BED file
ddt$chr <- ddt$seqnames
ddt$newname <- paste0(ddt$chr, "_", ddt$start, "_", ddt$end)
input <- ddt[, c(1:3, 15)]
write.table(input, "DMR.DMRcate.bed", sep = "\t", row.names = FALSE, quote = FALSE)

# Perform annotation from GREAT
krt <- read.csv("DMR.DMRcate.anno.txt", sep = "\t", header = TRUE)
colnames(krt) <- c("Region", "Gene")
krt$Gene <- gsub("\\s*\\([^\\)]+\\)", "", as.character(krt$Gene))
krt$Gene <- gsub(" ", "", as.character(krt$Gene))
krt$Gene <- gsub(",", ";", as.character(krt$Gene))

# Annotate DMR results
for (i in seq_along(ddt$newname)) {
  reg <- ddt$newname[i]
  indx <- which(krt$Region %in% reg)
  
  if (length(indx) == 0) {
    ddt$newname[i] <- "NoGene"
  } else {
    ddt$newname[i] <- krt$Gene[indx]
  }
}

colnames(ddt)[8] <- "Genesymbol"
write.xlsx(ddt, "DMR.DMRcate.Genesymbol.xlsx")