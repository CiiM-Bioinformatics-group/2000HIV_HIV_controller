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
library(tidyverse)
library(batchelor)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(qqman)
library(edgeR)
library(biomaRt)
library(org.Hs.eg.db)
library(bacon)

# Load datasets
beta2 <- readRDS("../ValidatedCPG.dis.MVal.rds")
result <- readRDS("../cor_info.rds")
dat1 <- readRDS("CombineRNAseq.rds")
samp <- readRDS("2000HIV_bulk_transcriptomics_sample_table.RDS")
dat11 <- readRDS("2000HIV_bulk_transcriptomics_normalized_counts.RDS")
pheno3 <- readRDS("../../2Pheno/NewPheno.Elite.April23.dis.rds")
pheno5 <- read.csv("../../2Pheno//covid_groups_selected_variables_may.csv")

# Annotation data
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotdata <- ann850k[c("chr", "pos", "UCSC_RefGene_Name", "Relation_to_Island")]

# Merge result with annotation data
Ml <- merge(result, annotdata, by = "row.names")
Ml <- data.frame(Ml)
Ml$BP <- as.integer(Ml$pos)
Ml$SNP <- as.character(Ml$Row.names)
Ml$CHR <- as.numeric(gsub("[^0-9]", "", Ml$chr))
Ml$P <- Ml$pval
Ml1 <- Ml %>% filter(!is.na(P))

# Update RNAseq data column names with DONOR_ID
for (i in seq_along(samp$ID)) {
  idx <- which(samp$ID == colnames(dat1)[i])
  if (length(idx) > 0) {
    colnames(dat1)[i] <- samp$DONOR_ID[idx]
  }
}

# Filter out sex chromosome genes
dat11 <- dat11 %>% filter(CHR %in% c("X", "Y"))
dat1 <- dat1 %>% filter(!SYMBOL %in% dat11$SYMBOL)
rownames(dat1) <- dat1$SYMBOL

# Filter phenotype data
pheno3 <- pheno3 %>%
  filter(ETHNICITY == "White") %>%
  filter(CONTROLLER_1B != "") %>%
  mutate(CONTROLLER_YESNO = case_when(
    CONTROLLER_1B == "EC_persistent" ~ 1,
    CONTROLLER_1B == "EC_transient" ~ 1,
    CONTROLLER_1B == "Non-EC" ~ 0,
    TRUE ~ 0
  ))

# Merge phenotype data with COVID data
newdata <- pheno3 %>%
  rowwise() %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .))) %>%
  left_join(pheno5, by = c("ID_idat" = "ID"))

pheno3 <- newdata

# Find common IDs
combs <- list(
  a = pheno3$ID_idat,
  b = colnames(beta2),
  c = colnames(dat1)
)
ItemsList <- venn(combs, show.plot = FALSE)
Common_up <- attributes(ItemsList)$intersections$`a:b:c`

# Subset phenotype data
pheno3 <- pheno3 %>% filter(ID_idat %in% Common_up)

# Align data matrices with phenotype data
beta2 <- beta2[, colnames(beta2) %in% pheno3$ID_idat]
dat1 <- dat1[, colnames(dat1) %in% pheno3$ID_idat]
pheno3 <- pheno3 %>% filter(ID_idat %in% colnames(beta2))

# Ensure the order of columns in data matrices match the phenotype data
beta2 <- beta2[, match(pheno3$ID_idat, colnames(beta2))]
dat1 <- dat1[, match(pheno3$ID_idat, colnames(dat1))]

# EQTM analysis
tx_annotation <- read.delim("ID2SYMBOL_gencode_v27_transcript.txt", header = FALSE, 
                            col.names = c("GENEID", "TXNAME", "SYMBOL", "GENETYPE"))

correction <- function(toptable) {
  bc <- bacon(effectsizes = as.matrix(toptable$logFC), standarderrors = as.matrix(toptable$logFC / toptable$t))
  toptable$logFC.cor <- es(bc)
  toptable$t.cor <- tstat(bc)
  toptable$P.Value.cor <- pval(bc)
  toptable$adj.P.Val.cor <- p.adjust(toptable$P.Value.cor, method = "BH")
  return(toptable)
}

do.twas <- function(cpg, covariates, sample_data, counts_data, filter = 0.8) {
  sample_data[, cpg] <- dnam[cpg, ]
  vars <- sample_data[, c(cpg, covariates)]
  vars <- na.omit(vars)
  counts_data <- counts_data[, match(rownames(vars), colnames(counts_data))]
  design <- model.matrix(~ ., vars)
  counts_data <- DGEList(counts = counts_data)
  counts_data <- counts_data[rowSums(counts_data$counts > 0) > filter * ncol(counts_data), ]
  counts_data <- calcNormFactors(counts_data)
  counts_data <- voom(counts_data, design)
  fit <- lmFit(counts_data, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = 2, n = Inf)
  results2 <- cbind(results, tx_annotation[match(rownames(results), tx_annotation$SYMBOL), ])
  return(results2)
}

cpgname <- rownames(beta2)
dnam <- beta2
counts <- dat1
ss <- pheno3
ss$EC <- ss$CONTROLLER_YESNO
cov1 <- c("AGE", "SEX_BIRTH", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu", "Sample_Plate")

for (i in seq_along(cpgname)) {
  twas <- do.twas(cpgname[i], cov1, ss, counts)
  twas <- correction(twas)
  tdata <- cbind(rownames(twas), twas[, "t"], twas[, "P.Value"], twas[, "SYMBOL"], nrow(ss))
  colnames(tdata) <- c("SNP", "t", "pvalue", "Gene", "num")
  saveRDS(tdata, paste0("Result//", cpgname[i], ".rds"))
}
