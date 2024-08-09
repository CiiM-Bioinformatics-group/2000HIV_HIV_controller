# Load required libraries
library(tidyverse)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(qqman)
library(batchelor)
library(edgeR)
library(biomaRt)
library(org.Hs.eg.db)
library(bacon)

# Read data
beta2 <- readRDS("../../../Discovery/European/ValidatedCPG.val.MVal.rds")
result <- readRDS("../cor_info.rds")
dat1 <- readRDS("../../../Discovery/European/EWAS_RNAseq/CombineRNAseq.rds")
samp <- readRDS("../../../Discovery/European/EWAS_RNAseq/2000HIV_bulk_transcriptomics_sample_table.RDS")
pheno3 <- readRDS("../../2Pheno/NewPheno.Elite.April23.val.rds")
pheno5 <- read.csv("../../../Discovery/2Pheno//covid_groups_selected_variables_may.csv")

# Process annotation data
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotdata <- ann850k[c("chr", "pos", "UCSC_RefGene_Name", "Relation_to_Island")]

# Merge result with annotation data
Ml <- merge(result, annotdata, by.x = "row.names", by.y = "row.names") %>%
  mutate(
    BP = as.integer(pos),
    SNP = Row.names,
    CHR = as.numeric(gsub("[^0-9]", "", chr)),
    P = pval
  ) %>%
  select(CHR, BP, P, SNP) %>%
  filter(!is.na(P))

# Process RNAseq data
for (i in seq_along(samp$ID)) {
  idx <- which(samp$ID == colnames(dat1)[i])
  if (length(idx) > 0) {
    colnames(dat1)[i] <- samp$DONOR_ID[idx]
  }
}

dat11 <- readRDS("../../../Discovery/European/EWAS_RNAseq/2000HIV_bulk_transcriptomics_normalized_counts.RDS")
dat1 <- dat1 %>%
  filter(!SYMBOL %in% dat11$SYMBOL[dat11$CHR %in% c("X", "Y")]) %>%
  as.data.frame()
rownames(dat1) <- dat1$SYMBOL

# Process phenotype data
pheno3 <- pheno3 %>%
  filter(ETHNICITY == "White", CONTROLLER_1B != "") %>%
  mutate(
    ID_idat = ID,
    CONTROLLER_YESNO = case_when(
      CONTROLLER_1B %in% c("EC_persistent", "EC_transient") ~ 1,
      CONTROLLER_1B == "Non-EC" ~ 0,
      TRUE ~ 0
    )
  )

# Merge additional phenotype data
newdata <- map_dfr(pheno3$ID_idat, function(id) {
  match <- which(pheno5$ID == id)
  if (length(match) > 0) {
    cbind(pheno3[pheno3$ID_idat == id, ], pheno5[match, ])
  } else {
    cbind(pheno3[pheno3$ID_idat == id, ], NA)
  }
})

pheno3 <- newdata

# Find common IDs across datasets
common_ids <- Reduce(intersect, list(pheno3$ID_idat, colnames(beta2), colnames(dat1)))

# Filter datasets to common IDs
pheno3 <- pheno3 %>% filter(ID_idat %in% common_ids)
beta3 <- beta2[, common_ids]
dat3 <- dat1[, common_ids]

# Ensure alignment
stopifnot(all(pheno3$ID_idat == colnames(beta3)))
stopifnot(all(pheno3$ID_idat == colnames(dat3)))

# Setup for EQTM analysis
annotation_filename <- "../../../Discovery/European/EWAS_RNAseq/ID2SYMBOL_gencode_v27_transcript.txt"
tx_annotation <- read.delim(file.path(annotation_filename),
                            header = FALSE,
                            stringsAsFactors = FALSE,
                            col.names = c("GENEID", "TXNAME", "SYMBOL", "GENETYPE"))

cov1 <- c("AGE")
dir.create("Result")

# Define correction function
correction <- function(toptable) {
  bc <- bacon(effectsizes = as.matrix(toptable$logFC), 
              standarderrors = as.matrix(toptable$logFC/toptable$t))
  toptable$logFC.cor <- es(bc)
  toptable$t.cor <- tstat(bc)
  toptable$P.Value.cor <- pval(bc)
  toptable$adj.P.Val.cor <- p.adjust(toptable$P.Value.cor, method = "BH")
  return(toptable)
}

# Define TWAS function
do.twas <- function(ss_variable = "EC", ss_covariates = cov1, ss. = ss, counts. = counts, filter = 0.8) {
  vars <- ss.[, c(ss_variable, ss_covariates)] %>% na.omit()
  rownames(vars) <- ss$ID_idat
  counts. <- counts.[, match(rownames(vars), colnames(counts.))]
  design <- model.matrix(~ ., vars)
  counts. <- DGEList(counts = counts.) %>%
    calcNormFactors() %>%
    voom(design = design)
  fit <- lmFit(counts., design) %>% eBayes()
  results <- topTable(fit, coef = 2, n = Inf)
  results2 <- left_join(results, tx_annotation, by = c("row.names" = "SYMBOL"))
  return(results2)
}

# Run EQTM analysis
cpgname <- rownames(beta3)
dnam <- beta3
counts <- dat3
ss <- pheno3 %>% mutate(EC = CONTROLLER_YESNO)

for (i in seq_along(cpgname)) {
  ss[, cpgname[i]] <- dnam[cpgname[i], ]
  twas <- do.twas(cpgname[i], cov1)
  twas <- correction(twas)
  tdata <- tibble(
    SNP = rownames(twas),
    t = twas$t,
    pvalue = twas$P.Value,
    Gene = twas$SYMBOL,
    num = nrow(ss)
  )
  output <- paste0("Result//", cpgname[i], ".rds")
  saveRDS(tdata, output)
}