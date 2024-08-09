# Load necessary libraries
library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(stringr)
library(ggplot2)
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
library(openxlsx)

# Load data
beta2 <- readRDS("../../../Discovery/European/ValidatedCPG.val.MVal.rds")
pheno3 <- read.xlsx("../../2Pheno/NewPheno.Elite.April23.val.xlsx")
metab <- read.xlsx("../../../Common_Data/UpdatedFlowcytometry/2000HIV_FLOW_PER_panel123merged_QCed_untransformed(raw)data_1434samples_356vars_Nov062023.xlsx")

# Create result directory
dir.create("Result")

# Filter phenotype data
pheno3$ID_idat <- pheno3$ID
pheno3 <- dplyr::filter(pheno3, ETHNICITY == "White", CONTROLLER_1B != "")
EuAgeHIV <- pheno3
EuAgeHIV$CONTROLLER_YESNO <- gsub("EC_persistent|EC_transient", 1, EuAgeHIV$CONTROLLER_1B)
EuAgeHIV$CONTROLLER_YESNO <- gsub("Non-EC", 0, EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO[is.na(EuAgeHIV$CONTROLLER_YESNO)] <- 0
pheno3 <- EuAgeHIV

# Find common elements
comm1 <- intersect(metab[, 1], pheno3$ID_idat)
comm <- intersect(comm1, colnames(beta2))

# Filter metabolite and phenotype data
metab <- metab[metab[, 1] %in% comm, ]
pheno3 <- pheno3[pheno3$ID_idat %in% comm, ]

# Match beta data with filtered phenotype data
idx <- intersect(colnames(beta2), pheno3$ID_idat)
beta3 <- beta2[, colnames(beta2) %in% idx]
pheno4 <- pheno3[pheno3$ID_idat %in% idx, ]
beta4 <- beta3[, match(pheno4$ID_idat, colnames(beta3))]

# Remove outliers function
removeOutliers <- function(probes) {
  require(matrixStats)
  rowIQR <- rowIQRs(probes, na.rm = TRUE)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = TRUE)
  maskL <- probes < row2575[, 1] - 3 * rowIQR
  maskU <- probes > row2575[, 2] + 3 * rowIQR
  initial_NAs <- rowSums(is.na(probes))
  probes[maskL] <- NA
  probes[maskU] <- NA
  Log <- data.frame(initial_NAs, removed_lower = rowSums(is.na(probes)) - initial_NAs, N_for_probe = rowSums(!is.na(probes)))
  return(list(probes, Log))
}

# Remove outliers from METH
system.time(OutlierResults <- removeOutliers(beta3))
METH.2 <- OutlierResults[[1]]
Log <- OutlierResults[[2]]
save(Log, file = "Outlier_log.Rdata")

# Prepare data for analysis
mVals.T0 <- t(METH.2)
pheno.T0 <- pheno4
TrImmRes.T0 <- pheno4

# Correlation test function
cor_test <- function(gene) {
  sub <- df[, gene, drop = FALSE]
  ind <- which(is.na(sub) == TRUE)
  colnames(sub) <- 'methylation'
  test.df <- cbind(pheno.T0, sub)
  
  if (length(ind) > 50) {
    return(rep(NA, 4))
  }
  
  tryCatch({
    ML <- rlm(as.numeric(TrImmRes.T0$SAM2) ~ methylation + AGE, data = test.df)
    cf <- coeftest(ML, vcov = vcovHC(ML, type = "HC0"))
    return(cf["methylation", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")])
  }, error = function(e) {
    return(rep(NA, 4))
  })
}

# Execute correlation tests
genes <- unique(colnames(mVals.T0))
Result <- vapply(genes, cor_test, numeric(4))
cor_info <- data.frame(t(Result))
colnames(cor_info) <- c("Estimate", "StdError", "z-score", "pval")
cor_info$FDR <- p.adjust(cor_info$pval, method = "fdr")
output <- paste0("Result//cor_info.protein", ikt, ".rds")
saveRDS(cor_info, output, compress = "xz")