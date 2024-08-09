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
dir.create("Result", showWarnings = FALSE)

# Load phenotype data and filter
pheno3 <- read.xlsx("../../2Pheno/NewPheno.Elite.April23.val.xlsx")
pheno3$ID_idat <- pheno3$ID
pheno3 <- pheno3 %>% filter(ETHNICITY == "White" & CONTROLLER_1B != "")
pheno3$CONTROLLER_YESNO <- recode(pheno3$CONTROLLER_1B, "EC_persistent" = 1, "EC_transient" = 1, "Non-EC" = 0)
pheno3$CONTROLLER_YESNO[is.na(pheno3$CONTROLLER_YESNO)] <- 0

# Load metabolic data
metab <- read.xlsx("../../../Common_Data/2000HIV_EX_VIVO_24H_52meas_1760samples_afterQC_RAW.xlsx")
colnames(metab)[1] <- "Row.names"

# Find common elements
comm <- Reduce(intersect, list(metab[,1], pheno3$ID_idat, colnames(beta2)))

# Subset data
metab <- metab %>% filter(Row.names %in% comm)
pheno3 <- pheno3 %>% filter(ID %in% comm)

# Combine data
result <- metab[metab$Row.names %in% comm, ]
result2 <- pheno3[pheno3$ID_idat %in% comm, ]

metab <- result
pheno3 <- result2

# Align and filter beta values
idx <- intersect(colnames(beta2), pheno3$ID_idat)
beta3 <- beta2[, idx]
pheno4 <- pheno3 %>% filter(ID_idat %in% idx)
beta4 <- beta3[, match(pheno4$ID_idat, colnames(beta3))]

# Remove outliers
removeOutliers <- function(probes) {
  require(matrixStats)
  rowIQR <- rowIQRs(probes, na.rm = TRUE)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = TRUE)
  maskL <- probes < row2575[,1] - 3 * rowIQR
  maskU <- probes > row2575[,2] + 3 * rowIQR
  probes[maskL] <- NA
  probes[maskU] <- NA
  list(probes, data.frame(
    initial_NAs = rowSums(is.na(probes)),
    removed_lower = rowSums(maskL, na.rm = TRUE),
    removed_upper = rowSums(maskU, na.rm = TRUE),
    N_for_probe = rowSums(!is.na(probes))
  ))
}

OutlierResults <- removeOutliers(beta4)
METH.2 <- OutlierResults[[1]]
Log <- OutlierResults[[2]]
save(Log, file = "Outlier_log.Rdata")

mVals.T0 <- t(METH.2)
pheno.T0 <- pheno4
TrImmRes.T0 <- pheno4

# Main analysis loop
for (ikt in seq_along(metab)) {
  result <- result2
  id1 <- ikt + 1
  result$SAM <- as.numeric(metab[, id1])
  result$SAM2 <- qnorm((rank(result$SAM, na.last = "keep") - 0.5) / sum(!is.na(result$SAM)))
  TrImmRes.T0 <- result
  
  cor_test <- function(gene) {
    sub <- data.frame(methylation = df[, gene])
    test.df <- cbind(pheno.T0, sub)
    cyto1 <- TrImmRes.T0$SAM2
    result <- tryCatch({
      ML <- rlm(as.numeric(cyto1) ~ methylation + AGE, data = test.df)
      cf <- coeftest(ML, vcov = vcovHC(ML, type = "HC0"))
      if (class(cf) == "try-error") rep(NA, 4) else cf["methylation", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
    }, error = function(e) rep(NA, 4))
    result
  }
  
  P_value_Index <- unique(colnames(mVals.T0))
  df <- as.data.frame(mVals.T0)
  Result <- vapply(P_value_Index, cor_test, numeric(4))
  cor_info <- as.data.frame(t(Result))
  colnames(cor_info) <- c("Estimate", "StdError", "z-score", "pval")
  cor_info$CpGsite <- rownames(cor_info)
  cor_info$FDR <- p.adjust(cor_info$pval, method = "fdr")
  cor_info$Protein <- colnames(metab)[id1]
  saveRDS(cor_info, file = sprintf("Result/cor_info.protein%03d.rds", ikt), compress = "xz")
  gc()
}
