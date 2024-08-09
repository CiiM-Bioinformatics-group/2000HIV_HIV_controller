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

# Load and filter phenotype data
pheno3 <- read.xlsx("../../2Pheno/NewPheno.Elite.April23.val.xlsx")
pheno3$ID_idat <- pheno3$ID
pheno3 <- pheno3 %>%
  filter(ETHNICITY == "White", CONTROLLER_1B != "")

# Process EuAgeHIV data
EuAgeHIV <- pheno3
EuAgeHIV$CONTROLLER_YESNO <- EuAgeHIV$CONTROLLER_1B
EuAgeHIV$CONTROLLER_YESNO <- gsub("EC_persistent|EC_transient", 1, EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO <- gsub("Non-EC", 0, EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO[is.na(EuAgeHIV$CONTROLLER_YESNO)] <- 0

pheno3 <- EuAgeHIV

# Load metabolite data
metab <- read.xlsx("../../../Common_Data/2000HIV_EX_VIVO_7DAY_38meas_1754samples_afterQC_RAW.xlsx")
colnames(metab)[1] <- "Row.names"

# Find common elements
comm1 <- intersect(metab[,1], pheno3$ID_idat)
comm <- intersect(comm1, colnames(beta2))

# Filter data based on common elements
metab <- metab[metab[,1] %in% comm, ]
pheno3 <- pheno3[pheno3$ID %in% comm, ]

# Ensure matching of IDs
result <- data.frame()
for (i in seq_along(comm)) {
  id <- which(metab[,1] == comm[i])
  if (length(id) > 0) {
    result <- rbind(result, metab[id, ])
  }
}
metab <- result

result2 <- data.frame()
for (i in seq_along(comm)) {
  id <- which(pheno3$ID_idat == comm[i])
  if (length(id) > 0) {
    result2 <- rbind(result2, pheno3[id, ])
  }
}
pheno3 <- result2

# Filter beta2 data
idx <- intersect(colnames(beta2), pheno3$ID_idat)
beta3 <- beta2[, idx]
pheno4 <- pheno3[pheno3$ID_idat %in% idx, ]

# Ensure order matching
beta4 <- beta3[, match(pheno4$ID_idat, colnames(beta3))]
beta3 <- beta4
pheno3 <- pheno4

# Function to remove outliers
removeOutliers <- function(probes) {
  if (nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = TRUE)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = TRUE)
  maskL <- probes < row2575[,1] - 3 * rowIQR
  maskU <- probes > row2575[,2] + 3 * rowIQR
  initial_NAs <- rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes)) - initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes)) - removed_lower - initial_NAs
  N_for_probe <- rowSums(!is.na(probes))
  Log <- data.frame(initial_NAs, removed_lower, removed_upper, N_for_probe)
  return(list(probes, Log))
}

# Remove outliers from METH
system.time(OutlierResults <- removeOutliers(beta3))
METH.2 <- OutlierResults[[1]]
Log <- OutlierResults[[2]]
save(Log, file = "Outlier_log.Rdata")

# Prepare data for correlation analysis
mVals.T0 <- t(METH.2)
pheno.T0 <- pheno3
TrImmRes.T0 <- pheno3
df <- as.data.frame(mVals.T0)
genes <- unique(colnames(mVals.T0))

# Iterate over  for correlation analysis
for (ikt in seq_len(ncol(metab) - 1)) {
  print(ikt)
  result <- pheno3
  result$SAM <- as.numeric(metab[, ikt + 1])
  result$SAM1 <- log(result$SAM, 2)
  result$SAM2 <- qnorm((rank(result$SAM, na.last = "keep") - 0.5) / sum(!is.na(result$SAM)))
  
  TrImmRes.T0 <- result
  
  cor_test <- function(gene) {
    sub <- df[, gene, drop = FALSE]
    colnames(sub) <- 'methylation'
    test.df <- cbind(pheno.T0, sub)
    cyto1 <- TrImmRes.T0$SAM2
    
    if (sum(is.na(sub)) > 50) {
      return(rep(NA, 4))
    }
    
    tryCatch({
      ML <- rlm(as.numeric(cyto1) ~ methylation + AGE, data = test.df)
      cf <- coeftest(ML, vcov = vcovHC(ML, type = "HC0"))
      if ("try-error" %in% class(cf)) {
        return(rep(NA, 4))
      }
      return(cf["methylation", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")])
    }, error = function(e) {
      return(rep(NA, 4))
    })
  }
  
  P_value_Index <- genes
  Result <- vapply(P_value_Index, cor_test, numeric(4))
  cor_info <- as.data.frame(t(Result))
  colnames(cor_info) <- c("Estimate", "StdError", "z-score", "pval")
  cor_info$CpGsite <- rownames(cor_info)
  cor_info$FDR <- p.adjust(cor_info$pval, method = "fdr")
  cor_info$Protein <- colnames(metab)[ikt + 1]
  
  output <- paste0("Result/cor_info.protein", ikt, ".rds")
  saveRDS(cor_info, output, compress = "xz")
  
  gc()
}
