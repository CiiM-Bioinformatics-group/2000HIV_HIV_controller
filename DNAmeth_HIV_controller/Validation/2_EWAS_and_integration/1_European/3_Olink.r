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

beta2 <- readRDS("../../../Discovery/European/ValidatedCPG.val.MVal.rds")
dir.create("Result")

pheno3 <- read.xlsx("../../2Pheno/NewPheno.Elite.April23.val.xlsx")
pheno3$ID_idat <- pheno3$ID
pheno3 <- dplyr::filter(pheno3, ETHNICITY == "White", CONTROLLER_1B != "")
EuAgeHIV <- pheno3
EuAgeHIV$CONTROLLER_YESNO <- gsub("EC_persistent|EC_transient", 1, EuAgeHIV$CONTROLLER_1B)
EuAgeHIV$CONTROLLER_YESNO <- gsub("Non-EC", 0, EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO[is.na(EuAgeHIV$CONTROLLER_YESNO)] <- 0
pheno3 <- EuAgeHIV

metab <- read.csv("../../../Common_Data/updated_OLINK_Explore3072_1910samples_2367proteins_after_bridging_normalization_QCed_21nov2022.tsv", sep = "\t")

comm1 <- intersect(metab[, 1], pheno3$ID_idat)
comm <- intersect(comm1, colnames(beta2))

metab <- metab[metab[, 1] %in% comm, ]
pheno3 <- pheno3[pheno3$ID %in% comm, ]

result <- data.frame()
for (i in seq_along(comm)) {
  id <- which(metab[, 1] == comm[i])
  if (length(id) > 0) {
    result <- rbind(result, metab[id, ])
  }
}

result2 <- data.frame()
for (i in seq_along(comm)) {
  id <- which(pheno3$ID_idat == comm[i])
  if (length(id) > 0) {
    result2 <- rbind(result2, pheno3[id, ])
  }
}

metab <- result
pheno3 <- result2

idx <- intersect(colnames(beta2), pheno3$ID_idat)
beta3 <- beta2[, which(colnames(beta2) %in% idx)]
pheno4 <- pheno3[pheno3$ID_idat %in% idx, ]
beta4 <- beta3[, match(pheno4$ID_idat, colnames(beta3))]

removeOutliers <- function(probes) {
  require(matrixStats)
  if (nrow(probes) < ncol(probes)) warning("Expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = TRUE)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = TRUE)
  maskL <- probes < row2575[, 1] - 3 * rowIQR
  maskU <- probes > row2575[, 2] + 3 * rowIQR
  initial_NAs <- rowSums(is.na(probes))
  probes[maskL] <- NA
  probes[maskU] <- NA
  removed_lower <- rowSums(is.na(probes)) - initial_NAs
  removed_upper <- rowSums(is.na(probes)) - removed_lower - initial_NAs
  N_for_probe <- rowSums(!is.na(probes))
  Log <- data.frame(initial_NAs, removed_lower, removed_upper, N_for_probe)
  return(list(probes, Log))
}

system.time(OutlierResults <- removeOutliers(beta3))
METH.2 <- OutlierResults[[1]]
Log <- OutlierResults[[2]]
save(Log, file = "Outlier_log.Rdata")

mVals.T0 <- t(METH.2)
pheno.T0 <- pheno4

df <- mVals.T0 %>% as.data.frame()
genes <- unique(colnames(mVals.T0))

for (ikt in seq_len(ncol(metab))) {
  result$SAM <- as.numeric(metab[, ikt + 1])
  result$SAM1 <- log(result$SAM, 2)
  result$SAM2 <- qnorm((rank(result$SAM, na.last = "keep") - 0.5) / sum(!is.na(result$SAM)))
  TrImmRes.T0 <- result
  cor_test <- function(gene) {
    sub <- df[, gene] %>% data.frame()
    ind <- which(is.na(sub) == TRUE)
    colnames(sub) <- 'methylation'
    test.df <- cbind(pheno.T0, sub)
    cyto1 <- TrImmRes.T0$SAM2
    bad <- as.numeric(rep(NA, 4))
    names(bad) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    result <- bad
    if (length(ind) <= 50) {
      tryCatch({
        ML <- rlm(as.numeric(cyto1) ~ methylation + AGE, data = test.df)
        cf <- try(coeftest(ML, vcov = vcovHC(ML, type = "HC0")))
        result <- cf["methylation", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
      }, error = function(error_message) {
        message("Error: ", error_message)
        return(NA)
      })
    }
    return(result)
  }
  
  P_value_Index <- genes
  Result <- vapply(P_value_Index, cor_test, numeric(4))
  cor_info <- as.data.frame(t(Result))
  colnames(cor_info) <- c("Estimate", "StdError", "z-score", "pval")
  cor_info$FDR <- p.adjust(cor_info$pval, method = "fdr")
  cor_info$Protein <- colnames(metab)[ikt + 1]
  output <- paste0("Result//cor_info.protein", ikt, ".rds")
  saveRDS(cor_info, output, compress = "xz")
  gc()
}