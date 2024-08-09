# Load necessary libraries
library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(stringr)
library(ggplot2)
# library(wateRmelon)
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

# Create a directory to store results
dir.create("Result")

# Read validated CpG data
beta2 <- readRDS("../ValidatedCPG.dis.MVal.rds")

# Load phenotype data
pheno3 <- read.xlsx("../../2Pheno/NewPheno.Elite.April23.dis.xlsx")
pheno3$ID_idat <- pheno3$ID
pheno3 <- dplyr::filter(pheno3, ETHNICITY == "White")
pheno3 <- dplyr::filter(pheno3, CONTROLLER_1B != "")
EuAgeHIV <- pheno3
EuAgeHIV$CONTROLLER_YESNO <- EuAgeHIV$CONTROLLER_1B
EuAgeHIV$CONTROLLER_YESNO <- gsub("EC_persistent", 1, EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO <- gsub("EC_transient", 1, EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO <- gsub("Non-EC", 0, EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO[is.na(EuAgeHIV$CONTROLLER_YESNO)] <- 0
pheno3 <- EuAgeHIV

# Load metabolite data
metab <- read.xlsx("../../../Common_Data/UpdatedFlowcytometry/2000HIV_FLOW_PER_panel123merged_QCed_untransformed(raw)data_1434samples_356vars_Nov062023.xlsx")

# Find common IDs between datasets
comm1 <- intersect(metab[, 1], pheno3$ID_idat)
comm <- intersect(comm1, colnames(beta2))

# Filter metabolite and phenotype data to keep only common IDs
metab <- metab[which(metab[, 1] %in% comm),]
pheno3 <- pheno3[which(pheno3$ID %in% comm),]

# Ensure the data order matches
result <- data.frame()
for (i in 1:length(comm)) {
  id <- which(metab[, 1] == comm[i])
  if (length(id) > 0) {
    result <- rbind(result, metab[id,])
  }
}
metab <- result

result2 <- data.frame()
for (i in 1:length(comm)) {
  id <- which(pheno3$ID_idat == comm[i])
  if (length(id) > 0) {
    result2 <- rbind(result2, pheno3[id,])
  }
}
pheno3 <- result2

# Align methylation and phenotype data
idx <- intersect(colnames(beta2), pheno3$ID_idat)
beta3 <- beta2[, which(colnames(beta2) %in% idx)]
pheno4 <- pheno3[which(pheno3$ID_idat %in% idx),]

beta4 <- beta3[, match(pheno4$ID_idat, colnames(beta3))]
pheno3 <- pheno4

# Remove outliers from methylation data
removeOutliers <- function(probes) {
  require(matrixStats)
  if (nrow(probes) < ncol(probes)) warning("Expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = TRUE)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = TRUE)
  maskL <- probes < row2575[, 1] - 3 * rowIQR
  maskU <- probes > row2575[, 2] + 3 * rowIQR
  initial_NAs <- rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes)) - initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes)) - removed_lower - initial_NAs
  N_for_probe <- rowSums(!is.na(probes))
  Log <- data.frame(initial_NAs, removed_lower, removed_upper, N_for_probe)
  return(list(probes, Log))
}

OutlierResults <- removeOutliers(beta3)
METH.2 <- OutlierResults[[1]]
Log <- OutlierResults[[2]]
save(Log, file = "Outlier_log.Rdata")

mVals.T0 <- t(METH.2)
pheno.T0 <- pheno6 <- pheno3

genes <- unique(colnames(mVals.T0))
df <- as.data.frame(mVals.T0)

# Correlation test function
cor_test <- function(gene) {
  sub <- df[, gene] %>% data.frame()
  colnames(sub) <- 'methylation'
  test.df <- cbind(pheno.T0, sub)
  cyto1 <- TrImmRes.T0$SAM2
  
  result <- tryCatch({
    ML <- rlm(as.numeric(cyto1) ~ methylation + AGE + as.factor(SEX_BIRTH) + as.numeric(CD8T) +
                as.numeric(CD4T) + as.numeric(NK) + as.numeric(Bcell) + as.numeric(Mono) +
                as.numeric(Neu) + as.factor(Sample_Plate), data = test.df)
    cf <- coeftest(ML, vcov = vcovHC(ML, type = "HC0"))
    if (class(cf) == "try-error") {
      rep(NA, 4)
    } else {
      cf["methylation", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
    }
  }, error = function(e) {
    rep(NA, 4)
  })
  
  names(result) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  return(result)
}

for (ikt in 1:ncol(metab)) {
  result <- result2
  id1 <- ikt + 1
  result$SAM <- as.numeric(metab[, id1])
  result$SAM1 <- log(result$SAM, 2)
  result$SAM2 <- qnorm((rank(result$SAM, na.last = "keep") - 0.5) / sum(!is.na(result$SAM)))
  TrImmRes.T0 <- result
  
  P_value_Index <- genes
  cor_info <- vapply(P_value_Index, cor_test, numeric(4))
  cor_info <- as.data.frame(t(cor_info))
  colnames(cor_info) <- c("Estimate", "StdError", "z-score", "pval")
  cor_info$CpGsite <- rownames(cor_info)
  cor_info$FDR <- p.adjust(cor_info$pval, method = "fdr")
  cor_info$Protein <- colnames(metab)[id1]
  
  output <- paste0("Result//cor_info.protein", ikt, ".rds")
  saveRDS(cor_info, output, compress = "xz")
  gc()
}
