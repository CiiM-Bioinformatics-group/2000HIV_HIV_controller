library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(stringr)
library(ggplot2)
#library(wateRmelon)
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

# Load the validated CpG methylation data
beta2 <- readRDS("../ValidatedCPG.dis.MVal.rds")

# Create a directory to save the results
dir.create("Result")

# Load the phenotype data
pheno3 <- read.xlsx("../../2Pheno/NewPheno.Elite.April23.dis.xlsx")

# Add a new column to the phenotype data
pheno3$ID_idat <- pheno3$ID

# Filter the phenotype data to include only White ethnicity and those with CONTROLLER_1B value
pheno3 <- dplyr::filter(pheno3, ETHNICITY == "White")
pheno3 <- dplyr::filter(pheno3, CONTROLLER_1B != "")

# Create a new dataframe EuAgeHIV and manipulate the CONTROLLER_YESNO column
EuAgeHIV <- pheno3
EuAgeHIV$CONTROLLER_YESNO <- EuAgeHIV$CONTROLLER_1B
EuAgeHIV$CONTROLLER_YESNO <- gsub("EC_persistent", 1, EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO <- gsub("EC_transient", 1, EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO <- gsub("Non-EC", 0, EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO[is.na(EuAgeHIV$CONTROLLER_YESNO)] <- 0

# Update pheno3 with the manipulated data
pheno3 <- EuAgeHIV

# Load the metabolite data
metab <- read.xlsx("../../../Common_Data/2000HIV_EX_VIVO_24H_52meas_1760samples_afterQC_RAW.xlsx")

# Find common IDs between metab and pheno3
comm1 <- intersect(metab[,1], pheno3$ID_idat)
comm <- intersect(comm1, colnames(beta2))

# Filter metab and pheno3 based on the common IDs
metab <- metab[which(metab[,1] %in% comm),]
pheno3 <- pheno3[which(pheno3$ID %in% comm),]

# Create a new data frame result with rows from metab that match comm
result <- data.frame()
for(i in 1:length(comm)) {
  id <- which(metab[,1] == comm[i])
  if(length(id) > 0) {
    result <- rbind(result, metab[id,])
  }
}
sum(result[,1] == comm)

# Create a new data frame result2 with rows from pheno3 that match comm
result2 <- data.frame()
for(i in 1:length(comm)) {
  id <- which(pheno3$ID_idat == comm[i])
  if(length(id) > 0) {
    result2 <- rbind(result2, pheno3[id,])
  }
}
sum(result2$ID_idat == comm)

# Update metab and pheno3 with the filtered data
metab <- result
pheno3 <- result2

# Find the intersection of IDs between beta2 columns and pheno3
idx <- intersect(colnames(beta2), pheno3$ID_idat)
beta3 <- beta2[,which(colnames(beta2) %in% idx)]
pheno4 <- pheno3[which(pheno3$ID_idat %in% idx),]

# Ensure the order of columns in beta4 matches the order of IDs in pheno4
beta4 <- beta3[,match(pheno4$ID_idat, colnames(beta3))]
sum(pheno4$ID_idat == colnames(beta4))

# Update beta3 and pheno3 with the ordered data
beta3 <- beta4
pheno3 <- pheno4

# Verify the consistency of IDs across different data frames
sum(pheno3$ID_idat == colnames(beta4))
sum(pheno3$ID == colnames(beta4))
sum(pheno3$ID == metab[,1])
sum(colnames(beta4) == metab[,1])

# Copy pheno3 to pheno6 for later use
pheno6 <- pheno3

# Function to remove outliers from methylation data
removeOutliers <- function(probes) {
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
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

# Remove outliers from methylation data and save the log
OutlierResults <- removeOutliers(beta3)
METH.2 <- OutlierResults[[1]]
Log <- OutlierResults[[2]]
save(Log, file = "Outlier_log.Rdata")

# Transpose the methylation data for further analysis
mVals.T0 <- t(METH.2)
pheno.T0 <- pheno6
TrImmRes.T0 <- pheno6

# Ensure dimensions are consistent
dim(mVals.T0)
dim(pheno.T0)

# Extract unique gene names
genes <- unique(colnames(mVals.T0))
df <- as.data.frame(mVals.T0)

# Rename the first column of metab to "Row.names"
colnames(metab)[1] <- "Row.names"

# Loop over each column in metab (starting from the second column)
for(ikt in 1:(ncol(metab) - 1)) {
  print(ikt)
  
  # Update result with the current column of metab
  result <- result2
  id1 <- ikt + 1
  result$SAM <- as.numeric(metab[,id1])
  
  # Apply log transformation and inverse rank transformation to SAM
  result$SAM1 <- log(result$SAM, 2)
  result$SAM2 <- qnorm((rank(result$SAM, na.last = "keep") - 0.5) / sum(!is.na(result$SAM)))
  
  # Update TrImmRes.T0 with the transformed SAM values
  TrImmRes.T0 <- result
  
  # Define a function for correlation testing
  cor_test <- function(gene) {
    sub <- df[, gene] %>% data.frame()
    ind <- which(is.na(sub))
    colnames(sub) <- 'methylation'
    test.df <- cbind(pheno.T0, sub)
    cyto1 <- TrImmRes.T0$SAM2
    bad <- as.numeric(rep(NA, 4))
    names(bad) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    result <- bad
    
    if(length(ind) > 50) {
      return(bad)
    }
    
    if(length(ind) <= 50) {
      tryCatch({
        ML <- rlm(as.numeric(cyto1) ~ methylation + AGE + as.factor(SEX_BIRTH) + as.numeric(CD8T) + as.numeric(CD4T) + as.numeric(NK) + as.numeric(Bcell) + as.numeric(Mono) + as.numeric(Neu) + as.factor(Sample_Plate), data = test.df)
        cf <- try(coeftest(ML, vcov = vcovHC(ML, type = "HC0")))
        if(class(cf) == "try-error" || length(cf) < 4) {
          return(bad)
        }
        result <- cf["methylation", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
      }, error = function(error_message) {
        message("Error encountered: ", error_message)
        return(NA)
      })
    }
    result
  }
  
  # Initialize an empty data frame for correlation results
  cor_info <- data.frame()
  
  # Apply the correlation test to each gene
  Result <- vapply(genes, cor_test, numeric(4))
  
  # Convert the results to a data frame
  cor_info <- as.data.frame(t(Result))
  colnames(cor_info) <- c("Estimate", "StdError", "z-score", "pval")
  cor_info$CpGsite <- rownames(cor_info)
  cor_info$FDR <- p.adjust(cor_info$pval, method = "fdr")
  cor_info$Protein <- colnames(metab)[id1]
  
  # Save the correlation results
  output <- paste0("Result/cor_info.protein", ikt, ".rds")
  saveRDS(cor_info, output, compress = "xz")
  
  gc()
}
