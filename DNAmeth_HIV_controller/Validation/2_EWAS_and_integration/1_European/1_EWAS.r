# Load required libraries
library(tidyverse)
library(matrixStats)
library(sandwich)
library(lmtest)

# Read data
beta2 <- readRDS("/vol/projects/CIIM/2000HIV/AnalysisDV/Validation_MethylationDataBased/1qc6/2000HIV.Mvalue.rds")
pheno3 <- readRDS("/vol/projects/CIIM/2000HIV/AnalysisDV/Validation_MethylationDataBased/4EWAS_newELite/2Pheno/NewPheno.Elite.April23.val.rds")
pheno5 <- read.csv("/vol/projects/CIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_newElite/PCA/covid_groups_selected_variables_may.csv")

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

# Align beta values and phenotype data
idx <- intersect(colnames(beta2), pheno3$ID_idat)
beta3 <- beta2[, idx]
pheno4 <- pheno3[match(idx, pheno3$ID_idat), ]

# Remove outliers function
removeOutliers <- function(probes) {
  rowIQR <- rowIQRs(probes, na.rm = TRUE)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = TRUE)
  maskL <- probes < row2575[,1] - 3 * rowIQR 
  maskU <- probes > row2575[,2] + 3 * rowIQR 
  probes[maskL | maskU] <- NA
  N_for_probe <- rowSums(!is.na(probes))
  Log <- data.frame(
    initial_NAs = rowSums(is.na(probes)),
    removed_lower = rowSums(maskL),
    removed_upper = rowSums(maskU),
    N_for_probe = N_for_probe
  )
  list(probes = probes, Log = Log)
}

# Remove outliers
OutlierResults <- removeOutliers(beta3)
METH.2 <- OutlierResults$probes
save(OutlierResults$Log, file = "Outlier_log.1.Rdata")

# Prepare data for analysis
mVals.T0 <- t(METH.2)
pheno.T0 <- TrImmRes.T0 <- pheno4
df <- as.data.frame(mVals.T0)

# Correlation test function
cor_test <- function(gene) {
  sub <- df[, gene, drop = FALSE]
  colnames(sub) <- 'methylation'
  test.df <- cbind(pheno.T0, sub)
  
  if (sum(is.na(sub)) > 50) {
    return(rep(NA, 4))
  }
  
  tryCatch({
    ML <- glm(as.factor(CONTROLLER_YESNO) ~ methylation + AGE, 
              data = test.df, family = "binomial")
    
    cf <- coeftest(ML, vcov = vcovHC(ML, type = "HC0"))
    cf["methylation", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
  }, error = function(e) {
    message("Error in gene: ", gene)
    message(e)
    rep(NA, 4)
  })
}

# Run correlation test
genes <- colnames(df)
x <- 0
Result <- vapply(genes, function(gene) {
  x <<- x + 1
  cat("Processing gene", x, "of", length(genes), "\n")
  cor_test(gene)
}, numeric(4))

# Process results
cor_info <- as.data.frame(t(Result))
colnames(cor_info) <- c("Estimate", "StdError", "z-score", "pval")
cor_info$CpGsite <- rownames(cor_info)
cor_info$FDR <- p.adjust(cor_info$pval, method = "fdr")

# Save results
saveRDS(cor_info, "cor_info.rds", compress = "xz")

