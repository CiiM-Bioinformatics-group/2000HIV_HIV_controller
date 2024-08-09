# Load required libraries
library(tidyverse)
library(matrixStats)
library(sandwich)
library(lmtest)

# Read data
beta2 <- readRDS("/vol/projects/CIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_newElite/1qcwithFamily/2000HIV.Mvalue.rds")
pheno3 <- readRDS("/vol/projects/CIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_newElite/2Pheno/NewPheno.Elite.April23.dis.rds")
pheno5 <- read.csv("../../PCA/covid_groups_selected_variables_may.csv")

# Remove duplicates
dup <- which(colnames(beta2) %in% c("OLV282","OLV283","OLV284"))[1:3]
beta2 <- beta2[, -dup]

# Process phenotype data
pheno3$ID_idat <- pheno3$ID
pheno3 <- pheno3 %>% 
  filter(ETHNICITY == "White", CONTROLLER_1B != "") %>%
  mutate(CONTROLLER_YESNO = case_when(
    CONTROLLER_1B %in% c("EC_persistent", "EC_transient") ~ 1,
    CONTROLLER_1B == "Non-EC" ~ 0,
    TRUE ~ 0
  ))

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
  list(probes = probes, 
       log = data.frame(
         initial_NAs = rowSums(is.na(probes)),
         removed_lower = rowSums(maskL),
         removed_upper = rowSums(maskU),
         N_for_probe = rowSums(!is.na(probes))
       ))
}

# Remove outliers
OutlierResults <- removeOutliers(beta3)
METH.2 <- OutlierResults$probes
save(OutlierResults$log, file = "Outlier_log.1.Rdata")

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
    ML <- glm(as.factor(CONTROLLER_YESNO) ~ methylation + AGE + as.factor(SEX_BIRTH) + 
                CD8T + CD4T + NK + Bcell + Mono + Neu + as.factor(Sample_Plate), 
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
Result <- vapply(genes, cor_test, numeric(4))

# Process results
cor_info <- as.data.frame(t(Result))
colnames(cor_info) <- c("Estimate", "StdError", "z-score", "pval")
cor_info$CpGsite <- rownames(cor_info)
cor_info$FDR <- p.adjust(cor_info$pval, method = "fdr")

# Save results
saveRDS(cor_info, "cor_info.rds", compress = "xz")

