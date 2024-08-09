# Load required libraries
library(tidyverse)
library(openxlsx)
library(minfi)
library(sandwich)
library(lmtest)
library(progress)

# Read data
cpg1 <- read.xlsx("../../DisCPGSannotation_pval.xlsx")$CpG
beta2 <- readRDS("beta2.rds")  # Assuming this file exists
pheno3 <- read.xlsx("../../../2Pheno/NewPheno.Elite.April23.dis.xlsx")

# Process phenotype data
process_phenotype <- function(pheno) {
  pheno %>%
    filter(PC1 != "", ETHNICITY == "White", CONTROLLER_1B != "") %>%
    mutate(
      ID_idat = ID,
      CONTROLLER_YESNO = case_when(
        CONTROLLER_1B %in% c("EC_persistent", "EC_transient") ~ 1,
        CONTROLLER_1B == "Non-EC" ~ 0,
        TRUE ~ 0
      )
    )
}

pheno3 <- process_phenotype(pheno3)

# Write population info
pop1 <- table(pheno3$CONTROLLER_YESNO) %>% as.data.frame()
write.xlsx(pop1, "popinfo.xlsx")

# Align data
idx <- intersect(colnames(beta2), pheno3$ID_idat)
beta3 <- beta2[, idx]
pheno4 <- pheno3 %>% filter(ID_idat %in% idx)

# Ensure alignment
beta4 <- beta3[, match(pheno4$ID_idat, colnames(beta3))]
stopifnot(all(pheno4$ID_idat == colnames(beta4)))

# Prepare data for analysis
mVals.T0 <- t(beta4)
pheno.T0 <- pheno4
TrImmRes.T0 <- pheno4

genes <- colnames(mVals.T0)
df <- as.data.frame(mVals.T0)

# Define correlation test function
cor_test <- function(gene) {
  sub <- df[, gene, drop = FALSE]
  colnames(sub) <- 'methylation'
  test.df <- cbind(pheno.T0, sub)
  
  cyto1 <- TrImmRes.T0$CONTROLLER_YESNO
  
  if (sum(is.na(sub)) > 50) {
    return(rep(NA, 4))
  }
  
  tryCatch({
    ML <- glm(as.factor(cyto1) ~ methylation, data = test.df, family = "binomial")
    cf <- coeftest(ML, vcov = vcovHC(ML, type = "HC0"))
    result <- cf["methylation", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
    
    if (length(result) != 4) {
      result <- rep(NA, 4)
    }
    
    result
  },
  error = function(e) {
    message("Error in gene: ", gene)
    message(e)
    rep(NA, 4)
  })
}

# Run analysis
pb <- progress_bar$new(total = length(genes))
Result <- map_dfr(genes, function(gene) {
  pb$tick()
  res <- cor_test(gene)
  tibble(
    CpGsite = gene,
    Estimate = res[1],
    StdError = res[2],
    `z-score` = res[3],
    pval = res[4]
  )
})

# Process results
cor_info <- Result %>%
  mutate(FDR = p.adjust(pval, method = "fdr"))

# Print and save results
print(head(cor_info))
saveRDS(cor_info, "cor_info.rds", compress = "xz")