# Load required libraries
library(tidyverse)
library(matrixStats)
library(sandwich)
library(lmtest)
library(knitr)
library(minfi)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(wateRmelon)
library(future)
library(gtools)
library(data.table)
library(MASS)
library(parallel)
library(R.utils)
library(readxl)
library(openxlsx)
library(dplyr)
library(reshape2)
library(ggsci)
library(scales)

# Read data
beta2 <- readRDS("/vol/projects/CIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_newElite/1qcwithFamily/2000HIV.Mvalue.rds")
pheno3 <- readRDS("/vol/projects/CIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_newElite/2Pheno/NewPheno.Elite.April23.dis.rds")
pheno5 <- read.csv("../../PCA/covid_groups_selected_variables_may.csv")

# Remove duplicates
dup <- which(colnames(beta2) %in% c("OLV282", "OLV283", "OLV284"))[1:3]
beta2 <- beta2[, -dup]

# Process phenotype data
pheno3$ID_idat <- pheno3$ID
pheno3 <- pheno3 %>% 
  filter(CONTROLLER_1B != "", PC1 != "") %>%
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

# Perform PCA
df <- beta3
n_princomps <- 30
df3 <- as.matrix(df)
ntop <- 5000
rv <- rowVars(df3)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pcs <- prcomp(t(df[select, ]))
saveRDS(pcs, "2000HIV.Mvalues.pcs.rds", compress = "xz")
saveRDS(pheno3, "pheno.pca.rds", compress = "xz")

# Calculate variance explained by principal components
pca_var <- pcs$sdev^2
pca_var_per <- (pca_var / sum(pca_var) * 100)
pca_var_per[1:30]