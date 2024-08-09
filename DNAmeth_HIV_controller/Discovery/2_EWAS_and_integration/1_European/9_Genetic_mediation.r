# Load necessary libraries
library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(wateRmelon)
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
library(mediation)

# Function to extract mediation summary
extract_mediation_summary <- function(x) {
  clp <- 100 * x$conf.level
  isLinear.y <- inherits(x$model.y, c("lm", "rq", "glm", "survreg")) && 
    (inherits(x$model.y, "glm") && x$model.y$family$family == "gaussian" && x$model.y$family$link == "identity") || 
    (inherits(x$model.y, "survreg") && x$model.y$dist == "gaussian")
  
  smat <- if (!x$INT && isLinear.y) {
    rbind(
      c(x$d1, x$d1.ci, x$d1.p),
      c(x$z0, x$z0.ci, x$z0.p),
      c(x$tau.coef, x$tau.ci, x$tau.p),
      c(x$n0, x$n0.ci, x$n0.p)
    )
  } else {
    rbind(
      c(x$d0, x$d0.ci, x$d0.p),
      c(x$d1, x$d1.ci, x$d1.p),
      c(x$z0, x$z0.ci, x$z0.p),
      c(x$z1, x$z1.ci, x$z1.p),
      c(x$tau.coef, x$tau.ci, x$tau.p),
      c(x$n0, x$n0.ci, x$n0.p),
      c(x$n1, x$n1.ci, x$n1.p),
      c(x$d.avg, x$d.avg.ci, x$d.avg.p),
      c(x$z.avg, x$z.avg.ci, x$z.avg.p),
      c(x$n.avg, x$n.avg.ci, x$n.avg.p)
    )
  }
  
  rownames(smat) <- if (!x$INT && isLinear.y) {
    c("ACME", "ADE", "Total Effect", "Prop. Mediated")
  } else {
    c("ACME (control)", "ACME (treated)", "ADE (control)", "ADE (treated)", "Total Effect",
      "Prop. Mediated (control)", "Prop. Mediated (treated)", "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
  }
  
  colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), paste(clp, "% CI Upper", sep = ""), "p-value")
  smat
}

# Read data
beta2 <- readRDS("CPGpassedMR.rds")

pheno3 <- readRDS("/vol/projects/CIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_newElite/2Pheno/NewPheno.Elite.April23.dis.rds")
pheno3 <- dplyr::filter(pheno3, ETHNICITY == "White" & CONTROLLER_1B != "")
pheno3$CONTROLLER_YESNO <- gsub(c("EC_persistent" = 1, "EC_transient" = 1, "Non-EC" = 0), pheno3$CONTROLLER_1B)
pheno3$CONTROLLER_YESNO[is.na(pheno3$CONTROLLER_YESNO)] <- 0

pheno5 <- read.csv("../../../PCA/covid_groups_selected_variables_may.csv")
newdata <- do.call(rbind, lapply(pheno4$ID_idat, function(id) {
  id_match <- which(pheno5$ID == id)
  if (length(id_match) > 0) {
    cbind(pheno4[pheno4$ID_idat == id, ], pheno5[id_match, ])
  } else {
    cbind(pheno4[pheno4$ID_idat == id, ], pheno5[1, ])
  }
}))
pheno3 <- newdata

# Filter beta2 and pheno3
idx <- intersect(colnames(beta2), pheno3$ID_idat)
beta3 <- beta2[, colnames(beta2) %in% idx]
pheno4 <- pheno3[pheno3$ID_idat %in% idx, ]
beta4 <- beta3[, match(pheno4$ID_idat, colnames(beta3))]

# Process metabolite data
metab <- read.csv("/vol/projects/CIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_newElite/GWASsummary/dosages_SNPs_of_interest_Manoj_discovery_QCd.txt", sep="\t", header=T)
metab <- t(metab)
colnames(metab) <- metab[1, ]
metab <- data.frame(ID = rownames(metab), apply(metab[-1, ], 2, as.numeric))
metab <- dplyr::filter(metab, ID %in% pheno3$ID_idat)
metab <- metab[, c(1, 10, 12)]

newresult1 <- newresult2 <- newresult3 <- data.frame()

for (ikt in 1:ncol(metab)) {
  if (ikt > 1) {
    pheno3 <- pheno3[, !colnames(pheno3) %in% c("SAM", "SAM1", "SAM2")]
  }
  
  result <- do.call(rbind, lapply(pheno3$ID_idat, function(id) {
    id_match <- which(metab[, 1] == id)
    if (length(id_match) > 0) {
      cbind(pheno3[pheno3$ID_idat == id, ], metab[id_match, ikt + 1])
    } else {
      pheno3[pheno3$ID_idat == id, ] %>% mutate(SAM = NA)
    }
  }))
  
  result <- result %>%
    mutate(SAM1 = log(SAM, 2),
           SAM2 = qnorm((rank(SAM, na.last = "keep") - 0.5) / sum(!is.na(SAM))))
  
  beta3 <- beta4
  pheno3 <- result
  pheno6 <- pheno3
  
  # Remove outliers
  removeOutliers <- function(probes) {
    require(matrixStats)
    rowIQR <- rowIQRs(probes, na.rm = TRUE)
    row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = TRUE)
    maskL <- probes < row2575[, 1] - 3 * rowIQR
    maskU <- probes > row2575[, 2] + 3 * rowIQR
    initial_NAs <- rowSums(is.na(probes))
    probes[maskL] <- NA
    probes[maskU] <- NA
    Log <- data.frame(
      initial_NAs = initial_NAs,
      removed_lower = rowSums(is.na(probes)) - initial_NAs,
      removed_upper = rowSums(is.na(probes)) - initial_NAs,
      N_for_probe = rowSums(!is.na(probes))
    )
    list(probes, Log)
  }
  
  system.time(OutlierResults <- removeOutliers(beta3))  
  METH.2 <- OutlierResults[[1]]
  Log <- OutlierResults[[2]]
  save(Log, file = "Outlier_log.Rdata") # save log
  
  mVals.T0 <- t(METH.2)
  pheno.T0 <- pheno6
  TrImmRes.T0 <- pheno6
  
  genes <- unique(colnames(mVals.T0))
  df <- as.data.frame(mVals.T0)
  
  for (in1 in 1:length(genes)) {
    gene1 <- genes[in1]
    sub <- data.frame(methylation = df[, gene1])
    test.df <- cbind(pheno.T0, sub)
    
    Medition.data.frame <- test.df %>%
      mutate(
        HIC = as.numeric(CONTROLLER_YESNO),
        cpg = as.numeric(methylation),
        snp = as.numeric(SAM2),
        Age = as.numeric(AGE),
        Gender = as.numeric(SEX_BIRTH),
        CD8T = as.numeric(CD8T),
        CD4T = as.numeric(CD4T),
        NK = as.numeric(NK),
        Bcell = as.numeric(Bcell),
        Mono = as.numeric(Mono),
        Neu = as.numeric(Neu),
        Sample_Plate = as.factor(Sample_Plate)
      ) %>%
      filter(complete.cases(.))
    
    med.fit <- lm(cpg ~ snp + Age + Gender + CD8T + CD4T + NK + Bcell + Mono + Neu + Sample_Plate, data = Medition.data.frame)
    out.fit <- glm(HIC ~ cpg + snp + Age + Gender + CD8T + CD4T + NK + Bcell + Mono + Neu + Sample_Plate, data = Medition.data.frame, family = "binomial")
    
    med <- mediate(med.fit, out.fit, treat = "snp", mediator = "cpg", sims = 1000)
    tmp <- extract_mediation_summary(med)
    newresult1 <- rbind(newresult1, cbind(gene1, "SNP", tmp))
  }
}

# Save the final result
write.csv(newresult1, "Final_Mediation_Results.csv", row.names = FALSE)
