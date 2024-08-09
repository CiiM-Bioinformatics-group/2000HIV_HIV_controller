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

pcs=readRDS("../../2000HIV.Mvalues.pcs.rds")
head(pcs$x[,1:2])
pheno4=readRDS("../../pheno.pca.rds")

# Prepare phenotype data for correlation analysis
pheno2 <- data.frame(
  Age = as.numeric(pheno4$AGE), 
  Gender = as.factor(pheno4$SEX_BIRTH),
  CD4T = as.numeric(pheno4$CD4T),
  CD8T = as.numeric(pheno4$CD8T),
  NK = as.numeric(pheno4$NK),
  Bcell = as.numeric(pheno4$Bcell),
  Mono = as.numeric(pheno4$Mono),
  Neu = as.numeric(pheno4$Neu),
  Sample_Plate = droplevels(as.factor(pheno4$Sample_Plate))
)


pheno=pheno2

library(dplyr)
dim(pheno)
dim(pcs$x)

  n_princomps=30
  princomps <- pcs$x[, 1:n_princomps]
  corr.pvalues <- matrix(nrow = ncol(princomps),
                         ncol = ncol(pheno),
                         dimnames = list(colnames(princomps), colnames(pheno)))
  
  for (var1 in colnames(pheno)) {
    for (pc in colnames(princomps)) {
      if (type_sum(pheno[, var1]) == "fct") {
        a <- aov(princomps[, pc] ~ pheno[, var1])
        pvali <- summary(a)[[1]][["Pr(>F)"]][1]
        corr.pvalues[pc, var1] <- pvali
      } else {
        a <- rlm(pheno2[, var1] ~ princomps[, pc])
        cf <- try(coeftest(a, vcov = vcovHC(a, type = "HC0")))
        pvali <- cf["princomps[, pc]", c("Pr(>|z|)")]
        corr.pvalues[pc, var1] <- pvali
      }
    }
  }
  # for (col in colnames(corr.pvalues)){
  #     corr.pvalues[, col] <- p.adjust(corr.pvalues[, col], method = "fdr")
  # }
  library(ggplot2)
  df.res <- melt(corr.pvalues)
  colnames(df.res)=c("PC","phenotype","pval")
  df.res$fdr <- p.adjust(df.res$pval, method="fdr") # bonferroni
  
  df.res$significance <- cut(df.res$pval, breaks=c(-Inf, 1e-10, 1e-5, 0.01, 0.05, Inf), label=c("P<1e-10","P<1e-5", "P<0.01", "P<0.05", "NS"))
  PC=colnames(pcs$x[, 1:n_princomps])
  cov=colnames(pheno)
  df.res$PC<-factor(df.res$PC, levels=PC)
  df.res$phenotype<-factor(df.res$phenotype)
  
  output <- paste("Pic.2.HIV.PCA_n3w.jpeg", sep = "")
  jpeg(output, width = 12, height = 4, units = 'in', res = 300)
  p1=ggplot(data = df.res, aes(x = PC, y = phenotype, fill=significance)) +
    geom_tile() + scale_fill_manual(values=c("#E31A1C","#FC4E2A","#FEB24C","#FFF3B2","white")) +#scale_fill_manual(values=c("#E31A1C","#FC4E2A","#FEB24C","#FFF3B2","white"))
    ylab("Covariates")+xlab("Principal component of DNA methylation in HIV")+
    theme_bw()+
    theme(axis.text.x=element_text(angle = 90, hjust = 0))+
    theme(legend.position="bottom")#+ggtitle(label)
  print(p1)
  dev.off()

# Calculate variance explained by principal components
pca_var <- pcs$sdev^2
pca_var_per <- (pca_var / sum(pca_var) * 100)
pca_var_per[1:30]
