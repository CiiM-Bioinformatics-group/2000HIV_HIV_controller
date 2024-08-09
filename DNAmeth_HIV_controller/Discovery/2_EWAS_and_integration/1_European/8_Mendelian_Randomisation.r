# Load necessary libraries
library(TwoSampleMR)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(openxlsx)

# Read region information from an Excel file
reg <- read.xlsx("../Regioninfo.xlsx")
reg <- reg[2,] # Select the second row of the region information

# Read the correlation results
result <- readRDS("../cor_info.rds")

# Load annotation data
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(qqman)
library(ggplot2)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotdata <- ann850k[c("chr", "pos", "UCSC_RefGene_Name", "Relation_to_Island")]

# Merge result data with annotation data
Ml <- data.frame(merge(result, annotdata, by="row.names"))

# Prepare data for manhattan or qqplot
Ml$BP <- as.integer(Ml$pos)
Ml$SNP <- as.character(Ml$Row.names)
Ml$CHR <- as.numeric(gsub("[^0-9]", "", (Ml$chr)))
Ml$P <- Ml$pval

# Clean the region information
reg$chr <- gsub(" ", "", reg$chr)
reg1 <- reg

# Filter methylation data based on the region of interest
Ml1 <- dplyr::filter(Ml, chr == reg1$chr)
Ml1 <- dplyr::filter(Ml1, pos >= reg1$start - 500000)
Ml1 <- dplyr::filter(Ml1, pos <= reg1$end + 500000)
Ml2 <- dplyr::filter(Ml1, FDR < 0.05)
Ml2

# Load GWAS and mQTL data
gwas1 <- readRDS("gwasinfo2.rds")
mQTLinfo21 <- readRDS("mQTLinfo2.rds")

# Filter mQTL data based on the markers present in Ml2
mQTLinfo21 <- dplyr::filter(mQTLinfo21, marker %in% Ml2$CpGsite) %>% na.omit()
intersect(Ml2$SNP, mQTLinfo21$marker)

# Further process GWAS data
gwas2 <- dplyr::filter(gwas1, trait == "HIV 1 control") %>%
  mutate(n = n[1]) %>%
  mutate_at(c('beta', 'se'), as.numeric)
gwas2 <- na.omit(gwas2)

# Filter mQTL data based on GWAS SNPs
mQTLinfo21 <- dplyr::filter(mQTLinfo21, SNP %in% gwas2$SNP)

# Initialize data frames to store results
mr_mqtl_gwas_res <- data.frame()
mr_gwas_mqtl_res <- data.frame()

# Loop through each marker in mQTL data and perform MR analysis
for (ik in 1:length(mQTLinfo21$marker)) {
  mQTLinfo2 <- mQTLinfo21[ik,] %>% mutate_at(c('beta', 'se', 'p'), as.numeric)
  
  # Format data for MR analysis
  exp_mqtl <- format_data(mQTLinfo2, type = "exposure", effect_allele_col = 'a1', other_allele_col = 'a2',
                          eaf_col = 'eur', pval_col = 'p', samplesize_col = 'n')
  out_mqtl <- format_data(mQTLinfo2, type = "outcome", effect_allele_col = 'a1', other_allele_col = 'a2',
                          eaf_col = 'eur', pval_col = 'p', samplesize_col = 'n')
  
  exp_gwas <- format_data(gwas2, type = "exposure", effect_allele_col = 'a1', other_allele_col = 'a2',
                          eaf_col = 'eur', pval_col = 'P', samplesize_col = 'n')
  out_gwas <- format_data(gwas2, type = "outcome", effect_allele_col = 'a1', other_allele_col = 'a2',
                          eaf_col = 'eur', pval_col = 'P', samplesize_col = 'n')
  
  # Harmonize data for MR analysis
  dat_mqtl_gwas <- harmonise_data(exp_mqtl, out_gwas)
  dat_gwas_mqtl <- harmonise_data(exp_gwas, out_mqtl)
  
  # Perform MR analysis
  mr_mqtl_gwas <- mr_singlesnp(dat_mqtl_gwas)
  mr_gwas_mqtl <- mr_singlesnp(dat_gwas_mqtl)
  
  # Add additional information to the results
  mr_mqtl_gwas$CPG <- mQTLinfo2$marker
  mr_mqtl_gwas$exposure <- "Methylation"
  mr_mqtl_gwas$outcome <- "Genetic"
  mr_mqtl_gwas$gene <- mQTLinfo2$hgnc
  mr_mqtl_gwas$study <- mQTLinfo2$study
  mr_mqtl_gwas$pmid <- mQTLinfo2$pmid
  mr_mqtl_gwas$ancestry <- mQTLinfo2$ancestry
  mr_mqtl_gwas$tissue <- mQTLinfo2$tissue
  
  mr_gwas_mqtl$CPG <- mQTLinfo2$marker
  mr_gwas_mqtl$exposure <- "Genetic"
  mr_gwas_mqtl$outcome <- "Methylation"
  mr_gwas_mqtl$gene <- mQTLinfo2$hgnc
  mr_gwas_mqtl$study <- mQTLinfo2$study
  mr_gwas_mqtl$pmid <- mQTLinfo2$pmid
  mr_gwas_mqtl$ancestry <- mQTLinfo2$ancestry
  mr_gwas_mqtl$tissue <- mQTLinfo2$tissue
  
  # Append results to the main result data frames
  mr_mqtl_gwas_res <- rbind(mr_mqtl_gwas_res, mr_mqtl_gwas[1,])
  mr_gwas_mqtl_res <- rbind(mr_gwas_mqtl_res, mr_gwas_mqtl[1,])
  
  # Print progress
  print(ik)
  print(mr_mqtl_gwas)
  print(mr_gwas_mqtl)
  
  # Clear variables to save memory
  exp_mqtl <- ""
  out_mqtl <- ""
  exp_gwas <- ""
  out_gwas <- ""
}

# Select relevant columns from the results
mr_mqtl_gwas_res <- mr_mqtl_gwas_res %>% 
  select(exposure, outcome, CPG, SNP, gene, colnames(mr_mqtl_gwas_res)[6:9], tissue, study, pmid, ancestry)
mr_gwas_mqtl_res <- mr_gwas_mqtl_res %>% 
  select(exposure, outcome, CPG, SNP, colnames(mr_gwas_mqtl_res)[6:9], tissue, study, pmid, ancestry)

# Save results to RDS files
saveRDS(mr_mqtl_gwas_res, "mr_mqtl_gwas_res.rds")
saveRDS(mr_gwas_mqtl_res, "mr_gwas_mqtl_res.rds")

# Load and filter correlation results based on significant CpG sites
hc <- readRDS("../cor_info.rds")
hc <- dplyr::filter(hc, CpGsite %in% mr_mqtl_gwas_res$CPG)

# Rename column for merging
colnames(hc)[5] <- "CPG"
hc1 <- merge(hc, mr_mqtl_gwas_res, by="CPG") %>% na.omit()

# Save merged results to RDS and Excel files
saveRDS(hc1, "MRresult.rds")
library(openxlsx)
hc1 <- hc1 %>% arrange(FDR)
write.xlsx(hc1, "MRresult.xlsx")

# Further filter results based on estimate direction
hc2 <- hc1
hc2$Estimate1 <- ifelse(hc2$Estimate > 0, "Up", "Down")
hc2$Estimate2 <- ifelse(hc2$b > 0, "Up", "Down")
id <- which(hc2$Estimate1 == hc2$Estimate2)
hc1 <- hc1[id,]
write.xlsx(hc1, "MRresult2.xlsx")
