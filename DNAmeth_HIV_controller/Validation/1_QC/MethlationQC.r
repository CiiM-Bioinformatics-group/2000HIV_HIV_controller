# Import required packages
library(preprocessCore)
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
library(FlowSorted.Blood.EPIC)
library(minfiData)

# Set file paths
loc.idat <- "/vol/projects/CIIM/2000HIV/2000HIV_DNAMethylation/MET2021-279-014_Data/idats/"
loc.sheet <- "/vol/projects/CIIM/2000HIV/2000HIV_DNAMethylation/MET2021-279-014_Data/STS/"

# Read sample sheet
targets <- read.metharray.sheet(loc.sheet, pattern="*.csv", recursive = TRUE)
targets$Basename <- apply(targets, 1, function(target) {paste(target[9], target[8], sep="_")})
targets$Sample_Group <- "HIV"

# Filter targets
a <- grep("ETZ",targets$Sample_Name,value=T)
b <- grep("EMC", targets$Sample_Name, value=T)
c <- grep("OLV", targets$Sample_Name, value=T)
d <- grep("RAD", targets$Sample_Name, value=T)
e <- grep("FAM", targets$Sample_Name, value=T)
targets <- dplyr::filter(targets, Sample_Name %in% c(a))

cat("Here is the size of the discovery cohort:", dim(targets), "\n")

# Read methylation data
RG.set <- read.metharray.exp(base=file.path(loc.idat, basename=targets$Slide), targets=targets, recursive=TRUE, extended = TRUE)
colnames(RG.set) <- targets$Sample_Name
Pheno.data <- pData(RG.set)
saveRDS(Pheno.data, "Pheno.data.rds", compress = "xz")

# Get the QC report for the raw data
qcReport(RG.set, sampNames = Pheno.data$Sample_Name, 
         sampGroups = Pheno.data$Sample_Group, 
         pdf = "Pic.1.Before filtering.minfi_qcReport.pdf", 
         maxSamplesPerPage = 50, controls = c())

# Preprocess raw data
MSet.raw <- preprocessRaw(RG.set)
qc <- getQC(MSet.raw)
saveRDS(qc, "qc1.rds", compress="xz")

# QC plot function
plotQC1 <- function(qc, badSampleCutoff = 10.5) {
  meds <- (qc$mMed + qc$uMed) / 2
  whichBad <- which((meds < badSampleCutoff))
  plot(qc$mMed, qc$uMed, xlim = c(8, 14), ylim = c(8, 14), xaxt = "n", yaxt = "n",
       xlab = "Meth median intensity (log2)", ylab = "Unmeth median intensity (log2)",
       col = ifelse(1:nrow(qc) %in% whichBad, "red", "black"))
  axis(side = 1, at = c(9, 11, 13))
  axis(side = 2, at = c(9, 11, 13))
  abline(badSampleCutoff * 2, -1, lty = 2)
  if (length(whichBad) > 0) {
    text(qc$mMed[whichBad], qc$uMed[whichBad] - 0.25, labels = rownames(qc)[whichBad], col = "red")
    legend("topleft", legend = c("good", "bad, with sample index"), pch = 1, col = c("black", "red"), bty = "n")
  }
}

# Generate QC plots
jpeg("Pic.2.Before filtering.minfi_qcplot.jpeg", width = 5, height = 5, units = 'in', res = 300)
plotQC1(qc)
dev.off()

# Control probes plot
pdf("Pic 3.Control probes plot qcReport.pdf")
controlStripPlot(RG.set, controls="BISULFITE CONVERSION II")
dev.off()

# Sample QC
det.p <- detectionP(RG.set)
failed.probes <- det.p > 0.01
failed.fraction <- colMeans(failed.probes)
failed.fraction <- data.frame(failed.fraction)
failed.fraction <- cbind(Pheno.data$Sample_Name, failed.fraction)

# Generate callrate plots and tables
ff <- as.data.frame(failed.fraction)
ff$callrate <- 1 - failed.fraction$failed.fraction
write.csv(ff, file="Table1.Callrate of all samples.csv")
ff <- ff[, c(1, 3)]
colnames(ff) <- c("SampleID", "callrate")

jpeg("Pic 4.barplot_callrate.all samples.jpeg", width = 15, height = 10, units = 'in', res = 300)
p <- ggplot(data=ff, aes(x=SampleID, y=callrate)) + 
  geom_bar(stat="identity") + 
  ylim(0, 1) + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = 'none', 
        axis.text=element_text(size=15, colour="black"), 
        axis.title=element_text(size=15, face="bold", colour="black"))
print(p)
dev.off()

# Check Gender Concordance
ratioSet <- ratioConvert(MSet.raw, what = "both", keepCN = TRUE)
gset <- mapToGenome(ratioSet)
rm(MSet.raw)
det.p <- detectionP(filtered_RG.set)
estSex <- getSex(gset, cutoff = -2)
gset <- addSex(gset, sex = estSex)

# Gender mismatch check
predictedSex = getSex(gset, cutoff = -2)$predictedSex
documentedSex <- Pheno.data$Gender
levels(documentedSex) <- c("F", "M")

gcheck <- function(docsex, predsex, samnam) {
  docsex <- as.character(docsex)
  sex.match <- identical(docsex, predsex)
  if (!sex.match) {
    sex.mismatch <- samnam[which(docsex != predsex)]
    write.csv(sex.mismatch, file="Table3.sexmismatch.csv")
  }
}

gcheck(documentedSex, predictedSex, Pheno.data$Sample_Name)

# Save gset
saveRDS(gset, "gset_sex.rds", compress = "xz")

# Plot sex concordance
jpeg("Pic 6.sex_concordance_plot.jpeg", width = 6, height = 6, units = 'in', res = 300)
plotSex(gset, id = gset$Sample_Name)
dev.off()

# Gender mismatch table
sex.mismatch = read.csv("Table3.sexmismatch.csv")
fail.detp = failed.fraction[, 1]

# Remove failed samples and gender mismatch samples
fail.sample <- unique(c(sex.mismatch$x, fail.detp))
fail.sample.index <- which(Pheno.data$Sample_Name %in% fail.sample)

if (length(fail.sample.index) == 0) {
  filtered_RG.set <- RG.set
} else {
  filtered_RG.set <- RG.set[, -fail.sample.index]
}

# Probes QC
m.set <- preprocessRaw(filtered_RG.set)
ratioSet <- ratioConvert(m.set, what = "both", keepCN = TRUE)
gm.set <- mapToGenome(ratioSet)
rm(m.set)
rm(ratioSet)

det.p <- detectionP(filtered_RG.set)
saveRDS(det.p, "AllprobesdetP.rds", compress="xz")

bad.probes <- rowMeans(det.p > 0.01) > 0.1
saveRDS(bad.probes, "bad.probe.rds", compress="xz")

bad.probe.names.detP <- rownames(det.p[bad.probes, ])
saveRDS(bad.probe.names.detP, "bad.probe.names.detP.rds", compress="xz")

jpeg("Pic 8.barplot.FilteredRGset.probes.only bad.jpeg", width = 20, height = 10, units = 'in', res = 300)
barplot(bad.probes, las=2, cex.names=0.8, ylim = c(0, 1), ylab="Mean detection p-values", cex.lab = 1)
dev.off()

# Remove cross-reactive and polymorphic probes
## "13059_2016_1066_MOESM4_ESM.csv", "13059_2016_1066_MOESM5_ESM.csv", and "13059_2016_1066_MOESM1_ESM.csv" were download from the reference (Pidsley et al., 2016).
snpprobes1 <- read.csv(file=paste("cros", "13059_2016_1066_MOESM4_ESM.csv", sep="/"), sep=",", stringsAsFactors=FALSE)
keep2 <- !(featureNames(gm.set) %in% snpprobes1$PROBE[snpprobes1$EUR_AF > 0.05])
bad.probe.names.pm1 <- rownames(gm.set[!keep2, ])

snpprobes2 <- read.csv(file=paste("cros", "13059_2016_1066_MOESM5_ESM.csv", sep="/"), sep=",", stringsAsFactors=FALSE)
keep3 <- !(featureNames(gm.set) %in% snpprobes2$PROBE[snpprobes2$EUR_AF > 0.05])
bad.probe.names.pm2 <- rownames(gm.set[!keep3, ])

reactive.probes1 <- read.csv(file=paste("cros", "13059_2016_1066_MOESM1_ESM.csv", sep="/"), sep=",", stringsAsFactors=FALSE)
reactive.probes2 <- unlist(reactive.probes1)
keep1 <- !(featureNames(gm.set) %in% reactive.probes2)
bad.probe.names.cr <- rownames(gm.set[!keep1, ])

crossprobes <- unique(c(bad.probe.names.pm1, bad.probe.names.pm2, bad.probe.names.cr))
saveRDS(crossprobes, "crossprobes.rds", compress="xz")

# Get probes on X and Y chromosomes
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
X <- featureNames(gm.set) %in% ann850k$Name[ann850k$chr %in% c("chrX")]
Y <- featureNames(gm.set) %in% ann850k$Name[ann850k$chr %in% c("chrY")]
probe.names.x <- rownames(gm.set[X, ])
probe.names.y <- rownames(gm.set[Y, ])
sexprobes <- unique(c(probe.names.x, probe.names.y))
saveRDS(sexprobes, "sexprobes.rds", compress="xz")

# Remove bad probes and probes on X and Y chromosomes
remove.probenoxy <- unique(c(bad.probe.names.detP, crossprobes))

MSetraw <- preprocessRaw(filtered_RG.set)
MSetrp <- MSetraw[!(rownames(MSetraw) %in% remove.probenoxy), ]

gmsetrp <- mapToGenome(MSetrp)
MSet.sq <- preprocessQuantile(gmsetrp)

# Final filtering
m.set.flt <- MSet.sq[!(rownames(MSet.sq) %in% sexprobes), ]

# Save beta values
saveRDS(getBeta(m.set.flt), "2000HIV.betavalues.rds", compress="xz")
# Save M-values,
saveRDS(getM(m.set.flt),"2000HIV.Mvalue.rds",compress="xz")