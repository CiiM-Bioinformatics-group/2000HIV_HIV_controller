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
library(data.table)
library(readxl)
library(matrixStats)
library(MASS)
library(dplyr)
library(tidyverse)
library(openxlsx)
###################################################
###########                                   #####

cpg1=read.xlsx("../../DisCPGSannotation_pval.xlsx")$CpG
cpg1

# 
# beta4=readRDS("../../Bval.dis.eur.rds")
# idrt=which(rownames(beta4)%in%cpg1)
#beta2 <-beta4[idrt,]

#saveRDS(beta2,"beta2.rds")

dim(beta2)
pheno3=read.xlsx("../../../2Pheno/NewPheno.Elite.April23.dis.xlsx")
pheno3$ID_idat=pheno3$ID
table(pheno3$ETHNICITY)


pheno3=dplyr::filter(pheno3,PC1!= "") 
pheno3=dplyr::filter(pheno3,ETHNICITY=="Asian") #1104
pheno3=dplyr::filter(pheno3,CONTROLLER_1B!= "") #1104 ### removing which are not in either of elite or non-elite
EuAgeHIV=pheno3
EuAgeHIV$CONTROLLER_YESNO=EuAgeHIV$CONTROLLER_1B
table(EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO=gsub("EC_persistent",1,EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO=gsub("EC_transient",1,EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO=gsub("Non-EC",0,EuAgeHIV$CONTROLLER_YESNO)
table(EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO[is.na(EuAgeHIV$CONTROLLER_YESNO)] = 0


pop1=as.data.frame(table(EuAgeHIV$CONTROLLER_YESNO))
library(openxlsx)
write.xlsx(pop1,"popinfo.xlsx")

pheno3=EuAgeHIV




which(pheno3$ID!=pheno3$ID_idat)
idx <- intersect(colnames(beta2), pheno3$ID_idat)
length(idx)


beta3 <- beta2[,which(colnames(beta2) %in% idx)]
pheno4=pheno3[which(pheno3$ID_idat %in% idx),]

sum(pheno4$ID_idat == colnames(beta3))

beta4 <- beta3[,match(pheno4$ID_idat, colnames(beta3))]
sum(pheno4$ID_idat == colnames(beta4))


idx2 <- match(colnames(beta3), pheno4$ID_idat)


for(i in 1:length(colnames(beta4)))
{
  if(colnames(beta4)[i]!=pheno4$ID_idat[i])
  {
    tt=paste0(i, colnames(beta4)[i],"to", pheno4$ID_idat[i])
    print(tt)
  }
}


beta3=beta4
pheno3=pheno4


table(pheno3$CONTROLLER_YESNO)

sum(pheno3$ID_idat == colnames(beta4))
sum(pheno3$ID == colnames(beta4))


####################################################################

### test example here
####
#beta3=beta4[c("cg26928153", "cg16269199"),]





#####
########################################################
pheno6 <- pheno3


# removeOutliers<-function(probes){
#   require(matrixStats)
#   if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
#   rowIQR <- rowIQRs(probes, na.rm = T)
#   row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
#   maskL <- probes < row2575[,1] - 3 * rowIQR 
#   maskU <- probes > row2575[,2] + 3 * rowIQR 
#   initial_NAs<-rowSums(is.na(probes))
#   probes[maskL] <- NA
#   removed_lower <- rowSums(is.na(probes))-initial_NAs
#   probes[maskU] <- NA
#   removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
#   N_for_probe<-rowSums(!is.na(probes))
#   Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
#   return(list(probes, Log))
# }
# 
# #Remove outliers from METH (methylation data where probes are rows and samples are columns)
# system.time(OutlierResults<-removeOutliers(beta3))  
# METH.2<-OutlierResults[[1]]
# Log<-OutlierResults[[2]]
# rm(beta2)
# rm(beta3)
# save(Log,file="Outlier_log.1.Rdata") #save log

#62225
#######
mVals.T0=t(beta3)
pheno.T0=pheno6
TrImmRes.T0=pheno6 ##62 elite controller

library("tidyr")
#https://sparkbyexamples.com/r-programming/remove-rows-with-na-in-r/

dim(mVals.T0)#[1] 1525  500
dim(pheno.T0)#1525   32
dim(mVals.T0)#1525  500   

genes <- unique(colnames(mVals.T0))#test remove [1:100]
df=mVals.T0 


df = df %>% as.data.frame()

gene= "cg26928153"

gene=genes


cor_test=function(gene)
{ 
  #print(gene)
  sub <- df[, gene] %>% data.frame()  #df=normalized gene count matrix
  ind=which(is.na(sub)==TRUE)
  colnames(sub) <- 'methylation'
  test.df <- cbind(pheno.T0, sub)
  
  # a = cor.test(test.df$gene_expression, as.numeric(test.df[,Condition_1_temp]), method = "spearman",exact = FALSE)
  # values = c(a$estimate, a$p.value)
  #return (values);
  
  cyto1 <- TrImmRes.T0$CONTROLLER_YESNO
  
  print(ind)
  
  #  test.df1=test.df[-ind,]
  # cyto2 =cyto1[-ind] 
  
  
  #
  
  bad <- as.numeric(rep(NA, 4))
  names(bad) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  result <- bad
  
  
  if(length(ind)>50)
  {
    bad <- as.numeric(rep(NA, 4))
    names(bad)<- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    result <- bad
  }  
  
  
  if(length(ind)<=50)
  {
    
    
    
    tryCatch(
      # This is what I want to do...
      {
        
        
        
        #ML=rlm(methylation ~ as.factor(cyto1) + AGE + as.factor(SEX_BIRTH) + CD8T + CD4T + NK + Bcell + Mono + Neu+as.factor(Sample_Plate)+as.factor(cov)+as.factor(iso2)+as.numeric(season_cos)+as.numeric(season_sin)+as.numeric(HIV_DURATION), data=test.df)
        #ML=rlm(methylation ~ as.factor(cyto1) +AGE + as.factor(SEX_BIRTH)  + CD8T + CD4T + NK + Bcell + Mono + Neu+
        #         as.factor(Sample_Plate)+HIV_DURATION+season_cos+season_sin, data=test.df)
        
        
        # ML=glm(as.factor(cyto1) ~ methylation +AGE + as.factor(SEX_BIRTH)  + CD8T + CD4T + NK + Bcell + Mono + Neu+
        #          as.factor(Sample_Plate)+HIV_DURATION+season_cos+season_sin, data=test.df, family = "binomial")
        # 
        
        ML=glm(as.factor(cyto1)  ~ methylation, data=test.df, family = "binomial")
        
        
        
        
        
        cf <- try(coeftest(ML, vcov=vcovHC(ML, type="HC0")))
        
        
        
        ##
        
        ##keep track
        x <<- x+1
        a <- sum(x:length(P_value_Index))
        cat("I am at ",x, " of total=",length(P_value_Index),".\n",sep="")
        cat("Left =",length(P_value_Index)-x,".\n",sep="")
        
        ####
        if (class(cf)=="try-error") {
          bad <- as.numeric(rep(NA, 4))
          names(bad)<- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
          result <- bad
        }
        else{
          result <- cf["methylation", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
        }
        
        if(length(result)<4)
        {
          bad <- as.numeric(rep(NA, 4))
          names(bad)<- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
          result <- bad
        }
        if(length(result)==4)
        {
          result <- cf["methylation", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
        }
        
        
      },
      # ... but if an error occurs, tell me what happened: 
      error=function(error_message) {
        message("This is my custom message.")
        message("And below is the error message from R:")
        message(error_message)
        return(NA)
      }
    )
  } 
  result
}




L = list()
P_value_Index=genes


cor_info<- data.frame()

x <- 0
Result=vapply(P_value_Index,cor_test,numeric(4))

cor_info=as.data.frame(t(Result))


colnames(cor_info)=c("Estimate", "StdError", "z-score", "pval")

cor_info[,5]=rownames(cor_info)#cor_info[,4]=genes
colnames(cor_info)=c("Estimate", "StdError", "z-score", "pval", "CpGsite")

cor_info$FDR=p.adjust(cor_info$pval, method = "fdr")

print(head(cor_info))
saveRDS(cor_info,"cor_info.rds",compress="xz")






