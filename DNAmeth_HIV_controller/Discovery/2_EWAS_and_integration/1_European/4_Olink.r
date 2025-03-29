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
library(data.table)
library(readxl)
library(matrixStats)
library(MASS)
library(dplyr)
library(tidyverse)
###################################################
###########                                   #####
########### 


#####
library(openxlsx)
beta2 <- readRDS("../ValidatedCPG.dis.MVal.rds")
dir.create("Result")

###########





pheno3=read.xlsx("../../2Pheno/NewPheno.Elite.April23.dis.xlsx")
pheno3$ID_idat=pheno3$ID
pheno3=dplyr::filter(pheno3,ETHNICITY=="White") #1104
pheno3=dplyr::filter(pheno3,CONTROLLER_1B!= "") #1104 ### removing which are not in either of elite or non-elite
EuAgeHIV=pheno3
EuAgeHIV$CONTROLLER_YESNO=EuAgeHIV$CONTROLLER_1B
table(EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO=gsub("EC_persistent",1,EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO=gsub("EC_transient",1,EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO=gsub("Non-EC",0,EuAgeHIV$CONTROLLER_YESNO)
table(EuAgeHIV$CONTROLLER_YESNO)
EuAgeHIV$CONTROLLER_YESNO[is.na(EuAgeHIV$CONTROLLER_YESNO)] = 0
table(EuAgeHIV$CONTROLLER_YESNO)


pheno3=EuAgeHIV





library(openxlsx)

metab=read.csv("../../../Common_Data/updated_OLINK_Explore3072_1910samples_2367proteins_after_bridging_normalization_QCed_21nov2022.tsv", sep = "\t")



#####################
# Assuming metab, pheno3, and beta2 are already defined

# Find the common elements between the first column of metab and pheno3$ID_idat
comm1 = intersect(metab[,1], pheno3$ID_idat)

# Find the common elements between the result and row names of beta4
comm = intersect(comm1, colnames(beta2))

# Print the common elements
print(comm)


###########

metab=metab[which(metab[,1]%in%comm),]

pheno3=pheno3[which(pheno3$ID%in%comm),]

result=data.frame()


for( i in 1:length(comm))
{
  print(i)
  
  id=which(metab[,1]==comm[i])
  if(length(id)>0)
  {
    result=rbind(result,metab[id,])
  }
}



sum(result[,1]==comm)


result2=data.frame()
for( i in 1:length(comm))
{
  print(i)
  
  id=which(pheno3$ID_idat==comm[i])
  if(length(id)>0)
  {
    result2=rbind(result2,pheno3[id,])
  }
}



sum(result2$ID_idat==comm)


metab=result
pheno3=result2

##############

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
sum(pheno3$ID == metab[,1])

sum(colnames(beta4) == metab[,1])


##################

####                                               #####
########################################################
pheno6 <- pheno3


removeOutliers<-function(probes){
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR 
  maskU <- probes > row2575[,2] + 3 * rowIQR 
  initial_NAs<-rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe<-rowSums(!is.na(probes))
  Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  return(list(probes, Log))
}

#Remove outliers from METH (methylation data where probes are rows and samples are columns)
system.time(OutlierResults<-removeOutliers(beta3))  
METH.2<-OutlierResults[[1]]
Log<-OutlierResults[[2]]
#rm(beta2)
#rm(beta3)
save(Log,file="Outlier_log.Rdata") #save log

#62225
#######
mVals.T0=t(METH.2)
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

gene= "cg04784635"

gene=genes


########################













######################
#head(metab)
colnames(metab)[1]="Row.names"






ikt=1


#metab=metab[,c(1,1230:dim(metab)[2])]

for( ikt in 1:dim(metab)[2])
{
  
  print(ikt )
  
  result=result2 ### THIS IS MATCH PHENO
  
  id1=ikt+1 ##since metabolites start from 2nd position
  
  result$SAM=as.numeric(metab[,id1])

  result$SAM1=log(result$SAM, 2)
  #shapiro.test(as.numeric(result$SAM1)) 
  
  #he data is normal if the p-value is above 0.05 
  #p-value=0.9
  #https://www.r-bloggers.com/2019/08/shapiro-wilk-test-for-normality-in-r/
  
  
  ##data normalisation only if p-value <0.05. since p value  not less than 0.05 we performeed
  # the inverse rank transform
  result$SAM2 <- qnorm ((rank(result$SAM, na.last="keep") - 0.5) / sum(!is.na(result$SAM)))
  head(result)
  
  TrImmRes.T0= result
  
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
    
    cyto1 <- TrImmRes.T0$SAM2
    
    print(ind)
    
    #  test.df1=test.df[-ind,]
    # cyto2 =cyto1[-ind] 
    
    
    #
    
    
    bad <- as.numeric(rep(NA, 4))
    names(bad)<- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
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
          
          
          ML=  rlm(as.numeric(cyto1) ~ methylation +AGE+ as.factor(SEX_BIRTH)+ as.numeric(CD8T) + as.numeric(CD4T)+ as.numeric(NK) + as.numeric(NK)+ as.numeric(Bcell)+ as.numeric(Mono)+ as.numeric(Neu)+as.factor(Sample_Plate), data=test.df)
          
          
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
  cor_info$Protein=colnames(metab)[id1]
 # print(head(cor_info))
  output=paste0("Result//cor_info.protein",ikt,".rds")
  
  saveRDS(cor_info,output,compress="xz")
  
  
  cor_info<- data.frame()
  gc()
}









