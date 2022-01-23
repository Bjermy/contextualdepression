#Load libraries
library(ggplot2)
library(data.table)

#Read in all 'main effect' results
MDDPRS <- read.csv('AllDefinitionsTable_MDDPRS_withselfreportbirth.csv')
MDDPRS <- MDDPRS[MDDPRS$Pt=='Pt < 0.05',-1]

BPDPRS <- read.csv('AllDefinitionsTable_BPDPRS_withselfreportbirth.csv')
BPDPRS <- BPDPRS[BPDPRS$Pt=='Pt < 0.01',-1] 

EnvMDDexc <- read.csv('./Heterogeneous/EnvironmentalAssociations_Heterogeneous_FullSample.csv')
EnvMDDexc <- EnvMDDexc[EnvMDDexc$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'),-1] 

EnvMDDthenDep <- read.csv('./Chronic_Disease/EnvironmentalAssociationsMDDthenDep_FullSample.csv')
EnvMDDthenDep <- EnvMDDthenDep[EnvMDDthenDep$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'),-1] 

EnvPPD <- read.csv('./Postpartum/EnvironmentalAssociationsPostPartum_FullSample_withselfreportbirth.csv')
EnvPPD <- EnvPPD[EnvPPD$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'),-1] 

#Subset to p-values only
allpvals <- c(MDDPRS$p, BPDPRS$p, EnvMDDexc$p, EnvMDDthenDep$p, EnvPPD$p)

## multiple testing correction based on:
## http://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/FDR

allpvals
sorted.pvalue<-sort(allpvals) 
sorted.pvalue 
j.alpha <- (1:30)*(.05/30) # input total number of pvalues
diff <- sorted.pvalue-j.alpha 
neg.diff <- diff[diff<0] 
pos.diff <- neg.diff[length(neg.diff)] 
index <- diff==pos.diff 
p.cutoff <-sorted.pvalue[index] 
p.cutoff # this value is the critical pvalue. Values below should be considered significant
p.sig <- allpvals[allpvals <= p.cutoff] 
p.sig # which pvalues are significant based on the critical pvalue

p.sig <- allpvals[allpvals > p.cutoff] 

#Histogram of p-values
ggplot(as.data.frame(sorted.pvalue), aes(x=sorted.pvalue)) + geom_histogram()

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

#Read in all 'case-case' results
MDDPRS_CC <- read.csv('AllDefinitionsTable_CaseCase_MDDPRS.csv')
MDDPRS_CC <- MDDPRS_CC[MDDPRS_CC$Pt=='Pt < 0.05',-1]

BPDPRS_CC <- read.csv('AllDefinitionsTable_CaseCase_BPDPRS.csv')
BPDPRS_CC$Pt <- factor(BPDPRS_CC$Pt, levels = c("BIPO02_Pt_0.01", "BIPO02_Pt_1"), 
                    labels = c("Pt < 0.01", "Pt < 1"))
BPDPRS_CC <- BPDPRS_CC[BPDPRS_CC$Pt=='Pt < 0.01' & BPDPRS_CC$Pheno!='PostPartumFemaleOnly',-1] 

EnvMDDthenDep_CC <- read.csv('./CaseCase/EnvironmentalAssociationsMDDthenDep_CaseCase_FullSample.csv')
EnvMDDthenDep_CC <- EnvMDDthenDep_CC[EnvMDDthenDep_CC$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'),-1] 

EnvPPD_CC <- read.csv('./CaseCase/EnvironmentalAssociationsPostPartum_CaseCase_FullSample.csv')
EnvPPD_CC <- EnvPPD_CC[EnvPPD_CC$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'),-1] 

#Subset to p-values only
allpvals_CC <- c(MDDPRS_CC$p, BPDPRS_CC$p, EnvMDDthenDep_CC$p, EnvPPD_CC$p)

## multiple testing correction based on:
## http://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/FDR

allpvals_CC
sorted.pvalue<-sort(allpvals_CC) 
sorted.pvalue 
j.alpha <- (1:20)*(.05/20) # input total number of pvalues
diff <- sorted.pvalue-j.alpha 
neg.diff <- diff[diff<0] 
pos.diff <- neg.diff[length(neg.diff)] 
index <- diff==pos.diff 
p.cutoff <-sorted.pvalue[index] 
p.cutoff # this value is the critical pvalue. Values below should be considered significant
p.sig <- allpvals_CC[allpvals_CC <= p.cutoff] 
p.sig # which pvalues are significant based on the critical pvalue

#Histogram of p-values 
ggplot(as.data.frame(sorted.pvalue), aes(x=sorted.pvalue)) + geom_histogram()

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

#Read in all 'case-case' results when the heterogeneous definition is subset to female only
MDDPRS_CC <- read.csv('./CaseCase/PostPartum_CaseCase_MDDPRSresults_FemaleOnly')
MDDPRS_CC$Pt <- factor(MDDPRS_CC$Pt, levels = c("DEPR06_Pt_0.05", "DEPR06_Pt_1"), 
                       labels = c("Pt < 0.05", "Pt < 1"))
MDDPRS_CC <- MDDPRS_CC[MDDPRS_CC$Pt=='Pt < 0.05',-1]

BPDPRS_CC <- read.csv('./CaseCase/PostPartum_CaseCase_BPDPRSresults_FemaleOnly')
BPDPRS_CC$Pt <- factor(BPDPRS_CC$Pt, levels = c("BIPO02_Pt_0.01", "BIPO02_Pt_1"), 
                       labels = c("Pt < 0.01", "Pt < 1"))
BPDPRS_CC <- BPDPRS_CC[BPDPRS_CC$Pt=='Pt < 0.01',-1] 

EnvPPD_CC <- read.csv('./CaseCase/EnvironmentalAssociationsPostPartum_CaseCase_FemaleOnly_FullSample.csv')
EnvPPD_CC <- EnvPPD_CC[EnvPPD_CC$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'),-1] 

#Subset to p-values only
allpvals_CC <- c(MDDPRS_CC$p, BPDPRS_CC$p, EnvPPD_CC$p)

## multiple testing correction based on:
## http://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/FDR

allpvals_CC
sorted.pvalue<-sort(allpvals_CC) 
sorted.pvalue 
j.alpha <- (1:10)*(.05/10) # input total number of pvalues
diff <- sorted.pvalue-j.alpha 
neg.diff <- diff[diff<0] 
pos.diff <- neg.diff[length(neg.diff)] 
index <- diff==pos.diff 
p.cutoff <-sorted.pvalue[index] 
p.cutoff # this value is the critical pvalue. Values below should be considered significant
p.sig <- allpvals_CC[allpvals_CC <= p.cutoff] 
p.sig # which pvalues are significant based on the critical pvalue

#Histogram of p-values 
ggplot(as.data.frame(sorted.pvalue), aes(x=sorted.pvalue)) + geom_histogram()

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

#Read in all 'control-control' results
MDDPRS_CC_PPD <- read.csv('PostPartum_ControlControl_MDDPRSresults_withselfreport.csv')
MDDPRS_CC_PPD <- MDDPRS_CC_PPD[MDDPRS_CC_PPD$Pt=='DEPR06_Pt_0.05',-c(1,2)]

MDDPRS_CC_Disease <- read.csv('MDDwithDisorder_ControlControl_MDDPRSresults.csv')
MDDPRS_CC_Disease <- MDDPRS_CC_Disease[MDDPRS_CC_Disease$Pt=='Pt < 0.05',-1]

BPDPRS_CC_PPD <- read.csv('PostPartum_ControlControl_BPDPRSresults_withselfreport.csv')
BPDPRS_CC_PPD <- BPDPRS_CC_PPD[BPDPRS_CC_PPD$Pt=='BIPO02_Pt_0.01',-c(1,2,3)]

BPDPRS_CC_Disease <- read.csv('MDDwithDisorder_ControlControl_BPDPRSresults.csv')
BPDPRS_CC_Disease <- BPDPRS_CC_Disease[BPDPRS_CC_Disease$Pt=='Pt < 0.01',-1]

EnvMDDthenDep_CC <- read.csv('EnvironmentalAssociationsMDDthenDep_ControlControl_FullSample.csv')
EnvMDDthenDep_CC <- EnvMDDthenDep_CC[EnvMDDthenDep_CC$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'),-1] 

EnvPPD_CC <- read.csv('EnvironmentalAssociationsPostPartum_ControlControl_FullSample_withselfreport.csv')
EnvPPD_CC <- EnvPPD_CC[EnvPPD_CC$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'),-1] 

#Subset to p-values only
allpvals_CC <- c(MDDPRS_CC_PPD$p, MDDPRS_CC_Disease$p, BPDPRS_CC_PPD$p, BPDPRS_CC_Disease$p, EnvMDDthenDep_CC$p, EnvPPD_CC$p)

## multiple testing correction based on:
## http://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/FDR

allpvals_CC
sorted.pvalue<-sort(allpvals_CC) 
sorted.pvalue 
j.alpha <- (1:16)*(.05/16) # input total number of pvalues
diff <- sorted.pvalue-j.alpha 
neg.diff <- diff[diff<0] 
pos.diff <- neg.diff[length(neg.diff)] 
index <- diff==pos.diff 
p.cutoff <-sorted.pvalue[index] 
p.cutoff # this value is the critical pvalue. Values below should be considered significant
p.sig <- allpvals_CC[allpvals_CC <= p.cutoff] 
p.sig # which pvalues are significant based on the critical pvalue

p.sig <- allpvals_CC[allpvals_CC > p.cutoff] 

#Histogram of p-values 
ggplot(as.data.frame(sorted.pvalue), aes(x=sorted.pvalue)) + geom_histogram()
