#MDD according to the CIDI-SF with PPD removed and never suffered a disorder associated with depression. 

#Load libraries
library(data.table)
library(dplyr)
library(sqldf)
library(tidyverse)
library(lubridate)
library(ukbkings)
library(fmsb)

#Load the lifetime depression cases
CIDIDep <- fread('LifetimeDepression.csv', data.table=FALSE)
CIDIDep <- CIDIDep[,-1]

#Remove any that have endorsed childbirth as a possible reason for the depressive episode.
ChildbirthRelated <- readRDS('./PostPartum/PostPartumDepression.rds')

CIDIDep <- inner_join(CIDIDep, ChildbirthRelated, by='eid')
colnames(CIDIDep)[3] <- 'Childbirth'

CIDIDep <- subset(CIDIDep, LifetimeDepression == 1 & (Childbirth==0 | is.na(Childbirth) | Childbirth==-313)) 

#Load in GP depression - any cases of depression
DepressionGP <- fread('DepressionGP', data.table=FALSE)
DepressionGP <- DepressionGP[,-1]

#Join cases to LifetimeMDD definition
DepGPandCIDI <- full_join(CIDIDep, DepressionGP, by='eid') 

#Remove any cases of postnatal depression identified through other means.
PPD <- fread('./PostPartum/PostPartumMHQandEHR', data.table=FALSE)
PPD <- PPD[,-1]

DepGPandCIDI <- subset(DepGPandCIDI, !(DepGPandCIDI$eid %in% PPD))

#Make sure the cases have not suffered from one of the disorders considered. 

##Removal Step 1: Hospitalisations - Anti-join on existing disorder list.
Disorders <- fread('./Diseases/Disorders', data.table=FALSE)

DepGPandCIDI <- anti_join(DepGPandCIDI, Disorders, by='eid')

##Removal Step 2: Self-reported any of the disorders at nurse interview. 
CancerSR <- readRDS('SelfReportMedicalIllnessCancer.rds')

#Only concerned with diagnosis not year
CancerSR <- CancerSR[,c(1,3:26)]

DepGPandCIDI <- left_join(DepGPandCIDI, CancerSR, by='eid')

#If the individual hasn't ever reported a cancer the total NAs for each person will equal the total number of columns in which the person had the opportunity to report which is 24
DepGPandCIDI$CancerNumber <- apply(DepGPandCIDI[,c(6:29)], 1, function(x) sum(is.na(x)))

DepGPandCIDI <- subset(DepGPandCIDI, CancerNumber==24)

DepGPandCIDI <- DepGPandCIDI[,c(1:3)]

NonCancerSR <- readRDS('SelfReportMedicalIllnessNonCancer.rds')

#Only concerned with diagnosis not year
NonCancerSR <- NonCancerSR[,c(1,3:138)]

DepGPandCIDI <- left_join(DepGPandCIDI, NonCancerSR, by='eid')

DepGPandCIDI$NonCancer <- apply(DepGPandCIDI, 1, function(r) any(r %in% c(1264,1222,1223,1220,1081,1425,1491,1086,1075,1093,1261,1259,1311,1456,1154,1459,1463,1331,1453,1464,1382,1226,1462,1531,1377,1428,1082,1437,1313,1381,1461,1376,1521,1076,1083,1260,1345,1583)))

DepGPandCIDI <- subset(DepGPandCIDI, NonCancer==FALSE)

#NonCancer values left that are relevant
values <- apply(DepGPandCIDI[,c(4:139)], 2, function(x) unique(x))
values <- unlist(values)
values <- unique(values)
values <- sorted(values)

write.csv(values, 'uniqueNonCancer.csv')

#Exclusion criteria - taken from script LifetimeDepressionMHQDefinition.R
exclusion <- fread('ExclusionCriteriaApplied.csv', data.table=FALSE)
exclusion <- exclusion[,-1]

DepGPandCIDIExc <- inner_join(x=DepGPandCIDI, y=exclusion, by="eid") 
DepGPandCIDIExc <- unique(DepGPandCIDIExc[,1]) #34700 cases

write.csv(DepGPandCIDIExc, 'MHQCaseExcPPDandDisorders')

#Europeans only.
europeans <- fread('ukb18177_glanville_post_qc_id_list.fam' , data.table=FALSE)
colnames(europeans)[1] <- 'eid'

DepGPandCIDIExc <- as.data.frame(DepGPandCIDIExc)
colnames(DepGPandCIDIExc) <- 'eid'

DepGPandCIDIEur <- inner_join(DepGPandCIDIExc, europeans)

DepGPandCIDIcase <- DepGPandCIDIEur[,1] 

write.csv(DepGPandCIDIcase, 'MHQCaseExcPPDandDisorders_QC')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Regression analysis 

#MDD PRS
cases <- fread('MHQCaseExcPPDandDisorders_QC', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

controls <- fread('controlsMDDExcPPDandDisorders_QC.csv', data.table=FALSE)
controls$MDD <- 0
colnames(controls)[2] <- 'eid'
controls <- controls[,-1]

MDD <- rbind(cases, controls)

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('ukb18177_glanville_covariates.txt' , data.table=FALSE)
colnames(covariates)[2] <- 'eid'
covariates <- covariates[,-1]

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/DEPR06_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c('eid','DEPR06_Pt_0.05','DEPR06_Pt_1')]

#Merge to the postpartum depression dataset
MDD <- left_join(MDD, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
MDD$DEPR06_Pt_0.05 <- scale(MDD$DEPR06_Pt_0.05)
MDD$DEPR06_Pt_1 <- scale(MDD$DEPR06_Pt_1)

#Make batch and assessment centre factors.
MDD$assessment_centre <- factor(MDD$assessment_centre, ordered=FALSE)
MDD$batch <- factor(MDD$batch, ordered=FALSE)

scores<-c("DEPR06_Pt_0.05","DEPR06_Pt_1")

MDDresult <- data.frame()

for (j in 1:length(scores)){
  MDD$scores=MDD[,scores[j]]
  
  # create model
  # full model
  model<-glm(MDD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=MDD,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(MDD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=MDD,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  MDDresult <- rbind(MDDresult, output)
} 

MDDresult$Pheno <- 'All' 

write.csv(MDDresult, 'Heterogeneous_MDDPRSresults')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#BPD PRS

cases <- fread('MHQCaseExcPPDandDisorders_QC', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

controls <- fread('controlsMDDExcPPDandDisorders_QC.csv', data.table=FALSE)
controls$MDD <- 0
colnames(controls)[2] <- 'eid'
controls <- controls[,-1]

MDD <- rbind(cases, controls)

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/BIPO02_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c('eid','BIPO02_Pt_0.01','BIPO02_Pt_1')]

#Merge to the postpartum depression dataset
MDD <- left_join(MDD, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
MDD$BIPO02_Pt_0.01 <- scale(MDD$BIPO02_Pt_0.01)
MDD$BIPO02_Pt_1 <- scale(MDD$BIPO02_Pt_1)

#Make batch and assessment centre factors.
MDD$assessment_centre <- factor(MDD$assessment_centre, ordered=FALSE)
MDD$batch <- factor(MDD$batch, ordered=FALSE)

scores<-c("BIPO02_Pt_0.01","BIPO02_Pt_1")

MDDresult <- data.frame()

for (j in 1:length(scores)){
  MDD$scores=MDD[,scores[j]]
  
  # create model
  # full model
  model<-glm(MDD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=MDD,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(MDD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=MDD,family=binomial,na.action=na.exclude)
  # get results
  #print(summary(model)$coefficients)
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  MDDresult <- rbind(MDDresult, output)
} 

MDDresult$Pheno <- 'All' 

write.csv(MDDresult, 'Heterogeneous_BPDPRSresults')

###########################################################################################################################################################################
###########################################################################################################################################################################

#Supplementary analysis 1

#Only take female cases and controls and repeat the analysis so we can remove the sex effect from this group and make results more comparable to PPD. 

#MDD PRS

#Read in cases
cases <- fread('MHQCaseExcPPDandDisorders_QC', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

#Read in controls
controls <- fread('controlsMDDExcPPDandDisorders_QC.csv', data.table=FALSE)
controls$MDD <- 0
colnames(controls)[2] <- 'eid'
controls <- controls[,-1]

#Merge
MDD <- rbind(cases, controls)

#Read in sex
sex <- readRDS("sex.rds")
colnames(sex)[2] <- 'sex'

#Merge
MDD <- inner_join(MDD, sex, by='eid')

#Subset to females only
MDDFem <- subset(MDD, sex==0) 

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('ukb18177_glanville_covariates.txt' , data.table=FALSE)
colnames(covariates)[2] <- 'eid'
covariates <- covariates[,-1]

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/DEPR06_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c('eid','DEPR06_Pt_0.05','DEPR06_Pt_1')]

#Merge to the postpartum depression dataset
MDDFem <- left_join(MDDFem, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
MDDFem$DEPR06_Pt_0.05 <- scale(MDDFem$DEPR06_Pt_0.05)
MDDFem$DEPR06_Pt_1 <- scale(MDDFem$DEPR06_Pt_1)

#Make batch and assessment centre factors.
MDDFem$assessment_centre <- factor(MDDFem$assessment_centre, ordered=FALSE)
MDDFem$batch <- factor(MDDFem$batch, ordered=FALSE)

scores<-c("DEPR06_Pt_0.05","DEPR06_Pt_1")

MDDFemresult <- data.frame()

for (j in 1:length(scores)){
  MDDFem$scores=MDDFem[,scores[j]]
  
  # create model
  # full model
  model<-glm(MDD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=MDDFem,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(MDD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=MDDFem,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  MDDFemresult <- rbind(MDDFemresult, output)
} 

MDDFemresult$Pheno <- 'All' 

write.csv(MDDFemresult, 'Heterogeneous_MDDPRSresults_FemaleOnly')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#BPD PRS

#Read in cases
cases <- fread('MHQCaseExcPPDandDisorders_QC', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

#Read in controls
controls <- fread('controlsMDDExcPPDandDisorders_QC.csv', data.table=FALSE)
controls$MDD <- 0
colnames(controls)[2] <- 'eid'
controls <- controls[,-1]

#Merge
MDD <- rbind(cases, controls)

#Read in sex
sex <- readRDS("sex.rds")
colnames(sex)[2] <- 'sex'

#Merge
MDD <- inner_join(MDD, sex, by='eid')

#Subset to females only
MDDFem <- subset(MDD, sex==0)

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('/mnt/lustre/groups/ukbiobank/ukb18177_glanville/2019_new_qc/ukb18177_glanville_covariates.txt' , data.table=FALSE)
colnames(covariates)[2] <- 'eid'
covariates <- covariates[,-1]

#Read in the PRS - thresholds of 0.01 and 1.
PRS <- fread('./PRS/BIPO02_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c('eid','BIPO02_Pt_0.01','BIPO02_Pt_1')]

#Merge to the postpartum depression dataset
MDDFem <- left_join(MDDFem, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
MDDFem$BIPO02_Pt_0.01 <- scale(MDDFem$BIPO02_Pt_0.01)
MDDFem$BIPO02_Pt_1 <- scale(MDDFem$BIPO02_Pt_1)

#Make batch and assessment centre factors.
MDDFem$assessment_centre <- factor(MDDFem$assessment_centre, ordered=FALSE)
MDDFem$batch <- factor(MDDFem$batch, ordered=FALSE)

scores<-c("BIPO02_Pt_0.01","BIPO02_Pt_1")

MDDFemresult <- data.frame()

for (j in 1:length(scores)){
  MDDFem$scores=MDDFem[,scores[j]]
  
  # create model
  # full model
  model<-glm(MDD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=MDDFem,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(MDD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=MDDFem,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  MDDFemresult <- rbind(MDDFemresult, output)
} 

MDDFemresult$Pheno <- 'All' 

write.csv(MDDFemresult, 'Heterogeneous_BPDPRSresults_FemaleOnly')

###########################################################################################################################################################################
###########################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Environmental Variables

environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

#Read in cases
cases <- fread('MHQCaseExcPPDandDisorders', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

#Read in controls
controls <- fread('controlsMDDExcPPDandDisorders.csv', data.table=FALSE)
controls$MDD <- 0
colnames(controls)[2] <- 'eid'
controls <- controls[,-1]

#Merge
MDD <- rbind(cases, controls)

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

MDD <- inner_join(MDD, covariates, by='eid') %>%
  inner_join(., environment)

#Calculate frequency of endorsement among trauma variables - are we losing disproportionately from cases or controls
trauma <- c('loved_as_childBINARY','doctor_as_childBINARY','confiding_relationshipBINARY','pay_mortgageBINARY','serious_accidentBINARY','exposed_combatBINARY','life_threatening_illnessBINARY','violent_crimeBINARY','violent_deathBINARY','family_member_hateBINARY','physical_abuseBINARY','sexual_molestBINARY','partner_belittleBINARY','partner_violenceBINARY','partner_sexualBINARY','sexual_assaultBINARY')

for(i in trauma){
  print(i)
  MDDcase <- subset(MDD, MDD==1)
  MDDcontrol <- subset(MDD, MDD==0)
  print(table(MDD[[i]]))
  print(sum(is.na(MDD[[i]])))
  print(table(MDDcase[[i]]))
  print(sum(is.na(MDDcase[[i]])))
  print(table(MDDcontrol[[i]]))
  print(sum(is.na(MDDcontrol[[i]])))
}

env_variables <- c('TDI','neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

for(i in env_variables){
  print(i)
  variable <- MDD[is.na(MDD[[i]]),]
  variabletwo <- subset(variable, MDD==1)
  print(length(variable$eid))
  print(length(variabletwo$eid))
}

#All
table(MDD$education_attainment)

#Cases
case <- subset(MDD, MDD==1)
table(case$education_attainment)

#Controls
control <- subset(MDD, MDD==0)
table(control$education_attainment)

MDD <- subset(MDD, !(education_attainment == "missing" | education_attainment == "pnta" | is.na(TDI)))

MDD$assessment_centre <- factor(MDD$assessment_centre, ordered=FALSE)
MDD$education_attainment <- factor(MDD$education_attainment, ordered=FALSE)
MDD$education_attainment <- relevel(MDD$education_attainment, ref='degree')

MDD$TDI <- scale(MDD$TDI)
MDD$neuroticism <- scale(MDD$neuroticism)

env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

for(i in env_variables){
  print(i)
  variable <- MDD[is.na(MDD[[i]]),]
  variabletwo <- subset(variable, MDD==1)
  print(length(variable$eid))
  print(length(variabletwo$eid))
}

#Regression
MDDEnvResults <- data.frame()

for(j in 1:length(env_variables)) {
  
  MDD$variable=MDD[,env_variables[j]]
  
  print(env_variables[j])
  
  # create model
  
  # full model
  model<-glm(MDD ~ variable + assessment_centre + YOB + TDI + education_attainment,
             data=MDD,family=binomial,na.action=na.exclude)
  
  # get results
  print(summary(model)$coefficients)
  EnvVariable=env_variables[j]
  beta=(summary(model)$coefficients["variable",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["variable",2])
  p=(summary(model)$coefficients["variable",4])
  output=data.frame(EnvVariable,beta,se,OR,p)
  print(output)
  MDDEnvResults <- rbind(MDDEnvResults, output)
} 

#Perform separate regressions with just educational attainment and SES
model <- glm(MDD ~ education_attainment + assessment_centre + TDI + YOB,
             data=MDD,family=binomial,na.action=na.exclude)

beta=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),1])
OR=exp(beta)
se=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),2])
p=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),4])
EnvVariable <- c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI")
output=data.frame(EnvVariable,beta,se,OR,p)
print(output)
MDDEnvResults <- rbind(MDDEnvResults, output)

write.csv(MDDEnvResults, 'EnvironmentalAssociations_Heterogeneous_FullSample.csv')

###########################################################################################################################################################################
###########################################################################################################################################################################

#Supplementary Analysis - Subset to females and re-run

#Read in cases and controls for PPD - get an idea of degree of missingness for each variable
environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

#Read in cases
cases <- fread('MHQCaseExcPPDandDisorders', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

#Read in controls
controls <- fread('controlsMDDExcPPDandDisorders.csv', data.table=FALSE)
controls$MDD <- 0
colnames(controls)[2] <- 'eid'
controls <- controls[,-1]

#Merge
MDD <- rbind(cases, controls)

#Read in sex
sex <- readRDS("sex.rds")
colnames(sex)[2] <- 'sex'

#Merge
MDD <- left_join(MDD, sex, by='eid')

#Subset to females only
MDDFem <- subset(MDD, sex==0)

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

MDDFem <- inner_join(MDDFem, covariates, by='eid') %>%
  inner_join(., environment)

#Remove missing cases of SES and education attainment and see level of missingness of environmental variable
MDD <- subset(MDD, !(education_attainment == "missing" | education_attainment == "pnta" | is.na(TDI)))

MDDFem$assessment_centre <- factor(MDDFem$assessment_centre, ordered=FALSE)
MDDFem$education_attainment <- factor(MDDFem$education_attainment, ordered=FALSE)
MDDFem$education_attainment <- relevel(MDDFem$education_attainment, ref='degree')

MDDFem$TDI <- scale(MDDFem$TDI)
MDDFem$neuroticism <- scale(MDDFem$neuroticism)

#Most of the missingness driven by controls - of the two case definitions used - is one contributing to missingness more than the other?
env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

#Regression
MDDEnvResults <- data.frame()

for(j in 1:length(env_variables)) {
  
  MDDFem$variable=MDDFem[,env_variables[j]]
  
  print(env_variables[j])
  
  # create model
  
  # full model
  model<-glm(MDD~ variable + assessment_centre + TDI + YOB + education_attainment,
             data=MDDFem,family=binomial,na.action=na.exclude)
  
  # get results
  EnvVariable=env_variables[j]
  beta=(summary(model)$coefficients["variable",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["variable",2])
  p=(summary(model)$coefficients["variable",4])
  output=data.frame(EnvVariable,beta,se,OR,p)
  print(output)
  MDDEnvResults <- rbind(MDDEnvResults, output)
} 

#Perform separate regressions with just educational attainment and SES
model <- glm(MDD ~ education_attainment + assessment_centre + TDI + YOB,
             data=MDDFem,family=binomial,na.action=na.exclude)

beta=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),1])
OR=exp(beta)
se=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),2])
p=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),4])
EnvVariable <- c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI")
output=data.frame(EnvVariable,beta,se,OR,p)
print(output)
MDDEnvResults <- rbind(MDDEnvResults, output)

MDDEnvResults <- rbind(MDDEnvResults, output)

write.csv(MDDEnvResults, 'EnvironmentalAssociations_Heterogeneous_FemaleOnly.csv')