#The depression group excluding the other two definitions (MDDexc) will act as the reference in this test

#TEST 1: MDDExc vs DiseaseThenDep

### Load dependencies
library(data.table)
library(tidyverse)
library(ukbkings)
library(dplyr)
library(ggplot2)
library(polycor)
library(lubridate)
library(qdapRegex)
library(fmsb)

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#MDD PRS

cases <- fread('depression_within_year_selfreport_ScreenedandQCd', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

controls <- fread('MHQCaseExcPPDandDisorders_QC', data.table=FALSE)
controls$MDD <- 0
controls <- controls[,-1]
colnames(controls)[1] <- 'eid'

MDD <- rbind(cases, controls)

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('covariates.txt' , data.table=FALSE)
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
  print(scores[j])
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

write.csv(MDDresult, 'MDDwithDisorder_CaseCase_MDDPRSresults')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#BPD PRS

cases <- fread('depression_within_year_selfreport_ScreenedandQCd', data.table=FALSE)
cases$BPD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

controls <- fread('MHQCaseExcPPDandDisorders_QC', data.table=FALSE)
controls$BPD <- 0
controls <- controls[,-1]
colnames(controls)[1] <- 'eid'

BPD <- rbind(cases, controls)

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('covariates.txt' , data.table=FALSE)
colnames(covariates)[2] <- 'eid'
covariates <- covariates[,-1]

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/BIPO02_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c("eid","BIPO02_Pt_0.01","BIPO02_Pt_1")]

#Merge to the postpartum depression dataset
BPD <- left_join(BPD, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
BPD$BIPO02_Pt_0.01 <- scale(BPD$BIPO02_Pt_0.01)
BPD$BIPO02_Pt_1 <- scale(BPD$BIPO02_Pt_1)

#Make batch and assessment centre factors.
BPD$assessment_centre <- factor(BPD$assessment_centre, ordered=FALSE)
BPD$batch <- factor(BPD$batch, ordered=FALSE)

scores<-c("BIPO02_Pt_0.01","BIPO02_Pt_1")

BPDresult <- data.frame()

for (j in 1:length(scores)){
  BPD$scores=BPD[,scores[j]]
  print(scores[j])
  # create model
  # full model
  model<-glm(BPD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=BPD,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(BPD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=BPD,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  BPDresult <- rbind(BPDresult, output)
} 

write.csv(BPDresult, 'MDDwithDisorder_CaseCase_BPDPRSresults')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Environmental Variables Analysis
environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

#European ancestry only
cases <- fread('depression_within_year_selfreport_Screened', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

controls <- fread('MHQCaseExcPPDandDisorders', data.table=FALSE)
controls$MDD <- 0
controls <- controls[,-1]
colnames(controls)[1] <- 'eid'

MDD <- rbind(cases, controls)

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

MDD <- inner_join(MDD, covariates, by='eid') %>%
  inner_join(., environment)

MDD <- subset(MDD, !(education_attainment == "missing" | education_attainment == "pnta" | is.na(TDI)))

MDD$assessment_centre <- factor(MDD$assessment_centre, ordered=FALSE)
MDD$education_attainment <- factor(MDD$education_attainment, ordered=FALSE)
MDD$education_attainment <- relevel(MDD$education_attainment, ref='degree')

MDD$TDI <- scale(MDD$TDI)
MDD$neuroticism <- scale(MDD$neuroticism)

env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

#Regression
MDDEnvResults <- data.frame()

for(j in 1:length(env_variables)) {
  
  MDD$variable=MDD[,env_variables[j]]
  
  print(env_variables[j])
  
  # create model
  
  # full model
  model<-glm(MDD ~ variable + assessment_centre + education_attainment + TDI + YOB,
             data=MDD,family=binomial,na.action=na.exclude)
  
  # get results
  #print(summary(model)$coefficients)
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


write.csv(MDDEnvResults, 'EnvironmentalAssociationsMDDthenDep_CaseCase_FullSample.csv')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#TEST 2: MDDExc vs PostPartum

#MDD PRS

cases <- fread('./PostPartum/PostPartumMHQandEHR_QC', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,c(3,8)]
colnames(cases)[1] <- 'eid'

controls <- fread('MHQCaseExcPPDandDisorders_QC', data.table=FALSE)
controls$MDD <- 0
controls <- controls[,-1]
colnames(controls)[1] <- 'eid'

MDD <- rbind(cases, controls)

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('/mnt/lustre/groups/ukbiobank/ukb18177_glanville/2019_new_qc/ukb18177_glanville_covariates.txt' , data.table=FALSE)
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
  print(scores[j])
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

write.csv(MDDresult, 'PostPartum_CaseCase_MDDPRSresults')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#BPD PRS

cases <- fread('./PostPartum/PostPartumMHQandEHR_QC', data.table=FALSE)
cases$BPD <- 1
cases <- cases[,c(3,8)]
colnames(cases)[1] <- 'eid'

controls <- fread('MHQCaseExcPPDandDisorders_QC', data.table=FALSE)
controls$BPD <- 0
controls <- controls[,-1]
colnames(controls)[1] <- 'eid'

BPD <- rbind(cases, controls)

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('covariates.txt' , data.table=FALSE)
colnames(covariates)[2] <- 'eid'
covariates <- covariates[,-1]

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/BIPO02_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c("eid","BIPO02_Pt_0.01","BIPO02_Pt_1")]

#Merge to the postpartum depression dataset
BPD <- left_join(BPD, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
BPD$BIPO02_Pt_0.01 <- scale(BPD$BIPO02_Pt_0.01)
BPD$BIPO02_Pt_1 <- scale(BPD$BIPO02_Pt_1)

#Make batch and assessment centre factors.
BPD$assessment_centre <- factor(BPD$assessment_centre, ordered=FALSE)
BPD$batch <- factor(BPD$batch, ordered=FALSE)

scores<-c("BIPO02_Pt_0.01","BIPO02_Pt_1")

BPDresult <- data.frame()

for (j in 1:length(scores)){
  BPD$scores=BPD[,scores[j]]
  print(scores[j])
  # create model
  # full model
  model<-glm(BPD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=BPD,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(BPD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=BPD,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  BPDresult <- rbind(BPDresult, output)
} 

write.csv(BPDresult, 'PostPartum_CaseCase_BPDPRSresults')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Environmental Variables Analysis
environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

#European ancestry only
cases <- fread('./PostPartum/PostPartumMHQandEHR', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

controls <- fread('MHQCaseExcPPDandDisorders', data.table=FALSE)
controls$MDD <- 0
controls <- controls[,-1]
colnames(controls)[1] <- 'eid'

MDD <- rbind(cases, controls)

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

MDD <- inner_join(MDD, covariates, by='eid') %>%
  inner_join(., environment)

MDD <- subset(MDD, !(education_attainment == "missing" | education_attainment == "pnta" | is.na(TDI)))

MDD$assessment_centre <- factor(MDD$assessment_centre, ordered=FALSE)
MDD$education_attainment <- factor(MDD$education_attainment, ordered=FALSE)
MDD$education_attainment <- relevel(MDD$education_attainment, ref='degree')

MDD$TDI <- scale(MDD$TDI)
MDD$neuroticism <- scale(MDD$neuroticism)

#Most of the missingness driven by controls - of the two case definitions used - is one contributing to missingness more than the other?
env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

#Regression
MDDEnvResults <- data.frame()

for(j in 1:length(env_variables)) {
  
  MDD$variable=MDD[,env_variables[j]]
  
  print(env_variables[j])
  
  # create model
  
  # full model
  model<-glm(MDD ~ variable + assessment_centre + education_attainment + TDI + YOB,
             data=MDD,family=binomial,na.action=na.exclude)
  
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
             data=MDD,family=binomial,na.action=na.exclude)

beta=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),1])
OR=exp(beta)
se=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),2])
p=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),4])
EnvVariable <- c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI")
output=data.frame(EnvVariable,beta,se,OR,p)
print(output)
MDDEnvResults <- rbind(MDDEnvResults, output)

write.csv(MDDEnvResults, 'EnvironmentalAssociationsPostPartum_CaseCase_FullSample.csv')

#Model everything jointly to see how it changes the association
model <- glm(MDD ~ neuroticism + child_trauma_PC1 + adult_trauma_PC1 + FamHist + education_attainment + assessment_centre + TDI + YOB,
             data=MDD,family=binomial,na.action=na.exclude)

beta=(summary(model)$coefficients[c(2:8,29),1])
OR=exp(beta)
se=(summary(model)$coefficients[c(2:8,29),2])
p=(summary(model)$coefficients[c(2:8,29),4])
EnvVariable <- c("neuroticism","child_trauma_PC1","adult_trauma_PC1","FamHist","education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI")
output=data.frame(EnvVariable,beta,se,OR,p)
print(output)

MDD$NAsum <- rowSums(is.na(MDD[,c("neuroticism","child_trauma_PC1","adult_trauma_PC1","FamHist","education_attainment","TDI")]))

sum(is.na(MDD$education_attainment))

reduced <- subset(MDD, NAsum == 0)

write.csv(output, 'PostPartum_CaseCase_FullModelAssociations.csv')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Supplementary Analysis: MDDExc vs PostPartum - Subsetting the MDDExc group to be female only

#MDD PRS

cases <- fread('./PostPartum/PostPartumMHQandEHR_QC', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,c(3,8)]
colnames(cases)[1] <- 'eid'

controls <- fread('MHQCaseExcPPDandDisorders_QC', data.table=FALSE)
controls$MDD <- 0
controls <- controls[,-1]
colnames(controls)[1] <- 'eid'

MDD <- rbind(cases, controls)

#Read in sex
sex <- readRDS("sex.rds")
colnames(sex)[2] <- 'sex'

#Merge
MDD <- inner_join(MDD, sex, by='eid')

#Subset to females only
MDD <- subset(MDD, sex==0)

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('covariates.txt' , data.table=FALSE)
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
  print(scores[j])
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

write.csv(MDDresult, 'PostPartum_CaseCase_MDDPRSresults_FemaleOnly')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#BPD PRS

cases <- fread('./PostPartum/PostPartumMHQandEHR_QC', data.table=FALSE)
cases$BPD <- 1
cases <- cases[,c(3,8)]
colnames(cases)[1] <- 'eid'

controls <- fread('MHQCaseExcPPDandDisorders_QC', data.table=FALSE)
controls$BPD <- 0
controls <- controls[,-1]
colnames(controls)[1] <- 'eid'

BPD <- rbind(cases, controls)

#Read in sex
sex <- readRDS("sex.rds")
colnames(sex)[2] <- 'sex'

#Merge
BPD <- inner_join(BPD, sex, by='eid')

#Subset to females only
BPD <- subset(BPD, sex==0)

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('covariates.txt' , data.table=FALSE)
colnames(covariates)[2] <- 'eid'
covariates <- covariates[,-1]

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/BIPO02_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c("eid","BIPO02_Pt_0.01","BIPO02_Pt_1")]

#Merge to the postpartum depression dataset
BPD <- left_join(BPD, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
BPD$BIPO02_Pt_0.01 <- scale(BPD$BIPO02_Pt_0.01)
BPD$BIPO02_Pt_1 <- scale(BPD$BIPO02_Pt_1)

#Make batch and assessment centre factors.
BPD$assessment_centre <- factor(BPD$assessment_centre, ordered=FALSE)
BPD$batch <- factor(BPD$batch, ordered=FALSE)

scores<-c("BIPO02_Pt_0.01","BIPO02_Pt_1")

BPDresult <- data.frame()

for (j in 1:length(scores)){
  BPD$scores=BPD[,scores[j]]
  print(scores[j])
  # create model
  # full model
  model<-glm(BPD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=BPD,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(BPD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=BPD,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  BPDresult <- rbind(BPDresult, output)
} 

write.csv(BPDresult, 'PostPartum_CaseCase_BPDPRSresults_FemaleOnly')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Environmental Variables Analysis
#Read in cases and controls for PPD - get an idea of degree of missingness for each variable
environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

#European ancestry only
cases <- fread('./PostPartum/PostPartumMHQandEHR', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

controls <- fread('MHQCaseExcPPDandDisorders', data.table=FALSE)
controls$MDD <- 0
controls <- controls[,-1]
colnames(controls)[1] <- 'eid'

MDD <- rbind(cases, controls)

#Read in sex
sex <- readRDS("sex.rds")
colnames(sex)[2] <- 'sex'

#Merge
MDD <- inner_join(MDD, sex, by='eid')

#Subset to females only
MDD <- subset(MDD, sex==0)

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

MDD <- inner_join(MDD, covariates, by='eid') %>%
  inner_join(., environment)

MDD <- subset(MDD, !(education_attainment == "missing" | education_attainment == "pnta" | is.na(TDI)))

MDD$assessment_centre <- factor(MDD$assessment_centre, ordered=FALSE)
MDD$education_attainment <- factor(MDD$education_attainment, ordered=FALSE)
MDD$education_attainment <- relevel(MDD$education_attainment, ref='degree')

MDD$TDI <- scale(MDD$TDI)
MDD$neuroticism <- scale(MDD$neuroticism)

#Most of the missingness driven by controls - of the two case definitions used - is one contributing to missingness more than the other?
env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

#Regression
MDDEnvResults <- data.frame()

for(j in 1:length(env_variables)) {
  
  MDD$variable=MDD[,env_variables[j]]
  
  print(env_variables[j])
  
  # create model
  
  # full model
  model<-glm(MDD ~ variable + assessment_centre + education_attainment + TDI + YOB,
             data=MDD,family=binomial,na.action=na.exclude)
  
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
             data=MDD,family=binomial,na.action=na.exclude)

beta=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),1])
OR=exp(beta)
se=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),2])
p=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),4])
EnvVariable <- c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI")
output=data.frame(EnvVariable,beta,se,OR,p)
print(output)
MDDEnvResults <- rbind(MDDEnvResults, output)

write.csv(MDDEnvResults, 'EnvironmentalAssociationsPostPartum_CaseCase_FemaleOnly_FullSample.csv')