#Control comparisons can show the importance of control selection dependent on risk exposure also. 

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

cases <- fread('controlsMedicalDisorder_QC.csv', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,c('eid','MDD')]

controls <- fread('controlsMDDExcPPDandDisorders_QC.csv', data.table=FALSE)
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
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  MDDresult <- rbind(MDDresult, output)
} 

write.csv(MDDresult, 'MDDwithDisorder_ControlControl_MDDPRSresults')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#BPD PRS

cases <- fread('controlsMedicalDisorder_QC.csv', data.table=FALSE)
cases$BPD <- 1
cases <- cases[,c('eid','BPD')]

controls <- fread('controlsMDDExcPPDandDisorders_QC.csv', data.table=FALSE)
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

write.csv(BPDresult, 'MDDwithDisorder_ControlControl_BPDPRSresults')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Environmental Variables Analysis
environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

#European ancestry only
cases <- fread('controlsMedicalDisorder.csv', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,c('eid','MDD')]

controls <- fread('controlsMDDExcPPDandDisorders.csv', data.table=FALSE)
controls$MDD <- 0
controls <- controls[,-1]
colnames(controls)[1] <- 'eid'

MDD <- rbind(cases, controls)

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

MDD <- inner_join(MDD, covariates, by='eid') %>%
  inner_join(., environment)

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

write.csv(MDDEnvResults, 'EnvironmentalAssociationsMDDthenDep_ControlControl_FullSample.csv')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#TEST 2: MDDExc vs PostPartum

#MDD PRS

cases <- fread('controlsPPD_with_selfreport_QC.csv', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,c('eid','MDD')]

controls <- fread('controlsMDDExcPPDandDisorders_QC.csv', data.table=FALSE)
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
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  MDDresult <- rbind(MDDresult, output)
} 

write.csv(MDDresult, 'PostPartum_ControlControl_MDDPRSresults_withselfreport')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#BPD PRS

cases <- fread('controlsPPD_with_selfreport_QC.csv', data.table=FALSE)
cases$BPD <- 1
cases <- cases[,c('eid','BPD')]

controls <- fread('controlsMDDExcPPDandDisorders_QC.csv', data.table=FALSE)
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

write.csv(BPDresult, 'PostPartum_ControlControl_BPDPRSresults_withselfreport')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Environmental Variables Analysis
environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

cases <- fread('controlsPPD_with_selfreport.csv', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,c('eid','MDD')]

controls <- fread('controlsMDDExcPPDandDisorders.csv', data.table=FALSE)
controls$MDD <- 0
controls <- controls[,-1]
colnames(controls)[1] <- 'eid'

MDD <- rbind(cases, controls)

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

MDD <- inner_join(MDD, covariates, by='eid') %>%
  inner_join(., environment)

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

write.csv(MDDEnvResults, 'EnvironmentalAssociationsPostPartum_ControlControl_FullSample_withselfreport.csv')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
