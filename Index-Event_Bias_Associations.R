#Depression following diagnosis of a chronic disease

#Association of the index event - diagnosis of a chronic disease

#Load libraries
library(data.table)
library(dplyr)
library(ukbkings)
library(tidyverse)

#Cases - any HES record of a disease
Disorders <- fread('./Diseases/Disorders', data.table=FALSE)
Disorders <- Disorders[,-1]

Disorders$disease <- 1 
Disorders <- Disorders[,c(1,42)]

#Controls - No HES record and no self-report of the disease. 

##Removal Step 2: Self-reported any of the disorders at nurse interview. 

##Cancer Diagnoses
CancerSR <- readRDS('SelfReportMedicalIllnessCancer.rds')

#Only concerned with diagnosis not year
CancerSR <- CancerSR[,c(1,3:26)]

#If the individual hasn't ever reported a cancer the total NAs for each person will equal the total number of columns in which the person had the opportunity to report which is 24
CancerSR$CancerNumber <- apply(CancerSR[,c(2:25)], 1, function(x) sum(is.na(x)))

CancerSR <- subset(CancerSR, CancerNumber==24)

CancerSR <- CancerSR[,c(1,26)]

##Non-cancer Diagnoses
NonCancerSR <- readRDS('SelfReportMedicalIllnessNonCancer.rds')

#Only concerned with diagnosis not year
NonCancerSR <- NonCancerSR[,c(1,3:138)]

NonCancerSR$NonCancer <- apply(NonCancerSR, 1, function(r) any(r %in% c(1264,1222,1223,1220,1081,1425,1491,1086,1075,1093,1261,1259,1311,1456,1154,1459,1463,1331,1453,1464,1382,1226,1462,1531,1377,1428,1082,1437,1313,1381,1461,1376,1521,1076,1083,1260,1345,1583)))

NonCancerSR <- subset(NonCancerSR, NonCancer==FALSE)

NonCancerSR <- NonCancerSR[,c(1,138)]

#Merge cancer and non-cancer to get controls - needs to be an inner join as cant have reported either a cancer or non-cancer medical illness
controls <- inner_join(CancerSR, NonCancerSR, by='eid')

#Anti join with the diseases above to make mutually exclusive
controls <- anti_join(controls, Disorders, by='eid')

controls$disease <- 0 
controls <- controls[,c(1,4)]

#Merge cases and controls
disease <- full_join(Disorders, controls)

#Apply exclusion critera - check it is the full sample
exclusion <- fread('ExclusionCriteriaApplied.csv', data.table=FALSE)
exclusion <- exclusion[,-1]

disease <- inner_join(x=disease, y=exclusion, by="eid") 

#Regression analysis

#Environmental Variables

environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

disease <- inner_join(disease, covariates, by='eid') %>%
  inner_join(., environment)

env_variables <- c('TDI','neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

for(i in env_variables){
  print(i)
  variable <- disease[is.na(disease[[i]]),]
  variabletwo <- subset(variable, disease==1)
  print(length(variable$eid))
  print(length(variabletwo$eid))
}

#All
table(disease$education_attainment)

#Remove missing cases of SES and education attainment and see level of missingness of environmental variable
disease <- subset(disease, !(education_attainment == "missing" | education_attainment == "pnta" | is.na(TDI)))

disease$assessment_centre <- factor(disease$assessment_centre, ordered=FALSE)
disease$education_attainment <- factor(disease$education_attainment, ordered=FALSE)
disease$education_attainment <- relevel(disease$education_attainment, ref='degree')

disease$TDI <- scale(disease$TDI)
disease$neuroticism <- scale(disease$neuroticism)

env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

for(i in env_variables){
  print(i)
  variable <- disease[is.na(disease[[i]]),]
  variabletwo <- subset(variable, disease==1)
  print(length(variable$eid))
  print(length(variabletwo$eid))
}

#Regression
DiseaseEnvResults <- data.frame()

for(j in 1:length(env_variables)) {
  
  disease$variable=disease[,env_variables[j]]
  
  print(env_variables[j])
  
  # create model
  
  # full model
  model<-glm(disease ~ variable + assessment_centre + YOB + TDI + education_attainment,
             data=disease,family=binomial,na.action=na.exclude)
  
  # get results
  print(summary(model)$coefficients)
  EnvVariable=env_variables[j]
  beta=(summary(model)$coefficients["variable",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["variable",2])
  p=(summary(model)$coefficients["variable",4])
  output=data.frame(EnvVariable,beta,se,OR,p)
  print(output)
  DiseaseEnvResults <- rbind(DiseaseEnvResults, output)
} 

#Perform separate regressions with just educational attainment and SES
model <- glm(disease ~ education_attainment + assessment_centre + TDI + YOB,
             data=disease,family=binomial,na.action=na.exclude)

beta=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),1])
OR=exp(beta)
se=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),2])
p=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),4])
EnvVariable <- c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI")
output=data.frame(EnvVariable,beta,se,OR,p)
print(output)
DiseaseEnvResults <- rbind(DiseaseEnvResults, output)

write.csv(DiseaseEnvResults, 'IndexDisease_EnvironmentalAssociations')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#MDD PRS

#Subset to european individuals
europeans <- fread('ukb18177_glanville_post_qc_id_list.fam' , data.table=FALSE)
colnames(europeans)[1] <- 'eid'

disease <- inner_join(disease, europeans)

disease <- disease[,c(1,2)]

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('ukb18177_glanville_covariates.txt' , data.table=FALSE)
colnames(covariates)[2] <- 'eid'
covariates <- covariates[,-1]

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/DEPR06_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c('eid','DEPR06_Pt_0.05','DEPR06_Pt_1')]

#Merge to the postpartum depression dataset
disease <- left_join(disease, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
disease$DEPR06_Pt_0.05 <- scale(disease$DEPR06_Pt_0.05)
disease$DEPR06_Pt_1 <- scale(disease$DEPR06_Pt_1)

#Make batch and assessment centre factors.
disease$assessment_centre <- factor(disease$assessment_centre, ordered=FALSE)
disease$batch <- factor(disease$batch, ordered=FALSE)

scores<-c("DEPR06_Pt_0.05","DEPR06_Pt_1")

Diseaseresult <- data.frame()

for (j in 1:length(scores)){
  disease$scores=disease[,scores[j]]
  
  # create model
  # full model
  model<-glm(disease~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=disease,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(disease~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=disease,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  Diseaseresult <- rbind(Diseaseresult, output)
} 

write.csv(Diseaseresult, 'IndexDisease_MDDPRS')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#BPD PRS

disease <- disease[,c(1,2)]

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/BIPO02_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c('eid','BIPO02_Pt_0.01','BIPO02_Pt_1')]

#Merge to the postpartum depression dataset
disease <- left_join(disease, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
disease$BIPO02_Pt_0.01 <- scale(disease$BIPO02_Pt_0.01)
disease$BIPO02_Pt_1 <- scale(disease$BIPO02_Pt_1)

#Make batch and assessment centre factors.
disease$assessment_centre <- factor(disease$assessment_centre, ordered=FALSE)
disease$batch <- factor(disease$batch, ordered=FALSE)

scores<-c("BIPO02_Pt_0.01","BIPO02_Pt_1")

Diseaseresult <- data.frame()

for (j in 1:length(scores)){
  disease$scores=disease[,scores[j]]
  
  # create model
  # full model
  model<-glm(disease~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=disease,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(disease~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=disease,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  Diseaseresult <- rbind(Diseaseresult, output)
} 

write.csv(Diseaseresult, 'IndexDisease_BPDPRS')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Associations with giving birth vs not giving birth

#Extract all cases from the maternity list where the ID corresponds to the IDs found in noncases. 
project_dir <- "/scratch/datasets/ukbiobank/ukb18177"

hesin_maternity <- bio_hesin(project_dir, "maternity")
hesin_maternity <- hesin_maternity[,c('eid','ins_index')] 

hesin_diag <- bio_hesin(project_dir, "diag")

#Extract codes for live births in the diagnosis dataset. 
hesin_birth <- hesin_diag[with(hesin_diag, grepl("Z370|Z371|Z372|Z375|Z38|Z390", paste(diag_icd9, diag_icd10), ignore.case=FALSE)),]
hesin_birth <- hesin_birth[,c('eid','ins_index')] 

#Full join so that complete cases from both variables are identified
hesin_maternity <- full_join(hesin_birth, hesin_maternity, by=c('eid','ins_index'))

#Make so only unique IDs. 
cases <- unique(hesin_maternity$eid)
cases <- as.data.frame(cases)
colnames(cases) <- 'eid'
cases$give_birth <- 1

#Self report Cases
childbirth <- readRDS('ChildbirthSelfReport.rds') 
colnames(childbirth) <- c('eid','children_one','children_two','children_three','children_four')
childbirth$children <- apply(childbirth[,2:5], 1, function(r) any(r > 0))

#subset to cases and merge with above
childbirthcases <- subset(childbirth, children==TRUE)
childbirthcases$give_birth <- 1
childbirthcases <- childbirthcases[,c('eid','give_birth')]

#Merge with current case set
childbirthcases <- full_join(childbirthcases, cases)

#Define control set
controls_wmales <- subset(childbirth, children==FALSE | is.na(children))

#Anti-join with cases to make exclusive
controls_wmales <- anti_join(controls_wmales, childbirthcases, by='eid')
controls_wmales$give_birth <-  0
controls_wmales <- controls_wmales[,c('eid','give_birth')]

#Control set two - female only - in line with the definition of heterogeneous depression initially
sex <- readRDS('sex.rds')
colnames(sex) <- c('eid','sex')

controls_womales <- inner_join(controls_wmales, sex)
controls_womales <- subset(controls_womales, sex==0)
controls_womales <- controls_womales[,c('eid','give_birth')]

#Merge the cases and controls

birth_wmales <- full_join(childbirthcases, controls_wmales)

birth_womales <- full_join(childbirthcases, controls_womales)

##############################################################################################################################################
##############################################################################################################################################

#Regression analysis - First set of controls (including males)

#Environmental Variables

environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

birth_wmales <- inner_join(birth_wmales, covariates, by='eid') %>%
  inner_join(., environment)

env_variables <- c('TDI','neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

for(i in env_variables){
  print(i)
  variable <- birth_wmales[is.na(birth_wmales[[i]]),]
  variabletwo <- subset(variable, give_birth==1)
  print(length(variable$eid))
  print(length(variabletwo$eid))
}

#All
table(birth_wmales$education_attainment)

birth_wmales <- subset(birth_wmales, !(education_attainment == "missing" | education_attainment == "pnta" | is.na(TDI)))

birth_wmales$assessment_centre <- factor(birth_wmales$assessment_centre, ordered=FALSE)
birth_wmales$education_attainment <- factor(birth_wmales$education_attainment, ordered=FALSE)
birth_wmales$education_attainment <- relevel(birth_wmales$education_attainment, ref='degree')

birth_wmales$TDI <- scale(birth_wmales$TDI)
birth_wmales$neuroticism <- scale(birth_wmales$neuroticism)

env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

for(i in env_variables){
  print(i)
  variable <- birth_wmales[is.na(birth_wmales[[i]]),]
  variabletwo <- subset(variable, give_birth==1)
  print(length(variable$eid))
  print(length(variabletwo$eid))
}

#Regression
BirthEnvResults <- data.frame()

for(j in 1:length(env_variables)) {
  
  birth_wmales$variable=birth_wmales[,env_variables[j]]
  
  print(env_variables[j])
  
  # create model
  
  # full model
  model<-glm(give_birth ~ variable + assessment_centre + YOB + TDI + education_attainment,
             data=birth_wmales,family=binomial,na.action=na.exclude)
  
  # get results
  print(summary(model)$coefficients)
  EnvVariable=env_variables[j]
  beta=(summary(model)$coefficients["variable",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["variable",2])
  p=(summary(model)$coefficients["variable",4])
  output=data.frame(EnvVariable,beta,se,OR,p)
  print(output)
  BirthEnvResults <- rbind(BirthEnvResults, output)
} 

#Perform separate regressions with just educational attainment and SES
model <- glm(give_birth ~ education_attainment + assessment_centre + TDI + YOB,
             data=birth_wmales,family=binomial,na.action=na.exclude)

beta=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),1])
OR=exp(beta)
se=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),2])
p=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),4])
EnvVariable <- c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI")
output=data.frame(EnvVariable,beta,se,OR,p)
print(output)
BirthEnvResults <- rbind(BirthEnvResults, output)

write.csv(BirthEnvResults, 'IndexGivingBirth_EnvironmentalAssociations_withselfreport')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#MDD PRS

#Subset to european individuals
europeans <- fread('ukb18177_glanville_post_qc_id_list.fam' , data.table=FALSE)
colnames(europeans)[1] <- 'eid'

birth_wmales <- inner_join(birth_wmales, europeans)

birth_wmales <- birth_wmales[,c(1,2)]

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('ukb18177_glanville_covariates.txt' , data.table=FALSE)
colnames(covariates)[2] <- 'eid'
covariates <- covariates[,-1]

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/DEPR06_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c('eid','DEPR06_Pt_0.05','DEPR06_Pt_1')]

#Merge to the postpartum depression dataset
birth_wmales <- left_join(birth_wmales, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
birth_wmales$DEPR06_Pt_0.05 <- scale(birth_wmales$DEPR06_Pt_0.05)
birth_wmales$DEPR06_Pt_1 <- scale(birth_wmales$DEPR06_Pt_1)

#Make batch and assessment centre factors.
birth_wmales$assessment_centre <- factor(birth_wmales$assessment_centre, ordered=FALSE)
birth_wmales$batch <- factor(birth_wmales$batch, ordered=FALSE)

scores<-c("DEPR06_Pt_0.05","DEPR06_Pt_1")

Birthresult <- data.frame()

for (j in 1:length(scores)){
  birth_wmales$scores=birth_wmales[,scores[j]]
  
  # create model
  # full model
  model<-glm(give_birth~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=birth_wmales,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(give_birth~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=birth_wmales,family=binomial,na.action=na.exclude)
  # get results
  print(summary(model)$coefficients)
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  Birthresult <- rbind(Birthresult, output)
} 

write.csv(Birthresult, 'IndexGivingBirth_MDDPRS_withselfreport')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#BPD PRS

birth_wmales <- birth_wmales[,c(1,2)]

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/BIPO02_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c('eid','BIPO02_Pt_0.01','BIPO02_Pt_1')]

#Merge to the postpartum depression dataset
birth_wmales <- left_join(birth_wmales, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
birth_wmales$BIPO02_Pt_0.01 <- scale(birth_wmales$BIPO02_Pt_0.01)
birth_wmales$BIPO02_Pt_1 <- scale(birth_wmales$BIPO02_Pt_1)

#Make batch and assessment centre factors.
birth_wmales$assessment_centre <- factor(birth_wmales$assessment_centre, ordered=FALSE)
birth_wmales$batch <- factor(birth_wmales$batch, ordered=FALSE)

scores <- c("BIPO02_Pt_0.01","BIPO02_Pt_1")

Birthresult <- data.frame()

for (j in 1:length(scores)){
  birth_wmales$scores=birth_wmales[,scores[j]]
  
  # create model
  # full model
  model<-glm(give_birth~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=birth_wmales,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(give_birth~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=birth_wmales,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  Birthresult <- rbind(Birthresult, output)
} 

write.csv(Birthresult, 'IndexGivingBirth_BPDPRS_withselfreport')

##############################################################################################################################################
##############################################################################################################################################

#Regression analysis - Second set of controls (excluding males)

#Environmental Variables
environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

birth_womales <- inner_join(birth_womales, covariates, by='eid') %>%
  inner_join(., environment)

env_variables <- c('TDI','neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

for(i in env_variables){
  print(i)
  variable <- birth_womales[is.na(birth_womales[[i]]),]
  variabletwo <- subset(variable, give_birth==1)
  print(length(variable$eid))
  print(length(variabletwo$eid))
}

#All
table(birth_womales$education_attainment)

birth_womales <- subset(birth_womales, !(education_attainment == "missing" | education_attainment == "pnta" | is.na(TDI)))

birth_womales$assessment_centre <- factor(birth_womales$assessment_centre, ordered=FALSE)
birth_womales$education_attainment <- factor(birth_womales$education_attainment, ordered=FALSE)
birth_womales$education_attainment <- relevel(birth_womales$education_attainment, ref='degree')

birth_womales$TDI <- scale(birth_womales$TDI)
birth_womales$neuroticism <- scale(birth_womales$neuroticism)

env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

for(i in env_variables){
  print(i)
  variable <- birth_womales[is.na(birth_womales[[i]]),]
  variabletwo <- subset(variable, give_birth==1)
  print(length(variable$eid))
  print(length(variabletwo$eid))
}

#Regression
BirthEnvResults <- data.frame()

for(j in 1:length(env_variables)) {
  
  birth_womales$variable=birth_womales[,env_variables[j]]
  
  print(env_variables[j])
  
  # create model
  
  # full model
  model<-glm(give_birth ~ variable + assessment_centre + YOB + TDI + education_attainment,
             data=birth_womales,family=binomial,na.action=na.exclude)
  
  # get results
  print(summary(model)$coefficients)
  EnvVariable=env_variables[j]
  beta=(summary(model)$coefficients["variable",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["variable",2])
  p=(summary(model)$coefficients["variable",4])
  output=data.frame(EnvVariable,beta,se,OR,p)
  print(output)
  BirthEnvResults <- rbind(BirthEnvResults, output)
} 

#Perform separate regressions with just educational attainment and SES
model <- glm(give_birth ~ education_attainment + assessment_centre + TDI + YOB,
             data=birth_womales,family=binomial,na.action=na.exclude)

beta=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),1])
OR=exp(beta)
se=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),2])
p=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),4])
EnvVariable <- c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI")
output=data.frame(EnvVariable,beta,se,OR,p)
print(output)
BirthEnvResults <- rbind(BirthEnvResults, output)

write.csv(BirthEnvResults, 'IndexGivingBirth_EnvironmentalAssociations_FemaleOnly_withselfreport')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#MDD PRS

#Subset to european individuals
europeans <- fread('/mnt/lustre/groups/ukbiobank/ukb18177_glanville/2019_new_qc/ukb18177_glanville_post_qc_id_list.fam' , data.table=FALSE)
colnames(europeans)[1] <- 'eid'

birth_womales <- inner_join(birth_womales, europeans)

birth_womales <- birth_womales[,c(1,2)]

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('/mnt/lustre/groups/ukbiobank/ukb18177_glanville/2019_new_qc/ukb18177_glanville_covariates.txt' , data.table=FALSE)
colnames(covariates)[2] <- 'eid'
covariates <- covariates[,-1]

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/DEPR06_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c('eid','DEPR06_Pt_0.05','DEPR06_Pt_1')]

#Merge to the postpartum depression dataset
birth_womales <- left_join(birth_womales, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
birth_womales$DEPR06_Pt_0.05 <- scale(birth_womales$DEPR06_Pt_0.05)
birth_womales$DEPR06_Pt_1 <- scale(birth_womales$DEPR06_Pt_1)

#Make batch and assessment centre factors.
birth_womales$assessment_centre <- factor(birth_womales$assessment_centre, ordered=FALSE)
birth_womales$batch <- factor(birth_womales$batch, ordered=FALSE)

scores<-c("DEPR06_Pt_0.05","DEPR06_Pt_1")

Birthresult <- data.frame()

for (j in 1:length(scores)){
  birth_womales$scores=birth_womales[,scores[j]]
  
  # create model
  # full model
  model<-glm(give_birth~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=birth_womales,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(give_birth~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=birth_womales,family=binomial,na.action=na.exclude)
  # get results
  print(summary(model)$coefficients)
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  Birthresult <- rbind(Birthresult, output)
} 

write.csv(Birthresult, 'IndexGivingBirth_MDDPRS_FemaleOnly_withselfreport')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#BPD PRS

birth_womales <- birth_womales[,c(1,2)]

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/BIPO02_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c('eid','BIPO02_Pt_0.01','BIPO02_Pt_1')]

#Merge to the postpartum depression dataset
birth_womales <- left_join(birth_womales, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
birth_womales$BIPO02_Pt_0.01 <- scale(birth_womales$BIPO02_Pt_0.01)
birth_womales$BIPO02_Pt_1 <- scale(birth_womales$BIPO02_Pt_1)

#Make batch and assessment centre factors.
birth_womales$assessment_centre <- factor(birth_womales$assessment_centre, ordered=FALSE)
birth_womales$batch <- factor(birth_womales$batch, ordered=FALSE)

scores <- c("BIPO02_Pt_0.01","BIPO02_Pt_1")

Birthresult <- data.frame()

for (j in 1:length(scores)){
  birth_womales$scores=birth_womales[,scores[j]]
  
  # create model
  # full model
  model<-glm(give_birth~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=birth_womales,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(give_birth~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=birth_womales,family=binomial,na.action=na.exclude)
  # get results
  #print(summary(model)$coefficients)
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,p)
  print(output)
  Birthresult <- rbind(Birthresult, output)
} 

write.csv(Birthresult, 'IndexGivingBirth_BPDPRS_FemaleOnly_withselfreport')

##############################################################################################################################################
##############################################################################################################################################