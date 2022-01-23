### Load dependencies
library(data.table)
library(tidyverse)
library(ukbkings)
library(dplyr)
library(ggplot2)
library(polycor)
library(lubridate)
library(rlist)
library(qdapRegex)
library(sqldf)
library(fmsb)
library(car)
library(psych)
library(factoextra)
library(ggplot2)

#Definition 1:
#Post-natal depression - self-report using MHQ

#Read in current lifetime depression - includes current exclusion criteria applied except for genetic quality control. 
depressed <- fread('LifetimeDepression.csv', data.table=FALSE)
depressed <- depressed[,-1]

#Read in the depression datasets 
ppdep <- readRDS("./PostPartum/PostPartumDepression.rds")
colnames(ppdep)[2] <- 'Due.to.childbirth' 

#Merge with lifetime depression dataset 
ppdep <- inner_join(ppdep,depressed)

ppdep$PPDep <- case_when(
  ppdep$Due.to.childbirth == 1 & ppdep$LifetimeDepression == 1 ~ 1,
  TRUE ~ 999
) #3109

#Remove any cases mentioning substance abuse, schizophrenia and bipolar disorder
exclusion <- fread('ExclusionCriteriaApplied.csv', data.table=FALSE)
exclusion <- exclusion[,-1]

#Filter by exclusion
ppdep.and.exclusion <- subset(ppdep, eid %in% exclusion$eid)  

#Remove columns relating to exclusion criteria
ppdep.and.exclusion <- ppdep.and.exclusion[,c('eid','PPDep')]
ppdep.and.exclusion <- subset(ppdep.and.exclusion, PPDep == 1)
write.csv(ppdep.and.exclusion, './PostPartum/PostPartumDepression_MHQ') 

#How many of these survive genetic QC
#Get the fam file and inner join to work out people of european ancestries
europeans <- fread('ukb18177_glanville_post_qc_id_list.fam' , data.table=FALSE)
colnames(europeans)[1] <- 'eid'

ppdepQC <- inner_join(ppdep.and.exclusion, europeans)
ppdepQC <- ppdepQC[,-c(5:9)] 

write.csv(ppdepQC, './PostPartum/PostPartumDepression_MHQandQC') 

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Definition 2:
#Postpartum depression as recorded by GP records. 
project_dir <- "ukb18177"
gp_clinical <- bio_gp(project_dir, "clinical")

#Extract codes using GP data. 
#gp_clinicalSub <- gp_clinical[with(gp_clinical, grepl("12K8|62T1|8ID|E204|Eu530|Xaeft|XaX2L", paste(read_2, read_3), ignore.case=FALSE)),]

#Codes E204, 62T1 and Eu530 were the only ones that were relevant so further subset by these codes.
gp_clinicalSub <- gp_clinical[with(gp_clinical, grepl("62T1|E204|Eu530", paste(read_2, read_3), ignore.case=FALSE)),]

#Case or needs more info status - as is the case for E204 as it may be reactive subtype as opposed to postnatal depression.
gp_clinicalSub$case <- case_when(
  (gp_clinicalSub$read_2 == '62T1.' | gp_clinicalSub$read_2 == 'Eu530') | (gp_clinicalSub$read_3 == '62T1.' | gp_clinicalSub$read_3 == 'Eu530') ~ 1,
  TRUE ~ 0
)

case <- subset(gp_clinicalSub, case == 1)

#Extract first episode of these cases - Number of unique participants - 604 unique participants
case <- case %>%
          group_by(eid) %>%
          arrange(event_dt) %>%
          slice(1L)

case <- as.data.frame(case)

colnames(case)[3] <- 'PPD_Date'

#We dont know if this is the first episode from the primary care records - load GP depression - and see if previous history of depression
DepressionGPCodes <- fread('DepressionGP', data.table=FALSE)
DepressionGPCodes <- DepressionGPCodes[,-1]

#Left join with cases 
case <- left_join(case, DepressionGPCodes, by='eid')
case$PPD_Date <- as.Date(case$PPD_Date,"%d/%m/%Y")
case$difference <- time_length(difftime(case$PPD_Date, case$event_dt), "years")

dim(case[!is.na(case$difference) & case$difference > 0,])

case <- subset(case, is.na(difference) | difference <= 0)

#Take IDs as these are cases regardless whether we can determine a pregnancy or not. 
cases <- unique(case$eid)

#Change so that the code is now where it is equal to E204
E204 <- subset(gp_clinicalSub, read_2=='E204.' | read_3=='E204.')

E204 <- E204 %>%
          group_by(eid) %>%
          arrange(event_dt) %>%
          slice(1L)

#Extract maternity hospital information.

#Extract all cases from the maternity list where the ID corresponds to the IDs found in noncases. 
hesin_maternity <- bio_hesin(project_dir, "maternity")
hesin_maternity <- hesin_maternity[,c('eid','ins_index')] 

hesin_diag <- bio_hesin(project_dir, "diag")

#Extract codes for live births in the diagnosis dataset. 
hesin_birth <- hesin_diag[with(hesin_diag, grepl("Z370|Z371|Z372|Z375|Z38|Z390", paste(diag_icd9, diag_icd10), ignore.case=FALSE)),]
hesin_birth <- hesin_birth[,c('eid','ins_index')] 

#Anti join to see how many more cases are identified by using this hesin_birth variable
maternity <- anti_join(hesin_birth, hesin_maternity, by=c('eid','ins_index')) 

#Full join so that complete cases from both variables are identified
hesin_maternity <- full_join(hesin_birth, hesin_maternity, by=c('eid','ins_index'))

#Subset so only considering the instances where GP data found a code for possible postnatal depression. 
hesin_maternity <- subset(hesin_maternity, hesin_maternity$eid %in% E204$eid) 

#What's the breakdown of these noncases? 
maternity_dates <- bio_hesin(project_dir, "hesin")

#Inner join the two tables
maternity <- inner_join(hesin_maternity, maternity_dates, by=c('eid','ins_index'))
maternity <- maternity[,c('eid','ins_index','dsource','epistart','epiend','admidate')]

#For each hospital admission, calculate the time between admission date and postpartum depression. 

#Subset maternity so that it only contains people we are interested in potential postnatal cases
maternityE204 <- subset(maternity, maternity$eid %in% E204$eid)
maternityE204 <- left_join(maternityE204, E204, by='eid')
maternityE204$epistart <- as.Date(maternityE204$epistart,"%d/%m/%Y")
maternityE204$event_dt <- as.Date(maternityE204$event_dt,"%d/%m/%Y")

#Calculate time difference between each pregnancy date and the earliest diagnosis of the postnatal depression case.
maternityE204$TimeDiff <- time_length(difftime(maternityE204$event_dt, maternityE204$epistart), "years")

confirmedcases <- subset(maternityE204, TimeDiff >= 0 & TimeDiff <=1) 

#Take the earliest date of pregnancy 
confirmedcases <- confirmedcases %>%
                    group_by(eid) %>%
                    arrange(epistart) %>%
                    slice(1L)

confirmedcases <- as.data.frame(confirmedcases)
colnames(confirmedcases)[8] <- 'PPD_Date'

#Is this there first episode of depression - if not remove. 
DepressionGPCodes <- fread('DepressionGP', data.table=FALSE)
DepressionGPCodes <- DepressionGPCodes[,-1]

#Left join with cases 
confirmedcases <- left_join(confirmedcases, DepressionGPCodes, by='eid')
confirmedcases$difference <- time_length(difftime(confirmedcases$PPD_Date, confirmedcases$event_dt), "years")
confirmedcases <- subset(confirmedcases, is.na(difference) | difference <= 0) 

confirmedcases <- confirmedcases$eid

caseupdate <- append(cases,confirmedcases)
caseupdate <- unique(caseupdate) 
#write.csv(caseupdate, 'PostPartumConfirmedIDs')

#Identify pregnancies from GP records - need to know if the pregnancy actually happened so postnatal codes only to be considered...
gp_pregnancies <- gp_clinical[with(gp_clinical, grepl("6331.|6333.|6336.", paste(read_2, read_3), ignore.case=FALSE)),]

#Merge with E204 cases
gp_pregnancies <- subset(gp_pregnancies, gp_pregnancies$eid %in% E204$eid) 
colnames(E204)[3] <- 'PPD_Date'

gp_pregE204 <- left_join(gp_pregnancies, E204, by='eid')
gp_pregE204$PPD_Date <- as.Date(gp_pregE204$PPD_Date,"%d/%m/%Y")
gp_pregE204$event_dt <- as.Date(gp_pregE204$event_dt,"%d/%m/%Y")

gp_pregE204$difference <- time_length(difftime(gp_pregE204$PPD_Date, gp_pregE204$event_dt), "years")

confirmedGPcases <- subset(gp_pregE204, difference >= 0 & difference <= 1)

#Extract relevant columns only 
confirmedGPcases <- confirmedGPcases[,c('eid','PPD_Date')]

#Was it their first episode as evidenced from primary care data?
DepressionGPCodes <- fread('DepressionGP', data.table=FALSE)
DepressionGPCodes <- DepressionGPCodes[,-1]

#Left join with cases 
confirmedGPcases <- left_join(confirmedGPcases, DepressionGPCodes, by='eid')
confirmedGPcases$difference <- time_length(difftime(confirmedGPcases$PPD_Date, confirmedGPcases$event_dt), "years")
confirmedGPcases <- subset(confirmedGPcases, is.na(difference) | difference <= 0) #5 cases
confirmedGPcases <- unique(confirmedGPcases$eid)

caseupdate <- append(caseupdate,confirmedGPcases)
caseupdate <- unique(caseupdate) 

#write.csv(caseupdate, 'PostPartumConfirmedIDs') 

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#What happens if you just use depression GP codes as opposed to postpartum codes. 

#Read in all depression codes - only interested in first episode so only pull these out
DepressionGPCodes <- fread('DepressionGP', data.table=FALSE)
DepressionGPCodes <- DepressionGPCodes[,-1]

#HES data maternity
#Extract all cases from the maternity list where the ID corresponds to the IDs found in noncases. 
project_dir <- "ukb18177"
hesin_maternity <- bio_hesin(project_dir, "maternity")
hesin_maternity <- hesin_maternity[,c('eid','ins_index')] 

hesin_diag <- bio_hesin(project_dir, "diag")

#Extract codes for live births in the diagnosis dataset. 
hesin_birth <- hesin_diag[with(hesin_diag, grepl("Z370|Z371|Z372|Z375|Z38|Z390", paste(diag_icd9, diag_icd10), ignore.case=FALSE)),]
hesin_birth <- hesin_birth[,c('eid','ins_index')] 

#Full join so that complete cases from both variables are identified
hesin_maternity <- full_join(hesin_birth, hesin_maternity, by=c('eid','ins_index'))

#Subset so only considering the instances from maternity where GP data found a code for depression. 
hesin_maternity <- subset(hesin_maternity, hesin_maternity$eid %in% DepressionGPCodes$eid)

#Inner join the two tables
maternity <- inner_join(hesin_maternity, maternity_dates, by=c('eid','ins_index'))
maternity <- maternity[,c('eid','ins_index','dsource','epistart','epiend','admidate')]

#Subset maternity so that it only contains people we are interested in potential postnatal cases
maternityGP <- left_join(maternity, DepressionGPCodes, by='eid')
maternityGP$epistart <- as.Date(maternityGP$epistart,"%d/%m/%Y")
maternityGP$event_dt <- as.Date(maternityGP$event_dt,"%d/%m/%Y")

#Calculate time difference between each pregnancy date and the earliest diagnosis of the postnatal depression case.
maternityGP$TimeDiff <- time_length(difftime(maternityGP$event_dt, maternityGP$epistart), "years")

confirmedcases <- subset(maternityGP, TimeDiff >= 0 & TimeDiff <=1) 

#Take the earliest date of pregnancy 
confirmedcases <- confirmedcases %>%
  group_by(eid) %>%
  arrange(epistart) %>%
  slice(1L)

confirmedcases <- as.data.frame(confirmedcases) 

caseupdate <- append(caseupdate,confirmedcases) 
caseupdate <- unique(caseupdate) 

#Gp data maternity

#Identify pregnancies from GP records 
gp_pregnancies <- gp_clinical[with(gp_clinical, grepl("6331.|6333.|6336.", paste(read_2, read_3), ignore.case=FALSE)),]

#How many of these gp_pregnancies have depression codes. Only interested in these people. 
gp_pregnancies <- subset(gp_pregnancies, gp_pregnancies$eid %in% DepressionGPCodes$eid) 

colnames(DepressionGPCodes)[2] <- 'PPD_Date'
  
#Subset maternity so that it only contains people we are interested in potential postnatal cases
GPpregnancies <- left_join(gp_pregnancies, DepressionGPCodes, by='eid')
GPpregnancies$event_dt <- as.Date(GPpregnancies$event_dt,"%d/%m/%Y")

#Calculate time difference between each pregnancy date and the earliest diagnosis of the postnatal depression case.
GPpregnancies$TimeDiff <- time_length(difftime(GPpregnancies$PPD_Date, GPpregnancies$event_dt), "years")

confirmedcases <- subset(GPpregnancies, TimeDiff >= 0 & TimeDiff <=1) 

caseupdate <- append(caseupdate,confirmedcases) 
caseupdate <- unique(caseupdate) 

exclusion <- fread('ExclusionCriteriaApplied.csv', data.table = FALSE)

#Innerjoin to phenotypes dataset
case_after_exc <- subset(exclusion, exclusion$eid %in% caseupdate) 
case_after_exc <- case_after_exc$eid

write.csv(case_after_exc, './PostPartum/PostPartumEHR') 

#Subset so european only.
europeans <- fread('ukb18177_glanville_post_qc_id_list.fam' , data.table=FALSE)
colnames(europeans)[1] <- 'eid'

case_after_exc_QC <- subset(europeans, europeans$eid %in% case_after_exc)

#Remove columns relating to exclusion criteria
case_after_exc_QC <- case_after_exc_QC$eid
write.csv(case_after_exc_QC, './PostPartum/PostPartumEHR_QC')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Aggregating both definitions - All ancestries

MHQ <- fread('./PostPartum/PostPartumDepression_MHQ', data.table=FALSE)
EHR <- fread('./PostPartum/PostPartumEHR', data.table=FALSE)
colnames(EHR)[2] <- 'eid'

allPP <- full_join(MHQ, EHR, by='eid') 
allPP <- allPP$eid
write.csv(allPP, './PostPartum/PostPartumMHQandEHR') 

#Including both definitions - Europeans only

MHQ <- fread('./PostPartum/PostPartumDepression_MHQandQC', data.table=FALSE)
EHR <- fread('./PostPartum/PostPartumEHR_QC', data.table=FALSE)
colnames(EHR)[2] <- 'eid'

allPP <- full_join(MHQ, EHR, by='eid') 
write.csv(allPP, './PostPartum/PostPartumMHQandEHR_QC') 

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Regression

#MDD PRS

#PPD full

cases <- fread('./PostPartum/PostPartumMHQandEHR_QC', data.table=FALSE)
cases$PPD <- 1
cases <- cases[,c(3,8)]
colnames(cases)[1] <- 'eid'

controls <- fread('controlsPPD_with_selfreport_QC.csv', data.table=FALSE)
controls <- controls[,-1]

PPD <- rbind(cases, controls)

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('/mnt/lustre/groups/ukbiobank/ukb18177_glanville/2019_new_qc/ukb18177_glanville_covariates.txt' , data.table=FALSE)
colnames(covariates)[2] <- 'eid'
covariates <- covariates[,-1]

#Read in the PRS - thresholds of 0.05 and 1.
PRS <- fread('./PRS/DEPR06_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c('eid','DEPR06_Pt_0.05','DEPR06_Pt_1')]

#Merge to the postpartum depression dataset
PPD <- left_join(PPD, covariates, by='eid') %>%
          left_join(., PRS, by = 'eid')

#Standardise PRS
PPD$DEPR06_Pt_0.05 <- scale(PPD$DEPR06_Pt_0.05)
PPD$DEPR06_Pt_1 <- scale(PPD$DEPR06_Pt_1)

#Make batch and assessment centre factors.
PPD$assessment_centre <- factor(PPD$assessment_centre, ordered=FALSE)
PPD$batch <- factor(PPD$batch, ordered=FALSE)

scores<-c("DEPR06_Pt_0.05","DEPR06_Pt_1")

PPDresult <- data.frame()

for (j in 1:length(scores)){
  PPD$scores=PPD[,scores[j]]
    
  # create model
  # full model
  model<-glm(PPD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=PPD,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(PPD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=PPD,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,beta,p)
  print(output)
  PPDresult <- rbind(PPDresult, output)
} 

PPDresult$Pheno <- 'All' 

######################################################################################################################################################################################################################################

#Break results into MHQ and HES/GP definitions and review results. 

MHQcase <- fread('./PostPartum/PostPartumDepression_MHQandQC', data.table=FALSE)
MHQcase$PPD <- 1
MHQcase <- MHQcase[,c(2,6)]

#MHQ Self-report definition
MHQ <- rbind(MHQcase, controls)

#Merge to the postpartum depression dataset
MHQ <- left_join(MHQ, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
MHQ$DEPR06_Pt_0.05 <- scale(MHQ$DEPR06_Pt_0.05)
MHQ$DEPR06_Pt_1 <- scale(MHQ$DEPR06_Pt_1)

#Make batch and assessment centre factors.
MHQ$assessment_centre <- factor(MHQ$assessment_centre, ordered=FALSE)
MHQ$batch <- factor(MHQ$batch, ordered=FALSE)

#Compute association statistics on observation scale
scores<-c("DEPR06_Pt_0.05","DEPR06_Pt_1")

MHQresult <- data.frame()

for (j in 1:length(scores)){
  MHQ$scores=MHQ[,scores[j]]
  
  # create model
  # full model
  model<-glm(PPD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=MHQ,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(PPD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=MHQ,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,beta,p)
  print(output)
  MHQresult <- rbind(MHQresult, output)
} 

MHQresult$Pheno <- 'MHQ' 

######################################################################################################################################################################################################################################

#PostPartum cases from HES and GP definition

HESGPcase <- fread('./PostPartum/PostPartumEHR_QC', data.table=FALSE)
HESGPcase$PPD <- 1
HESGPcase <- HESGPcase[,c(2,3)]
colnames(HESGPcase)[1] <- 'eid'

HESGP <- rbind(HESGPcase, controls)

#Merge to the postpartum depression dataset
HESGP <- left_join(HESGP, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
HESGP$DEPR06_Pt_0.05 <- scale(HESGP$DEPR06_Pt_0.05)
HESGP$DEPR06_Pt_1 <- scale(HESGP$DEPR06_Pt_1)

#Make batch and assessment centre factors.
HESGP$assessment_centre <- factor(HESGP$assessment_centre, ordered=FALSE)
HESGP$batch <- factor(HESGP$batch, ordered=FALSE)

#Compute association statistics on the observation scale
scores<-c("DEPR06_Pt_0.05","DEPR06_Pt_1")

HESGPresult <- data.frame()

for (j in 1:length(scores)){
  HESGP$scores=HESGP[,scores[j]]
  
  # create model
  # full model
  model<-glm(PPD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=HESGP,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(PPD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=HESGP,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,beta,p)
  print(output)
  HESGPresult <- rbind(HESGPresult, output)
} 

HESGPresult$Pheno <- 'EHR' 

######################################################################################################################################################################################################################################

#Merge all results together
result <- bind_rows(PPDresult, MHQresult, HESGPresult)
write.csv(result, 'PPD_MDDPRSresults_withselfreportbirth')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Bipolar PRS
#PPD full

cases <- fread('./PostPartum/PostPartumMHQandEHR_QC', data.table=FALSE)
cases$PPD <- 1
cases <- cases[,c(3,8)]
colnames(cases)[1] <- 'eid'

controls <- fread('controlsPPD_with_selfreport_QC.csv', data.table=FALSE)
controls <- controls[,-1]

PPD <- rbind(cases, controls)

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('ukb18177_glanville_covariates.txt' , data.table=FALSE)
colnames(covariates)[2] <- 'eid'
covariates <- covariates[,-1]

#Read in the PRS - thresholds of 0.01 and 1.
PRS <- fread('./PRS/BIPO02_header.all.score', data.table=FALSE)
colnames(PRS)[2] <- 'eid'
PRS <- PRS[,c('eid','BIPO02_Pt_0.01','BIPO02_Pt_1')]

#Merge to the postpartum depression dataset
PPD <- left_join(PPD, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
PPD$BIPO02_Pt_0.01 <- scale(PPD$BIPO02_Pt_0.01)
PPD$BIPO02_Pt_1 <- scale(PPD$BIPO02_Pt_1)

#Make batch and assessment centre factors.
PPD$assessment_centre <- factor(PPD$assessment_centre, ordered=FALSE)
PPD$batch <- factor(PPD$batch, ordered=FALSE)

scores<-c("BIPO02_Pt_0.01","BIPO02_Pt_1")

PPDresult <- data.frame()

for (j in 1:length(scores)){
  PPD$scores=PPD[,scores[j]]
  
  # create model
  # full model
  model<-glm(PPD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=PPD,family=binomial,na.action=na.exclude)
  # null model without PRS  
  modelNULL<-glm(PPD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=PPD,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,beta,p)
  print(output)
  PPDresult <- rbind(PPDresult, output)
} 

PPDresult$Pheno <- 'All' 
######################################################################################################################################################################################################################################

#Break results into MHQ and HES/GP definitions and review results. 

MHQcase <- fread('./PostPartum/PostPartumDepression_MHQandQC', data.table=FALSE)
MHQcase$PPD <- 1
MHQcase <- MHQcase[,c(2,6)]

#MHQ Self-report definition
MHQ <- rbind(MHQcase, controls)

#Merge to the postpartum depression dataset
MHQ <- left_join(MHQ, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
MHQ$BIPO02_Pt_0.01 <- scale(MHQ$BIPO02_Pt_0.01)
MHQ$BIPO02_Pt_1 <- scale(MHQ$BIPO02_Pt_1)

#Make batch and assessment centre factors.
MHQ$assessment_centre <- factor(MHQ$assessment_centre, ordered=FALSE)
MHQ$batch <- factor(MHQ$batch, ordered=FALSE)

#Compute association statistics on observation scale
scores<-c("BIPO02_Pt_0.01","BIPO02_Pt_1")

MHQresult <- data.frame()

for (j in 1:length(scores)){
  MHQ$scores=MHQ[,scores[j]]
  
  # create model
  # full model
  model<-glm(PPD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=MHQ,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(PPD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=MHQ,family=binomial,na.action=na.exclude)
  # get results
  #print(summary(model)$coefficients)
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,beta,p)
  print(output)
  MHQresult <- rbind(MHQresult, output)
} 

MHQresult$Pheno <- 'MHQ' 

######################################################################################################################################################################################################################################

#PostPartum cases from HES and GP definition
HESGPcase <- fread('./PostPartum/PostPartumEHR_QC', data.table=FALSE)
HESGPcase$PPD <- 1
colnames(HESGPcase)[2] <- 'eid'
HESGPcase <- HESGPcase[,c(2,3)]

HESGP <- rbind(HESGPcase, controls)

#Merge to the postpartum depression dataset
HESGP <- left_join(HESGP, covariates, by='eid') %>%
  left_join(., PRS, by = 'eid')

#Standardise PRS
HESGP$BIPO02_Pt_0.01 <- scale(HESGP$BIPO02_Pt_0.01)
HESGP$BIPO02_Pt_1 <- scale(HESGP$BIPO02_Pt_1)

#Make batch and assessment centre factors.
HESGP$assessment_centre <- factor(HESGP$assessment_centre, ordered=FALSE)
HESGP$batch <- factor(HESGP$batch, ordered=FALSE)

#Compute association statistics on the observation scale
scores<-c("BIPO02_Pt_0.01","BIPO02_Pt_1")

HESGPresult <- data.frame()

for (j in 1:length(scores)){
  HESGP$scores=HESGP[,scores[j]]
  
  # create model
  # full model
  model<-glm(PPD~scores+PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
             data=HESGP,family=binomial,na.action=na.exclude)
  # null model without PRS
  modelNULL<-glm(PPD~PC1+PC2+PC3+PC4+PC5+PC6+batch+assessment_centre,
                 data=HESGP,family=binomial,na.action=na.exclude)
  # get results
  Pt=scores[j]
  beta=(summary(model)$coefficients["scores",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["scores",2])
  r2=NagelkerkeR2(model)$R2-NagelkerkeR2(modelNULL)$R2
  p=(summary(model)$coefficients["scores",4])
  output=data.frame(Pt,OR,se,beta,p)
  print(output)
  HESGPresult <- rbind(HESGPresult, output)
} 

HESGPresult$Pheno <- 'EHR' 

######################################################################################################################################################################################################################################

#Merge all results together

result <- bind_rows(PPDresult,MHQresult,HESGPresult)
write.csv(result, 'PPD_BPDPRSresults_withselfreport')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Environmental variable analysis
#UKB Kings to extract the post-partum question
#Extract Years of education, neuroticism, townsend deprivation index, trauma variables
project_dir <- "/scratch/datasets/ukbiobank/ukb18177"
f <- bio_field(project_dir)

#Select data
f %>%
  select(field, name) %>%
  filter(str_detect(name, "f6138_0_0|f6138_0_1|f6138_0_2|f6138_0_3|f6138_0_4|f6138_0_5|f189_0_0|f20127|f20489|f20488|f20487|f20490|f20491|f20522|f20523|f20521|f20524|f20525|f20531|f20529|f20526|f20530|f20528|f20527|f20107|f20110")) %>%
  bio_field_add("EnvironmentalVariables.txt")

system("cat EnvironmentalVariables.txt")

#Read data into file
bio_phen(
  project_dir,
  field = "EnvironmentalVariables.txt",
  out = "EnvironmentalVariables"
)

system("ls -lh EnvironmentalVariables.rds")

#Get raw data for family history

######################################################################################################################################################################################################################################

#Read in Environmental Variables - already includes PCA for trauma from prior script. 
environment <- readRDS('EnvironmentalVariables.rds') 

#Family History Variables - extract from environmental variables
fh <- environment[,c(1,20:103)] 

fh$FamHist <- as.numeric(apply(fh[,-1], 1, function(r) any(r %in% 12)))

#Assign missingness such that an individual responds don't know or PNTA and provides no other information
fh$MissingCheck <- apply(fh[,-c(1,86)], 1, function(r) any(r %in% c(-11,-13,-21,-23)) & !any(r %in% c(1:13,-17,-27)))

fh$FamHist <- case_when(
    fh$FamHist==1 ~ 1,
    fh$FamHist==0 & fh$MissingCheck==1 ~ NA_real_,
    TRUE ~ 0
  )

fh <- fh[,c('eid','FamHist')]

environment <- full_join(environment, fh, by='eid')

#Education Attainment
education <- environment[,c('eid','educ_one','educ_two','educ_three','educ_four','educ_five','educ_six')]
educ <- c('educ_one','educ_two','educ_three','educ_four','educ_five','educ_six')

for(i in educ){
  education[,paste(i)] <- as.numeric(education[,paste(i)])
  education[,paste(i,'_factored',sep="")] <- factor(
         ifelse(is.na(education[,paste(i)]), "missing",
                ifelse(!is.na(education[,paste(i)]) & education[,paste(i)] == -7, "none",
                       ifelse(!is.na(education[,paste(i)]) & education[,paste(i)] == 1, "degree",
                              ifelse(!is.na(education[,paste(i)]) & (education[,paste(i)] == 2 | education[,paste(i)] == 5 | education[,paste(i)] == 6), "alevels",
                                     ifelse(!is.na(education[,paste(i)]) & education[,paste(i)] > 2 & education[,paste(i)] < 5, "gcse",
                                            ifelse(!is.na(education[,paste(i)]) & education[,paste(i)] == -3, "pnta", "missing")))))),
    levels=c("none", "gcse", "alevels", "degree", "pnta", "missing"))
}

education$best_educ <- case_when(
  education$educ_one_factored == 'degree' | education$educ_two_factored == 'degree' | education$educ_three_factored == 'degree' | education$educ_four_factored == 'degree' | education$educ_five_factored == 'degree' | education$educ_six_factored == 'degree' ~ 'degree',
  education$educ_one_factored == 'alevels' | education$educ_two_factored == 'alevels' | education$educ_three_factored == 'alevels' | education$educ_four_factored == 'alevels' | education$educ_five_factored == 'alevels' | education$educ_six_factored == 'alevels' ~ 'alevels',
  education$educ_one_factored == 'gcse' | education$educ_two_factored == 'gcse' | education$educ_three_factored == 'gcse' | education$educ_four_factored == 'gcse' | education$educ_five_factored == 'gcse' | education$educ_six_factored == 'gcse' ~ 'gcse',
  education$educ_one_factored == 'none' | education$educ_two_factored == 'none' | education$educ_three_factored == 'none' | education$educ_four_factored == 'none' | education$educ_five_factored == 'none' | education$educ_six_factored == 'none' ~ 'none',
  education$educ_one_factored == 'pnta' | education$educ_two_factored == 'pnta' | education$educ_three_factored == 'pnta' | education$educ_four_factored == 'pnta' | education$educ_five_factored == 'pnta' | education$educ_six_factored == 'pnta' ~ 'pnta',
  education$educ_one_factored == 'missing' | education$educ_two_factored == 'missing' | education$educ_three_factored == 'missing' | education$educ_four_factored == 'missing' | education$educ_five_factored == 'missing' | education$educ_six_factored == 'missing' ~ 'missing'
)

degree <- sum(table(education$educ_one_factored)['degree'], table(education$educ_two_factored)['degree'], table(education$educ_three_factored)['degree'], table(education$educ_four_factored)['degree'], table(education$educ_five_factored)['degree'], table(education$educ_six_factored)['degree'])
alevels <- sum(table(education$educ_one_factored)['alevels'], table(education$educ_two_factored)['alevels'], table(education$educ_three_factored)['alevels'], table(education$educ_four_factored)['alevels'], table(education$educ_five_factored)['alevels'], table(education$educ_six_factored)['alevels'])
gcse <- sum(table(education$educ_one_factored)['gcse'], table(education$educ_two_factored)['gcse'], table(education$educ_three_factored)['gcse'], table(education$educ_four_factored)['gcse'], table(education$educ_five_factored)['gcse'], table(education$educ_six_factored)['gcse'])
none <- sum(table(education$educ_one_factored)['none'], table(education$educ_two_factored)['none'], table(education$educ_three_factored)['none'], table(education$educ_four_factored)['none'], table(education$educ_five_factored)['none'], table(education$educ_six_factored)['none'])
pnta <- sum(table(education$educ_one_factored)['pnta'], table(education$educ_two_factored)['pnta'], table(education$educ_three_factored)['pnta'], table(education$educ_four_factored)['pnta'], table(education$educ_five_factored)['pnta'], table(education$educ_six_factored)['pnta'])

total <- sum(degree,alevels,gcse,none,pnta)

education$best_educ <- as.factor(education$best_educ)
education$best_educ <- relevel(education$best_educ, ref='degree')
education <- education[,c('eid','best_educ')]

environment <- full_join(environment, education, by='eid')

#Include year of birth for cohort effect
yob <- fread('SelfReportAAO_Depression', data.table=FALSE)
yob <- yob[,c(2,4)]

environment <- full_join(environment, yob, by='eid')
environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-c(1,2,127)]
write.csv(environment, 'EnvironmentalVariables.csv')

##########################################################################################################
##########################################################################################################

#Analysis
#Read in cases and controls for PPD - get an idea of degree of missingness for each variable
environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

cases <- fread('./PostPartum/PostPartumMHQandEHR', data.table=FALSE)
cases$PPD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

controls <- fread('controlsPPD_with_selfreport.csv', data.table=FALSE)
controls <- controls[,-1]

PPD <- rbind(cases, controls)

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

PPD <- inner_join(PPD, covariates, by='eid') %>%
  inner_join(., environment)

#Most of the missingness driven by controls - of the two case definitions used - is one contributing to missingness more than the other?
env_variables <- c('TDI','neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

for(i in env_variables){
  print(i)
  full <- PPD[is.na(PPD[[i]]),]
  case <- subset(full, PPD==1)
  print(length(full$eid))
  print(length(case$eid))
}

#All
table(PPD$education_attainment)

#Cases
case <- subset(PPD, PPD==1)
table(case$education_attainment)

#Controls
control <- subset(PPD, PPD==0)
table(control$education_attainment)

#Remove missing cases of SES and education attainment and see level of missingness of environmental variable
PPD <- subset(PPD, !(education_attainment == "missing" | education_attainment == "pnta" | is.na(TDI)))

PPD$assessment_centre <- factor(PPD$assessment_centre, ordered=FALSE)
PPD$education_attainment <- factor(PPD$education_attainment, ordered=FALSE)
PPD$education_attainment <- relevel(PPD$education_attainment, ref='degree')

#PPD$TDI <- scale(PPD$TDI)
#PPD$neuroticism <- scale(PPD$neuroticism) 

env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

for(i in env_variables){
  print(i)
  variable <- PPD[is.na(PPD[[i]]),]
  variabletwo <- subset(variable, PPD==1)
  print(length(variable$eid))
  print(length(variabletwo$eid))
}

#Regression of other environmental variables not SES or education attainment - do this separately as slightly different model

env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

PPDEnvResults <- data.frame()

for(j in 1:length(env_variables)) {
  
  PPD$variable=PPD[,env_variables[j]]
  
  print(env_variables[j])
  
  # create model
  
  # full model
  model<-glm(PPD ~ variable + assessment_centre + education_attainment + TDI + YOB,
             data=PPD,family=binomial,na.action=na.exclude)
  
  # get results
  print(summary(model)$coefficients)
  EnvVariable=env_variables[j]
  beta=(summary(model)$coefficients["variable",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["variable",2])
  p=(summary(model)$coefficients["variable",4])
  output=data.frame(EnvVariable,beta,se,OR,p)
  print(output)
  PPDEnvResults <- rbind(PPDEnvResults, output)
} 

#Perform separate regressions with just educational attainment and SES
model <- glm(PPD ~ education_attainment + assessment_centre + TDI + YOB,
             data=PPD,family=binomial,na.action=na.exclude)

beta=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),1])
OR=exp(beta)
se=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),2])
p=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),4])
EnvVariable <- c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI")
output=data.frame(EnvVariable,beta,se,OR,p)
print(output)
PPDEnvResults <- rbind(PPDEnvResults, output)

write.csv(PPDEnvResults, 'EnvironmentalAssociationsPostPartum_FullSample_withselfreport.csv')

##########################################################################################################
##########################################################################################################

#Repeat the analysis with just MHQ

#Read in cases and controls for PPD - get an idea of degree of missingness for each variable
environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

cases <- fread('./PostPartum/PostPartumDepression_MHQ', data.table=FALSE)
cases$PPD <- 1
cases <- cases[,c(2,4)]
colnames(cases)[1] <- 'eid'

controls <- fread('controlsPPD_with_selfreport.csv', data.table=FALSE)
controls <- controls[,-1]

PPD <- rbind(cases, controls)

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

PPD <- inner_join(PPD, covariates, by='eid') %>%
  inner_join(., environment)

#Remove missing cases of SES and education attainment and see level of missingness of environmental variable
PPD <- subset(PPD, !(education_attainment == "missing" | education_attainment == "pnta" | is.na(TDI)))

PPD$assessment_centre <- factor(PPD$assessment_centre, ordered=FALSE)
PPD$education_attainment <- factor(PPD$education_attainment, ordered=FALSE)
PPD$education_attainment <- relevel(PPD$education_attainment, ref='degree')

PPD$TDI <- scale(PPD$TDI)
PPD$neuroticism <- scale(PPD$neuroticism) 

#Regression of other environmental variables not SES or education attainment - do this separately as slightly different model

env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

PPDEnvResults <- data.frame()

for(j in 1:length(env_variables)) {
  
  PPD$variable=PPD[,env_variables[j]]
  
  print(env_variables[j])
  
  # create model
  
  # full model
  model<-glm(PPD ~ variable + assessment_centre + education_attainment + TDI + YOB,
             data=PPD,family=binomial,na.action=na.exclude)
  
  # get results
  print(summary(model)$coefficients)
  EnvVariable=env_variables[j]
  beta=(summary(model)$coefficients["variable",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["variable",2])
  p=(summary(model)$coefficients["variable",4])
  output=data.frame(EnvVariable,beta,se,OR,p)
  print(output)
  PPDEnvResults <- rbind(PPDEnvResults, output)
} 

#Perform separate regressions with just educational attainment and SES
model <- glm(PPD ~ education_attainment + assessment_centre + TDI + YOB,
             data=PPD,family=binomial,na.action=na.exclude)

beta=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),1])
OR=exp(beta)
se=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),2])
p=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),4])
EnvVariable <- c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI")
output=data.frame(EnvVariable,beta,se,OR,p)
print(output)
PPDEnvResults <- rbind(PPDEnvResults, output)

write.csv(PPDEnvResults, 'EnvironmentalAssociationsPostPartum_MHQOnly_withselfreport.csv')

##########################################################################################################
##########################################################################################################

#Repeat the environmental analysis with just EHR definitions

#Read in cases and controls for PPD - get an idea of degree of missingness for each variable
environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

cases <- fread('./PostPartum/PostPartumEHR', data.table=FALSE)
cases$PPD <- 1
cases <- cases[,c(2,3)]
colnames(cases)[1] <- 'eid'

controls <- fread('controlsPPD_with_selfreport.csv', data.table=FALSE)
controls <- controls[,-1]

PPD <- rbind(cases, controls)

covariates <- readRDS('assessment_centre.rds')
colnames(covariates)[2] <- 'assessment_centre'

PPD <- inner_join(PPD, covariates, by='eid') %>%
  inner_join(., environment)

#Remove missing cases of SES and education attainment and see level of missingness of environmental variable
PPD <- subset(PPD, !(education_attainment == "missing" | education_attainment == "pnta" | is.na(TDI)))

PPD$assessment_centre <- factor(PPD$assessment_centre, ordered=FALSE)
PPD$education_attainment <- factor(PPD$education_attainment, ordered=FALSE)
PPD$education_attainment <- relevel(PPD$education_attainment, ref='degree')

PPD$TDI <- scale(PPD$TDI)
PPD$neuroticism <- scale(PPD$neuroticism) 

#Regression of other environmental variables not SES or education attainment - do this separately as slightly different model

env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

PPDEnvResults <- data.frame()

for(j in 1:length(env_variables)) {
  
  PPD$variable=PPD[,env_variables[j]]
  
  print(env_variables[j])
  
  # create model
  
  # full model
  model<-glm(PPD ~ variable + assessment_centre + education_attainment + TDI + YOB,
             data=PPD,family=binomial,na.action=na.exclude)
  
  # get results
  print(summary(model)$coefficients)
  EnvVariable=env_variables[j]
  beta=(summary(model)$coefficients["variable",1])
  OR=exp(beta)
  se=(summary(model)$coefficients["variable",2])
  p=(summary(model)$coefficients["variable",4])
  output=data.frame(EnvVariable,beta,se,OR,p)
  print(output)
  PPDEnvResults <- rbind(PPDEnvResults, output)
} 

#Perform separate regressions with just educational attainment and SES
model <- glm(PPD ~ education_attainment + assessment_centre + TDI + YOB,
             data=PPD,family=binomial,na.action=na.exclude)

beta=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),1])
OR=exp(beta)
se=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),2])
p=(summary(model)$coefficients[c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI"),4])
EnvVariable <- c("education_attainmentalevels","education_attainmentgcse","education_attainmentnone","TDI")
output=data.frame(EnvVariable,beta,se,OR,p)
print(output)
PPDEnvResults <- rbind(PPDEnvResults, output)

write.csv(PPDEnvResults, 'EnvironmentalAssociationsPostPartum_EHROnly_withselfreport.csv')
