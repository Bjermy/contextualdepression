#Controls
library(data.table)
library(dplyr)
library(sqldf)
library(tidyverse)
library(lubridate)
library(ukbkings)

##Requirement - Must have experienced the disorder - match the proportion within the case definition. Random selection from within the subgroup will be ok. 

#Take the healthy at baseline people - these haven't reported depression but we also need to remove people who have ever reported ICD or GP depression. The inverse join will do well here. 
HealthyAtBaseline <- fread('HealthyAtBaseline', data.table=FALSE)
HealthyAtBaseline <- HealthyAtBaseline[,-1]
colnames(HealthyAtBaseline)[1] <- 'eid'

#Exclusion criteria - taken from script LifetimeDepressionMHQDefinition.R
exclusion <- fread('ExclusionCriteriaApplied.csv', data.table=FALSE)
exclusion <- exclusion[,-1]

HealthyAtBaseline <- inner_join(x=HealthyAtBaseline, y=exclusion, by="eid") 

NonCancerSR <- readRDS('SelfReportMedicalIllnessNonCancer.rds')

#Only concerned with diagnosis not year
NonCancerSR <- NonCancerSR[,c(1,3:138)]

HealthyAtBaseline <- left_join(HealthyAtBaseline, NonCancerSR, by='eid')

#Instances of depression will arrive after baseline - remove these cases
HealthyAtBaseline$InterviewDep <- apply(HealthyAtBaseline, 1, function(r) any(r %in% 1286))
HealthyAtBaseline$InterviewDep <- as.integer(as.logical(HealthyAtBaseline$InterviewDep))

Healthy_noDep <- subset(HealthyAtBaseline, InterviewDep==0)

DepressionGP <- fread('DepressionGP', data.table=FALSE)
DepressionGP <- DepressionGP[,-1]

DepressionICD <- fread('DepressionICD', data.table=FALSE)
DepressionICD <- DepressionICD[,-1]

#Remove people who have reported an ICD code or GP code of depression
Healthy_noDep <- anti_join(Healthy_noDep, DepressionGP, by='eid') %>%
                  anti_join(., DepressionICD, by='eid')

#Remove people who have not been identified according to the Smith et al definition

#Read in the datasets
exclusion <- fread('icd10_treatment_srinterview_Depressed_Ever_Phenotype_File.txt', data.table=FALSE)
colnames(exclusion)[1] <- 'eid'
exclusion$SmithMood <- with(exclusion, ifelse(is.na(Bipolar.and.major.depression.status.0.0), NA,
                                            ifelse(!is.na(Bipolar.and.major.depression.status.0.0) & (Bipolar.and.major.depression.status.0.0 == 3) |
                                                     !is.na(Bipolar.and.major.depression.status.0.0) & (Bipolar.and.major.depression.status.0.0 == 4) |
                                                       !is.na(Bipolar.and.major.depression.status.0.0) & (Bipolar.and.major.depression.status.0.0 == 5), 1, 0)))

exclusion_depressed <- subset(exclusion, SmithMood == 1)

Healthy_noDep <- anti_join(Healthy_noDep, exclusion_depressed, by='eid')

#Remove people who do not report lifetime depression - if you're not a case you're a control. 
LifetimeDepression <- fread('LifetimeDepression.csv', data.table=FALSE)
LifetimeDepression <- subset(LifetimeDepression, LifetimeDepression == 1)

Healthy_noDep <- anti_join(Healthy_noDep, LifetimeDepression, by='eid')

#Remove any postnatal depression cases
PostPartum <- fread('./PostPartum/PostPartumMHQandEHR', data.table=FALSE)
colnames(PostPartum)[2] <- 'eid'

Healthy_noDep <- anti_join(Healthy_noDep, PostPartum, by='eid')

#Make sure the cases have not suffered from one of the disorders considered. 

##Removal Step 1: Hospitalisations - Anti-join on existing disorder list.
Disorders <- fread('./Diseases/Disorders', data.table=FALSE)

MDDexcControls <- anti_join(Healthy_noDep, Disorders, by='eid')

##Removal Step 2: Self-reported any of the disorders at nurse interview. 
CancerSR <- readRDS('SelfReportMedicalIllnessCancer.rds')

#Only concerned with diagnosis not year
CancerSR <- CancerSR[,c(1,3:26)]

MDDexcControls <- left_join(MDDexcControls, CancerSR, by='eid')

#If the individual hasn't ever reported a cancer the total NAs for each person will equal the total number of columns in which the person had the opportunity to report which is 24
MDDexcControls$CancerNumber <- apply(MDDexcControls[,c(11:34)], 1, function(x) sum(is.na(x)))

MDDexcControls <- subset(MDDexcControls, CancerNumber==24)

MDDexcControls <- MDDexcControls[,c(1:2)]

NonCancerSR <- readRDS('SelfReportMedicalIllnessNonCancer.rds')

#Only concerned with diagnosis not year
NonCancerSR <- NonCancerSR[,c(1,3:138)]

MDDexcControls <- left_join(MDDexcControls, NonCancerSR, by='eid')

MDDexcControls$NonCancer <- apply(MDDexcControls, 1, function(r) any(r %in% c(1286,1264,1222,1223,1220,1081,1425,1491,1086,1075,1093,1261,1259,1311,1456,1154,1459,1463,1331,1453,1464,1382,1226,1462,1531,1377,1428,1082,1437,1313,1381,1461,1376,1521,1076,1083,1260,1345,1583)))

MDDexcControls <- subset(MDDexcControls, NonCancer==FALSE)
MDDexcControls <- unique(MDDexcControls[,1])
write.csv(MDDexcControls, 'controlsMDDExcPPDandDisorders.csv') #76097

#Get the fam file and inner join to work out people of european ancestries
europeans <- fread('ukb18177_glanville_post_qc_id_list.fam' , data.table=FALSE)
colnames(europeans)[1] <- 'eid'

MDDexcControls <- as.data.frame(MDDexcControls)
colnames(MDDexcControls) <- 'eid'

MDDexcControls <- inner_join(MDDexcControls, europeans, by='eid')
controls <- MDDexcControls[,1] 

write.csv(controls, 'controlsMDDExcPPDandDisorders_QC.csv')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Medical Disorder Controls

#Read the Disorders dataset
Disorders <- fread('Disorders', data.table=FALSE)
Disorders <- Disorders[,-1]

disorderduplicates <- Disorders %>%
                        count(eid) %>%
                          filter(n >= 2)

#Merge the two datasets - inner_join as must have the disorder
control <- inner_join(Healthy_noDep, Disorders, by='eid')

#Must not have died as a result of the disorder 
DeathDates <- read.csv.sql(file="/scratch/groups/ukbiobank/KCL_Data/Phenotypes/Full_Dataset_August_2017/24072020/death.txt", 
                            sql="SELECT eid,date_of_death
                            FROM file
                            WHERE eid in (
                            SELECT eid
                            FROM control
                            )",
                            header=TRUE,
                            sep = "\t")

deathduplicates <- DeathDates %>%
                       count(eid) %>%
                           filter(n >= 2)

Duplicates <- subset(DeathDates, eid %in% deathduplicates$eid)

DeathDates <- DeathDates %>% distinct(eid, .keep_all = TRUE)

controlWithDeath <- left_join(control, DeathDates, by='eid')
controlWithDeath$date_of_death <- as.Date(controlWithDeath$date_of_death,"%d/%m/%Y")

#Calculate the years after diagnosis that somebody died. If it is less than a year after then do not consider this individual a control for this particular disorder. 
disorder <- c("Stroke","Diabetes","Cancer","Autoimmune","MyocardialInfarction","MultipleSclerosis","MotorNeuroneDisease","Epilepsy")
for(i in disorder) {
  controlWithDeath[,paste('Death.After.', i, sep="")] <- time_length(difftime(controlWithDeath$date_of_death,controlWithDeath[,paste('Date.of.',i, sep="")]), "years")
}

#Things required to be considered a control - You have had the disorder for over a year and did not die within the year of having the diagnosis. 
disorder <- c("Stroke","Diabetes","Cancer","Autoimmune","MyocardialInfarction","MultipleSclerosis","MotorNeuroneDisease","Epilepsy")
for(i in disorder) {
  controlWithDeath[,paste('Control.For.', i, sep="")] <- case_when(
    (!is.na(controlWithDeath[,paste('Date.of.', i, sep="")]) & controlWithDeath[,paste('Date.of.', i, sep="")] < '2019-03-31') & is.na(controlWithDeath[,paste('Death.After.', i, sep="")]) | controlWithDeath[,paste('Death.After.', i, sep="")] > 1 ~ 1,
    TRUE ~ 0
  )
}

DeathDates$date_of_death <- as.Date(DeathDates$date_of_death,"%d/%m/%Y")
DeathDates %>%
  arrange(desc(date_of_death)) %>%
  slice(1L)

#Keep only the columns where it tells us if they can be considered a control for a particular disorder
controlWithDeath <- controlWithDeath[,c("eid","Date.of.Stroke","Date.of.Diabetes","Date.of.Cancer","Date.of.Autoimmune","Date.of.MyocardialInfarction","Date.of.MultipleSclerosis","Date.of.MotorNeuroneDisease","Date.of.Epilepsy",
                                        "Control.For.Stroke","Control.For.Diabetes","Control.For.Cancer","Control.For.Autoimmune","Control.For.MyocardialInfarction","Control.For.MultipleSclerosis","Control.For.MotorNeuroneDisease","Control.For.Epilepsy")]

controls <- subset(controlWithDeath, Control.For.Stroke == 1 | Control.For.Diabetes == 1 | Control.For.Cancer == 1 | Control.For.Autoimmune == 1 | Control.For.MyocardialInfarction == 1 | Control.For.MultipleSclerosis == 1 | Control.For.MotorNeuroneDisease == 1 | Control.For.Epilepsy == 1)

sum(duplicated(controls$eid))
             
write.csv(controls, 'controlsMedicalDisorder.csv') 

#How many pass Genetic QC
controlWithDeath <- fread('controlsMedicalDisorder.csv', data.table=FALSE)

#Get the fam file and inner join to work out people of european ancestries
europeans <- fread('ukb18177_glanville_post_qc_id_list.fam' , data.table=FALSE)
colnames(europeans)[1] <- 'eid'

controlWithDeath.euro <- inner_join(controlWithDeath, europeans)
controlWithDeath.euro <- controlWithDeath.euro[,-c(1,19:23)] 
write.csv(controlWithDeath.euro, 'controlsMedicalDisorder_QC.csv')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#PostPartum Depression Controls

#Have given birth - live birth
project_dir <- "/scratch/datasets/ukbiobank/ukb18177"
hesin_maternity <- bio_hesin(project_dir, "maternity")
hesin_maternity <- hesin_maternity[,c('eid','ins_index')] #30907

hesin_diag <- bio_hesin(project_dir, "diag")

#Extract codes for live births in the diagnosis dataset. 
hesin_birth <- hesin_diag[with(hesin_diag, grepl("Z370|Z371|Z372|Z375|Z38|Z390", paste(diag_icd9, diag_icd10), ignore.case=FALSE)),]
hesin_birth <- hesin_birth[,c('eid','ins_index')] 

#Full join so that complete cases from both variables are identified
hesin_maternity <- full_join(hesin_birth, hesin_maternity, by=c('eid','ins_index'))

DatesLabour <- read.csv.sql(file="/scratch/groups/ukbiobank/KCL_Data/Phenotypes/Full_Dataset_August_2017/24072020/hesin.txt", 
                              sql="SELECT eid,ins_index,epistart
                              FROM file
                              WHERE eid in (
                              SELECT eid
                              FROM hesin_maternity
                              )",
                            header=TRUE,
                            sep = "\t")

hesin_maternity <- inner_join(hesin_maternity, DatesLabour, by=c('eid','ins_index'))

colnames(hesin_maternity)[3] <- 'Date.of.Labour'

sex <- readRDS('sex.rds')

hesin_maternity <- left_join(hesin_maternity, sex, by='eid')
colnames(hesin_maternity)[4] <- 'sex'

hesin_maternity <- subset(hesin_maternity, sex==0)

#Didn't die within the year of having birth or as a result of the birth. 
DeathDates <- read.csv.sql(file="/scratch/groups/ukbiobank/KCL_Data/Phenotypes/Full_Dataset_August_2017/24072020/death.txt", 
                           sql="SELECT eid,date_of_death
                           FROM file
                           WHERE eid in (
                           SELECT eid
                           FROM hesin_maternity
                           )",
                            header=TRUE,
                           sep = "\t")

deathduplicates <- DeathDates %>%
  count(eid) %>%
  filter(n >= 2)

#Remove duplicates
DeathDates <- DeathDates %>% distinct(eid, .keep_all = TRUE)

controlWithDeath <- left_join(hesin_maternity, DeathDates, by='eid')
controlWithDeath$date_of_death <- as.Date(controlWithDeath$date_of_death,"%d/%m/%Y")
controlWithDeath$Date.of.Labour <- as.Date(controlWithDeath$Date.of.Labour,"%d/%m/%Y")

#Calculate the years after diagnosis that somebody died. If it is less than a year after then do not consider this individual a control for this particular disorder. 
controlWithDeath$Death.After.Labour <- time_length(difftime(controlWithDeath$date_of_death,controlWithDeath$Date.of.Labour), "years")

#Things required to be considered a control - You have had the disorder for over a year and did not die within the year of having the diagnosis. 
controlWithDeath$Control.For.Labour <- case_when(
    (!is.na(controlWithDeath$Date.of.Labour) & controlWithDeath$Date.of.Labour < '2019-03-31' & is.na(controlWithDeath$Death.After.Labour)) | controlWithDeath$Death.After.Labour > 1 ~ 1,
    TRUE ~ 0
  )

controlPPD <- subset(controlWithDeath, Control.For.Labour == 1)

controlPPD <- controlPPD[,c(1,2)] 

#Merge with the set of participants who have been found to not report major depressive disorder previously. 
controlPPD <- inner_join(Healthy_noDep, controlPPD, by='eid') 
controlPPD <- controlPPD %>% distinct(eid, .keep_all = TRUE) 

DeathDates$date_of_death <- as.Date(DeathDates$date_of_death,"%d/%m/%Y")
DeathDates %>%
  arrange(desc(date_of_death)) %>%
  slice(1L)

write.csv(controlPPD[,c(1,2)], 'controlsPPD.csv') 

#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################

#Also add in people who self-report giving birth
childbirth <- readRDS('ChildbirthSelfReport.rds')
colnames(childbirth) <- c('eid','numofchildren_one','numofchildren_two','numofchildren_three','numofchildren_four')

childbirth$children <- apply(childbirth[,2:5], 1, function(r) any(r > 0))

sex <- readRDS('sex.rds')

childbirth <- left_join(childbirth, sex, by='eid')
colnames(childbirth)[7] <- 'sex'

childbirth <- subset(childbirth, children==TRUE)
childbirthcontrol <- inner_join(Healthy_noDep, childbirth, by='eid') #116747

#Merge with current control set
controlPPD <- fread('controlsPPD.csv', data.table=FALSE)

childbirthcontrol <- full_join(childbirthcontrol, controlPPD, by='eid')
childbirthcontrol$PPD <- 0
childbirthcontrol <- childbirthcontrol[,c('eid','PPD')]

write.csv(childbirthcontrol, 'controlsPPD_with_selfreport.csv')

#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################

#How many pass Genetic QC
controlPPD <- fread('controlsPPD.csv', data.table=FALSE)

#Get the fam file and inner join to work out people of european ancestries
europeans <- fread('ukb18177_glanville_post_qc_id_list.fam' , data.table=FALSE)
colnames(europeans)[1] <- 'eid'

controlPPD.euro <- inner_join(controlPPD, europeans, by='eid')
controlPPD.euro <- controlPPD.euro[,c(2,3)] 
controlPPD.euro <- controlPPD.euro %>% distinct(eid, .keep_all = TRUE)

write.csv(controlPPD.euro, 'controlsPPD_QC.csv')

#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################

#How many pass Genetic QC when including self-report controls
controlPPD <- fread('controlsPPD_with_selfreport.csv', data.table=FALSE)

#Get the fam file and inner join to work out people of european ancestries
europeans <- fread('ukb18177_glanville_post_qc_id_list.fam' , data.table=FALSE)
colnames(europeans)[1] <- 'eid'

controlPPD.euro <- inner_join(controlPPD, europeans, by='eid')
controlPPD.euro <- controlPPD.euro[,c(2,3)] # controls using QC on full sample.
controlPPD.euro <- controlPPD.euro %>% distinct(eid, .keep_all = TRUE) 

write.csv(controlPPD.euro, 'controlsPPD_with_selfreport_QC.csv')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
