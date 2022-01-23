###################################################################################################################
################  GLOSSARY OF CODES - ALL FIELDS REQUIRED FOR exclusion CRITERIA      #############################
###################################################################################################################
### UK Biobank Code     Descriptive                                                                             ###
### f.eid	              ID										                                           col number = 1
### f.20446.0.0		      Ever.had.prolonged.feelings.of.sadness.or.depression			   		 col number = 12963
### f.20499.0.0		      Ever.sought.or.received.professional.help.for.mental.distress	   col_number = 12998
### f.20126.0		        Bipolar.and.major.depression.status                              col_number = 7684
### f.20544.0		        Mental.health.problems.ever.diagnosed.by.a.professional 	       col_number = 13039 - 13054
### f.20002.0		        Non.cancer.illness.code.self.reported							               col_number = 5793 - 5821
### f.41202.0		        Diagnoses.main.ICD10								   	                         col_number = 271 - 650
### f.41204.0		        Diagnoses.secondary.ICD10								                         col_number = 679 - 1113
###################################################################################################################

#Bash

#Get the individuals who have self-reported use of anti-psychotics
#awk '{print $1, $13, $15}' ./softlinks/identified_classes.onlyclass > drugs

#R

library(sqldf)
library(dplyr)
library(ukbkings)
library(tidyverse)
library(lubridate)
library(tidyverse)
library(data.table)
library(ggpubr)
library(png)
library(grid)
library(qdapRegex)

#Read in the datasets
exclusion <- fread('icd10_treatment_srinterview_Depressed_Ever_Phenotype_File.txt', data.table=FALSE)

antipsychotics <- fread('drugs', data.table=FALSE)
colnames(antipsychotics)[1] <- 'ID'
antipsychotics <- antipsychotics[,c(1,3)]

#Create all necessary variables

#####Physician Diagnosed Disorders

# To the question: "Have you been diagnosed with one or more of the following", there are 16 possible responses (i.e. array=16,  see MHQ pdf for details). Therefore, an individual can report up to 16 diagnoses. So for each individual, there are 16 columns of data for this question, such that each column contains an integer representing a possible response (e.g. Social Phobia =1) or NA. 

exclusion$SRSchizophrenia<-with(exclusion, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
                                                  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 2) |
                                                           (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 2), 1, 0)))

exclusion$SRPsychosisOther<-with(exclusion, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
                                                   ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 3) |
                                                            (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 3), 1, 0)))

# The next line of code is creating a column for SRPsychosisAny based on the last two columns created, i.e. if participant has scz OR psychosis(other), they will be a case in this column. 

exclusion$SRPsychosisAny<-with(exclusion, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
                                                 ifelse((!is.na(SRSchizophrenia) & SRSchizophrenia == 1) | (!is.na(SRPsychosisOther) & SRPsychosisOther == 1), 1, 0)))

exclusion$SRManiaBIP<-with(exclusion, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
                                             ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 10) |
                                                      (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 10), 1, 0)))


#####Drug-derived Diagnoses######

###Antipsychotics
exclusion <- full_join(x=exclusion, y=antipsychotics, by="ID")

###Mania and mood stabilisers
bpdcodes <- c('1140867504','1140867494','1140867498','1140867520','1140917270','1140867490')

exclusion$Lithium <- apply(exclusion[,grep("Treatment.medication.code.0.", colnames(exclusion))], 1, function(row) any(row %in% bpdcodes))
exclusion$Lithium <- as.integer(as.logical(exclusion$Lithium))

###Antidepressants at baseline only
antidepcodes <- c('1140879616', '1140921600', '1140879540', '1140867878', '1140916282', '1140909806', '1140867888', '1141152732', '1141180212', '1140879634', '1140867876', '1140882236', '1141190158',
                  '1141200564', '1140867726', '1140879620', '1140867818', '1140879630', '1140879628', '1141151946', '1140867948', '1140867624', '1140867756', '1140867884', '1141151978', '1141152736',
                  '1141201834', '1140867690', '1140867640', '1140867920', '1140867850', '1140879544', '1141200570', '1140867934', '1140867758', '1140867914', '1140867820', '1141151982', '1140882244',
                  '1140879556', '1140867852', '1140867860', '1140917460', '1140867938', '1140867856', '1140867922', '1140910820', '1140882312', '1140867944', '1140867784', '1140867812', '1140867668',
                  '1140867940')

exclusion$Antidepressants <- apply(exclusion[,grep("Treatment.medication.code.0.", colnames(exclusion))], 1, function(row) any(row %in% antidepcodes))
exclusion$Antidepressants <- as.integer(as.logical(exclusion$Antidepressants))

exclusion$anyDrug <- case_when(
      exclusion$Antidepressants == 1 | exclusion$Antipsychotics == 1 | exclusion$Lithium == 1 ~ 1,
      TRUE ~ 0 
      )

###Smith Depression Diagnosis 2013 - As this is not at baseline - do not remove any depression cases. 

# if column for Bipolar.and.major.depression.status is NA, print NA. Else, if entry is not NA and the score is either 1 or 2, mark individual as 1 for case, else 0 for control. The UKB coding for this variable is:
# 0	No Bipolar or Depression
# 1	Bipolar I Disorder
# 2 Bipolar II Disorder
# 3	Probable Recurrent major depression (severe)
# 4	Probable Recurrent major depression (moderate)
# 5	Single Probable major depression episode

exclusion$SmithMood<-with(exclusion, ifelse(is.na(Bipolar.and.major.depression.status.0.0), NA,
                                            ifelse(!is.na(Bipolar.and.major.depression.status.0.0) & (Bipolar.and.major.depression.status.0.0 == 1) |
                                                     !is.na(Bipolar.and.major.depression.status.0.0) & (Bipolar.and.major.depression.status.0.0 == 2), 1, 0)))

#Extract relevant columns only and remove rest
relevant <- c('ID', 'SRPsychosisAny','SRManiaBIP', 'SmithMood', 'anyDrug')
exclusion <- exclusion[,relevant]

#####Interview Diagnoses At Baseline######

#Breakdown of exclusion by disorder (note there is going to be overlap between cases as they are not mutually exclusive). 
#1289 = Schizophrenia 
#1291 = Bipolar 
#1408 = Alcohol dependency 
#1409 = Opioid dependency 
#1410 = Other dependency 
#1286 = Depression 

NonCancerSR <- readRDS('SelfReportMedicalIllnessNonCancer.rds')
NonCancerSR <- NonCancerSR[,c(1,3:138)]
colnames(NonCancerSR)[1] <- 'ID' 

#All apart from depression can be determined throughout all timepoints. 
NonCancerSR$InterviewExclusionOther <- apply(NonCancerSR, 1, function(row) any(row %in% c("1289","1291","1408","1409","1410")))
NonCancerSR$InterviewExclusionOther <- as.integer(as.logical(NonCancerSR$InterviewExclusionOther))

#Depression at initial baseline assessment
NonCancerSR$InterviewExclusionDep <-apply(NonCancerSR[,grep("20002-0.", colnames(NonCancerSR))], 1, function(row) any(row %in% "1286"))
NonCancerSR$InterviewExclusionDep <- as.integer(as.logical(NonCancerSR$InterviewExclusionDep))

NonCancerSR <- NonCancerSR[,c('ID','InterviewExclusionDep','InterviewExclusionOther')]

exclusion <- full_join(exclusion, NonCancerSR, by='ID')

####HES Inpatient Diagnoses######

# similarly to the above, applies function to search relevant ICD10 codes. Cells are TRUE if a participant has reported the ICD10 code for the column, else FALSE. TRUE = 598, FALSE = 5821907.
#Only removing schizophrenia, bipolar and substance abuse as ICD codes for depression will be explored later. Will form part of the diagnosis.

ICDFields<-c("F10","F100","F101","F102","F103","F104","F105","F106","F107","F108","F109",
             "F11","F110","F111","F112","F113","F114","F115","F116","F117","F118","F119",
             "F12","F120","F121","F122","F123","F124","F125","F126","F127","F128","F129",
             "F13","F130","F131","F132","F133","F134","F135","F136","F137","F138","F139",
             "F14","F140","F141","F142","F143","F144","F145","F146","F147","F148","F149",
             "F15","F150","F151","F152","F153","F154","F155","F156","F157","F158","F159",
             "F16","F160","F161","F162","F163","F164","F165","F166","F167","F168","F169",
             "F18","F180","F181","F182","F183","F184","F185","F186","F187","F188","F189",
             "F19","F190","F191","F192","F193","F194","F195","F196","F197","F198","F199",
             "F20","F200","F201","F202","F203","F204","F205","F206","F207","F208","F209",
             "F21",
             "F22","F228","F229",
             "F23","F230","F231","F230","F231","F232","F233","F238","F239",
             "F24",
             "F250","F251","F252","F258","F259",
             "F28",
             "F29",
             "F300","F301","F302","F308","F309",
             "F310","F311","F312","F313","F314","F315","F316","F317","F318","F319")

project_dir <- "ukb18177"

hesin_diag <- bio_hesin(project_dir, "diag")

#Extract codes for live births in the diagnosis dataset. 
HES <- subset(hesin_diag, diag_icd10 %in% ICDFields)

#Get a list of unique IDs. 
HES.IDs <- unique(HES$eid) #12765

#Make a field for the HES
exclusion$HES <- case_when(
  exclusion$ID %in% HES.IDs ~ 1,
  TRUE ~ 0
)

#Extract treatment seeking for nerves, anxiety or depression (psychiatrist/GP) and age at recruitment from Kens Package (ALL AT BASELINE!)
project_dir <- "ukb18177"
f <- bio_field(project_dir)

#Select data
f %>%
  select(field, name) %>%
  filter(str_detect(field, "2100-0.0|2090-0.0|21022-0.0")) %>%
  bio_field_add("treatment_seeking_and_age.txt")

#Read data into file
bio_phen(
  project_dir,
  field = "treatment_seeking_and_age.txt",
  out = "treatment_seeking_and_age"
)

#Read in the dataset 
treatment_seeking <- readRDS("treatment_seeking_and_age.rds")
treatment_seeking <- treatment_seeking[,-5]
colnames(treatment_seeking) <- c('ID','GPTreatment','psychiatristTreatment','AgeatRecruitment')

#Take a union of psychiatrist/GP treatment_seeking to condense the two columns
treatment_seeking$TreatmentSeeking <- case_when(
  treatment_seeking$psychiatristTreatment == 1 | treatment_seeking$GPTreatment == 1 ~ 1,
  TRUE ~ 0
)

#Merge with current exclusion criteria dataset
exclusion <- full_join(exclusion, treatment_seeking, by='ID')

#Extract relevant columns only and remove rest
relevant <- c('ID', 'InterviewExclusionDep', 'InterviewExclusionOther', 'SRPsychosisAny','SRManiaBIP', 'HES', 'SmithMood', 'anyDrug', 'TreatmentSeeking', 'AgeatRecruitment')
exclusion <- exclusion[,relevant]

#Make a field such that it is the total of all exclusion criteria - anything above 0 equates to a removal - ignore NAs in this function, i.e. treated as 0.
exclusion$total <- rowSums(exclusion[,2:9], na.rm=TRUE)

exclusion <- exclusion[exclusion$total==0,c('ID','AgeatRecruitment')]

write.csv(exclusion, 'HealthyAtBaseline')

#Merge this dataset with the disorders dataset without depression. 

#Search for depression diagnoses within HES data
#Depression ICD CODE EXTRACTION

DepressionCodes <- read.csv.sql(file="/scratch/groups/ukbiobank/KCL_Data/Phenotypes/Full_Dataset_August_2017/24072020/hesin_diag.txt",
                            sql="SELECT eid,ins_index,arr_index,level,diag_icd10
                            FROM file
                            WHERE diag_icd10 = 'F32'
                            OR diag_icd10 = 'F320'
                            OR diag_icd10 = 'F321'
                            OR diag_icd10 = 'F322'
                            OR diag_icd10 = 'F323'
                            OR diag_icd10 = 'F328'
                            OR diag_icd10 = 'F329'
                            OR diag_icd10 = 'F33'
                            OR diag_icd10 = 'F330'
                            OR diag_icd10 = 'F331'
                            OR diag_icd10 = 'F332'
                            OR diag_icd10 = 'F333'
                            OR diag_icd10 = 'F334'
                            OR diag_icd10 = 'F338'
                            OR diag_icd10 = 'F339'",
                            header=TRUE,
                            sep = "\t")

#Extract date of depression as well as treatment and main speciality of consultant making the diagnosis
DepressionDates <- read.csv.sql(file="hesin.txt", 
                            sql="SELECT eid,ins_index,dsource,epistart,admidate,mainspef_uni,mainspef,tretspef_uni,tretspef
                            FROM file
                            WHERE eid in (
                            SELECT eid
                            FROM DepressionCodes
                            )",
                            header=TRUE,
                            sep = "\t")

#Inner join the two tables
Depression <- inner_join(DepressionCodes, DepressionDates, by=c('eid','ins_index'))

#Create a new variable which is the epistart normally and admidate if epistart is missing
Depression$Date.of.Depression <- ifelse(is.na(Depression$epistart)|Depression$epistart=="", Depression$admidate, Depression$epistart)
Depression$Disorder <- 'Depression'

#As a test - where epistart is NA see if it has taken the date to be admidate
head(Depression[Depression$epistart=="",],10)
head(Depression[is.na(Depression$Date.of.Depression),],10)

Depression$epistart <- as.Date(Depression$epistart,"%d/%m/%Y")
Depression$admidate <- as.Date(Depression$admidate,"%d/%m/%Y")
Depression$Date.of.Depression <- as.Date(Depression$Date.of.Depression,"%d/%m/%Y")

Depression %>%
  arrange(desc(Date.of.Depression)) %>%
  slice(1L)

#Select the earliest date from the Depression admissions
Depression <- Depression %>%
  group_by(eid) %>%
  arrange(Date.of.Depression) %>%
  slice(1L)

#Select relevant variables only
Depression <- Depression[,c("eid","level","diag_icd10","Date.of.Depression","Disorder","mainspef_uni","tretspef_uni")]
Depression <- as.data.frame(Depression)
colnames(Depression)[2] <- c("icd10_Depression")

write.csv(Depression, "DepressionICD")
print("DepressionICD Finished")

Depression[is.na(Depression$Date.of.Depression),]

#Search for depresion diagnosis within GP data
#Cleaned GP data found here. 
DepressionGPCodes <- read.csv.sql(file="GP_psychiatric_codes_extracted_18022020_SPH.txt",
                                  sql="SELECT eid,event_dt,type
                                  FROM file
                                  WHERE type = 'depression'",
                                  header=TRUE,
                                  sep = "\t")

#Reformat the GP Dates as date type. 
DepressionGPCodes$event_dt <- as.Date(DepressionGPCodes$event_dt,"%Y-%m-%d")

#Some reports don't have an event date attached. In this case remove the entry as will mess up the earliest date otherwise. 
DepressionGPCodes <- subset(DepressionGPCodes, !is.na(event_dt))

#Select the earliest date from the Depression admissions
DepressionGPCodes <- DepressionGPCodes %>%
  group_by(eid) %>%
  arrange(event_dt) %>%
  slice(1L)

write.csv(DepressionGPCodes, "DepressionGP")
print("DepressionGP Finished")

#Merge the two datasets (full_join) and take the earliest diagnosis
Depression <- fread('DepressionICD', data.table=FALSE)
Depression <- Depression[,-1]
DepressionGPCodes <- fread('DepressionGP', data.table=FALSE)
DepressionGPCodes <- DepressionGPCodes[,-1]

AllDepression <- full_join(DepressionGPCodes, Depression, by='eid') 

AllDepression$First.Depression.Date <- with(AllDepression,
                                            case_when(!is.na(event_dt) & event_dt < Date.of.Depression |
                                                      !is.na(event_dt) & is.na(Date.of.Depression) ~ event_dt, 
                                                      TRUE ~ Date.of.Depression)
                                            )

AllDepression$First.Depression.Date <- as.Date(AllDepression$First.Depression.Date, origin="1970-01-01")

#Remove cases which have an ICD diagnosis before the GP diagnosis as implies may have received a diagnosis prior without us knowing. 
AllDepression <- subset(AllDepression, event_dt==First.Depression.Date | mainspef_uni == 1020 | mainspef_uni == 1270 | mainspef_uni == 1600 | mainspef_uni == 1100 | 
                        tretspef_uni == 1060 | tretspef_uni == 1300 | tretspef_uni == 1470 | tretspef_uni == 1720 | tretspef_uni == 1940 | 
                          tretspef_uni == 2370) 

NAs <- subset(AllDepression, is.na(First.Depression.Date))

zerochar <- subset(AllDepression, First.Depression.Date == "")

#Merge to dataset that contains date of diagnosis for every disorder. 
Disorders <- fread('./Diseases/Disorders', data.table=FALSE)
Disorders <- Disorders[,-1]

Disorder_incDep <- inner_join(AllDepression, Disorders, by='eid')

#Calculate years between diagnosis of physical disorder and depression. 
disorder <- c("Stroke","Diabetes","Cancer","Autoimmune","MyocardialInfarction","MultipleSclerosis","MotorNeuroneDisease","Epilepsy")
for(i in disorder) {
  Disorder_incDep[,paste(i,'years_dep', sep="_")] <- time_length(difftime(Disorder_incDep$First.Depression.Date,Disorder_incDep[,paste('Date.of.',i, sep="")]), "years")
}

#Merge with people who were healthy at baseline
HealthyatBaseline <- fread('HealthyAtBaseline', data.table=FALSE)
HealthyatBaseline <- HealthyatBaseline[,-1]
colnames(HealthyatBaseline)[1] <- 'eid'

Healthy.Disorder.Dep <- inner_join(Disorder_incDep, HealthyatBaseline, by='eid') 

#Add the date of birth variable 
dep <- fread('SelfReportAAO_Depression', data.table=FALSE, select=c('eid','YOB','MOB'))

#Make the month of birth and year of birth into one column - assume the birthdate is the start of the month
dep$MOB <- as.character(dep$MOB)
dep$YOB <- as.character(dep$YOB)

dep <- dep %>% 
  unite(DOB, MOB, YOB, sep = "-", remove = FALSE)
  dep$DOB <- paste("01", dep$DOB, sep="-")
  dep$DOB <- as.Date(dep$DOB,"%d-%m-%Y")

dep <- dep[,c('eid','DOB')]

Healthy.Disorder.Dep <- inner_join(Healthy.Disorder.Dep,dep, by='eid')

#Age at onset of depression to the nearest integer
Healthy.Disorder.Dep$DepressionAAO <- floor(time_length(difftime(Healthy.Disorder.Dep$First.Depression.Date, Healthy.Disorder.Dep$DOB), "years"))

#If depression diagnosis AAO is greater than or equal to age at recruitment then can assume this is first episode of depression. 
Healthy.Disorder.Dep <- subset(Healthy.Disorder.Dep, AgeatRecruitment < DepressionAAO) 

#Calculate sample sizes based on time from 1 to 5 years. 

#5 year time difference

#How many people become depression between 0 and 5 years for each disorder - then take the union of this set
year_diff <- c('Stroke','Diabetes','Cancer','Autoimmune','MyocardialInfarction','MultipleSclerosis','MotorNeuroneDisease','Epilepsy')
for(i in year_diff){
  Healthy.Disorder.Dep[,paste(i,"ltfive",sep="_")] <- case_when(
    Healthy.Disorder.Dep[,paste(i,"_years_dep",sep="")] <= 5 & Healthy.Disorder.Dep[,paste(i,"_years_dep",sep="")] >= 0 ~ 1,
    TRUE ~ 0
  )
}

#Calculate sample sizes for each disorder
for(i in year_diff){
  print(table(Healthy.Disorder.Dep[,paste(i,"ltfive",sep="_")]))
}

#Take the union of the set for total case number
Healthy.Disorder.Dep$Case_ltfive <- case_when(
  Healthy.Disorder.Dep$Stroke_ltfive == 1 | Healthy.Disorder.Dep$Diabetes_ltfive == 1 | Healthy.Disorder.Dep$Cancer_ltfive == 1 | Healthy.Disorder.Dep$Autoimmune_ltfive == 1 |
    Healthy.Disorder.Dep$MyocardialInfarction_ltfive == 1 | Healthy.Disorder.Dep$MultipleSclerosis_ltfive == 1 | Healthy.Disorder.Dep$MotorNeuroneDisease_ltfive == 1 | Healthy.Disorder.Dep$Epilepsy_ltfive == 1 ~ 1,
  TRUE ~ 0
)

table(Healthy.Disorder.Dep$Case_ltfive)
head(subset(Healthy.Disorder.Dep, Case_ltfive==1),10)

#4 year time difference

#How many people become depression between 0 and 4 years for each disorder - then take the union of this set
year_diff <- c('Stroke','Diabetes','Cancer','Autoimmune','MyocardialInfarction','MultipleSclerosis','MotorNeuroneDisease','Epilepsy')
for(i in year_diff){
  Healthy.Disorder.Dep[,paste(i,"ltfour",sep="_")] <- case_when(
    Healthy.Disorder.Dep[,paste(i,"_years_dep",sep="")] <= 4 & Healthy.Disorder.Dep[,paste(i,"_years_dep",sep="")] >= 0 ~ 1,
    TRUE ~ 0
  )
}

#Calculate sample sizes for each disorder
for(i in year_diff){
  print(table(Healthy.Disorder.Dep[,paste(i,"ltfour",sep="_")]))
}

#Take the union of the set for total case number
Healthy.Disorder.Dep$Case_ltfour <- case_when(
  Healthy.Disorder.Dep$Stroke_ltfour == 1 | Healthy.Disorder.Dep$Diabetes_ltfour == 1 | Healthy.Disorder.Dep$Cancer_ltfour == 1 | Healthy.Disorder.Dep$Autoimmune_ltfour == 1 |
    Healthy.Disorder.Dep$MyocardialInfarction_ltfour == 1 | Healthy.Disorder.Dep$MultipleSclerosis_ltfour == 1 | Healthy.Disorder.Dep$MotorNeuroneDisease_ltfour == 1 | Healthy.Disorder.Dep$Epilepsy_ltfour == 1 ~ 1,
  TRUE ~ 0
)

table(Healthy.Disorder.Dep$Case_ltfour)

#3 year time difference

#How many people become depression between 0 and 4 years for each disorder - then take the union of this set
year_diff <- c('Stroke','Diabetes','Cancer','Autoimmune','MyocardialInfarction','MultipleSclerosis','MotorNeuroneDisease','Epilepsy')
for(i in year_diff){
  Healthy.Disorder.Dep[,paste(i,"ltthree",sep="_")] <- case_when(
    Healthy.Disorder.Dep[,paste(i,"_years_dep",sep="")] <= 3 & Healthy.Disorder.Dep[,paste(i,"_years_dep",sep="")] >= 0 ~ 1,
    TRUE ~ 0
  )
}

#Calculate sample sizes for each disorder
for(i in year_diff){
  print(table(Healthy.Disorder.Dep[,paste(i,"ltthree",sep="_")]))
}

#Take the union of the set for total case number
Healthy.Disorder.Dep$Case_ltthree <- case_when(
  Healthy.Disorder.Dep$Stroke_ltthree == 1 | Healthy.Disorder.Dep$Diabetes_ltthree == 1 | Healthy.Disorder.Dep$Cancer_ltthree == 1 | Healthy.Disorder.Dep$Autoimmune_ltthree == 1 |
    Healthy.Disorder.Dep$MyocardialInfarction_ltthree == 1 | Healthy.Disorder.Dep$MultipleSclerosis_ltthree == 1 | Healthy.Disorder.Dep$MotorNeuroneDisease_ltthree == 1 | Healthy.Disorder.Dep$Epilepsy_ltthree == 1 ~ 1,
  TRUE ~ 0
)

table(Healthy.Disorder.Dep$Case_ltthree)

#2 year time difference

#How many people become depression between 0 and 4 years for each disorder - then take the union of this set
year_diff <- c('Stroke','Diabetes','Cancer','Autoimmune','MyocardialInfarction','MultipleSclerosis','MotorNeuroneDisease','Epilepsy')
for(i in year_diff){
  Healthy.Disorder.Dep[,paste(i,"lttwo",sep="_")] <- case_when(
    Healthy.Disorder.Dep[,paste(i,"_years_dep",sep="")] <= 2 & Healthy.Disorder.Dep[,paste(i,"_years_dep",sep="")] >= 0 ~ 1,
    TRUE ~ 0
  )
}

#Calculate sample sizes for each disorder
for(i in year_diff){
  print(table(Healthy.Disorder.Dep[,paste(i,"lttwo",sep="_")]))
}

#Take the union of the set for total case number
Healthy.Disorder.Dep$Case_lttwo <- case_when(
  Healthy.Disorder.Dep$Stroke_lttwo == 1 | Healthy.Disorder.Dep$Diabetes_lttwo == 1 | Healthy.Disorder.Dep$Cancer_lttwo == 1 | Healthy.Disorder.Dep$Autoimmune_lttwo == 1 |
    Healthy.Disorder.Dep$MyocardialInfarction_lttwo == 1 | Healthy.Disorder.Dep$MultipleSclerosis_lttwo == 1 | Healthy.Disorder.Dep$MotorNeuroneDisease_lttwo == 1 | Healthy.Disorder.Dep$Epilepsy_lttwo == 1 ~ 1,
  TRUE ~ 0
)

table(Healthy.Disorder.Dep$Case_lttwo)

#1 year time difference

#How many people become depression between 0 and 4 years for each disorder - then take the union of this set
year_diff <- c('Stroke','Diabetes','Cancer','Autoimmune','MyocardialInfarction','MultipleSclerosis','MotorNeuroneDisease','Epilepsy')
for(i in year_diff){
  Healthy.Disorder.Dep[,paste(i,"ltone",sep="_")] <- case_when(
    Healthy.Disorder.Dep[,paste(i,"_years_dep",sep="")] <= 1 & Healthy.Disorder.Dep[,paste(i,"_years_dep",sep="")] >= 0 ~ 1,
    TRUE ~ 0
  )
}

#Calculate sample sizes for each disorder
for(i in year_diff){
  print(table(Healthy.Disorder.Dep[,paste(i,"ltone",sep="_")]))
}

#Take the union of the set for total case number
Healthy.Disorder.Dep$Case_ltone <- case_when(
  Healthy.Disorder.Dep$Stroke_ltone == 1 | Healthy.Disorder.Dep$Diabetes_ltone == 1 | Healthy.Disorder.Dep$Cancer_ltone == 1 | Healthy.Disorder.Dep$Autoimmune_ltone == 1 |
    Healthy.Disorder.Dep$MyocardialInfarction_ltone == 1 | Healthy.Disorder.Dep$MultipleSclerosis_ltone == 1 | Healthy.Disorder.Dep$MotorNeuroneDisease_ltone == 1 | Healthy.Disorder.Dep$Epilepsy_ltone == 1 ~ 1,
  TRUE ~ 0
)

table(Healthy.Disorder.Dep$Case_ltone)

#Extract cases within the year
oneyearcasesICDandGP <- subset(Healthy.Disorder.Dep, Case_ltone==1)

#Assign all temporal triggers of major depression as if one cannot be identified by 
cause <- apply(oneyearcasesICDandGP[,c('Stroke_ltone','Diabetes_ltone','Cancer_ltone','Autoimmune_ltone','MyocardialInfarction_ltone','MultipleSclerosis_ltone','MotorNeuroneDisease_ltone','Epilepsy_ltone')] == 1 , 1, function(x) names(which(x)))

#Split into multiple trigger columns
oneyearcasesICDandGP$cause <- cause

#Turn dataframe into a tibble
oneyearcasesICDandGP <- oneyearcasesICDandGP %>%
  separate(col = cause, into = LETTERS[1:8], sep = ", ")

#Rename columns
colnames(oneyearcasesICDandGP)[107:114] <- c('Trigger1','Trigger2','Trigger3','Trigger4','Trigger5','Trigger6','Trigger7','Trigger8')

for(i in colnames(oneyearcasesICDandGP)[107:114]){
  oneyearcasesICDandGP[,paste(i)] <- gsub('_ltone', '', oneyearcasesICDandGP[,paste(i)] )
}

#Remove quotations and brackets
for(i in 1:8) {
  oneyearcasesICDandGP[[paste('Trigger',i,'intermediate',sep="")]] <- rm_between(oneyearcasesICDandGP[[paste('Trigger',i,sep="")]], '"', '"', extract=TRUE)
  oneyearcasesICDandGP[[paste('Trigger',i,'intermediate',sep="")]]  <- as.character(oneyearcasesICDandGP[[paste('Trigger',i,'intermediate',sep="")]])
  oneyearcasesICDandGP[[paste('Trigger',i,'intermediate',sep="")]] <- ifelse(is.na(oneyearcasesICDandGP[[paste('Trigger',i,'intermediate',sep="")]]), oneyearcasesICDandGP[[paste('Trigger',i,sep="")]], oneyearcasesICDandGP[[paste('Trigger',i,'intermediate',sep="")]])
  oneyearcasesICDandGP[[paste('Trigger',i,sep="")]] <- oneyearcasesICDandGP[[paste('Trigger',i,'intermediate',sep="")]]
}

write.csv(oneyearcasesICDandGP, 'incident_depression')