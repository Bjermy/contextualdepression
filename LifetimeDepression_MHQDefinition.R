###################################################################################################################
##################  GLOSSARY OF CODES - ALL FIELDS REQUIRED TO DEFINE PHENOTYPES       ############################
###################################################################################################################
### UK Biobank Code     Descriptive                                                                             ###
### eid	        ID										                                            col number = 1
### 20435.0.0		Difficulty.concentrating.during.worst.depression					        col number = 12954
### 20436.0.0		Fraction.of.day.affected.during.worst.episode.of.depression				col number = 12955
### 20437.0.0		Thoughts.of.death.during.worst.depression						              col number = 12956
### 20438.0.0		Duration.of.worst.depression								                      col number = 12957
### 20439.0.0		Frequency.of.depressed.days.during.worst.episode.of.depression		col number = 12958		
### 20440.0.0		Impact.on.normal.roles.during.worst.period.of.depression				  col number = 12959
### 20441.0.0		Ever.had.prolonged.loss.of.interest.in.normal.activities				  col number = 12960
### 20442.0.0		Lifetime.number.of.depressed.periods							                col number = 12961
### 20446.0.0		Ever.had.prolonged.feelings.of.sadness.or.depression			   		  col number = 12963
### 20449.0.0		Feelings.of.tiredness.during.worst.episode.of.depression				  col number = 12966
### 20450.0.0		Feelings.of.worthlessness.during.worst.period.of.depression				col number = 12967
### 20532.0.0		Did.your.sleep.change									                            col number = 13027
### 20536.0.0		Weight.change.during.worst.episode.of.depression					        col number = 13031

###################################################################################################################

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

#Extract these variables using UKB Kings

#Extract Year of Birth and Age of Depression Onset from Kens Package
project_dir <- "/scratch/datasets/ukbiobank/ukb18177_glanville"
f <- bio_field(project_dir)

#Select data
f %>%
  select(field, name) %>%
  filter(str_detect(field, "20446.0.0|20441.0.0|20532.0.0|20435.0.0|20536.0.0|20437.0.0|20449.0.0|20450.0.0|20438.0.0|20439.0.0|20440.0.0|20442.0.0|20436.0.0")) %>%
  bio_field_add("LifetimeDepression.txt")

system("cat LifetimeDepression.txt")

#Read data into file
bio_phen(
  project_dir,
  field = "LifetimeDepression.txt",
  out = "LifetimeDepression"
)

system("ls -lh LifetimeDepression.rds")

#Pull down depressed ever phenotype - COPIED OVER WITHOUT THE USE OF R

#Read in the depression datasets and 
MHQ <- readRDS("LifetimeDepression.rds")

colnames(MHQ) <- c('eid', 'Difficulty.concentrating.during.worst.depression', 'Fraction.of.day.affected.during.worst.episode.of.depression', 'Thoughts.of.death.during.worst.depression', 'Duration.of.worst.depression',
                   'Frequency.of.depressed.days.during.worst.episode.of.depression','Impact.on.normal.roles.during.worst.period.of.depression','Ever.had.prolonged.loss.of.interest.in.normal.activities',
                   'Lifetime.number.of.depressed.periods','Ever.had.prolonged.feelings.of.sadness.or.depression','Feelings.of.tiredness.during.worst.episode.of.depression','Feelings.of.worthlessness.during.worst.period.of.depression',
                   'Did.your.sleep.change','Weight.change.during.worst.episode.of.depression')

#Remove NAs such that only individuals who have responded to the MHQ are retained. Sample size = 157358  
MHQ<-MHQ[!is.na(MHQ$Ever.had.prolonged.feelings.of.sadness.or.depression), ]

#Check there are no duplicated IDs 
sum(duplicated(MHQ$eid)) #=0 so that's good. 

#NOTE: When creating these phenotypes - the number of the phenotype in the field correspond to the phenotypes specified within the appendix of the pre-registration document

#Phenotype 1: Cardinal symptoms assessed only - All participants either endorse depressed mood or anhedonia get a 1, otherwise classed as 0. Prefer not to answers where the response to the second question in 
#the set is a 0 or another prefer not to answer are reclassified as NA. Checked against data showcase this should return a phenotype of 89,040 endorsements. It did.  
MHQ$Pheno1 <- case_when(
  MHQ$Ever.had.prolonged.feelings.of.sadness.or.depression == 0 & MHQ$Ever.had.prolonged.loss.of.interest.in.normal.activities == 0 ~ 0,
  MHQ$Ever.had.prolonged.feelings.of.sadness.or.depression == 1 | MHQ$Ever.had.prolonged.loss.of.interest.in.normal.activities == 1 ~ 1,
  TRUE ~ 999
  )

#Recode components such that NAs are changed to 999

#Recode the variables so they are suitable for the subsequent phenotypes, i.e. make them binary - Rules followed in accordance with the pre-registration
MHQ$BINARY_Impact.on.normal.roles.during.worst.period.of.depression <- case_when(
  MHQ$Impact.on.normal.roles.during.worst.period.of.depression == 2 | MHQ$Impact.on.normal.roles.during.worst.period.of.depression == 3 ~ 1,
  MHQ$Impact.on.normal.roles.during.worst.period.of.depression == 0 | MHQ$Impact.on.normal.roles.during.worst.period.of.depression == 1 ~ 0,
  TRUE ~ 999
  )

MHQ$BINARY_Fraction.of.day.affected.during.worst.episode.of.depression <- case_when(
  MHQ$Fraction.of.day.affected.during.worst.episode.of.depression == 3 | MHQ$Fraction.of.day.affected.during.worst.episode.of.depression == 4 ~ 1,
  MHQ$Fraction.of.day.affected.during.worst.episode.of.depression == 1 | MHQ$Fraction.of.day.affected.during.worst.episode.of.depression == 2 ~ 0,
  TRUE ~ 999
)

MHQ$BINARY_Frequency.of.depressed.days.during.worst.episode.of.depression <- case_when(
  MHQ$Frequency.of.depressed.days.during.worst.episode.of.depression == 2 | MHQ$Frequency.of.depressed.days.during.worst.episode.of.depression == 3 ~ 1,
  MHQ$Frequency.of.depressed.days.during.worst.episode.of.depression == 1 ~ 0,
  TRUE ~ 999
)

#Also create a new depressed.ever phenotype to act as the control for the heritability

### Depressed ever is equivalent to screening for caseness on the CIDI

MHQ$CIDI.MDD.No.Info<-with(MHQ,ifelse(Ever.had.prolonged.feelings.of.sadness.or.depression ==-818 & Ever.had.prolonged.loss.of.interest.in.normal.activities ==-818, 1, 0))

# if recent sadness or loss of interest > 2 weeks response is NOT NA and equal to 1 (yes) AND fraction of day > 2 (about half the day) AND frequency of depressed days > 1 (less often) AND impact on normal roles > 1 (a little), score as 1, else 0 (no info or lower scores for lifelong depression)
#If people have responded don't know or prefer not to answer, recode this to 0 for the sum of items. 
MHQ$Ever.had.prolonged.feelings.of.sadness.or.depression[MHQ$Ever.had.prolonged.feelings.of.sadness.or.depression==-818 | MHQ$Ever.had.prolonged.feelings.of.sadness.or.depression==-121] <- NA

MHQ$Ever.had.prolonged.loss.of.interest.in.normal.activities[MHQ$Ever.had.prolonged.loss.of.interest.in.normal.activities==-818 | MHQ$Ever.had.prolonged.loss.of.interest.in.normal.activities==-121] <- NA

MHQ$Did.your.sleep.change[MHQ$Did.your.sleep.change==-818 | MHQ$Did.your.sleep.change==-121] <- NA 

MHQ$Difficulty.concentrating.during.worst.depression[MHQ$Difficulty.concentrating.during.worst.depression==-818 | MHQ$Difficulty.concentrating.during.worst.depression==-121] <- NA 

MHQ$Weight.change.during.worst.episode.of.depression_Updated <- case_when(
  MHQ$Weight.change.during.worst.episode.of.depression == 1 | MHQ$Weight.change.during.worst.episode.of.depression == 2 | MHQ$Weight.change.during.worst.episode.of.depression == 3 ~ 1,
  MHQ$Weight.change.during.worst.episode.of.depression == 0 ~ 0,
  TRUE ~ 999
  )

MHQ$Weight.change.during.worst.episode.of.depression_Updated[MHQ$Weight.change.during.worst.episode.of.depression_Updated==999] <- NA

MHQ$Thoughts.of.death.during.worst.depression[MHQ$Thoughts.of.death.during.worst.depression==-818 | MHQ$Thoughts.of.death.during.worst.depression==-121] <- NA 

MHQ$Feelings.of.tiredness.during.worst.episode.of.depression[MHQ$Feelings.of.tiredness.during.worst.episode.of.depression==-818 | MHQ$Feelings.of.tiredness.during.worst.episode.of.depression==-121] <- NA 

MHQ$Feelings.of.worthlessness.during.worst.period.of.depression[MHQ$Feelings.of.worthlessness.during.worst.period.of.depression==-818 | MHQ$Feelings.of.worthlessness.during.worst.period.of.depression==-121] <- NA 

MHQ$SymptomNAs <- apply(MHQ[,c("Ever.had.prolonged.feelings.of.sadness.or.depression", "Ever.had.prolonged.loss.of.interest.in.normal.activities", "Did.your.sleep.change", "Difficulty.concentrating.during.worst.depression", 
                                   "Weight.change.during.worst.episode.of.depression_Updated", "Thoughts.of.death.during.worst.depression", "Feelings.of.tiredness.during.worst.episode.of.depression", "Feelings.of.worthlessness.during.worst.period.of.depression")], MARGIN = 1, function(x) sum(is.na(x)))

MHQ$SymptomCount <- rowSums(MHQ[,c("Ever.had.prolonged.feelings.of.sadness.or.depression", "Ever.had.prolonged.loss.of.interest.in.normal.activities", "Did.your.sleep.change", "Difficulty.concentrating.during.worst.depression", 
"Weight.change.during.worst.episode.of.depression_Updated", "Thoughts.of.death.during.worst.depression", "Feelings.of.tiredness.during.worst.episode.of.depression", "Feelings.of.worthlessness.during.worst.period.of.depression")], na.rm=TRUE)

#Defining MHQ depression cases - need to check 
MHQ$LifetimeDepression <- case_when(
  MHQ$Pheno1 == 1 & MHQ$BINARY_Impact.on.normal.roles.during.worst.period.of.depression == 1 & MHQ$BINARY_Frequency.of.depressed.days.during.worst.episode.of.depression  == 1 & MHQ$BINARY_Fraction.of.day.affected.during.worst.episode.of.depression == 1 & MHQ$SymptomCount > 4 ~ 1, 
  MHQ$Pheno1 == 1 & (MHQ$SymptomCount + MHQ$SymptomNAs >= 5) ~ 999,
  MHQ$Pheno1 == 0 | MHQ$BINARY_Impact.on.normal.roles.during.worst.period.of.depression == 0 | MHQ$BINARY_Frequency.of.depressed.days.during.worst.episode.of.depression  == 0 | MHQ$BINARY_Fraction.of.day.affected.during.worst.episode.of.depression == 0 ~ 0,
  TRUE ~ 999
  )

#Select the phenotypes only
MHQ <- MHQ[,c('eid','LifetimeDepression')]
MHQ <- MHQ[MHQ$LifetimeDepression!=999,]

write.csv(MHQ, 'LifetimeDepression.csv')

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
#Extract the necessary variables 

#Read in the datasets
exclusion <- fread('/mnt/lustre/groups/ukbiobank/usr/kylie/hla_analysis/phenotype/icd10_treatment_srinterview_Depressed_Ever_Phenotype_File.txt', data.table=FALSE)

antipsychotics <- fread('/scratch/groups/ukbiobank/usr/Bradley/Study3_Decomposition/Phenotype/antipsychotics', data.table=FALSE)

phenos <- fread('LifetimeDepression.csv', data.table=FALSE)
phenos <- phenos[,-1]

#Create all necessary variables

#####Physician Diagnosed Disorders

# To the question: "Have you been diagnosed with one or more of the following", there are 16 possible responses (i.e. array=16,  see MHQ pdf for details). Therefore, an individual can report up to 16 diagnoses. So for each individual, there are 16 columns of data for this question, such that each column contains an integer representing a possible response (e.g. Social Phobia =1) or NA. 
# The r code is creating a column for social anxiety or social phobia by first looking at the column Ever.sought.or.received.professional.help.for.mental.distress and adding an NA if the response was NA, then looking at each array column for Mental.health.problems.ever.diagnosed.by.a.professional and looking for the integer associated with a phenotype (e.g. 1 = social anxiety). If the column does not contain NA and does contain the relevant integer, participant is coded as 1 for case, else 0. 

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

#####Interview Diagnoses At Baseline######

# notes on code
	# apply function used as follows: apply(X, MARGIN, FUN, ...)
		# X is an array or a matrix if the dimension of the array is 2 (below, X = exclusion[,grep("Non.cancer.illness.code.self.reported", colnames(exclusion))] - this produces a dataframe (157365 rows and 29 columns [equal to the array for "Non.cancer.illness.code.self.reported"]))
		# MARGIN is a variable defining how the function is applied: when MARGIN=1, it applies over rows, whereas with MARGIN=2, it works over columns. Note that when you use the construct MARGIN=c(1,2), it applies to both rows and columns (below, MARGIN = 1, i.e applied over rows)
		# FUN, which is the function that you want to apply to the data. Below, FUN = function(row) "1286" %in% row - breakdown follows:
			# function(row) applied function across rows
			# %in%: A logical vector, indicating if a match was located for each element of x: thus the values are TRUE or FALSE and never NA
				# therefore, FUN searches each row and returns TRUE if codes relating to Bipolar, Schizophrenia or substance abuse are found, else FALSE
				
#Breakdown of exclusion by disorder (note there is going to be overlap between cases as they are not mutually exclusive). 
#1289 = Schizophrenia 
#1291 = Bipolar 
#1408 = Alcohol dependency 
#1409 = Opioid dependency 
#1410 = Other dependency 

NonCancerSR <- readRDS('SelfReportMedicalIllnessNonCancer.rds')
NonCancerSR <- NonCancerSR[,c(1,3:138)]
colnames(NonCancerSR)[1] <- 'ID' 

#All apart from depression can be determined throughout all timepoints. 
NonCancerSR$InterviewExclusion <- apply(NonCancerSR, 1, function(row) any(row %in% c("1289","1291","1408","1409","1410")))
NonCancerSR$InterviewExclusion <- as.integer(as.logical(NonCancerSR$InterviewExclusion))

NonCancerSR <- NonCancerSR[,c('ID','InterviewExclusion')]

exclusion <- full_join(exclusion, NonCancerSR, by='ID')

####HES Inpatient Diagnoses######

# similarly to the above, applies function to search relevant ICD10 codes - returns object with dimensions equal to the number of ICD10 codes by 157365 (participants). Cells are TRUE if a participant has reported the ICD10 code for the column, else FALSE. TRUE = 598, FALSE = 5821907.
#Make a note of frequencies for the R Markdown document

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

project_dir <- "/scratch/datasets/ukbiobank/ukb18177"

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

#####Drug-derived Diagnoses######

###Antipsychotics
colnames(antipsychotics)[1] <- "ID"

## List of codes for antipsychotics
exclusion <- full_join(x=exclusion, y=antipsychotics, by="ID")

###Antidepressants 
antidepcodes <- c('1140879616', '1140921600', '1140879540', '1140867878', '1140916282', '1140909806', '1140867888', '1141152732', '1141180212', '1140879634', '1140867876', '1140882236', '1141190158',
                  '1141200564', '1140867726', '1140879620', '1140867818', '1140879630', '1140879628', '1141151946', '1140867948', '1140867624', '1140867756', '1140867884', '1141151978', '1141152736',
                  '1141201834', '1140867690', '1140867640', '1140867920', '1140867850', '1140879544', '1141200570', '1140867934', '1140867758', '1140867914', '1140867820', '1141151982', '1140882244',
                  '1140879556', '1140867852', '1140867860', '1140917460', '1140867938', '1140867856', '1140867922', '1140910820', '1140882312', '1140867944', '1140867784', '1140867812', '1140867668',
                  '1140867940')

exclusion$Antidepressants <- apply(exclusion[,grep("Treatment.medication.code.0", colnames(exclusion))], 1, function(row) any(row %in% antidepcodes))
exclusion$Antidepressants <- as.integer(as.logical(exclusion$Antidepressants))

###Mania and mood stabilisers
bpdcodes <-c('1140867504','1140867494','1140867498','1140867520','1140917270','1140867490')

exclusion$Lithium <- apply(exclusion[,grep("Treatment.medication.code.0", colnames(exclusion))], 1, function(row) any(row %in% bpdcodes))
exclusion$Lithium <- as.integer(as.logical(exclusion$Lithium)) 

exclusion$BPDrug <- case_when(
exclusion$Lithium == 1 & (exclusion$Antidepressants == 0 | is.na(exclusion$Antidepressants)) ~ 1,
TRUE ~ 0
) 

###Smith Depression Diagnosis 2013

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
relevant <- c('ID', 'InterviewExclusion','SRPsychosisAny','SRManiaBIP', 'Antipsychotics', 'HES','SmithMood', 'BPDrug')
exclusion <- exclusion[,relevant]

#Make a field such that it is the total of all exclusion criteria - anything above 0 equates to a removal - ignore NAs in this function, i.e. treated as 0.
exclusion$total <- rowSums(exclusion[,2:8], na.rm=TRUE)
colnames(exclusion)[1] <- 'eid'

exclusion <- subset(exclusion, total==0)

write.csv(exclusion, 'ExclusionCriteriaApplied.csv')

#Innerjoin to phenotypes dataset
all <- inner_join(x=phenos, y=exclusion, by="eid")
all <- all[,c(1,2)]

write.csv(all, "LifetimeDepression.csv")

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

depression <- fread('LifetimeDepression.csv', data.table=FALSE)
depression <- depression[,-1]

disorders <- fread('./Diseases/Disorders', data.table=FALSE)
disorders <- disorders[,-1]

dep.and.dis <- inner_join(depression,disorders) #25030 of the depression sample have a disorder

#Total numbers reporting each of the disorders 
names <- c('Stroke','Diabetes','Cancer','Autoimmune','MyocardialInfarction','MultipleSclerosis','MotorNeuroneDisease','Epilepsy')

for(i in names){
  print(i)
  a <- dep.and.dis[which(dep.and.dis[[paste(i,".level", sep="")]] == 1 | dep.and.dis[[paste(i,".level", sep="")]] == 2),]
  print(nrow(a))
}

#Breakdown in depressed cases only
dep.and.dis.caseonly <- subset(dep.and.dis, LifetimeDepression == 1)

for(i in names){
  print(i)
  a <- dep.and.dis.caseonly[which(dep.and.dis.caseonly[[paste(i,".level", sep="")]] == 1 | dep.and.dis.caseonly[[paste(i,".level", sep="")]] == 2),]
  print(nrow(a))
}

#Read in Age at onset for case definitions
AAO <- fread('SelfReportAAO_Depression', data.table=FALSE)
AAO <- AAO[,-1]

dep.and.dis.caseonly <- inner_join(dep.and.dis.caseonly,AAO)

dep.and.dis.caseonly$MOB <- as.character(dep.and.dis.caseonly$MOB)
dep.and.dis.caseonly$YOB <- as.character(dep.and.dis.caseonly$YOB)

dep.and.dis.caseonly <- dep.and.dis.caseonly %>% 
  unite(DOB, MOB, YOB, sep = "-", remove = FALSE)
dep.and.dis.caseonly$DOB <- paste("01", dep.and.dis.caseonly$DOB, sep="-")
dep.and.dis.caseonly$DOB <- as.Date(dep.and.dis.caseonly$DOB,"%d-%m-%Y")

#Remove unnecessary columns
dep.and.dis.caseonly <- subset(dep.and.dis.caseonly, select=-c(YOB,MOB,Depressed.Ever))

#Calculate age at onset for each disorder 
date_columns <- c("Date.of.Stroke","Date.of.Diabetes","Date.of.Cancer","Date.of.Autoimmune","Date.of.MyocardialInfarction","Date.of.MultipleSclerosis","Date.of.MotorNeuroneDisease","Date.of.Epilepsy") 
for(i in date_columns){ 
  dep.and.dis.caseonly[,paste(i,'AAO', sep="_")] <- time_length(difftime(dep.and.dis.caseonly[,i], dep.and.dis.caseonly$DOB), "years")
  #Round each age to nearest integer - take floor so the current age is reflected
  dep.and.dis.caseonly[,paste(i,'AAO', sep="_")] <- floor(dep.and.dis.caseonly[,paste(i,'AAO', sep="_")])
}

#Rename columns as still includes reference to date
colnames(dep.and.dis.caseonly)[c(45:52)] <- c("Stroke.AAO","Diabetes.AAO","Cancer.AAO","Autoimmune.AAO","MyocardialInfarction.AAO","MultipleSclerosis.AAO","MotorNeuroneDisease.AAO","Epilepsy.AAO") 

#Subset so that only interested in MDD cases
dep.and.dis.caseonly$AAO_Dep[dep.and.dis.caseonly$AAO_Dep==-818|dep.and.dis.caseonly$AAO_Dep==-121] <- NA
dep.and.dis.caseonly <- subset(dep.and.dis.caseonly, LifetimeDepression==1 & !is.na(AAO_Dep))

#Age difference between each disorder and depression - positive values correspond to depression following disorder
AAO_Cols <- colnames(dep.and.dis.caseonly)[c(45:52)]
for(i in AAO_Cols){
  dep.and.dis.caseonly[,paste(i,"thenDep",sep="")] <- dep.and.dis.caseonly$AAO_Dep - dep.and.dis.caseonly[,i]   
}

colnames(dep.and.dis.caseonly)[c(53:60)] <- c("StrokethenDep","DiabetesthenDep","CancerthenDep","AutoimmunethenDep","MyocardialInfarctionthenDep","MultipleSclerosisthenDep","MotorNeuroneDiseasethenDep","EpilepsythenDep") 

#5 year time difference

#How many people become depression between 0 and 5 years for each disorder - then take the union of this set
year_diff <- colnames(dep.and.dis.caseonly)[c(53:60)]
for(i in year_diff){
  dep.and.dis.caseonly[,paste(i,"ltfive",sep="_")] <- case_when(
    dep.and.dis.caseonly[,paste(i)] <= 5 & dep.and.dis.caseonly[,paste(i)] >= 0 ~ 1,
    TRUE ~ 0
  )
}

#Calculate sample sizes for each disorder
fiveYearCase <- colnames(dep.and.dis.caseonly)[c(61:68)]
for(i in fiveYearCase){
  print(table(dep.and.dis.caseonly[,paste(i)]))
}

#Take the union of the set for total case number
dep.and.dis.caseonly$Case_ltfive <- case_when(
  dep.and.dis.caseonly$StrokethenDep_ltfive == 1 | dep.and.dis.caseonly$DiabetesthenDep_ltfive == 1 | dep.and.dis.caseonly$CancerthenDep_ltfive == 1 | dep.and.dis.caseonly$AutoimmunethenDep_ltfive == 1 |
    dep.and.dis.caseonly$MyocardialInfarctionthenDep_ltfive == 1 | dep.and.dis.caseonly$MultipleSclerosisthenDep_ltfive == 1 | dep.and.dis.caseonly$MotorNeuroneDiseasethenDep_ltfive == 1 | dep.and.dis.caseonly$EpilepsythenDep_ltfive == 1 ~ 1,
  TRUE ~ 0
)

table(dep.and.dis.caseonly$Case_ltfive)

#4 year time difference

#How many people become depression between 0 and 4 years for each disorder - then take the union of this set
year_diff <- colnames(dep.and.dis.caseonly)[c(53:60)]
for(i in year_diff){
  dep.and.dis.caseonly[,paste(i,"ltfour",sep="_")] <- case_when(
    dep.and.dis.caseonly[,paste(i)] <= 4 & dep.and.dis.caseonly[,paste(i)] >= 0 ~ 1,
    TRUE ~ 0
  )
}

#Calculate sample sizes for each disorder
fourYearCase <- colnames(dep.and.dis.caseonly)[c(70:77)]
for(i in fourYearCase){
  print(table(dep.and.dis.caseonly[,paste(i)]))
}

#Take the union of the set for total case number
dep.and.dis.caseonly$Case_ltfour <- case_when(
  dep.and.dis.caseonly$StrokethenDep_ltfour == 1 | dep.and.dis.caseonly$DiabetesthenDep_ltfour == 1 | dep.and.dis.caseonly$CancerthenDep_ltfour == 1 | dep.and.dis.caseonly$AutoimmunethenDep_ltfour == 1 |
    dep.and.dis.caseonly$MyocardialInfarctionthenDep_ltfour == 1 | dep.and.dis.caseonly$MultipleSclerosisthenDep_ltfour == 1 | dep.and.dis.caseonly$MotorNeuroneDiseasethenDep_ltfour == 1 | dep.and.dis.caseonly$EpilepsythenDep_ltfour == 1 ~ 1,
  TRUE ~ 0
)

table(dep.and.dis.caseonly$Case_ltfour)

#3 year time difference

#How many people become depression between 0 and 4 years for each disorder - then take the union of this set
year_diff <- colnames(dep.and.dis.caseonly)[c(53:60)]
for(i in year_diff){
  dep.and.dis.caseonly[,paste(i,"ltthree",sep="_")] <- case_when(
    dep.and.dis.caseonly[,paste(i)] <= 3 & dep.and.dis.caseonly[,paste(i)] >= 0 ~ 1,
    TRUE ~ 0
  )
}

#Calculate sample sizes for each disorder
threeYearCase <- colnames(dep.and.dis.caseonly)[c(79:86)]
for(i in threeYearCase){
  print(table(dep.and.dis.caseonly[,paste(i)]))
}

#Take the union of the set for total case number
dep.and.dis.caseonly$Case_ltthree <- case_when(
  dep.and.dis.caseonly$StrokethenDep_ltthree == 1 | dep.and.dis.caseonly$DiabetesthenDep_ltthree == 1 | dep.and.dis.caseonly$CancerthenDep_ltthree == 1 | dep.and.dis.caseonly$AutoimmunethenDep_ltthree == 1 |
    dep.and.dis.caseonly$MyocardialInfarctionthenDep_ltthree == 1 | dep.and.dis.caseonly$MultipleSclerosisthenDep_ltthree == 1 | dep.and.dis.caseonly$MotorNeuroneDiseasethenDep_ltthree == 1 | dep.and.dis.caseonly$EpilepsythenDep_ltthree == 1 ~ 1,
  TRUE ~ 0
)

table(dep.and.dis.caseonly$Case_ltthree)

#2 year time difference

#How many people become depression between 0 and 4 years for each disorder - then take the union of this set
year_diff <- colnames(dep.and.dis.caseonly)[c(53:60)]
for(i in year_diff){
  dep.and.dis.caseonly[,paste(i,"lttwo",sep="_")] <- case_when(
    dep.and.dis.caseonly[,paste(i)] <= 2 & dep.and.dis.caseonly[,paste(i)] >= 0 ~ 1,
    TRUE ~ 0
  )
}

#Calculate sample sizes for each disorder
twoYearCase <- colnames(dep.and.dis.caseonly)[c(88:95)]
for(i in twoYearCase){
  print(table(dep.and.dis.caseonly[,paste(i)]))
}

#Take the union of the set for total case number
dep.and.dis.caseonly$Case_lttwo <- case_when(
  dep.and.dis.caseonly$StrokethenDep_lttwo == 1 | dep.and.dis.caseonly$DiabetesthenDep_lttwo == 1 | dep.and.dis.caseonly$CancerthenDep_lttwo == 1 | dep.and.dis.caseonly$AutoimmunethenDep_lttwo == 1 |
    dep.and.dis.caseonly$MyocardialInfarctionthenDep_lttwo == 1 | dep.and.dis.caseonly$MultipleSclerosisthenDep_lttwo == 1 | dep.and.dis.caseonly$MotorNeuroneDiseasethenDep_lttwo == 1 | dep.and.dis.caseonly$EpilepsythenDep_lttwo == 1 ~ 1,
  TRUE ~ 0
)

table(dep.and.dis.caseonly$Case_lttwo)

#1 year time difference

#How many people become depression between 0 and 4 years for each disorder - then take the union of this set
year_diff <- colnames(dep.and.dis.caseonly)[c(53:60)]
for(i in year_diff){
  dep.and.dis.caseonly[,paste(i,"ltone",sep="_")] <- case_when(
    dep.and.dis.caseonly[,paste(i)] <= 1 & dep.and.dis.caseonly[,paste(i)] >= 0 ~ 1,
    TRUE ~ 0
  )
}

#Calculate sample sizes for each disorder
oneYearCase <- colnames(dep.and.dis.caseonly)[c(97:104)]
for(i in oneYearCase){
  print(table(dep.and.dis.caseonly[,paste(i)]))
}

#Take the union of the set for total case number
dep.and.dis.caseonly$Case_ltone <- case_when(
  dep.and.dis.caseonly$StrokethenDep_ltone == 1 | dep.and.dis.caseonly$DiabetesthenDep_ltone == 1 | dep.and.dis.caseonly$CancerthenDep_ltone == 1 | dep.and.dis.caseonly$AutoimmunethenDep_ltone == 1 |
    dep.and.dis.caseonly$MyocardialInfarctionthenDep_ltone == 1 | dep.and.dis.caseonly$MultipleSclerosisthenDep_ltone == 1 | dep.and.dis.caseonly$MotorNeuroneDiseasethenDep_ltone == 1 | dep.and.dis.caseonly$EpilepsythenDep_ltone == 1 ~ 1,
  TRUE ~ 0
)

table(dep.and.dis.caseonly$Case_ltone)


######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Now figure out the causes for each depressive episode as in ICD Incidence depression

#Assign temporal causes of major depression

oneyearcases <- subset(dep.and.dis.caseonly, Case_ltone==1)

#Assign all temporal triggers of major depression as if one cannot be identified by 
cause <- apply(oneyearcases[,c('StrokethenDep_ltone','DiabetesthenDep_ltone','CancerthenDep_ltone','AutoimmunethenDep_ltone','MyocardialInfarctionthenDep_ltone','MultipleSclerosisthenDep_ltone','MotorNeuroneDiseasethenDep_ltone','EpilepsythenDep_ltone')] == 1 , 1, function(x) names(which(x)))

#Split into multiple trigger columns
oneyearcases$cause <- cause

#Turn dataframe into a tibble
oneyearcases <- oneyearcases %>%
  separate(col = cause, into = LETTERS[1:8], sep = ", ")

#Rename columns
colnames(oneyearcases)[106:113] <- c('Trigger1','Trigger2','Trigger3','Trigger4','Trigger5','Trigger6','Trigger7','Trigger8')

for(i in colnames(oneyearcases)[106:113]){
  oneyearcases[,paste(i)] <- gsub('thenDep_ltone', '', oneyearcases[,paste(i)] )
}

#Remove quotations and brackets
for(i in 1:8) {
  oneyearcases[[paste('Trigger',i,'intermediate',sep="")]] <- rm_between(oneyearcases[[paste('Trigger',i,sep="")]], '"', '"', extract=TRUE)
  oneyearcases[[paste('Trigger',i,'intermediate',sep="")]]  <- as.character(oneyearcases[[paste('Trigger',i,'intermediate',sep="")]])
  oneyearcases[[paste('Trigger',i,'intermediate',sep="")]] <- ifelse(is.na(oneyearcases[[paste('Trigger',i,'intermediate',sep="")]]), oneyearcases[[paste('Trigger',i,sep="")]], oneyearcases[[paste('Trigger',i,'intermediate',sep="")]])
  oneyearcases[[paste('Trigger',i,sep="")]] <- oneyearcases[[paste('Trigger',i,'intermediate',sep="")]]
}

oneyearcases <- oneyearcases[,-c(114:121)]

write.csv(oneyearcases, 'self_report_depression')

#Merge cases from two methods of identification
selfreport <- fread('self_report_depression', data.table=FALSE)
selfreport <- selfreport[,-1]

#Extract only relevant columns that are consistent between the two types - makes working with it later easier. 
selfreport <- selfreport[,c(1,3:42,106:113)]

incident <- fread('incident_depression', data.table=FALSE)
incident <- incident[,-1]

#Extract only relevant columns that are consistent between the two types - makes working with it later easier. 
incident <- incident[,c(1,11:50,107:114)]

#inner join to see which cases have been identified in both areas
overlap.depressed <- inner_join(selfreport, incident, by='eid')

#Identify rows where the causes are not equal to one another. All causes agree which means a simple full_join combining on both ID and cause
alldepressed <- full_join(selfreport, incident)

exclusion <- fread('ExclusionCriteriaApplied.csv', data.table=FALSE)
exclusion <- exclusion[,-1]

alldepressed <- inner_join(alldepressed, exclusion) 

alldepressed <- alldepressed[,-c(50:57)] #370 cases using QC on full sample. 
write.csv(alldepressed, 'allcases_within_one_year')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Regression analysis 

#MDD PRS

cases <- fread('depression_within_year_selfreport_ScreenedandQCd', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

controls <- fread('controlsMedicalDisorder_QC.csv', data.table=FALSE)
controls$MDD <- 0
colnames(controls)[2] <- 'eid'
controls <- controls[,c(2,19)]

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

MDDresult$Pheno <- 'All' 

write.csv(MDDresult, 'MDDwithDisorder_MDDPRSresults')

######################################################################################################################################################################################################################################

#BPD PRS

cases <- fread('depression_within_year_selfreport_ScreenedandQCd', data.table=FALSE)
cases$BPD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

controls <- fread('controlsMedicalDisorder_QC.csv', data.table=FALSE)
controls$BPD <- 0
colnames(controls)[2] <- 'eid'
controls <- controls[,c(2,19)]

BPD <- rbind(cases, controls)

#Get necessary genetic data and covariates (PCs, assessment centre and batch)
covariates <- fread('ukb18177_glanville_covariates.txt' , data.table=FALSE)
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

BPDresult$Pheno <- 'All' 

write.csv(BPDresult, 'MDDwithDisorder_BPDPRSresults')

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################

#Environmental Variables Analysis
#Read in cases and controls for PPD - get an idea of degree of missingness for each variable
environment <- fread('EnvironmentalVariables.csv', data.table=FALSE)
environment <- environment[,-1]

#European ancestry only
cases <- fread('depression_within_year_selfreport_Screened', data.table=FALSE)
cases$MDD <- 1
cases <- cases[,-1]
colnames(cases)[1] <- 'eid'

controls <- fread('controlsMedicalDisorder.csv', data.table=FALSE)
controls$MDD <- 0
controls <- controls[,-c(1,3:18)]

MDD <- rbind(cases, controls)

#Get necessary genetic data and covariates assessment centre
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
  print(table(MDDcase[[i]]))
  print(table(MDDcontrol[[i]]))
}

#Most of the missingness driven by controls - of the two case definitions used - is one contributing to missingness more than the other?
env_variables <- c('TDI','neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

for(i in env_variables){
  print(i)
  variable <- MDD[is.na(MDD[[i]]),]
  variabletwo <- subset(variable, MDD==1)
  #IncidentCases <- subset(variabletwo, variabletwo$eid %in% incidentnumber$eid) 
  #SelfReportCases <- subset(variabletwo, variabletwo$eid %in% selfreportnumber$eid)
  print(length(variable$eid))
  print(length(variabletwo$eid))
  #print(length(IncidentCases$eid)) 
  #print(length(SelfReportCases$eid))
}

#All
table(MDD$education_attainment)

#Cases
case <- subset(MDD, MDD==1)
table(case$education_attainment)

#Controls
control <- subset(MDD, MDD==0)
table(control$education_attainment)

#Remove missing cases of SES and education attainment and see level of missingness of environmental variable
MDD <- subset(MDD, !(education_attainment == "missing" | education_attainment == "pnta" | is.na(TDI)))

env_variables <- c('neuroticism','child_trauma_PC1','adult_trauma_PC1','adult_trauma_nolti_PC1','FamHist')

for(i in env_variables){
  print(i)
  variable <- MDD[is.na(MDD[[i]]),]
  variabletwo <- subset(variable, MDD==1)
  print(length(variable$eid))
  print(length(variabletwo$eid))
}

MDD$assessment_centre <- factor(MDD$assessment_centre, ordered=FALSE)
MDD$education_attainment <- factor(MDD$education_attainment, ordered=FALSE)
MDD$education_attainment <- relevel(MDD$education_attainment, ref='degree')

MDD$TDI <- scale(MDD$TDI)
MDD$neuroticism <- scale(MDD$neuroticism) 

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

write.csv(MDDEnvResults, 'EnvironmentalAssociationsMDDthenDep_FullSample.csv')

