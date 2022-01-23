library(data.table)
library(ggplot2)
library(tidyverse)
library(ukbkings)

#UKB Kings to extract the trauma items

project_dir <- "ukb18177"
f <- bio_field(project_dir)

#Select data
f %>%
  select(field, name) %>%
  filter(str_detect(name, "f20489|f20488|f20487|f20490|f20491|f20522|f20523|f20521|f20524|f20525|f20531|f20529|f20526|f20530|f20528|f20527")) %>%
  bio_field_add("TraumaVariables.txt")

#Read data into file
bio_phen(
  project_dir,
  field = "TraumaVariables.txt",
  out = "TraumaVariables"
)

#Read in Environmental Variables
trauma <- readRDS('TraumaVariables.rds') 

#Name columns appropriately - correct for assessment centre still in regression.
colnames(trauma)[c(1:24, 109)] <- c('eid','loved_as_child','doctor_as_child','confiding_relationship','pay_mortgage','serious_accident','exposed_combat','life_threatening_illness','violent_crime','violent_death','family_member_hate','physical_abuse','sexual_molest','partner_belittle','partner_violence','partner_sexual','sexual_assault')

#Reverse code variables so that each variable has the reference of 0 and anything higher represents more trauma
reversed <- c('loved_as_child','doctor_as_child','confiding_relationship','pay_mortgage')

for(i in reversed){
  print(table(trauma[[i]]))
  trauma[[i]] <- recode(trauma[[i]], "4=0; 3=1; 2=2; 1=3; 0=4")
  print(table(trauma[[i]]))
}

#Convert PNTA/Don't knows (-818) to NA.
indtrau <- c('loved_as_child','doctor_as_child','confiding_relationship','pay_mortgage','serious_accident','exposed_combat','life_threatening_illness','violent_crime','violent_death','family_member_hate','physical_abuse','sexual_molest','partner_belittle','partner_violence','partner_sexual','sexual_assault')
for(i in indtrau){
  print(table(trauma[[i]]))
  trauma[[i]][trauma[[i]]==-818] <- NA
  print(table(trauma[[i]]))
}

#Make each ordinal variable binary according to pre-defined cut-offs

#Often/sometimes true as cut-offs for these items - refer to supplementary material
oftentimes <- c('loved_as_child','family_member_hate','doctor_as_child','confiding_relationship','pay_mortgage')

for(i in oftentimes){
  print(i)
  print(table(trauma[[i]]))
  trauma[[paste(i,'BINARY',sep='')]] <- case_when(trauma[[i]] == 0 | trauma[[i]] == 1 ~ 0,
                                                       is.na(trauma[[i]]) ~ 999,
                                                       TRUE ~ 1)
  trauma[[paste(i,'BINARY',sep='')]][trauma[[paste(i,'BINARY',sep='')]]==999] <- NA
  print(table(trauma[[paste(i,'BINARY',sep='')]]))
}

#Never/any evidence used as cut-offs for these items - refer to supplementary materials
nevertimes <- c('serious_accident','exposed_combat','life_threatening_illness','violent_crime','violent_death','physical_abuse','sexual_molest','partner_belittle','partner_violence','partner_sexual','sexual_assault')

for(i in nevertimes){
  print(i)
  print(table(trauma[[i]]))
  trauma[[paste(i,'BINARY',sep='')]] <- case_when(trauma[[i]] == 0 ~ 0,
                                                       is.na(trauma[[i]]) ~ 999,
                                                       TRUE ~ 1)
  trauma[[paste(i,'BINARY',sep='')]][trauma[[paste(i,'BINARY',sep='')]]==999] <- NA
  print(table(trauma[[paste(i,'BINARY',sep='')]]))
}

write.csv(trauma, 'TraumaVariables.csv')


#Read in dataset and extract relevant variables
trauma <- fread('TraumaVariables.csv', data.table=FALSE)
trauma <- trauma[,-1]

child_trauma <- trauma[,c('eid','loved_as_childBINARY','doctor_as_childBINARY','physical_abuseBINARY','family_member_hateBINARY','sexual_molestBINARY')]
rownames(child_trauma) <- child_trauma$eid
child_trauma <- child_trauma[,-1]

adult_trauma <- trauma[,c('eid', 'confiding_relationshipBINARY','pay_mortgageBINARY','serious_accidentBINARY','exposed_combatBINARY','life_threatening_illnessBINARY','violent_crimeBINARY','violent_deathBINARY','partner_belittleBINARY','partner_violenceBINARY','partner_sexualBINARY','sexual_assaultBINARY')]
rownames(adult_trauma) <- adult_trauma$eid
adult_trauma <- adult_trauma[,-1]

#Convert variables to ordered factors. 

#No imputation - take complete cases only - N=153637
child_trauma <- child_trauma[complete.cases(child_trauma),]

child_trauma_corrmatrix <- tetrachoric(child_trauma)

#Manual to check getting similar results
s.eigen <- eigen(child_trauma_corrmatrix[[1]])

#Proportion of variance explained
for (s in s.eigen$values) {
  print(s / sum(s.eigen$values))
}

#Projection test - i.e. principal components for each person
child_trauma$eid <- rownames(child_trauma)

child_trauma <- child_trauma %>%
  mutate_if(is.ordered,as.numeric) 

Prin <- as.matrix(child_trauma[,-6]) %*% child_trauma_pca$loadings[1:5,1:5]
Prin <- as.data.frame(Prin)

child_trauma <- cbind(child_trauma, Prin)

#Plot density of childhood trauma
bitmap('Child_Trauma_PC_Density.png', res=1200)
  print(ggplot() + geom_density(data = child_trauma, aes(x=Comp.1)))
dev.off()

#Loadings of childhood trauma pca
bitmap('child_trauma_pca_loadings.png', res=1200)
  print(fviz_pca_var(child_trauma_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
  ))
dev.off()

#Plot people by first two dimensions - pca
bitmap('child_trauma_pca_individuals.png', res=1200)
  print(ggplot(child_trauma, aes(Comp.1, Comp.2)) +
          geom_hex() +
          scale_fill_viridis_c() +
          geom_point(shape = '.', col = 'white'))
dev.off()

#Select first principal component for childhood
child_trauma <- child_trauma[,c('eid','Comp.1')]
child_trauma$eid <- as.numeric(child_trauma$eid)
colnames(child_trauma) <- c('eid','child_trauma_PC1')

#############################################################################################
#############################################################################################

#Adult trauma PCA

#No imputation - take complete cases only - N=148441
adult_trauma <- adult_trauma[complete.cases(adult_trauma),]

adult_trauma_corrmatrix <- tetrachoric(adult_trauma)

#Manual to check getting similar results
s.eigen <- eigen(adult_trauma_corrmatrix[[1]])

#Proportion of variance explained
for (s in s.eigen$values) {
  print(s / sum(s.eigen$values))
}

#PCA using just the correlation matrix
adult_trauma_pca <- princomp(covmat=adult_trauma_corrmatrix[[1]])

#Projection test - i.e. principal components for each person
adult_trauma$eid <- rownames(adult_trauma)

adult_trauma <- adult_trauma %>%
  mutate_if(is.ordered,as.numeric) 

Prin <- as.matrix(adult_trauma[,-12]) %*% adult_trauma_pca$loadings[1:11,1:11]
Prin <- as.data.frame(Prin)

adult_trauma <- cbind(adult_trauma, Prin)

#Plot density of childhood trauma
bitmap('Adult_Trauma_PC_Density.png', res=1200)
print(ggplot() + geom_density(data = adult_trauma, aes(x=Comp.1)))
dev.off()

#Loadings of childhood trauma pca
bitmap('adult_trauma_pca_loadings.png', res=1200)
print(fviz_pca_var(adult_trauma_pca,
                   col.var = "contrib", # Color by contributions to the PC
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE     # Avoid text overlapping
))
dev.off()

#Plot people by first two dimensions - pca
bitmap('adult_trauma_pca_individuals.png', res=1200)
print(ggplot(adult_trauma, aes(Comp.1, Comp.2)) +
        geom_hex() +
        scale_fill_viridis_c() +
        geom_point(shape = '.', col = 'white'))
dev.off()

#Do the same plot with zeros removed as they are the dominant group. 
adult_trauma_nozero <- subset(adult_trauma, Comp.1!=0 & Comp.2!=0)

bitmap('adult_trauma_nozero_pca_individuals.png', res=1200)
print(ggplot(adult_trauma_nozero, aes(Comp.1, Comp.2)) +
        geom_hex() +
        scale_fill_viridis_c() +
        geom_point(shape = '.', col = 'white'))
dev.off()

#Select first principal component for childhood
adult_trauma <- adult_trauma[,c('eid','Comp.1')]
adult_trauma$eid <- as.numeric(adult_trauma$eid)
colnames(adult_trauma) <- c('eid','adult_trauma_PC1')