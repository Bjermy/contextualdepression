#Libraries
library(data.table)

#Step 1: Read in all Odds Ratios for each group and combine
EnvMDDexc <- read.csv('./Heterogeneous/EnvironmentalAssociations_Heterogeneous_FullSample.csv')
EnvMDDexc <- EnvMDDexc[EnvMDDexc$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'), c('EnvVariable','beta','se')]
EnvMDDexc$Pheno <- 'MDDexc'

EnvPPD <- read.csv('./Postpartum/EnvironmentalAssociationsPostPartum_FullSample_withselfreportbirth.csv')
EnvPPD <- EnvPPD[EnvPPD$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'), c('EnvVariable','beta','se')]
EnvPPD$Pheno <- 'PostPartum'

EnvDisease <- read.csv('./Chronic_Disease/EnvironmentalAssociationsMDDthenDep_FullSample.csv')
EnvDisease <- EnvDisease[EnvDisease$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'), c('EnvVariable','beta','se')]
EnvDisease$Pheno <- 'DiseaseThenMDD'

MDDPRS <- read.csv('AllDefinitionsTable_MDDPRS_withselfreportbirth.csv')
MDDPRS <- MDDPRS[MDDPRS$Pt == 'Pt < 0.05', c('Pheno','beta','se')]
MDDPRS$EnvVariable <- 'MDDPRS' 

BPDPRS <- read.csv('AllDefinitionsTable_BPDPRS_withselfreportbirth.csv')
BPDPRS <- BPDPRS[BPDPRS$Pt == 'Pt < 0.01', c('Pheno','beta','se')]
BPDPRS$EnvVariable <- 'BPDPRS' 

AllVars <- rbind(EnvMDDexc, EnvPPD, EnvDisease)

TDI <- AllVars[AllVars$EnvVariable == 'TDI',]

alevel <- AllVars[AllVars$EnvVariable == 'education_attainmentalevels',]
gcse <- AllVars[AllVars$EnvVariable == 'education_attainmentgcse',]
none <- AllVars[AllVars$EnvVariable == 'education_attainmentnone',]
  
neuroticism <- AllVars[AllVars$EnvVariable == 'neuroticism',]
  
child_trauma <- AllVars[AllVars$EnvVariable == 'child_trauma_PC1',]
  
adult_trauma <- AllVars[AllVars$EnvVariable == 'adult_trauma_PC1',]
  
FamHist <- AllVars[AllVars$EnvVariable == 'FamHist',]

#Step 2: Perform OR comparisons

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#TDI

#Subtract betas to get differences.
TDIdelta <- c(TDI[TDI$Pheno=='MDDexc','beta'] - TDI[TDI$Pheno=='PostPartum','beta'], TDI[TDI$Pheno=='MDDexc','beta'] - TDI[TDI$Pheno=='DiseaseThenMDD','beta'])

#Calculate Standard errors
TDIse <- c(sqrt((TDI[TDI$Pheno=='MDDexc','se']^2) + (TDI[TDI$Pheno=='PostPartum','se']^2)), sqrt((TDI[TDI$Pheno=='MDDexc','se']^2) + (TDI[TDI$Pheno=='DiseaseThenMDD','se']^2)))

#Calculate z-score
TDIz <- TDIdelta/TDIse

#Calculate p-value
TDIpval <- 2*pnorm(-abs(TDIz))
names(TDIpval) <- c('PostPartum','DiseaseThenMDD')

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#Education

#A-levels

#Subtract betas to get differences.
Aleveldelta <- c(alevel[alevel$Pheno=='MDDexc','beta'] - alevel[alevel$Pheno=='PostPartum','beta'], alevel[alevel$Pheno=='MDDexc','beta'] - alevel[alevel$Pheno=='DiseaseThenMDD','beta'])

#Calculate Standard errors
Alevelse <- c(sqrt((alevel[alevel$Pheno=='MDDexc','se']^2) + (alevel[alevel$Pheno=='PostPartum','se']^2)), sqrt((alevel[alevel$Pheno=='MDDexc','se']^2) + (alevel[alevel$Pheno=='DiseaseThenMDD','se']^2)))

#Calculate z-score
Alevelz <- Aleveldelta/Alevelse

#Calculate p-value
Alevelpval <- 2*pnorm(-abs(Alevelz))
names(Alevelpval) <- c('PostPartum','DiseaseThenMDD')

#GCSE

#Subtract betas to get differences.
gcsedelta <- c(gcse[gcse$Pheno=='MDDexc','beta'] - gcse[gcse$Pheno=='PostPartum','beta'], gcse[gcse$Pheno=='MDDexc','beta'] - gcse[gcse$Pheno=='DiseaseThenMDD','beta'])

#Calculate Standard errors
gcsese <- c(sqrt((gcse[gcse$Pheno=='MDDexc','se']^2) + (gcse[gcse$Pheno=='PostPartum','se']^2)), sqrt((gcse[gcse$Pheno=='MDDexc','se']^2) + (gcse[gcse$Pheno=='DiseaseThenMDD','se']^2)))

#Calculate z-score
gcsez <- gcsedelta/gcsese

#Calculate p-value
gcsepval <- 2*pnorm(-abs(gcsez))
names(gcsepval) <- c('PostPartum','DiseaseThenMDD')

#None

#Subtract betas to get differences.
nonedelta <- c(none[none$Pheno=='MDDexc','beta'] - none[none$Pheno=='PostPartum','beta'], none[none$Pheno=='MDDexc','beta'] - none[none$Pheno=='DiseaseThenMDD','beta'])

#Calculate Standard errors
nonese <- c(sqrt((none[none$Pheno=='MDDexc','se']^2) + (none[none$Pheno=='PostPartum','se']^2)), sqrt((none[none$Pheno=='MDDexc','se']^2) + (none[none$Pheno=='DiseaseThenMDD','se']^2)))

#Calculate z-score
nonez <- nonedelta/nonese

#Calculate p-value
nonepval <- 2*pnorm(-abs(nonez))
names(nonepval) <- c('PostPartum','DiseaseThenMDD')

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#Neuroticism

#Subtract betas to get differences.
Neuroticismdelta <- c(neuroticism[neuroticism$Pheno=='MDDexc','beta'] - neuroticism[neuroticism$Pheno=='PostPartum','beta'], neuroticism[neuroticism$Pheno=='MDDexc','beta'] - neuroticism[neuroticism$Pheno=='DiseaseThenMDD','beta'])

#Calculate Standard errors
Neuroticismse <- c(sqrt((neuroticism[neuroticism$Pheno=='MDDexc','se']^2) + (neuroticism[neuroticism$Pheno=='PostPartum','se']^2)), sqrt((neuroticism[neuroticism$Pheno=='MDDexc','se']^2) + (neuroticism[neuroticism$Pheno=='DiseaseThenMDD','se']^2)))

#Calculate z-score
Neuroticismz <- Neuroticismdelta/Neuroticismse

#Calculate p-value
Neuroticismpval <- 2*pnorm(-abs(Neuroticismz))
names(Neuroticismpval) <- c('PostPartum','DiseaseThenMDD')

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#Child Trauma

#Subtract betas to get differences.
Childtraumadelta <- c(child_trauma[child_trauma$Pheno=='MDDexc','beta'] - child_trauma[child_trauma$Pheno=='PostPartum','beta'], child_trauma[child_trauma$Pheno=='MDDexc','beta'] - child_trauma[child_trauma$Pheno=='DiseaseThenMDD','beta'])

#Calculate Standard errors
Childtraumase <- c(sqrt((child_trauma[child_trauma$Pheno=='MDDexc','se']^2) + (child_trauma[child_trauma$Pheno=='PostPartum','se']^2)), sqrt((child_trauma[child_trauma$Pheno=='MDDexc','se']^2) + (child_trauma[child_trauma$Pheno=='DiseaseThenMDD','se']^2)))

#Calculate z-score
Childtraumaz <- Childtraumadelta/Childtraumase

#Calculate p-value
Childtraumapval <- 2*pnorm(-abs(Childtraumaz))
names(Childtraumapval) <- c('PostPartum','DiseaseThenMDD')

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#Adult Trauma

#Subtract betas to get differences.
Adulttraumadelta <- c(adult_trauma[adult_trauma$Pheno=='MDDexc','beta'] - adult_trauma[adult_trauma$Pheno=='PostPartum','beta'], adult_trauma[adult_trauma$Pheno=='MDDexc','beta'] - adult_trauma[adult_trauma$Pheno=='DiseaseThenMDD','beta'])

#Calculate Standard errors
Adulttraumase <- c(sqrt((adult_trauma[adult_trauma$Pheno=='MDDexc','se']^2) + (adult_trauma[adult_trauma$Pheno=='PostPartum','se']^2)), sqrt((adult_trauma[adult_trauma$Pheno=='MDDexc','se']^2) + (adult_trauma[adult_trauma$Pheno=='DiseaseThenMDD','se']^2)))

#Calculate z-score
Adulttraumaz <- Adulttraumadelta/Adulttraumase

#Calculate p-value
Adulttraumapval <- 2*pnorm(-abs(Adulttraumaz))
names(Adulttraumapval) <- c('PostPartum','DiseaseThenMDD')

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#Family History

#Subtract betas to get differences.
FamHistdelta <- c(FamHist[FamHist$Pheno=='MDDexc','beta'] - FamHist[FamHist$Pheno=='PostPartum','beta'], FamHist[FamHist$Pheno=='MDDexc','beta'] - FamHist[FamHist$Pheno=='DiseaseThenMDD','beta'])

#Calculate Standard errors
FamHistse <- c(sqrt((FamHist[FamHist$Pheno=='MDDexc','se']^2) + (FamHist[FamHist$Pheno=='PostPartum','se']^2)), sqrt((FamHist[FamHist$Pheno=='MDDexc','se']^2) + (FamHist[FamHist$Pheno=='DiseaseThenMDD','se']^2)))

#Calculate z-score
FamHistz <- FamHistdelta/FamHistse

#Calculate p-value
FamHistpval <- 2*pnorm(-abs(FamHistz))
names(FamHistpval) <- c('PostPartum','DiseaseThenMDD')

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#MDD PRS

#Subtract betas to get differences.
MDDPRSdelta <- c(MDDPRS[MDDPRS$Pheno=='MDDexc','beta'] - MDDPRS[MDDPRS$Pheno=='PostPartum','beta'], MDDPRS[MDDPRS$Pheno=='MDDexc','beta'] - MDDPRS[MDDPRS$Pheno=='DiseaseThenMDD','beta'])

#Calculate Standard errors
MDDPRSse <- c(sqrt((MDDPRS[MDDPRS$Pheno=='MDDexc','se']^2) + (MDDPRS[MDDPRS$Pheno=='PostPartum','se']^2)), sqrt((MDDPRS[MDDPRS$Pheno=='MDDexc','se']^2) + (MDDPRS[MDDPRS$Pheno=='DiseaseThenMDD','se']^2)))

#Calculate z-score
MDDPRSz <- MDDPRSdelta/MDDPRSse

#Calculate p-value
MDDPRSpval <- 2*pnorm(-abs(MDDPRSz))
names(MDDPRSpval) <- c('PostPartum','DiseaseThenMDD')

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#BPD PRS

#Subtract betas to get differences.
BPDPRSdelta <- c(BPDPRS[MDDPRS$Pheno=='MDDexc','beta'] - BPDPRS[MDDPRS$Pheno=='PostPartum','beta'], BPDPRS[MDDPRS$Pheno=='MDDexc','beta'] - BPDPRS[BPDPRS$Pheno=='DiseaseThenMDD','beta'])

#Calculate Standard errors
BPDPRSse <- c(sqrt((BPDPRS[BPDPRS$Pheno=='MDDexc','se']^2) + (BPDPRS[BPDPRS$Pheno=='PostPartum','se']^2)), sqrt((BPDPRS[BPDPRS$Pheno=='MDDexc','se']^2) + (BPDPRS[BPDPRS$Pheno=='DiseaseThenMDD','se']^2)))

#Calculate z-score
BPDPRSz <- BPDPRSdelta/BPDPRSse

#Calculate p-value
BPDPRSpval <- 2*pnorm(-abs(BPDPRSz))
names(BPDPRSpval) <- c('PostPartum','DiseaseThenMDD')

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#Step 3: FDR correction

#Collate p-values
p_vals <- rbind(TDIpval, Alevelpval, gcsepval, nonepval, Neuroticismpval, Childtraumapval, Adulttraumapval, FamHistpval, MDDPRSpval, BPDPRSpval)
p_valsFDR <- c(TDIpval, Alevelpval, gcsepval, nonepval, Neuroticismpval, Childtraumapval, Adulttraumapval, FamHistpval, MDDPRSpval, BPDPRSpval)

sorted.pvalue<-sort(p_valsFDR) 
sorted.pvalue 
j.alpha <- (1:20)*(.05/20) # input total number of pvalues
diff <- sorted.pvalue-j.alpha 
neg.diff <- diff[diff<0] 
pos.diff <- neg.diff[length(neg.diff)] 
index <- diff==pos.diff 
p.cutoff <-sorted.pvalue[index] 
p.cutoff # this value is the critical pvalue. Values below should be considered significant
p.sig <- p_valsFDR[p_valsFDR <= p.cutoff] 
p.sig # which pvalues are significant based on the critical pvalue

#########################################################################################################################################################################################
#########################################################################################################################################################################################
#########################################################################################################################################################################################
#########################################################################################################################################################################################

#Repeat the above for heterogeneous definition of depression that excludes males

#Step 1: Read in all Odds Ratios for each group and combine
EnvMDDexc <- read.csv('./Heterogeneous/EnvironmentalAssociations_Heterogeneous_FemaleOnly.csv')
EnvMDDexc <- EnvMDDexc[EnvMDDexc$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'), c('EnvVariable','beta','se')]
EnvMDDexc$Pheno <- 'MDDexc'
EnvMDDexc <- EnvMDDexc[1:8,]

EnvPPD <- read.csv('./Postpartum/EnvironmentalAssociationsPostPartum_FullSample_withselfreportbirth.csv')
EnvPPD <- EnvPPD[EnvPPD$EnvVariable %in% c('TDI','education_attainmentalevels','education_attainmentgcse','education_attainmentnone','neuroticism','child_trauma_PC1','adult_trauma_PC1','FamHist'), c('EnvVariable','beta','se')]
EnvPPD$Pheno <- 'PostPartum'

MDDPRShet <- read.csv('./Heterogeneous/Heterogeneous_MDDPRSresults_FemaleOnly.csv')
MDDPRShet$beta <- log(MDDPRShet$OR)
MDDPRShet <- MDDPRShet[,-c(1,3,8,9)]
MDDPRShet$Pheno <- 'MDDexc'

MDDPRSppd <- read.csv('./Postpartum/PPD_MDDPRSresults_withselfreportbirth.csv')
MDDPRSppd <- MDDPRSppd[MDDPRSppd$Pheno=='All',c('Pt','beta','se','p','Pheno')]
MDDPRSppd$Pheno <- 'PostPartum'

MDDPRS <- rbind(MDDPRShet,MDDPRSppd)
MDDPRS$Pt <- factor(MDDPRS$Pt, levels = c("DEPR06_Pt_0.05", "DEPR06_Pt_1"), 
                    labels = c("Pt < 0.05", "Pt < 1"))

MDDPRS <- MDDPRS[MDDPRS$Pt == 'Pt < 0.05', c('Pheno','beta','se')]
MDDPRS$EnvVariable <- 'MDDPRS' 

BPDPRShet <- read.csv('./Heterogeneous/Heterogeneous_BPDPRSresults_FemaleOnly.csv')
BPDPRShet$beta <- log(BPDPRShet$OR)
BPDPRShet <- BPDPRShet[,-c(1,3,8,9)]
BPDPRShet$Pheno <- 'MDDexc'

BPDPRSppd <- read.csv('./Postpartum/PPD_BPDPRSresults_withselfreportbirth.csv')
BPDPRSppd <- BPDPRSppd[BPDPRSppd$Pheno=='All',c('Pt','beta','se','p','Pheno')]
BPDPRSppd$Pheno <- 'PostPartum'

BPDPRS <- rbind(BPDPRShet,BPDPRSppd)
BPDPRS$Pt <- factor(BPDPRS$Pt, levels = c("BIPO02_Pt_0.01", "BIPO02_Pt_1"), 
                    labels = c("Pt < 0.01", "Pt < 1"))

BPDPRS <- BPDPRS[BPDPRS$Pt == 'Pt < 0.01', c('Pheno','beta','se')]
BPDPRS$EnvVariable <- 'BPDPRS' 

AllVars <- rbind(EnvMDDexc, EnvPPD)

TDI <- AllVars[AllVars$EnvVariable == 'TDI',]

alevel <- AllVars[AllVars$EnvVariable == 'education_attainmentalevels',]
gcse <- AllVars[AllVars$EnvVariable == 'education_attainmentgcse',]
none <- AllVars[AllVars$EnvVariable == 'education_attainmentnone',]

neuroticism <- AllVars[AllVars$EnvVariable == 'neuroticism',]

child_trauma <- AllVars[AllVars$EnvVariable == 'child_trauma_PC1',]

adult_trauma <- AllVars[AllVars$EnvVariable == 'adult_trauma_PC1',]

FamHist <- AllVars[AllVars$EnvVariable == 'FamHist',]

#Step 2: Perform OR comparisons

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#TDI

#Subtract betas to get differences.
TDIdelta <- TDI[TDI$Pheno=='MDDexc','beta'] - TDI[TDI$Pheno=='PostPartum','beta']

#Calculate Standard errors
TDIse <- sqrt((TDI[TDI$Pheno=='MDDexc','se']^2) + (TDI[TDI$Pheno=='PostPartum','se']^2))

#Calculate z-score
TDIz <- TDIdelta/TDIse

#Calculate p-value
TDIpval <- 2*pnorm(-abs(TDIz))

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#A-levels

#Subtract betas to get differences.
aleveldelta <- alevel[alevel$Pheno=='MDDexc','beta'] - alevel[alevel$Pheno=='PostPartum','beta']

#Calculate Standard errors
alevelse <- sqrt((alevel[alevel$Pheno=='MDDexc','se']^2) + (alevel[alevel$Pheno=='PostPartum','se']^2))

#Calculate z-score
alevelz <- aleveldelta/alevelse

#Calculate p-value
alevelpval <- 2*pnorm(-abs(alevelz))

#GCSE

#Subtract betas to get differences.
gcsedelta <- gcse[gcse$Pheno=='MDDexc','beta'] - gcse[gcse$Pheno=='PostPartum','beta']

#Calculate Standard errors
gcsese <- sqrt((gcse[gcse$Pheno=='MDDexc','se']^2) + (gcse[gcse$Pheno=='PostPartum','se']^2))

#Calculate z-score
gcsez <- gcsedelta/gcsese

#Calculate p-value
gcsepval <- 2*pnorm(-abs(gcsez))

#None

#Subtract betas to get differences.
nonedelta <- none[none$Pheno=='MDDexc','beta'] - none[none$Pheno=='PostPartum','beta']

#Calculate Standard errors
nonese <- sqrt((none[none$Pheno=='MDDexc','se']^2) + (none[none$Pheno=='PostPartum','se']^2))

#Calculate z-score
nonez <- nonedelta/nonese

#Calculate p-value
nonepval <- 2*pnorm(-abs(nonez))

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#Neuroticism

#Subtract betas to get differences.
Neuroticismdelta <- neuroticism[neuroticism$Pheno=='MDDexc','beta'] - neuroticism[neuroticism$Pheno=='PostPartum','beta']

#Calculate Standard errors
Neuroticismse <- sqrt((neuroticism[neuroticism$Pheno=='MDDexc','se']^2) + (neuroticism[neuroticism$Pheno=='PostPartum','se']^2))

#Calculate z-score
Neuroticismz <- Neuroticismdelta/Neuroticismse

#Calculate p-value
Neuroticismpval <- 2*pnorm(-abs(Neuroticismz))

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#Child Trauma

#Subtract betas to get differences.
Childtraumadelta <- child_trauma[child_trauma$Pheno=='MDDexc','beta'] - child_trauma[child_trauma$Pheno=='PostPartum','beta']

#Calculate Standard errors
Childtraumase <- sqrt((child_trauma[child_trauma$Pheno=='MDDexc','se']^2) + (child_trauma[child_trauma$Pheno=='PostPartum','se']^2))

#Calculate z-score
Childtraumaz <- Childtraumadelta/Childtraumase

#Calculate p-value
Childtraumapval <- 2*pnorm(-abs(Childtraumaz))

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#Adult Trauma

#Subtract betas to get differences.
Adulttraumadelta <- adult_trauma[adult_trauma$Pheno=='MDDexc','beta'] - adult_trauma[adult_trauma$Pheno=='PostPartum','beta']

#Calculate Standard errors
Adulttraumase <- sqrt((adult_trauma[adult_trauma$Pheno=='MDDexc','se']^2) + (adult_trauma[adult_trauma$Pheno=='PostPartum','se']^2))

#Calculate z-score
Adulttraumaz <- Adulttraumadelta/Adulttraumase

#Calculate p-value
Adulttraumapval <- 2*pnorm(-abs(Adulttraumaz))

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#Family History

#Subtract betas to get differences.
FamHistdelta <- FamHist[FamHist$Pheno=='MDDexc','beta'] - FamHist[FamHist$Pheno=='PostPartum','beta']

#Calculate Standard errors
FamHistse <- sqrt((FamHist[FamHist$Pheno=='MDDexc','se']^2) + (FamHist[FamHist$Pheno=='PostPartum','se']^2))

#Calculate z-score
FamHistz <- FamHistdelta/FamHistse

#Calculate p-value
FamHistpval <- 2*pnorm(-abs(FamHistz))

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#MDD PRS

#Subtract betas to get differences.
MDDPRSdelta <- MDDPRS[MDDPRS$Pheno=='MDDexc','beta'] - MDDPRS[MDDPRS$Pheno=='PostPartum','beta']

#Calculate Standard errors
MDDPRSse <- sqrt((MDDPRS[MDDPRS$Pheno=='MDDexc','se']^2) + (MDDPRS[MDDPRS$Pheno=='PostPartum','se']^2))

#Calculate z-score
MDDPRSz <- MDDPRSdelta/MDDPRSse

#Calculate p-value
MDDPRSpval <- 2*pnorm(-abs(MDDPRSz))

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#BPD PRS

#Subtract betas to get differences.
BPDPRSdelta <- BPDPRS[MDDPRS$Pheno=='MDDexc','beta'] - BPDPRS[MDDPRS$Pheno=='PostPartum','beta']

#Calculate Standard errors
BPDPRSse <- sqrt((BPDPRS[BPDPRS$Pheno=='MDDexc','se']^2) + (BPDPRS[BPDPRS$Pheno=='PostPartum','se']^2))

#Calculate z-score
BPDPRSz <- BPDPRSdelta/BPDPRSse

#Calculate p-value
BPDPRSpval <- 2*pnorm(-abs(BPDPRSz))

#########################################################################################################################################################################################
#########################################################################################################################################################################################

#Step 3: FDR correction

#Collate p-values
p_vals <- rbind(TDIpval, alevelpval, gcsepval, nonepval, Neuroticismpval, Childtraumapval, Adulttraumapval, FamHistpval, MDDPRSpval, BPDPRSpval)
p_valsFDR <- c(TDIpval, alevelpval, gcsepval, nonepval, Neuroticismpval, Childtraumapval, Adulttraumapval, FamHistpval, MDDPRSpval, BPDPRSpval)

sorted.pvalue<-sort(p_valsFDR) 
sorted.pvalue 
j.alpha <- (1:10)*(.05/10) # input total number of pvalues
diff <- sorted.pvalue-j.alpha 
neg.diff <- diff[diff<0] 
pos.diff <- neg.diff[length(neg.diff)] 
index <- diff==pos.diff 
p.cutoff <-sorted.pvalue[index] 
p.cutoff # this value is the critical pvalue. Values below should be considered significant
p.sig <- p_valsFDR[p_valsFDR <= p.cutoff] 
p.sig # which pvalues are significant based on the critical pvalue