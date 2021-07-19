#------------SURVIVAL GWAS IN PROBAND------------#
#Using new 2020 version of clinical data

#---Load packages---####
library(ggplot2)
library(readstata13)
library(gridExtra)
library(gtable)
library(grid)
library(reshape2)
library(pROC)
library(ROCR)
library(readxl)
library(survival)
library(survminer)
library(lubridate)
library(lme4)
library(car)
library(MASS)
library(factoextra)
library(data.table)
library(tidyverse)
library(ggpubr)
library(Hmisc)
library(corrplot)

#---Load functions---####

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}


## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

#Normalise function
normFunc <- function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
}

#---Load clinical data---####

#Read in STATA dataset
#This is the new 2020 version (sent by Sofia Kanavou 17/06/2020)
all_visits <- as_tibble(read_csv("../clinical_dataset_2020/P3_TPD_V2_temp_06_17_2020/P3_TPD_Version2_temp_06_16_2020.csv"))

#Get variable names
colnames <- colnames(all_visits)
View(colnames)

#---Load education data---####

#Read in raw education data
education <- read_excel("../clinical_dataset_new_version2019/P3_raw_education_data_unenc.xlsx")

#Remove years_education_bin variable as this is already in the main dataset
education <- education %>% 
  select(-years_education_bin)

#Merge with main clinical dataset
all_visits <- all_visits %>% 
  left_join(education, by = "ID")

#Convert higher education binary variable into factor
all_visits$higher_education <- as.factor(all_visits$higher_education)

#---Load codebreak---####

codebreak <- read_excel("../genetic_data/Proband.Geno.Pheno.xlsx")

#---Merge clinical data with codebreak and fam file from genetic data---####

all_visits <- all_visits %>% 
  left_join(codebreak, by = c("ID" = "Kate.Patient.ID"))

#Read in fam file
fam <- read.table("../genetic_data/HTS_HD_Chips_MAF_GENO_HWE_MIND-updated.fam")

#Export list of cases in clinical dataset -merge with fam file, get FID and IID, then export as text file 
export_PDcases <- all_visits %>% 
  select(IID) %>% 
  filter(!is.na(IID))

export_PDcases_fam <- export_PDcases %>% 
  left_join(fam, by = c("IID" = "V2")) %>% 
  select(V1, IID) %>% 
  filter(!is.na(V1))

#---Filter out non-PD cases---####

#Look at different diagnoses
changeDiagnosis <- all_visits %>% 
  dplyr::group_by(change_diagnosis) %>% 
  dplyr::summarise(count = n())

#Count number of cases who received a different diagnosis
all_visits %>% 
  filter(!is.na(change_diagnosis))

#Filter out non-PD cases
PD_only <- all_visits %>% 
  filter(is.na(change_diagnosis))

#---Calculate time from baseline for each visit---####

PD_only %>% 
  select(contains("visit_date"), contains("_age"))

#Calculate time from baseline to each visit - using age variables
PD_only <- PD_only %>% 
  mutate(time_V1_to_V4 = age_V4 - age_V1,
         time_V1_to_V7 = age_V7 - age_V1,
         time_V1_to_V9 = age_V9 - age_V1,
         time_V1_to_V10 = age_V10 - age_V1,
         time_V1_to_V11 = age_V11 - age_V1)

#---Calculate time from PD diagnosis for each visit---####

#This is already calculated in the disease_duration_diag_V* variables

#---Calculate time from PD onset for each visit---####

PD_only <- PD_only %>% 
  mutate(V1_disease_duration_onset = age_V1 - age_onset,
         V4_disease_duration_onset = age_V4 - age_onset,
         V7_disease_duration_onset = age_V7 - age_onset,
         V10_disease_duration_onset = age_V10 - age_onset,
         V11_disease_duration_onset = age_V11 - age_onset)

#---Create annual change scores for MoCA and MDS-UPDRSIII---####

PD_only <- PD_only %>% 
  mutate(UPDRS_III_BLto18_annual = (V4_UPDRS_III_total - V1_UPDRS_III_total)/time_V1_to_V4,
         UPDRS_III_BLto36_annual = (V7_UPDRS_III_total - V1_UPDRS_III_total)/time_V1_to_V7,
         UPDRS_III_BLto54_annual = (V9_UPDRS_III_total - V1_UPDRS_III_total)/time_V1_to_V9,
         UPDRS_III_BLto72_annual = (V10_UPDRS_III_total - V1_UPDRS_III_total)/time_V1_to_V10,
         MOCA_BLto18_annual = (V4_MOCA_total_adj - V1_MOCA_total_adj)/time_V1_to_V4,
         MOCA_BLto36_annual = (V7_MOCA_total_adj - V1_MOCA_total_adj)/time_V1_to_V7,
         MOCA_BLto54_annual = (V9_MOCA_total_adj - V1_MOCA_total_adj)/time_V1_to_V9,
         MOCA_BLto72_annual = (V10_MOCA_total_adj - V1_MOCA_total_adj)/time_V1_to_V10)

#---Data tidying: Create imputed variables for disease duration, age onset, disease duration at follow-up visits---####

#Estimate missing age at onset from age at diagnosis using mean time from onset to diagnosis
PD_only %>% 
  summarise(mean_time_to_diagnosis = mean(time_onset_to_diag, na.rm = TRUE))
#Mean time from onset to diagnosis is 1.81 years

#Count how many people are missing age at onset
PD_only %>% 
  filter(is.na(age_onset)) %>% 
  summarise(count = n())

#Make new variable for age at diagnosis
PD_only <- PD_only %>% 
  mutate(age_diag = age_V1 - disease_duration_diag_V1)

#Create imputed age of onset variable using the mean time from onset to diagnosis from non-missing cases
PD_only <- PD_only %>% 
  mutate(age_onset_imput = ifelse(!is.na(age_onset), age_onset,
                                  ifelse(is.na(age_onset), age_diag - mean(time_onset_to_diag, na.rm = TRUE), NA)),
         disease_duration_onset_imput = ifelse(!is.na(V1_disease_duration_onset), V1_disease_duration_onset,
                                               ifelse(is.na(V1_disease_duration_onset), age_V1 - age_onset_imput, NA)))

PD_only <- PD_only %>% 
  mutate(V4_disease_duration_onset_imput = age_V4 - age_onset_imput,
         V7_disease_duration_onset_imput = age_V7 - age_onset_imput,
         V9_disease_duration_onset_imput = age_V9 - age_onset_imput,
         V10_disease_duration_onset_imput = age_V10 - age_onset_imput,
         V11_disease_duration_onset_imput = age_V11 - age_onset_imput)

#---Data tidying: Create new variable for grouped Hoehn and Yahr stage---####

#Create new variable for Hoehn and Yahr stage grouped - 0 to 1.5 vs 2 to 2.5 vs. 3+
PD_only <- PD_only %>% 
  mutate(V1_hoehn_and_yahr_grouped = ifelse(V1_hoehn_and_yahr_stage == 0 | V1_hoehn_and_yahr_stage == 1 | V1_hoehn_and_yahr_stage == 1.5, "0 to 1.5",
                                            ifelse(V1_hoehn_and_yahr_stage == 2 | V1_hoehn_and_yahr_stage == 2.5, "2 to 2.5",
                                                   ifelse(V1_hoehn_and_yahr_stage == 3 | V1_hoehn_and_yahr_stage == 4 | V1_hoehn_and_yahr_stage == 5, "3+", NA)))) %>% 
  mutate(V2_hoehn_and_yahr_grouped = ifelse(V4_hoehn_and_yahr_stage == 0 | V4_hoehn_and_yahr_stage == 1 | V4_hoehn_and_yahr_stage == 1.5, "0 to 1.5",
                                            ifelse(V4_hoehn_and_yahr_stage == 2 | V4_hoehn_and_yahr_stage == 2.5, "2 to 2.5",
                                                   ifelse(V4_hoehn_and_yahr_stage == 3 | V4_hoehn_and_yahr_stage == 4 | V4_hoehn_and_yahr_stage == 5, "3+", NA)))) %>% 
  mutate(V3_hoehn_and_yahr_grouped = ifelse(V7_hoehn_and_yahr_stage == 0 | V7_hoehn_and_yahr_stage == 1 | V7_hoehn_and_yahr_stage == 1.5, "0 to 1.5",
                                            ifelse(V7_hoehn_and_yahr_stage == 2 | V7_hoehn_and_yahr_stage == 2.5, "2 to 2.5",
                                                   ifelse(V7_hoehn_and_yahr_stage == 3 | V7_hoehn_and_yahr_stage == 4 | V7_hoehn_and_yahr_stage == 5, "3+", NA)))) %>% 
  mutate(V4_hoehn_and_yahr_grouped = ifelse(V9_hoehn_and_yahr_stage == 0 | V9_hoehn_and_yahr_stage == 1 | V9_hoehn_and_yahr_stage == 1.5, "0 to 1.5",
                                            ifelse(V9_hoehn_and_yahr_stage == 2 | V9_hoehn_and_yahr_stage == 2.5, "2 to 2.5",
                                                   ifelse(V9_hoehn_and_yahr_stage == 3 | V9_hoehn_and_yahr_stage == 4 | V9_hoehn_and_yahr_stage == 5, "3+", NA)))) %>% 
  mutate(V5_hoehn_and_yahr_grouped = ifelse(V10_hoehn_and_yahr_stage == 0 | V10_hoehn_and_yahr_stage == 1 | V10_hoehn_and_yahr_stage == 1.5, "0 to 1.5",
                                            ifelse(V10_hoehn_and_yahr_stage == 2 | V10_hoehn_and_yahr_stage == 2.5, "2 to 2.5",
                                                   ifelse(V10_hoehn_and_yahr_stage == 3 | V10_hoehn_and_yahr_stage == 4 | V10_hoehn_and_yahr_stage == 5, "3+", NA))))

#---TABLE 1: Count of how many people have completed each visit---####

#Count of the number of PD patients who have completed each visit - using age variable
PD_only %>% 
  summarise(V1_completed = sum(!is.na(age_V1)),
            V2_completed = sum(!is.na(age_V4)),
            V3_completed = sum(!is.na(age_V7)),
            V4_completed = sum(!is.na(age_V9)),
            V5_completed = sum(!is.na(age_V10)),
            V6_completed = sum(!is.na(age_V11)))

#---TABLE 1: Summary stats for MoCA scores and demographics at baseline and follow-up---####

#Summary stats for age at onset and age at entry at baseline
PD_only %>% 
  dplyr::summarise(AAO_mean = mean(age_onset_imput, na.rm = TRUE),
                   AAO_sd = sd(age_onset_imput, na.rm = TRUE),
                   agedx_mean = mean(age_diag, na.rm = TRUE),
                   agedx_sd = sd(age_diag, na.rm = TRUE),
                   age_mean = mean(age_V1, na.rm = TRUE),
                   age_sd = sd(age_V1, na.rm = TRUE))

#Mean disease duration at diagnosis
PD_only %>% 
  dplyr::summarise(disdur_dx_mean = mean(disease_duration_diag_V1, na.rm = TRUE),
                   disdur_dx_sd = sd(disease_duration_diag_V1, na.rm = TRUE),
                   disdur_onset_mean = mean(disease_duration_onset_imput, na.rm = TRUE),
                   disdur_onset_sd = sd(disease_duration_onset_imput, na.rm = TRUE))

#Mean age at visit
age_summary <- PD_only %>% 
  summarise(V1_age_mean = mean(age_V1, na.rm = TRUE),
            V1_age_sd = sd(age_V1, na.rm = TRUE),
            V2_age_mean = mean(age_V4, na.rm = TRUE),
            V2_age_sd = sd(age_V4, na.rm = TRUE),
            V3_age_mean = mean(age_V7, na.rm = TRUE),
            V3_age_sd = sd(age_V7, na.rm = TRUE),
            V4_age_mean = mean(age_V9, na.rm = TRUE),
            V4_age_sd = sd(age_V9, na.rm = TRUE),
            V5_age_mean = mean(age_V10, na.rm = TRUE),
            V5_age_sd = sd(age_V10, na.rm = TRUE))


#Summary stats for time from baseline to follow-up visits
PD_only %>% 
  dplyr::summarise(meanV2 = mean(time_V1_to_V4, na.rm = TRUE),
                   sdV2 = sd(time_V1_to_V4, na.rm = TRUE),
                   meanV3 = mean(time_V1_to_V7, na.rm = TRUE),
                   sdV3 = sd(time_V1_to_V7, na.rm = TRUE),
                   meanV4 = mean(time_V1_to_V9, na.rm = TRUE),
                   sdV4 = sd(time_V1_to_V9, na.rm = TRUE),
                   meanV5 = mean(time_V1_to_V10, na.rm = TRUE),
                   sdV5 = sd(time_V1_to_V10, na.rm = TRUE))

#Mean and SD MDS-UPDRS Part III at each timepoint
UPDRSIII_summary <- PD_only %>% 
  dplyr::summarise(UPDRSIII_mean = mean(V1_UPDRS_III_total, na.rm = TRUE),
                   UPDRSIII_sd = sd(V1_UPDRS_III_total, na.rm = TRUE),
                   V2_UPDRSIII_mean = mean(V4_UPDRS_III_total, na.rm = TRUE),
                   V2_UPDRSIII_sd = sd(V4_UPDRS_III_total, na.rm = TRUE),
                   V3_UPDRSIII_mean = mean(V7_UPDRS_III_total, na.rm = TRUE),
                   V3_UPDRSIII_sd = sd(V7_UPDRS_III_total, na.rm = TRUE),
                   V4_UPDRSIII_mean = mean(V9_UPDRS_III_total, na.rm = TRUE),
                   V4_UPDRSIII_sd = sd(V9_UPDRS_III_total, na.rm = TRUE),
                   V5_UPDRSIII_mean = mean(V10_UPDRS_III_total, na.rm = TRUE),
                   V5_UPDRSIII_sd = sd(V10_UPDRS_III_total, na.rm = TRUE))

#Mean MoCA at each timepoint
MOCA_summary <- PD_only %>% 
  dplyr::summarise(MOCA_mean = mean(V1_MOCA_total_adj, na.rm = TRUE),
                   MOCA_sd = sd(V1_MOCA_total_adj, na.rm = TRUE),
                   V2_MOCA_mean = mean(V4_MOCA_total_adj, na.rm = TRUE),
                   V2_MOCA_sd = sd(V4_MOCA_total_adj, na.rm = TRUE),
                   V3_MOCA_mean = mean(V7_MOCA_total_adj, na.rm = TRUE),
                   V3_MOCA_sd = sd(V7_MOCA_total_adj, na.rm = TRUE),
                   V4_MOCA_mean = mean(V9_MOCA_total_adj, na.rm = TRUE),
                   V4_MOCA_sd = sd(V9_MOCA_total_adj, na.rm = TRUE),
                   V5_MOCA_mean = mean(V10_MOCA_total_adj, na.rm = TRUE),
                   V5_MOCA_sd = sd(V10_MOCA_total_adj, na.rm = TRUE))


#Mean H&Y stage at baseline
HY_summary <- PD_only %>% 
  dplyr::summarise(HY_mean = mean(V1_hoehn_and_yahr_stage, na.rm = TRUE),
                   HY_sd =sd(V1_hoehn_and_yahr_stage, na.rm = TRUE),
                   V2_HY_mean = mean(V4_hoehn_and_yahr_stage, na.rm = TRUE),
                   V2_HY_sd = sd(V4_hoehn_and_yahr_stage, na.rm = TRUE),
                   V3_HY_mean = mean(V7_hoehn_and_yahr_stage, na.rm = TRUE),
                   V3_HY_sd = sd(V7_hoehn_and_yahr_stage, na.rm = TRUE),
                   V4_HY_mean = mean(V9_hoehn_and_yahr_stage, na.rm = TRUE),
                   V4_HY_sd = sd(V9_hoehn_and_yahr_stage, na.rm = TRUE),
                   V5_HY_mean = mean(V10_hoehn_and_yahr_stage, na.rm = TRUE),
                   V5_HY_sd = sd(V10_hoehn_and_yahr_stage, na.rm = TRUE))


#H&Y stage proportions at baseline
PD_only %>% 
  group_by(V1_hoehn_and_yahr_grouped) %>% 
  summarise(count = n())

#H&Y stage proportions at V2
PD_only %>% 
  group_by(V2_hoehn_and_yahr_grouped) %>% 
  summarise(count = n()) %>% 
  mutate(prop = count/ sum(!is.na(PD_only$V2_hoehn_and_yahr_grouped))*100)

#H&Y stage proportions at V3
PD_only %>% 
  group_by(V3_hoehn_and_yahr_grouped) %>% 
  summarise(count = n()) %>% 
  mutate(prop = count/ sum(!is.na(PD_only$V3_hoehn_and_yahr_grouped))*100)

#H&Y stage proportions at V4
PD_only %>% 
  group_by(V4_hoehn_and_yahr_grouped) %>% 
  summarise(count = n()) %>% 
  mutate(prop = count/ sum(!is.na(PD_only$V4_hoehn_and_yahr_grouped))*100)


#H&Y stage proportions at V5
PD_only %>% 
  group_by(V5_hoehn_and_yahr_grouped) %>% 
  summarise(count = n()) %>% 
  mutate(prop = count/ sum(!is.na(PD_only$V5_hoehn_and_yahr_grouped))*100)


#Semantic fluency scores
seman_flu_summary  <- PD_only %>% 
  dplyr::summarise(V1_sf_mean = mean(V1_seman_flu_score, na.rm = TRUE),
                   V1_sf_sd =sd(V1_seman_flu_score, na.rm = TRUE),
                   V2_sf_mean = mean(V4_seman_flu_score, na.rm = TRUE),
                   V2_sf_sd = sd(V4_seman_flu_score, na.rm = TRUE),
                   V3_sf_mean = mean(V7_seman_flu_score, na.rm = TRUE),
                   V3_sf_sd = sd(V7_seman_flu_score, na.rm = TRUE),
                   V4_sf_mean = mean(V9_seman_flu_score, na.rm = TRUE),
                   V4_sf_sd = sd(V9_seman_flu_score, na.rm = TRUE),
                   V5_sf_mean = mean(V10_seman_flu_score, na.rm = TRUE),
                   V5_sf_sd = sd(V10_seman_flu_score, na.rm = TRUE))


#---Summary of education variables---####

#Age left school summary
PD_only %>% 
  group_by(ed_years_before_further_ed) %>% 
  summarise(count = n())

#Age left school summary
PD_only %>% 
  group_by(higher_education) %>% 
  summarise(count = n())

#Composite variable
education_summary <- PD_only %>% 
  select(ed_years_before_further_ed, higher_education, years_education_bin)

#Code years_education_bin as a category
#12 years or less, >12 years or higher education
PD_only <- PD_only %>% 
  mutate(years_education_bin = ifelse(years_education_bin == 0, "more than 12 years or higher education",
                                      ifelse(years_education_bin == 1, "12 years or less", NA)))

PD_only %>% 
  group_by(years_education_bin, higher_education) %>% 
  summarise(count = n())

#---Gender summaries---####

PD_only %>% 
  group_by(gender) %>% 
  summarise(count = n())

#Code gender as text
PD_only <- PD_only %>% 
  mutate(gender = ifelse(gender == 0, "male",
                         ifelse(gender == 1, "female", NA)))

PD_only$gender <- as.factor(PD_only$gender)

#Change reference level
PD_only$gender <- relevel(PD_only$gender, ref = "male")

#Age of onset by gender
PD_only %>%
  group_by(gender) %>% 
  summarise(mean = mean(age_onset_imput, na.rm = TRUE))

#Create categories for EOPD vs LOPD using AAO 50
PD_only <- PD_only %>% 
  mutate(earlyonset_50 = ifelse(age_onset_imput <=50, "early onset",
                                ifelse(age_onset_imput > 50, "late onset", NA)))

PD_only %>% 
  group_by(gender, earlyonset_50) %>% 
  summarise(count = n())

#---Create age at onset groups---####

#Create groups for age at onset
PD_only <- PD_only %>% 
  mutate(aao_group = ifelse(age_onset_imput <= 45, "<=45",
                            ifelse(age_onset_imput > 45 & age_onset_imput <= 50, "45-50",
                                   ifelse(age_onset_imput > 50 & age_onset_imput <=55, "50-55",
                                          ifelse(age_onset_imput > 55 & age_onset_imput <=60, "55-60",
                                                 ifelse(age_onset_imput > 60 & age_onset_imput <=65, "60-65",
                                                        ifelse(age_onset_imput > 65 & age_onset_imput <=70, "65-70",
                                                               ifelse(age_onset_imput > 70 & age_onset_imput <=75, "70-75",
                                                                      ifelse(age_onset_imput > 75, "75+", NA)))))))))


PD_only %>% 
  group_by(aao_group) %>% 
  summarise(count = n())

#---Summary of disease duration at baseline---####

#Mean disease duration at baseline (time from onset to study entry)
PD_only %>% 
  summarise(mean_dd = mean(disease_duration_onset_imput, na.rm = TRUE))

#Plot histogram to see distribution
ggplot(data = PD_only, mapping = aes(disease_duration_onset_imput)) +
  geom_histogram(color="black", fill="white") +
  theme_bw()

#Create groups for disease duration at baseline
PD_only <- PD_only %>% 
  mutate(diseaseDuration_group = ifelse(disease_duration_onset_imput < 1, "0 to 1 year",
                                        ifelse(disease_duration_onset_imput >=1 & disease_duration_onset_imput < 2, "1 to 2 years",
                                               ifelse(disease_duration_onset_imput >=2 & disease_duration_onset_imput < 3, "2 to 3 years",
                                                      ifelse(disease_duration_onset_imput >=3, "3+ years", NA)))))

#Summary of groups
PD_only %>% 
  group_by(diseaseDuration_group) %>% 
  summarise(count = n())

#Count how many people have disease duration > 5 years
PD_only %>% 
  filter(disease_duration_onset_imput > 5) %>% 
  summarise(count = n())

#Count how many people have time from diagnosis > 3.5 years
PD_only %>% 
  filter(disease_duration_diag_V1 > 3.5) %>% 
  summarise(count = n())

#---Summary of Hoehn and Yahr stage at baseline---####

PD_only %>% 
  group_by(V1_hoehn_and_yahr_stage) %>% 
  summarise(count = n())

#Make grouped variable for HY2+ at baseline
PD_only <- PD_only %>% 
  mutate(BL_HY2 = ifelse(V1_hoehn_and_yahr_stage >=2, "BL HY2+",
                         ifelse(V1_hoehn_and_yahr_stage < 2, "BL HY<2", NA)))

#Check groups
PD_only %>% 
  group_by(V1_hoehn_and_yahr_stage, BL_HY2) %>% 
  summarise(count = n())

#---Calculate motor subtypes---####

#From Movement Disorders paper
#Tremor: mean of 2.10, 3.15a. 3.15b, 3.16a, 3.16b, 3.17a, 3.17b, 3.17c, 3.17d, 3.17e, 3.18
#PIGD: mean of 2.12, 2.13, 3.10, 3.11, 3.12
#Tremor score divided by PIGD score
#If ratio is >=1.15 = TD
#If ratio is <=0.90 = PIGD
#If ratio is between 0.90 and 1.15 = indeterminate

PD_only %>% 
  select(V1_UPDRS_III_15a)

#Calculate tremor dominant score
#Note that at the moment there is no leeway for missing scores - so if there is one item missing, cannot calculate TD score
TDscore <- PD_only %>% 
  select(ID, V1_UPDRS_II_10, V1_UPDRS_III_15a, V1_UPDRS_III_15b,
         V1_UPDRS_III_16a, V1_UPDRS_III_16b, 
         V1_UPDRS_III_17a, V1_UPDRS_III_17b, V1_UPDRS_III_17c, 
         V1_UPDRS_III_17d, V1_UPDRS_III_17e,
         V1_UPDRS_III_18) %>% 
  mutate(TDsum = rowSums(select(., contains("UPDRS")))) %>% 
  mutate(TDscore = TDsum/11)

#Write output to check with IPDGC spreadsheet
#write.csv(PIGDscore, "outputs/TDscore.csv", quote = F, row.names = F)

TDscore_merge <- TDscore %>% 
  select(ID, TDsum, TDscore)

#Calculate PIGD score
#PIGD: mean of 2.12, 2.13, 3.10, 3.11, 3.12
PIGDscore <- PD_only %>% 
  select(ID, V1_UPDRS_II_12, V1_UPDRS_II_13, V1_UPDRS_III_10, 
         V1_UPDRS_III_11, V1_UPDRS_III_12) %>% 
  mutate(PIGDsum = rowSums(select(., contains("UPDRS")))) %>% 
  mutate(PIGDscore = PIGDsum/5)

#Write output to check with IPDGC spreadsheet
#write.csv(PIGDscore, "outputs/PIGDscore.csv", quote = F, row.names = F)

PIGDscore_merge <- PIGDscore %>% 
  select(ID, PIGDsum, PIGDscore)


#Merge with main dataset
PD_only <- PD_only %>% 
  left_join(TDscore_merge, by = "ID") %>% 
  left_join(PIGDscore_merge, by = "ID")

#Calculate traditional motor subtype ratio which is TDscore/PIGDscore
#Assign motor subtype categories - if the PIGD score is 0, then subtype should be TD even though the ratio cannot be calculated
#However if both the TD and PIGD score is 0, then subtype should be indeterminate
PD_only <- PD_only %>% 
  mutate(subtype_ratio = TDscore/PIGDscore) %>% 
  mutate(subtype_cat = ifelse(PIGDscore == 0 & TDscore > 0, "TD",
                              ifelse(PIGDscore == 0 & TDscore == 0, "indeterminate",
                                     ifelse(subtype_ratio >=1.15, "TD",
                                            ifelse(subtype_ratio <= 0.90, "PIGD",
                                                   ifelse(subtype_ratio > 0.9 & subtype_ratio < 1.15, "indeterminate", NA))))))


#Summarise motor subtype categories
PD_only %>% 
  group_by(subtype_cat) %>% 
  summarise(count = n())

#Check against existing motor subtype variable
PD_only %>% 
  group_by(V1_UPDRS_phenotype_cat, subtype_cat) %>% 
  summarise(count = n())

#---Create grouping variable for MoCA score at baseline---####

PD_only %>% 
  group_by(V1_MOCA_total_adj) %>% 
  summarise(count = n())

#Create grouped variable - categories are arbitrary (this is just for illustration in the KM curve)
#Loosely based on Hu 2014 paper
PD_only <- PD_only %>% 
  mutate(V1_MOCA_total_adj_grouped = ifelse(V1_MOCA_total_adj <=21, "21 or less",
                                            ifelse(V1_MOCA_total_adj > 21 & V1_MOCA_total_adj <= 25, "22 to 25",
                                                   ifelse(V1_MOCA_total_adj >= 26, "26 to 30", NA))))

#Summary of categories
PD_only %>% 
  group_by(V1_MOCA_total_adj_grouped) %>% 
  summarise(count = n())

########## SUMMARISE WITHDRAWAL AND DEATH DATA ########## 

#---Summarise withdrawal data---####

#Count number of cases that have withdrawn
PD_only %>% 
  group_by(Withdrawn) %>% 
  summarise(count = n())

#Summarise by categories of withdrawal reason
withdrawals <- PD_only %>% 
  group_by(withd_res, withd_other) %>% 
  summarise(count = n())

#Look at the detail withdrawal reasons (freetext)
PD_only %>% 
  group_by(withd_other) %>% 
  summarise(count = n())

#---Encode withdrawal reasons that include dementia/cognitive impairment---####

#Text patterns to search for in withdrawal reason
to_match <- c("dementia", "PDD", "cognit", "capacity", "memory", "confus")

#Create new variable 'dementia' which codes for dementia in withdrawal reason
PD_only <- PD_only %>% 
  mutate(dementia = ifelse(is.na(withd_res), "FALSE",
                           ifelse((grepl(paste(to_match, collapse = "|"),
                                         withd_other)), "TRUE",
                                  ifelse(withd_res == "Developed dementia", "TRUE", "FALSE"))))

check <- PD_only %>% 
  select(withd_res, withd_other, dementia)

#Search for other text patterns that might indicate dementia
PD_only %>% 
  filter(grepl("dement", withd_other)) %>% 
  select(withd_other)


#---Check that withdrawal dates are within a sensible range---####

test <- PD_only %>% 
  select(Withdrawn, withdrawal_date)

#Convert withdrawal date into an actual date
PD_only <- PD_only %>% 
  mutate(withdrawal_date = as.Date(withdrawal_date, format = "%d%b%Y"))

#Correct dates out of range
#Sensible range is anything from 2012 to today
sensible_range <- interval(ymd("2012-01-01"), today())

PD_only %>% 
  filter(!withdrawal_date %within% sensible_range) %>% 
  select(ID, withdrawal_date)

#Correct dates that are outside sensible range (add 100 years to the one incorrect date)

#---Calculate time to withdrawal---####

#Convert dates to date format and calculate time in years from study entry to study withdrawal
PD_only <- PD_only %>% 
  mutate(age_withdrawal = age_V1 + time_baseline_to_withdrawal,
         time_onset_to_withdrawal = age_withdrawal - age_onset_imput) 

check <- PD_only %>% 
  select(age_V1, age_onset_imput, age_withdrawal, time_baseline_to_withdrawal, time_onset_to_withdrawal)

#---Check last visit date vs. withdrawal date---####

#Select time to each visit variables
visittimes <- PD_only %>% 
  select(ID, contains("time_V1_to_V"))

#Find the last time to visit value (the most recent follow-up visit)
lastvisittime <- visittimes %>% 
  gather(key = Key, value = Value, -ID) %>% 
  group_by(ID) %>% 
  filter(!is.na(Value)) %>% 
  slice(n()) %>% 
  select(-Key)

#Create time from study entry to most recent follow-up visit
PD_only <- PD_only %>% 
  left_join(lastvisittime, by = "ID") %>% 
  rename(time_V1_to_lastvisit = "Value")
#Note that if individuals withdrew after the first visit, there is no time to last visit value

#Sensibility check - all withdrawal dates should be after the last visit date (mostly)
check_lastDate <- PD_only %>% 
  filter(time_V1_to_lastvisit > time_baseline_to_withdrawal) %>% 
  select(ID, time_V1_to_lastvisit, withdrawal_date, time_baseline_to_withdrawal)
#These are fine - the last visit time is about the same as withdrawal time
#If withdrawal date is earlier than the last visit date, use the last visit date (the most recent date) as the time to event/censoring

#Create variable for time from PD onset to last visit date
PD_only <- PD_only %>% 
  mutate(time_onset_to_lastvisit = ifelse(is.na(time_V1_to_lastvisit), disease_duration_onset_imput,
                                          time_V1_to_lastvisit + disease_duration_onset_imput))

#There are some individuals who are missing withdrawal dates and time to last visit
missing_dates <- PD_only %>% 
  filter(is.na(time_V1_to_lastvisit) & is.na(time_baseline_to_withdrawal)) %>% 
  select(ID, time_baseline_to_withdrawal, time_V1_to_lastvisit, time_onset_to_lastvisit, withdrawal_date, Withdrawn, age_V1, age_V4, age_V7, age_V9, age_V10)
#There are only 2 individuals so just exclude these?
#These are individuals who have not withdrawn but the last visit is V1

########## DEMENTIA OUTCOMES ########## 
#---Calculate how many individuals meet MDS-UPDRS 1.1 and MOCA dementia criteria at baseline---####

#MoCA <=21 vs. all MDS-UPDRS 1.1 options
PD_only %>% 
  mutate(V1_MOCA_21cut = ifelse(V1_MOCA_total_adj <= 21, "MOCA 21 or less", 
                                ifelse(V1_MOCA_total_adj > 21, "MOCA greater than 21", NA))) %>% 
  group_by(V1_UPDRS_I_1, V1_MOCA_21cut) %>% 
  summarise(count = n())


#---Create variables for event and time to event dementia MoCA <=21 or withdrawal due to dementia---####
#Event is the first point to reach dementia MoCA <= 21 (first visit where this is met)
#Time to event is time from PD onset

#First select the Hoehn and Yahr variables at all visits
moca <- PD_only %>% 
  select(ID, ends_with("MOCA_total_adj"))

#Create new variable for event dementia if MOCA was <=21 at any time
moca <- moca %>% 
  mutate(count_NA = rowSums(is.na(moca))) %>% 
  mutate(any_MOCA21 = ifelse(V1_MOCA_total_adj <=21 | V4_MOCA_total_adj <=21 | V7_MOCA_total_adj <=21 | V9_MOCA_total_adj <=21 | V10_MOCA_total_adj <=21 | V11_MOCA_total_adj <=21, 1, 0)) %>% 
  mutate(event_MOCA21 = ifelse(count_NA == 6 & is.na(any_MOCA21), NA, #If all MoCAs are missing, MA
                               ifelse(count_NA < 6 & is.na(any_MOCA21), 0, #If at least one MoCA is done
                                      ifelse(count_NA < 6 & any_MOCA21 == 0, 0,
                                             ifelse(any_MOCA21 == 1, 1, "check"))))) %>%
  select(ID, event_MOCA21)

moca %>% 
  group_by(event_MOCA21) %>% 
  summarise(count = n())

#Merge with main dataset
PD_only <- PD_only %>% 
  left_join(moca, by = "ID")

#Code any individuals who have withdrawn because of cognitive impairment/dementia as meeting outcome
PD_only <- PD_only %>% 
  mutate(event_dementia = ifelse(dementia == "TRUE" & event_MOCA21 == 1, 1,
                                 ifelse(dementia == "TRUE" & event_MOCA21 == 0, 1,
                                        ifelse(dementia == "FALSE" & event_MOCA21 == 1, 1,
                                               ifelse(dementia == "FALSE" & event_MOCA21 == 0, 0, NA)))))

#Count how many individuals have withdrawn because of dememtia
PD_only %>% 
  group_by(dementia, event_MOCA21) %>% 
  summarise(count = n())

#Count how many individuals met dementia outcome (both withdrawal and MoCA<=21)
PD_only %>% 
  group_by(event_dementia) %>% 
  summarise(count = n())

#Create time to event for dementia
#For people who meet the outcome MoCA <=21, it is the first visit that this was met
PD_only_MOCA21 <- as.data.frame(PD_only %>% 
                                  filter(event_MOCA21 == 1) %>% 
                                  select(ID, event_MOCA21, disease_duration_onset_imput, contains("disease_duration_onset_imput"), ends_with("MOCA_total_adj")))

MOCA_names <- colnames(PD_only_MOCA21 %>% 
                         select(contains("MOCA_total_adj")))

disease_duration_onset_imput_names <- colnames(PD_only_MOCA21 %>% 
                                                 select(contains("disease_duration_onset_imput")))

#Reshape to long format
PD_only_MOCA21_long <- reshape(PD_only_MOCA21, idvar = "ID", direction = "long",
                               varying = list(Start = disease_duration_onset_imput_names, Value = MOCA_names),
                               v.names = c("time_from_onset", "MOCA_total_adj"))

#Arrange by ID and time from onset
PD_only_MOCA21_long <- PD_only_MOCA21_long %>% 
  arrange(ID, time_from_onset)

#Keep just the first instance of when MOCA<=21 was met
PD_only_yesMOCA21 <- PD_only_MOCA21_long %>% 
  filter(MOCA_total_adj <=21) %>% 
  distinct(ID, .keep_all = TRUE) %>% 
  mutate(timeToEvent_MOCA21 = time_from_onset) %>% 
  select(ID, timeToEvent_MOCA21)

#For people who do not meet the outcome, the time to event is either the last visit or the date of withdrawal
PD_only_noMOCA21 <- PD_only %>% 
  filter(event_MOCA21 == 0) %>% 
  #Use the most recent time (the longest duration) as the time to event - either the time to last visit, or time to withdrawal
  mutate(timeToEvent_MOCA21= ifelse(is.na(time_onset_to_withdrawal) & !is.na(time_onset_to_lastvisit), time_onset_to_lastvisit,
                                    ifelse(is.na(time_onset_to_lastvisit) & !is.na(time_onset_to_withdrawal), time_onset_to_withdrawal,
                                           ifelse(time_onset_to_lastvisit >= time_onset_to_withdrawal, time_onset_to_lastvisit,
                                                  ifelse(time_onset_to_lastvisit < time_onset_to_withdrawal, time_onset_to_withdrawal, NA))))) %>% 
  select(ID, timeToEvent_MOCA21)

#Merge the time to event data for people who have and have not met the outcome HY3
PD_only_MOCA21_timeToEvent <- rbind(PD_only_yesMOCA21, PD_only_noMOCA21)

#Merge with main dataset
#There are a handful of people who are missing time to event data - these are people who do not have any MOCA data
PD_only <- PD_only %>% 
  left_join(PD_only_MOCA21_timeToEvent, by = "ID")

#Code time to event for dementia (including people who withdrew due to dementia)
#time to MoCA <=21 should always be less than the time to withdrawal
#So time to dementia is the same as time to MoCA <=21
#For individuals who do not meet MoCA <=21 but withdrew due to dementia, the time to event is time to withdrawal which is already coded in timeToEvent_MOCA21
PD_only <- PD_only %>% 
  mutate(timeToEvent_dementia = timeToEvent_MOCA21)


#---Summary of dementia outcomes---####

PD_only %>% 
  group_by(event_dementia) %>% 
  summarise(count = n(),
            mean_time = mean(timeToEvent_dementia),
            sd_time = sd(timeToEvent_dementia))

#---Create survival object for dementia---####

surv_object_dementia <- Surv(time = PD_only$timeToEvent_dementia, event = PD_only$event_dementia)

#---Plot Kaplan-Meier curve for gender vs. dementia---####

fit_gender_dementia <- survfit(surv_object_dementia ~ gender, data = PD_only)

summary(fit_gender_dementia)

#Plot Kaplan Meier curve for gender vs. dementia
ggsurvplot(fit_gender_dementia, data = PD_only, pval = TRUE) +
  ggsave("plots/PROBAND_dementia_gender.png")

#Fit Cox Proportional Hazards model for dementia vs. gender
fit_gender_dementia_coxph = coxph(surv_object_dementia ~ gender, data = PD_only)
summary(fit_gender_dementia_coxph)

#Fit Cox Proportional Hazards model for dementia vs. gender with covariates
fit_gender_dementia_coxph_covars = coxph(surv_object_dementia ~ gender + age_onset_imput, data = PD_only)
summary(fit_gender_dementia_coxph_covars)

#---Merge with PROBAND imputed data fam file---####

#Read in fam file from final imputed dataset
imputed_fam <- read.table("../genetic_data/new_allchromosomes.converted.R2_0.8.MAF_0.01.fam")

imputed_fam <- imputed_fam %>% 
  separate(V1, c("FID", "IID", "pos")) %>% 
  mutate(IID_new = paste(IID, pos, sep = "_"))

#Merge with clinical dataset
#This is to get the correct IID for the PROBAND samples which are FID_IID
PD_only <- PD_only %>% 
  left_join(imputed_fam, by = c("IID" = "IID_new"))


#---Export clinical data for GWAS survival analysis of dementia---####

#Need to export: event, time to event, age at onset, gender
PROBAND_export_dementia <- PD_only %>% 
  select(ID, IID, event_dementia, timeToEvent_dementia, age_onset_imput, gender, V2,
         disease_duration_onset_imput, BL_HY2, subtype_cat, V1_MOCA_total_adj, years_education_bin) %>% 
  filter(!is.na(IID)) %>% #Remove individuals who are missing IID
  filter(!is.na(V2)) %>% #Remove individuals who are missing final IID in the fam file
  select(-ID, -IID) %>% 
  mutate(FID = V2) %>% 
  rename(IID = V2) %>% 
  filter(!is.na(timeToEvent_dementia)) #Remove individuals who are missing time to event

#Write as output
filename_dememtia <- paste("outputs/PROBAND_survival_dementia_", today(), ".txt", sep ="")

write.table(PROBAND_export_dementia, filename_dememtia,
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")


