#------------SURVIVAL GWAS IN OXFORD------------#


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

#Read in clinical dataset
#This is the  2019 version (sent by Michael Lawton 22/05/2019)
all_visits <- read_csv("../clinical_data/Discovery_all_FU_data_2019_07_29_v6_optout.csv")

#Read in list of individuals who want their genetic data withdrawn
remove_individuals <- read.table("../clinical_data/individuals_to_remove.txt")

#Note that this most recent dataset is missing time to withdrawal/ withdrawal dates
#Get variable names
#colnames <- colnames(all_visits)
#View(colnames)

#---Remove individuals who have withdrawn and want their data excluded---####

all_visits_removedIndividuals <- all_visits %>% 
  anti_join(remove_individuals, by = c("subjid" = "V1"))


#---Convert subject IDs to match format of genotype data (e.g. AH001)---####

#Split subject ID into 3 separate columsn then combine last two columns without / 
all_visits_removedIndividuals <- all_visits_removedIndividuals %>% 
  separate(subjid, c("subjid_1", "subjid_2", "subjid_3"), remove = FALSE) %>% 
  mutate(ID = paste(subjid_2, subjid_3, sep = ""))

#---PD probabilities---####

all_visits_removedIndividuals %>% 
  select(ID, pd_probability, pd_probability_bin)
#This dataset contains all individuals regardless of PD probability

#---Filter out non-PD cases---####

#Grouped variable (PD vs. non-PD)
all_visits_removedIndividuals %>% 
  group_by(alternative_diagnosis) %>% 
  summarise(count = n())

#Freetext variable
all_visits_removedIndividuals %>% 
  group_by(AltDiagnosis_type) %>% 
  summarise(count = n())

#Check for discrepancies between 2 variables
all_visits_removedIndividuals %>% 
  filter(is.na(AltDiagnosis_type)) %>% 
  filter(alternative_diagnosis == 1) %>% 
  select(ID, alternative_diagnosis, AltDiagnosis_type)

#Filter out non-PD cases using the alternative_diagnosis column
PD_only <- all_visits_removedIndividuals %>% 
  filter(is.na(alternative_diagnosis))


#---Calculate time from baseline for each visit---####

PD_only %>% 
  select(contains("visit_date"), contains("_age"))

#Calculate time from baseline to each visit - using age variables
PD_only <- PD_only %>% 
  mutate(time_V1_to_V2 = V2_age - V1_age,
         time_V1_to_V3 = V3_age - V1_age,
         time_V1_to_V4 = V4_age - V1_age,
         time_V1_to_V5 = V5_age - V1_age,
         time_V1_to_V6 = V6_age - V1_age)


#---Calculate time from PD diagnosis for each visit---####

#Calculate time from diagnosis to each visit
PD_only <- PD_only %>% 
  mutate(time_diag_to_V1 = V1_age - age_diag,
         time_diag_to_V2 = V2_age - age_diag,
         time_diag_to_V3 = V3_age - age_diag,
         time_diag_to_V4 = V4_age - age_diag,
         time_diag_to_V5 = V5_age - age_diag,
         time_diag_to_V6 = V6_age - age_diag)

#Summarise time from diagnosis
PD_only %>% 
  summarise(time_diag_to_V1_mean = mean(time_diag_to_V1),
            time_diag_to_V2_mean = mean(time_diag_to_V2, na.rm = TRUE),
            time_diag_to_V3_mean = mean(time_diag_to_V3, na.rm = TRUE))

#---Summarise age at onset and disease duration---####

#Mean age at onset
PD_only %>% 
  summarise(mean_AAO = mean(age_onset, na.rm = TRUE))

#Count how many individuals are missing AAO
PD_only %>% 
  filter(is.na(age_onset)) %>% 
  summarise(count = n())

#Count how many individuals are missing age at diagnosis
PD_only %>% 
  filter(is.na(age_diag)) %>% 
  summarise(count = n())

#Summarise disease duration variable
#Is this disease duration from diagnosis or onset?
PD_only %>% 
  summarise(mean_dd = mean(V1_disease_duration),
            V2_mean_dd = mean(V2_disease_duration, na.rm = TRUE),
            V3_mean_dd = mean(V3_disease_duration, na.rm = TRUE))
#Disease duration variable is from PD diagnosis NOT onset


#---Data tidying: Create imputed variables for age onset at follow-up visits---####

#Create new variable time from PD onset to diagnosis
PD_only <- PD_only %>% 
  mutate(time_onset_to_diag = age_diag - age_onset)

#Estimate missing age at onset from age at diagnosis using mean time from onset to diagnosis
PD_only %>% 
  summarise(mean_time_to_diagnosis = mean(time_onset_to_diag, na.rm = TRUE))
#Mean time from onset to diagnosis is 1.64 years

#Create new imputed age of onset variable which uses the mean time from onset to diagnosis for the other cases
PD_only <- PD_only %>% 
  mutate(age_onset_imput = ifelse(!is.na(age_onset), age_onset,
                                  ifelse(is.na(age_onset), age_diag - mean(time_onset_to_diag, na.rm = TRUE), NA)))

#Count if there are any more individuals missing age at onset with new imputed variable
PD_only %>% 
  filter(is.na(age_onset_imput)) %>% 
  summarise(count = n())

#---Data tidying: Create new variables for disease duration from PD symptom onset at each visit---####

#Create new variable for disease duration from PD symptom onset using each visit age minus the imputed age at onset
PD_only <- PD_only %>% 
  mutate(V1_disease_duration_onset = V1_age - age_onset_imput,
         V2_disease_duration_onset = V2_age - age_onset_imput,
         V3_disease_duration_onset = V3_age - age_onset_imput,
         V4_disease_duration_onset = V4_age - age_onset_imput,
         V5_disease_duration_onset = V5_age - age_onset_imput,
         V6_disease_duration_onset = V6_age - age_onset_imput)

#---Data tidying: Create new variable for grouped Hoehn and Yahr stage---####

#Create new variable for Hoehn and Yahr stage grouped - 0 to 1.5 vs 2 to 2.5 vs. 3+
PD_only <- PD_only %>% 
  mutate(V1_hoehn_and_yahr_grouped = ifelse(V1_af_hoehn_stage == 0 | V1_af_hoehn_stage == 1, "0 to 1",
                                            ifelse(V1_af_hoehn_stage == 2, "2",
                                                   ifelse(V1_af_hoehn_stage >=3, "3+", NA))),
         V2_hoehn_and_yahr_grouped = ifelse(V2_af_hoehn_stage == 0 | V2_af_hoehn_stage == 1, "0 to 1",
                                            ifelse(V2_af_hoehn_stage == 2, "2",
                                                   ifelse(V2_af_hoehn_stage >=3, "3+", NA))),
         V3_hoehn_and_yahr_grouped = ifelse(V3_af_hoehn_stage == 0 | V3_af_hoehn_stage == 1, "0 to 1",
                                            ifelse(V3_af_hoehn_stage == 2, "2",
                                                   ifelse(V3_af_hoehn_stage >=3, "3+", NA))),
         V4_hoehn_and_yahr_grouped = ifelse(V4_af_hoehn_stage == 0 | V4_af_hoehn_stage == 1, "0 to 1",
                                            ifelse(V4_af_hoehn_stage == 2, "2",
                                                   ifelse(V4_af_hoehn_stage >=3, "3+", NA))),
         V5_hoehn_and_yahr_grouped = ifelse(V5_af_hoehn_stage == 0 | V5_af_hoehn_stage == 1, "0 to 1",
                                            ifelse(V5_af_hoehn_stage == 2, "2",
                                                   ifelse(V5_af_hoehn_stage >=3, "3+", NA))),
         V6_hoehn_and_yahr_grouped = ifelse(V6_af_hoehn_stage == 0 | V6_af_hoehn_stage == 1, "0 to 1",
                                            ifelse(V6_af_hoehn_stage == 2, "2",
                                                   ifelse(V6_af_hoehn_stage >=3, "3+", NA))))


#---TABLE 1: Count of how many people have completed each visit---####

#Count of the number of PD patients who have completed each visit - using visit date variable
PD_only %>% 
  summarise(V1_completed = sum(!is.na(V1_visit_date)),
            V2_completed = sum(!is.na(V2_visit_date)),
            V3_completed = sum(!is.na(V3_visit_date)),
            V4_completed = sum(!is.na(V4_visit_date)),
            V5_completed = sum(!is.na(V5_visit_date)),
            V6_completed = sum(!is.na(V6_visit_date)),)

#---TABLE 1: Summary stats for MoCA scores and demographics at baseline and follow-up---####

#Summary stats for age at onset and age at entry at baseline
PD_only %>% 
  dplyr::summarise(AAO_mean = mean(age_onset_imput, na.rm = TRUE),
                   AAO_sd = sd(age_onset_imput, na.rm = TRUE),
                   agedx_mean = mean(age_diag, na.rm = TRUE),
                   agedx_sd = sd(age_diag, na.rm = TRUE),
                   age_mean = mean(V1_age, na.rm = TRUE),
                   age_sd = sd(V1_age, na.rm = TRUE))

#Mean disease duration at diagnosis
PD_only %>% 
  dplyr::summarise(disdur_dx_mean = mean(V1_disease_duration, na.rm = TRUE),
                   disdur_dx_sd = sd(V1_disease_duration, na.rm = TRUE),
                   disdur_onset_mean = mean(V1_disease_duration_onset, na.rm = TRUE),
                   disdur_onset_sd = sd(V1_disease_duration_onset, na.rm = TRUE))

#Mean age at visit
age_summary <- PD_only %>% 
  summarise(V1_age_mean = mean(V1_age, na.rm = TRUE),
            V1_age_sd = sd(V1_age, na.rm = TRUE),
            V2_age_mean = mean(V2_age, na.rm = TRUE),
            V2_age_sd = sd(V2_age, na.rm = TRUE),
            V3_age_mean = mean(V3_age, na.rm = TRUE),
            V3_age_sd = sd(V3_age, na.rm = TRUE),
            V4_age_mean = mean(V4_age, na.rm = TRUE),
            V4_age_sd = sd(V4_age, na.rm = TRUE),
            V5_age_mean = mean(V5_age, na.rm = TRUE),
            V5_age_sd = sd(V5_age, na.rm = TRUE),
            V6_age_mean = mean(V6_age, na.rm = TRUE),
            V6_age_sd = sd(V6_age, na.rm = TRUE))



#Mean disease duration from onset at each follow up visit
disduration_summary <- PD_only %>% 
  dplyr::summarise(disduration_mean = mean(V1_disease_duration_onset, na.rm = TRUE),
                   disduration_sd = sd(V1_disease_duration_onset, na.rm = TRUE),
                   V2_disduration_mean = mean(V2_disease_duration_onset, na.rm = TRUE),
                   V2_disduration_sd = sd(V2_disease_duration_onset, na.rm = TRUE),
                   V3_disduration_mean = mean(V3_disease_duration_onset, na.rm = TRUE),
                   V3_disduration_sd = sd(V3_disease_duration_onset, na.rm = TRUE),
                   V4_disduration_mean = mean(V4_disease_duration_onset, na.rm = TRUE),
                   V4_disduration_sd = sd(V4_disease_duration_onset, na.rm = TRUE),
                   V5_disduration_mean = mean(V5_disease_duration_onset, na.rm = TRUE),
                   V5_disduration_sd = sd(V5_disease_duration_onset, na.rm = TRUE),
                   V6_disduration_mean = mean(V6_disease_duration_onset, na.rm = TRUE),
                   V6_disduration_sd = sd(V6_disease_duration_onset, na.rm = TRUE))

#Mean time from diagnosis at each visit
time_from_dx_summary <- PD_only %>% 
  dplyr::summarise(time_from_dx_mean = mean(V1_disease_duration, na.rm = TRUE),
                   time_from_dx_sd = sd(V1_disease_duration, na.rm = TRUE),
                   V2_time_from_dx_mean = mean(V2_disease_duration, na.rm = TRUE),
                   V2_time_from_dx_sd = sd(V2_disease_duration, na.rm = TRUE),
                   V3_time_from_dx_mean = mean(V3_disease_duration, na.rm = TRUE),
                   V3_time_from_dx_sd = sd(V3_disease_duration, na.rm = TRUE),
                   V4_time_from_dx_mean = mean(V4_disease_duration, na.rm = TRUE),
                   V4_time_from_dx_sd = sd(V4_disease_duration, na.rm = TRUE),
                   V5_time_from_dx_mean = mean(V5_disease_duration, na.rm = TRUE),
                   V5_time_from_dx_sd = sd(V5_disease_duration, na.rm = TRUE),
                   V6_time_from_dx_mean = mean(V6_disease_duration, na.rm = TRUE),
                   V6_time_from_dx_sd = sd(V6_disease_duration, na.rm = TRUE))

#Summary stats for time from baseline to follow-up visits
PD_only %>% 
  dplyr::summarise(meanV2 = mean(time_V1_to_V2, na.rm = TRUE),
                   sdV2 = sd(time_V1_to_V2, na.rm = TRUE),
                   meanV3 = mean(time_V1_to_V3, na.rm = TRUE),
                   sdV3 = sd(time_V1_to_V3, na.rm = TRUE),
                   meanV4 = mean(time_V1_to_V4, na.rm = TRUE),
                   sdV4 = sd(time_V1_to_V4, na.rm = TRUE),
                   meanV5 = mean(time_V1_to_V5, na.rm = TRUE),
                   sdV5 = sd(time_V1_to_V5, na.rm = TRUE),
                   meanV6 = mean(time_V1_to_V6, na.rm = TRUE),
                   sdV6 = sd(time_V1_to_V6, na.rm = TRUE))

#Mean and SD MDS-UPDRS Part III at each timepoint
UPDRSIII_summary <- PD_only %>% 
  dplyr::summarise(UPDRSIII_mean = mean(V1_UPDRS_III, na.rm = TRUE),
                   UPDRSIII_sd = sd(V1_UPDRS_III, na.rm = TRUE),
                   V2_UPDRSIII_mean = mean(V2_UPDRS_III, na.rm = TRUE),
                   V2_UPDRSIII_sd = sd(V2_UPDRS_III, na.rm = TRUE),
                   V3_UPDRSIII_mean = mean(V3_UPDRS_III, na.rm = TRUE),
                   V3_UPDRSIII_sd = sd(V3_UPDRS_III, na.rm = TRUE),
                   V4_UPDRSIII_mean = mean(V4_UPDRS_III, na.rm = TRUE),
                   V4_UPDRSIII_sd = sd(V4_UPDRS_III, na.rm = TRUE),
                   V5_UPDRSIII_mean = mean(V5_UPDRS_III, na.rm = TRUE),
                   V5_UPDRSIII_sd = sd(V5_UPDRS_III, na.rm = TRUE),
                   V6_UPDRSIII_mean = mean(V6_UPDRS_III, na.rm = TRUE),
                   V6_UPDRSIII_sd = sd(V6_UPDRS_III, na.rm = TRUE))

#Mean MoCA at each timepoint
MOCA_summary <- PD_only %>% 
  dplyr::summarise(MOCA_mean = mean(V1_MOCA_total_adj, na.rm = TRUE),
                   MOCA_sd = sd(V1_MOCA_total_adj, na.rm = TRUE),
                   V2_MOCA_mean = mean(V2_MOCA_total_adj, na.rm = TRUE),
                   V2_MOCA_sd = sd(V2_MOCA_total_adj, na.rm = TRUE),
                   V3_MOCA_mean = mean(V3_MOCA_total_adj, na.rm = TRUE),
                   V3_MOCA_sd = sd(V3_MOCA_total_adj, na.rm = TRUE),
                   V4_MOCA_mean = mean(V4_MOCA_total_adj, na.rm = TRUE),
                   V4_MOCA_sd = sd(V4_MOCA_total_adj, na.rm = TRUE),
                   V5_MOCA_mean = mean(V5_MOCA_total_adj, na.rm = TRUE),
                   V5_MOCA_sd = sd(V5_MOCA_total_adj, na.rm = TRUE),
                   V6_MOCA_mean = mean(V6_MOCA_total_adj, na.rm = TRUE),
                   V6_MOCA_sd = sd(V6_MOCA_total_adj, na.rm = TRUE))


#Mean H&Y stage
HY_summary <- PD_only %>% 
  dplyr::summarise(HY_mean = mean(V1_af_hoehn_stage, na.rm = TRUE),
                   HY_sd =sd(V1_af_hoehn_stage, na.rm = TRUE),
                   V2_HY_mean = mean(V2_af_hoehn_stage, na.rm = TRUE),
                   V2_HY_sd = sd(V2_af_hoehn_stage, na.rm = TRUE),
                   V3_HY_mean = mean(V3_af_hoehn_stage, na.rm = TRUE),
                   V3_HY_sd = sd(V3_af_hoehn_stage, na.rm = TRUE),
                   V4_HY_mean = mean(V4_af_hoehn_stage, na.rm = TRUE),
                   V4_HY_sd = sd(V4_af_hoehn_stage, na.rm = TRUE),
                   V5_HY_mean = mean(V5_af_hoehn_stage, na.rm = TRUE),
                   V5_HY_sd = sd(V5_af_hoehn_stage, na.rm = TRUE),
                   V6_HY_mean = mean(V6_af_hoehn_stage, na.rm = TRUE),
                   V6_HY_sd = sd(V6_af_hoehn_stage, na.rm = TRUE))


#H&Y stage proportions at baseline
PD_only %>% 
  group_by(V1_hoehn_and_yahr_grouped) %>% 
  summarise(count = n()) %>% 
  mutate(prop = count/ sum(!is.na(PD_only$V1_hoehn_and_yahr_grouped))*100)

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

#H&Y stage proportions at V6
PD_only %>% 
  group_by(V6_hoehn_and_yahr_grouped) %>% 
  summarise(count = n()) %>% 
  mutate(prop = count/ sum(!is.na(PD_only$V6_hoehn_and_yahr_grouped))*100)


#Semantic fluency scores
seman_flu_summary  <- PD_only %>% 
  dplyr::summarise(V1_sf_mean = mean(V1_seman_flu_score, na.rm = TRUE),
                   V1_sf_sd =sd(V1_seman_flu_score, na.rm = TRUE),
                   V2_sf_mean = mean(V2_seman_flu_score, na.rm = TRUE),
                   V2_sf_sd = sd(V2_seman_flu_score, na.rm = TRUE),
                   V3_sf_mean = mean(V3_seman_flu_score, na.rm = TRUE),
                   V3_sf_sd = sd(V3_seman_flu_score, na.rm = TRUE),
                   V4_sf_mean = mean(V4_seman_flu_score, na.rm = TRUE),
                   V4_sf_sd = sd(V4_seman_flu_score, na.rm = TRUE),
                   V5_sf_mean = mean(V5_seman_flu_score, na.rm = TRUE),
                   V5_sf_sd = sd(V5_seman_flu_score, na.rm = TRUE),
                   V6_sf_mean = mean(V6_seman_flu_score, na.rm = TRUE),
                   V6_sf_sd = sd(V6_seman_flu_score, na.rm = TRUE))



#---Summary of education variables---####

PD_only %>% 
  select(i_2_further_ed, i_2_1_further_ed_years, ed_years_before_further_ed,
         education_years_bin)

#---Gender summaries---####

PD_only %>% 
  group_by(gender) %>% 
  summarise(count = n())


########## SUMMARISE WITHDRAWAL AND DEATH DATA ########## 
#---Summarise withdrawal data---####

#Count number of cases that have withdrawn
PD_only %>% 
  group_by(withdrawn) %>% 
  summarise(count = n())

#Look at freetext withdrawal reasons
withdrawals <- PD_only %>% 
  group_by(Withdrawn_incl_reason) %>% 
  select(ID, dead, Withdrawn_incl_reason, withdrawal_date)

#Unsure what D2 participant status indicates
PD_only %>% 
  select(D2participantstatus)

#Check number of individuals who have died
PD_only %>% 
  group_by(dead) %>% 
  summarise(count = n())

#---Encode new category for death as there are some participants not correctly coded---####

#Text patterns to match
to_match_death <- c("passed")
#Most deaths are coded correctly in the variable 'dead'
#But there were a few cases where the withdrawal reason said the participant passed away and this was not coded
#Have not included the word 'died' to match in the withdrawal reason as this may refer to a relative/spouse dying rather than the participant

#Code new variable for death
PD_only <- PD_only %>% 
  mutate(death = ifelse(grepl(paste(to_match_death, collapse = "|"),
                              Withdrawn_incl_reason), 1,
                        ifelse(is.na(dead), 0,
                               ifelse(dead == 1, 1, "check"))))

test <- PD_only %>% 
  select(Withdrawn_incl_reason, dead, death)

#Check numbers - there were 2 cases that were recoded
PD_only %>% 
  group_by(death, dead) %>% 
  summarise(count = n())

#---Encode withdrawal reasons that include dementia/cognitive impairment---####

#Text patterns to search for in withdrawal reason
to_match <- c("dementia", "PDD", "cognit", "capacity", "memory", "confus")

#Create new variable 'dementia' which codes for dementia in withdrawal reason
PD_only <- PD_only %>% 
  mutate(dementia = ifelse((grepl(paste(to_match, collapse = "|"),
                                  Withdrawn_incl_reason)), "TRUE", "FALSE"))

check <- PD_only %>% 
  select(ID, Withdrawn_incl_reason, dementia)

#There are one or two individuals where the withdrawal reason does not indicate dementia but there are pattern matches
#E.g. 'no evidence of dementia' is being coded as dementia
#E.g. 'found followups very confusing' - again this is not really dementia
#Need to identify these individuals and reocde on individual basis
PD_only <- PD_only %>% 
  mutate(dementia = ifelse(ID == "HH059" | ID == "MK062", "FALSE", dementia))


#Search for other text patterns that might indicate dementia
PD_only %>% 
  filter(grepl("confus", Withdrawn_incl_reason)) %>% 
  select(Withdrawn_incl_reason)

#---Convert all dates to date format---####

#Convert all dates to correct date format
datecols <- PD_only %>% 
  select(contains("date"))


#For loop
#This converts all the visit dates into actual dates from text format
for (nm in colnames(datecols)) {
  PD_only[[nm]] <- as.Date(PD_only[[nm]], format = "%d%b%Y")
}

#Check that this is all formatted correctly
PD_only %>% 
  select(contains("date"))


#---Check that withdrawal dates are within a sensible range---####

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
  mutate(time_V1_to_withdrawal = as.duration(withdrawal_date - V1_visit_date),
         age_withdrawal = dyears(V1_age) + time_V1_to_withdrawal,
         time_onset_to_withdrawal = age_withdrawal - dyears(age_onset_imput)) 

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
  filter(dyears(time_V1_to_lastvisit) > time_V1_to_withdrawal) %>% 
  select(ID, time_V1_to_lastvisit, withdrawal_date, time_V1_to_withdrawal)
#If withdrawal date is earlier than the last visit date, use the last visit date (the most recent date) as the time to event/censoring

#Can check individuals who have withdrawal dates earlier than the last visit date
#Sometimes this is because the withdrawal date has been incorrectly entered
# check <- PD_only %>% 
#   filter(ID == "HH045") %>% 
#   select(contains("date"), Withdrawn_incl_reason, withdrawn)

#Create variable for time from PD onset to last visit date
PD_only <- PD_only %>% 
  mutate(time_onset_to_lastvisit = ifelse(is.na(time_V1_to_lastvisit), V1_disease_duration_onset,
                                          time_V1_to_lastvisit + V1_disease_duration_onset))

#There are some individuals who are missing withdrawal dates and time to last visit
#These are individuals who have withdrawn but there is no withdrawal date
#Or individuals who have just completed Visit 1 but have not withdrawn
missing_dates <- PD_only %>% 
  filter(is.na(time_V1_to_lastvisit) & is.na(time_V1_to_withdrawal)) %>% 
  select(ID, time_V1_to_withdrawal, time_V1_to_lastvisit, time_onset_to_lastvisit, withdrawal_date, withdrawn,
         V1_age, V2_age, V3_age, V4_age, V5_age, V6_age)
#There are only 22 individuals so just exclude these?

########## DEMENTIA OUTCOMES ########## 
#---MoCA expected scores have been calculated in Discovery---####
#The expected score variables are calculated by using the mean of the answered questions where at least 80% of the questions are answered.  

#Check expected scores
PD_only %>% 
  select(V1_MOCA_total_adj, V1_MOCA_total_adj_exp) %>% 
  mutate(check = ifelse(V1_MOCA_total_adj_exp == V1_MOCA_total_adj, "ok", "check")) %>% 
  group_by(check) %>% 
  filter(is.na(check)) %>% 
  select(V1_MOCA_total_adj_exp, V1_MOCA_total_adj)

#---Calculate how many individuals meet MDS-UPDRS 1.1 and MOCA dementia criteria at baseline---####
#Using MoCA raw scores

#MoCA <=21 vs. all MDS-UPDRS 1.1 options
PD_only %>% 
  mutate(V1_MOCA_21cut = ifelse(V1_MOCA_total_adj <= 21, "MOCA 21 or less", 
                                ifelse(V1_MOCA_total_adj > 21, "MOCA greater than 21", NA))) %>% 
  group_by(V1_ae_1_1_cognitive, V1_MOCA_21cut) %>% 
  summarise(count = n())


#---Create variables for event and time to event dementia MoCA <=21 or withdrawal due to dementia---####
#Event is the first point to reach dementia MoCA <= 21 (first visit where this is met)
#Time to event is time from PD onset

#First select the Hoehn and Yahr variables at all visits
moca <- PD_only %>% 
  select(ID, ends_with("MOCA_total_adj"))

#Create new variable for event H&Y stage 3 if HY stage was >=3 at any time
moca <- moca %>% 
  mutate(count_NA = rowSums(is.na(moca))) %>% 
  mutate(any_MOCA21 = ifelse(V1_MOCA_total_adj <=21 | V2_MOCA_total_adj <=21 | V3_MOCA_total_adj <=21 | V4_MOCA_total_adj <=21 | V5_MOCA_total_adj <=21 | V6_MOCA_total_adj <=21, 1, 0)) %>% 
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
                                  select(ID, event_MOCA21, contains("disease_duration_onset"), ends_with("MOCA_total_adj")))

disease_duration_onset_names <- colnames(PD_only_MOCA21 %>% 
                                           select(contains("disease_duration_onset")))

MOCA_names <- colnames(PD_only_MOCA21 %>% 
                         select(contains("MOCA_total_adj")))

#Reshape to long format
PD_only_MOCA21_long <- reshape(PD_only_MOCA21, idvar = "ID", direction = "long",
                               varying = list(Start = disease_duration_onset_names, Value = MOCA_names),
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
  mutate(time_onset_to_lastvisit_duration = dyears(time_onset_to_lastvisit)) %>% #Convert time to last visit to a duration
  mutate(timeToEvent_MOCA21= ifelse(is.na(time_onset_to_withdrawal) & !is.na(time_onset_to_lastvisit_duration), time_onset_to_lastvisit_duration,
                                    ifelse(is.na(time_onset_to_lastvisit_duration) & !is.na(time_onset_to_withdrawal), time_onset_to_withdrawal,
                                           ifelse(time_onset_to_lastvisit_duration >= time_onset_to_withdrawal, time_onset_to_lastvisit_duration,
                                                  ifelse(time_onset_to_lastvisit_duration < time_onset_to_withdrawal, time_onset_to_withdrawal, NA))))) %>% 
  mutate(timeToEvent_MOCA21 = as.duration(timeToEvent_MOCA21)/dyears(1)) %>% 
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

#---Summary of dementia outcome---####

PD_only %>% 
  group_by(event_dementia) %>% 
  summarise(count = n(),
            mean_time = mean(timeToEvent_dementia),
            sd_time = sd(timeToEvent_dementia))

#---Create survival object for time to dementia---####

surv_object_dementia <- Surv(time = PD_only$timeToEvent_dementia, event = PD_only$event_dementia)

#---Plot Kaplan-Meier curve for gender vs. dementia---####

fit_gender_dementia <- survfit(surv_object_dementia ~ gender, data = PD_only)

summary(fit_gender_dementia)

#Plot Kaplan Meier curve for gender vs. mortality
ggsurvplot(fit_gender_dementia, data = PD_only, pval = TRUE)

#Fit Cox Proportional Hazards model for dementia vs. gender
fit_gender_dementia_coxph = coxph(surv_object_dementia ~ gender, data = PD_only)
summary(fit_gender_dementia_coxph)

#Fit Cox Proportional Hazards model for mortality vs. gender with covariates
fit_gender_dementia_coxph_covars = coxph(surv_object_dementia ~ gender + age_onset_imput, data = PD_only)
summary(fit_gender_dementia_coxph_covars)

#---Export clinical data for GWAS survival analysis of dementia---####

#Need to export: event, time to event, age at onset, gender
Oxford_export_dementia <- PD_only %>% 
  mutate(IID = paste(ID, "_", ID, sep = ""),
         FID = IID) %>% 
  select(IID, FID, event_dementia, timeToEvent_dementia, age_onset_imput, gender) %>% 
  filter(!is.na(timeToEvent_dementia)) #Remove individuals who are missing time to event

#Write as output
filename_dememtia <- paste("outputs/Oxford_survival_dementia_", today(), ".txt", sep ="")

write.table(Oxford_export_dementia, filename_dememtia,
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")


