### META-ANALYSIS SCRIPT ###

#---Load libraries---####
library(dplyr)
library(meta)
library(grid)
library(meta)
library(data.table)
library(readxl)


#---Read in data---####

#Read in PROBAND selected areas file - this has the full data, not all SNPs are in the excel tables
PROBAND_selected_areas <- fread("PROBAND_Selected_areas.txt")

#Read in Oxford data
Oxford_selected_areas <- fread("Oxford_Selected_areas.txt")

#---Read in list of rsIDs and positions---####

rsids <- read_excel("final_table_positions.xlsx", col_names = F)
colnames(rsids) <- c("SNP", "chrbp", "comments")

snp_rsids <- rsids$SNP
snp_positions <- rsids$chrbp

#---Meta-analysis---####

#Make results dataframe
re_metaanalysis_results <- as.data.frame(matrix(ncol = 11))
colnames(re_metaanalysis_results) <- c("rsid", "SNP", "beta", "se", "HR", "z", "pval", "HetIsq", "CochransQ_pval",
                                       "lower_CI", "upper_CI")

#From Table in paper

for (i in 1:length(snp_positions)) {
  
  snp <- snp_positions[i]
  
  PROBAND <- PROBAND_selected_areas %>% 
    filter(grepl(snp, SNP)) %>% 
    mutate(cohort = "PROBAND")
  
  Oxford <- Oxford_selected_areas %>% 
    filter(grepl(snp, SNP)) %>% 
    mutate(cohort = "Oxford")
  
  merged <- rbind(PROBAND, Oxford)
  
  #Do meta-analysis using random effects model
  meta <- metagen(Coeff,
                  se,
                  data = merged,
                  studlab = cohort,
                  comb.fixed = TRUE,
                  comb.random = TRUE,
                  method.tau = "DL",
                  prediction = TRUE,
                  sm = "HR")

  re_metaanalysis_results[i,1] <- snp_rsids[i]
  re_metaanalysis_results[i,2] <- snp
  re_metaanalysis_results[i,3] <- meta$TE.random
  re_metaanalysis_results[i,4] <- meta$seTE.random
  re_metaanalysis_results[i,5] <- exp(meta$TE.random)
  re_metaanalysis_results[i,6] <- meta$zval.random
  re_metaanalysis_results[i,7] <- meta$pval.random
  re_metaanalysis_results[i,8] <- meta$I2
  re_metaanalysis_results[i,9] <- meta$pval.Q
  re_metaanalysis_results[i,10] <- meta$lower.random
  re_metaanalysis_results[i,11] <- meta$upper.random
  
  
    
  
}

#If the SNP is only present in one study, then remove the meta-analysis HR and p-value
#Round to 2 dp
#Calculate 95% CI for Hazard Ratio
re_metaanalysis_results <- re_metaanalysis_results %>% 
  mutate(HR_final = round(ifelse(is.na(HetIsq), NA, HR),2),
         pval_final = round(ifelse(is.na(HetIsq), NA, pval),2),
         CI_HR_lower = round(exp(lower_CI),2),
         CI_HR_upper = round(exp(upper_CI),2))

#Save as table
write.table(re_metaanalysis_results, "./Manuela_outputs/meta_analysis_results.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")


