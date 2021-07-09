
for n in {0..898}
do
	echo "#!/usr/bin/Rscript
	library(data.table)
	library(survival)
	library(dplyr)

	#Load genetic data
	subset<- fread(\"/data/kronos/kronos/for_Maryam/PROBAND/file$n.raw\")
	#genetic data starts from column 7 of the raw file. 
	subset<- subset[,-c(3,4,5,6)]
	#Load clinical data
	clinical <- fread(\"/data/kronos/kronos/for_Maryam/PROBAND/PROBAND_survival_dementia_2020-07-01.txt\")

	#Load Principal Components
	PCs <- fread(\"/data/kronos/kronos/mtan/survival_GWAS/PROBAND/PCA.eigenvec\")

	#Select just the first 5 principal components
	PCs <- PCs %>%
		select(V1:V12) %>%
		rename(FID = V1,
			IID = V2,
			PC1 = V3,
			PC2 = V4,
			PC3 = V5,
			PC4 = V6,
			PC5 = V7,
			PC6 = V8,
			PC7 = V9,
			PC8 = V10,
			PC9 = V11,
			PC10 = V12)
	
	#Inner join all datasets (only keep individuals with values in all
	#Individuals missing the clinical data have been removed from the clinical dataset so they should be removed after the join
	TABLE <- clinical %>%
		inner_join(PCs, by = c(\"IID\", \"FID\")) %>%
		inner_join(subset, by = c(\"IID\" = \"FID\"))


	TABLE<- TABLE[, -22]
	TABLE<- as.data.frame(TABLE)
	
	coefficients<-as.data.frame(matrix(ncol= 9))
 
	names(coefficients) <- c(\"SNP\",\"Coeff\", \"se\", \"Pvalue\", \"Cox.zphPVal\", \"N\", \"ov.lik.ratio\",\"logrank\", \"r2\" )
 
	for (i in 22:ncol(TABLE)) {
	print(colnames(TABLE)[i])
  	snp <- TABLE[,c(i,1:21)]
  	model.cox<- coxph(Surv(snp\$timeToEvent_dementia, snp\$event_dementia) ~ snp[,1] + snp\$PC1+ snp\$PC2 + snp\$PC3 + snp\$PC4 + snp\$PC5+ snp\$PC6+ snp\$PC7+ snp\$PC8+ snp\$PC9+ snp\$PC10 + snp\$years_education_bin, data=snp)
  	kmz<- cox.zph(model.cox, transform = \"km\")
  	j= i-21
  	coefficients[j,1]<- paste(colnames(TABLE)[i])
  	coefficients[j,2]<- summary(model.cox)\$coefficients[1,1]
  	coefficients[j,3]<- summary(model.cox)\$coefficients[1,3]
  	coefficients[j,4]<- summary(model.cox)\$coefficients[1,5]
  	coefficients[j,5]<- kmz\$table[1,3]
  	coefficients[j,6]<- model.cox\$n
  	coefficients[j,7]<- summary(model.cox)\$logtest[[1]]
  	coefficients[j,8]<- summary(model.cox)\$sctest[[1]]
  	coefficients[j,9]<- summary(model.cox)\$rsq[[1]] # nagelkerke r square
  	}
 
	fwrite(coefficients, \"/data/kronos/kronos/for_Maryam/PROBAND/DEMENTIA/PROBAND_DEMENTIA_cxphres_$n.txt\", row.names=FALSE, sep=\"\t\", quote= FALSE)" > Dementia_PROBAND_$n.r
done