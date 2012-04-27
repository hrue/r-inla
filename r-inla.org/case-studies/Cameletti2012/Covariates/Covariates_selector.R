covariates_selector_funct = function(i_day, mean_covariates, sd_covariates){

	#--- Select the covariates for the chosen day and scale them using the mean and the sd computed on the station covariates values --> you get a 56x72 matrix
	
	#--- NB: the order of the covariates in mean/sd_covariates is: A, UTMX, UTMY, WS, TEMP, HMIX, PREC, EMI 
	
	WindSpeedGRID_i	= (WindSpeedGRID[,,i_day] - mean_covariates[4]) /sd_covariates[4]
	HMixMaxGRID_i	= (HMixMaxGRID[,,i_day] - mean_covariates[6]) / sd_covariates[6]	
	EmiGRID_i		= (EmiGRID[,,i_day] - mean_covariates[8]) / sd_covariates[8]
	Mean_Temp_i		= (Mean_Temp[,,i_day] - mean_covariates[5]) / sd_covariates[5]
	Prec_i			= (Prec[,,i_day] - mean_covariates[7]) / sd_covariates[7]
	Prec_i[is.na(Prec_i)]=0 #when all prec for day i are 0 the scale proc generates NA

	UTMX_std = (unique(Piemonte_grid[,1]) - mean_covariates[2]) / sd_covariates[2]
	UTMX_GRID = matrix(rep(UTMX_std,72),56,72) #constant for each row

	UTMY_std = (unique(Piemonte_grid[,2]) - mean_covariates[3]) / sd_covariates[3]
	UTMY_GRID = t(matrix(rep(UTMY_std,56),72,56)) #constant for each column

	AltitudeGRID = (AltitudeGRID - mean_covariates[1]) / sd_covariates[1]

	#--- Create the array (56x72x8) of standardized covariates for day i
	covariate_array_std = abind(AltitudeGRID,
		UTMX_GRID,
		UTMY_GRID,	
		WindSpeedGRID_i,
		Mean_Temp_i,
		HMixMaxGRID_i,
		Prec_i,
		EmiGRID_i,
		along=3)
	#dim(covariate_array_std)
	return(covariate_array_std)
	
	#filename = paste("covariate_array_std","_day",i,".Rdata",sep="")
	#save(covariate_array_std,file = filename)

}
