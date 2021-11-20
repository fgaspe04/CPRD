# #print(commandArgs(trailingOnly = TRUE))

args=commandArgs(trailingOnly = TRUE)
j      = as.integer(args[1])
gender = as.character(args[2]) #g can be male or female


LOCAL = FALSE

library(data.table)
library(nlme)
library(survival)
library(tidyverse)
library(haven)
library(feather)
library(dplyr)
library(miceadds)
library(dplyr)
library(survival)
library(pec)
library(rms)
library(riskRegression)
library(risksetROC)
library(survAUC)
library(dynpred)


source("mixoutsamp_v3.R")
source("pred_error_check.R")


if( LOCAL ){
  j = 60
  load.Rdata( paste0( "../../data_created/male/data_lm_",j,"_male.RData"), paste0( "data_ext") )
  path_to_save = "../../results/results_male/"
  mod_type = "RIRS_"
  load.Rdata( paste0( "../../results/results_male/test_fit_rirs",j,"_m.RData"), paste0( "fit_rirs") )
}else{
  load.Rdata( paste0( "../../data_created/", gender, "/data_lm_",j, "_", gender,".RData"), paste0( "data_ext") )
  mod_type = "RIRS_"
  load.Rdata( paste0( "../../results/results_", gender, "/fit_LMEM_rirs_",j,"_", gender,".RData"), paste0( "fit_rirs") )
  path_to_save = paste0("../../results/results_", gender, "/" )
}



#######################
#
# Type of analysis
#
#######################

#Number of investigated outcomes
outcomes_names = c('bmi','hdl','sbp', 'smokbin','tchol')
N_outcomes = length( outcomes_names )

#GLOBAL Horizon
GLOB_horizon = 10


#N subjects
N_patid_all = length( unique( data_ext$patid ) ) 
N_patid_der = length( unique( data_ext$patid[ which( data_ext$derivation == "derivation" ) ] ) ) 
N_patid_val = length( unique( data_ext$patid[ which( data_ext$derivation == "validation" ) ] ) ) 

patid_all = unique( data_ext$patid ) 
patid_der = unique( data_ext$patid[ which( data_ext$derivation == "derivation" ) ] ) 
patid_val = unique( data_ext$patid[ which( data_ext$derivation == "validation" ) ] ) 

#setting covariates for the cox model
#no more SLE ind and family ind nor Dementia
cox_cov_less60 = c( "diab_ind + bp_bin + hdl_blup_FE4 + sbp_blup_FE4 + smoke_blup_FE4 + tchol_blup_FE4 + bmi_blup_FE4 + townsend_20  +  depression_ind + Severe_mental_illness_ind + Migraine_ind ")
cox_cov_more60 = c( "diab_ind + bp_bin  + hdl_blup_FE4 + sbp_blup_FE4 + smoke_blup_FE4 + tchol_blup_FE4 + bmi_blup_FE4 + townsend_20 + renal_disease + atrial_fibrillation + rheumatoid_arthritis +  depression_ind + Severe_mental_illness_ind + Migraine_ind")#Dementia_ind +


if( j >= 60 )
{
  cox_cov = cox_cov_more60
  cox_cov_list = c("diab_ind", "bp_bin", "hdl_blup_FE4", "sbp_blup_FE4", "smoke_blup_FE4", "tchol_blup_FE4", "bmi_blup_FE4", "townsend_20", "renal_disease", "atrial_fibrillation", "rheumatoid_arthritis",   "depression_ind", "Severe_mental_illness_ind",  "Migraine_ind" ) #"Dementia_ind",
}else{
  cox_cov = cox_cov_less60
  cox_cov_list = c("diab_ind", "bp_bin", "hdl_blup_FE4", "sbp_blup_FE4", "smoke_blup_FE4", "tchol_blup_FE4", "bmi_blup_FE4", "townsend_20",  "depression_ind", "Severe_mental_illness_ind", "Migraine_ind"  )
}

N_cox_cov = length( cox_cov_list )

#######################
#
# Risk prediction est
#
#######################

risk_prediction =  data_ext[ exp_count == 1, .( patid, exp_age, statin_age, cvd_ind, cvd_age, death_ind, death_age, end_age ) ]

#at the baseline nobody was still prescribed statins


#######################
#
# Cumulative hazard est
#
#######################

cum_haz_est = matrix( 0, nrow = length( unique( data_ext$exp_age ) ), ncol = 6 )
cum_haz_est = as.data.table(cum_haz_est)

#######################
#
# Covariates matrix
#
#######################

covariates_matrix = as.data.table( unique( data_ext$patid ) )
#colnames( covariates_matrix ) = "patid"

cov_rec_before_lm_name = c("bp_bin", "statin_bin", "renal_disease",
                           "atrial_fibrillation", "rheumatoid_arthritis",
                           "Severe_mental_illness_ind", "Migraine_ind", "Dementia_ind", "depression_ind",
                           "SLE_ind")

cov_rec_before_lm = c( "bp_med_age", "statin_age", "renal_age",
                       "Atrial_fibrillation_age", "rheumatoid_arthritis_age",
                       "Severe_mental_illness_age", "Migraine_age", "Dementia_age", "depression_age",
                       "SLE_age")

covariates_matrix                = data_ext[ exp_count == 1, lapply(.SD, function(x) ifelse( x <= j & !is.na(x), 1, 0 )), .SDcols = cov_rec_before_lm  ]
colnames( covariates_matrix )    = cov_rec_before_lm_name
covariates_matrix                = cbind( covariates_matrix, data_ext[ exp_count == 1, .( patid, townsend_20, gen_ethnicity, family_history, derivation ) ] )
covariates_matrix$diab_ind       = data_ext[ exp_count == 1, ifelse( diab_age <= j  & !is.na(diab_age) & diab_type == 2 & !is.na(diab_type), 1, 0 ) ]

setcolorder( covariates_matrix, c(  "patid", "derivation", "diab_ind", "bp_bin",
                                    "statin_bin", "townsend_20", "gen_ethnicity",
                                    "renal_disease", "atrial_fibrillation",
                                    "rheumatoid_arthritis", "Severe_mental_illness_ind",
                                    "Migraine_ind", "Dementia_ind", 
                                    "depression_ind", "SLE_ind", "family_history" ))


#adding space for predicted time-varying covariates
covariates_matrix = cbind( covariates_matrix, as.data.table( matrix( 0, nrow = length( unique( data_ext$patid ) ), ncol = N_outcomes*( GLOB_horizon + 1 ) ) ) )
colnames( covariates_matrix ) = c( "patid", "derivation", "diab_ind", "bp_med", "statin_bin",
                                   'townsend_20', 'ethnia','renal_disease','atrial_fibrillation','rheumatoid_arthritis',
                                   'Severe_mental_illness_ind','Migraine_ind', 'Dementia_ind', 'depression_ind', 
                                   'SLE_ind', 'family_history',
                                   'bmi_blup_0', 'hdl_blup_0',  'sbp_blup_0', 'smoke_blup_0', 'tchol_blup_0',
                                   'bmi_blup_1', 'hdl_blup_1',   'sbp_blup_1',  'smoke_blup_1', 'tchol_blup_1',
                                   'bmi_blup_2', 'hdl_blup_2',   'sbp_blup_2',  'smoke_blup_2', 'tchol_blup_2',
                                   'bmi_blup_3', 'hdl_blup_3',   'sbp_blup_3',  'smoke_blup_3', 'tchol_blup_3',
                                   'bmi_blup_4', 'hdl_blup_4',   'sbp_blup_4',  'smoke_blup_4', 'tchol_blup_4',
                                   'bmi_blup_5', 'hdl_blup_5',   'sbp_blup_5',  'smoke_blup_5', 'tchol_blup_5',
                                   'bmi_blup_6', 'hdl_blup_6',   'sbp_blup_6',  'smoke_blup_6', 'tchol_blup_6',
                                   'bmi_blup_7', 'hdl_blup_7',   'sbp_blup_7',  'smoke_blup_7', 'tchol_blup_7',
                                   'bmi_blup_8', 'hdl_blup_8',   'sbp_blup_8',  'smoke_blup_8', 'tchol_blup_8',
                                   'bmi_blup_9', 'hdl_blup_9',   'sbp_blup_9',  'smoke_blup_9', 'tchol_blup_9',
                                   'bmi_blup_10', 'hdl_blup_10',  'sbp_blup_10', 'smoke_blup_10', 'tchol_blup_10'
)




pt_analysed = matrix( 0, nrow = GLOB_horizon + 1, ncol = 11 )


#computing fixed part of prediction
length( unique(data_ext$patid[ which( data_ext$exp_age_corr <= 0 ) ] ))
#SOME PEOPLE HAVE OBSERVATION ONLY IN THE FUTURE!!!! (25% at 60 yo)
length( unique(data_ext$patid[ which( data_ext$exp_age_corr <= 0 ) ] ))/ length(unique(data_ext$patid))

#All the predictions will be carried out from t_0 = Landamark age
dati_baseline = data.table( 
  patid = unique( data_ext$patid ),
  #bp_bin_base = as.numeric( ifelse( data_ext[ exp_count == 1, !is.na( bp_med_age ) & bp_med_age  <= j ], 1, 0 ) ),
  #statin_bin_base = as.numeric( ifelse( data_ext[ exp_count == 1, !is.na( statin_age ) & statin_age  <= j ], 1, 0 ) ),
  diab_ind = data_ext[ exp_count == 1, ifelse( diab_age <= j  & !is.na(diab_age) & diab_type == 2 & !is.na(diab_type), 1, 0 ) ],
  bp_bin = data_ext[ exp_count == 1, ifelse( bp_med_age <= j & !is.na(bp_med_age), 1, 0 ) ],
  statin_ind = data_ext[ exp_count == 1, ifelse( statin_age <= j  & !is.na(statin_age), 1, 0 ) ],
  townsend_20 = data_ext[ exp_count == 1, townsend_20 ],
  ethnia = data_ext[ exp_count == 1, gen_ethnicity ],
  renal_disease = data_ext[ exp_count == 1, ifelse( renal_age <= j  & !is.na(renal_age), 1, 0 )  ],
  atrial_fibrillation = data_ext[ exp_count == 1, ifelse( Atrial_fibrillation_age <= j  & !is.na(Atrial_fibrillation_age), 1, 0 )  ],
  rheumatoid_arthritis = data_ext[ exp_count == 1, ifelse( rheumatoid_arthritis_age <= j  & !is.na(rheumatoid_arthritis_age), 1, 0 )  ],
  Severe_mental_illness_ind = data_ext[ exp_count == 1, ifelse( Severe_mental_illness_age <= j  & !is.na(Severe_mental_illness_age), 1, 0 )  ],
  Dementia_ind = data_ext[ exp_count == 1, ifelse( Dementia_age <= j  & !is.na(Dementia_age), 1, 0 )  ],
  Migraine_ind = data_ext[ exp_count == 1, ifelse( Migraine_age <= j  & !is.na(Migraine_age), 1, 0 )  ],
  depression_ind = data_ext[ exp_count == 1, ifelse( depression_age <= j  & !is.na(depression_age), 1, 0 )  ],
  SLE_ind = data_ext[ exp_count == 1, ifelse( SLE_age <= j  & !is.na(SLE_age), 1, 0 )  ],
  family_history = data_ext[ exp_count == 1, "family_history" ],
  status_cvd_death = -1000,
  status_death = -1000,
  status_cvd = -1000,
  status_statin = -1000,
  time_cvd_death = -1000,
  time_death = -1000,
  time_cvd = -1000,
  time_statin = -1000, 
  cvd_age_glob_hor = -1000,
  cvd_status_glob_hor = -1000,
  death_age_glob_hor = -1000,
  death_status_glob_hor = -1000,
  composite_age_glob_hor = -1000,
  composite_status_glob_hor = -1000,
  derivation = data_ext[ exp_count ==1, derivation ],
  lm_subcohort_index = 1)

#ref_val = mixoutsamp(fit_rirs, as.data.frame( data_ext[ exp_age <= j & derivation == "validation", c("patid","exp_age_corr","exp_age","exposure", "sbp_ind","tchol_ind", "bp_bin", "statin_bin", "scaled_corr")]))$random
#ref_der = mixoutsamp(fit_rirs, as.data.frame( data_ext[ exp_age <= j & derivation == "derivation", c("patid","exp_age_corr","exp_age","exposure", "sbp_ind","tchol_ind", "bp_bin", "statin_bin", "scaled_corr")]))$random
#ref_all is composed of RE related only to those people that have measures 
# before landmark age j
ref_all = mixoutsamp(fit_rirs, as.data.frame( data_ext[ exp_age <= j, c("patid","exp_age_corr","exp_age","exposure", "sbp_ind","tchol_ind", "bp_bin", "statin_bin", "scaled_corr")]))$random

#setting as 0 al RE that are associated to those people that
# have no measurements before landmark age j
ref_all_complete = data_ext %>%
  distinct(patid, derivation) %>%
  left_join( ref_all, by = "patid") %>%
  mutate_if( is.numeric, funs(replace_na(.,0)))

for( l in 0:10 )
{
  blup_BMI_tmp = fit_rirs$coefficients$fixed["exposurebmi"] + ref_all_complete["exposurebmi"] +
    (fit_rirs$coefficients$fixed["exposurebmi:exp_age_corr"] + ref_all_complete["exposurebmi:exp_age_corr"] )  * l 

  blup_HDL_tmp =  fit_rirs$coefficients$fixed["exposurehdl"] + ref_all_complete["exposurehdl"] +
    (fit_rirs$coefficients$fixed["exposurehdl:exp_age_corr"] + ref_all_complete["exposurehdl:exp_age_corr"] )  * l 

  blup_SBP_tmp = fit_rirs$coefficients$fixed["exposuresbp"] + fit_rirs$coefficients$fixed["bp_bin:sbp_ind"] * dati_baseline[,bp_bin] + ref_all_complete["exposuresbp"] +
    (fit_rirs$coefficients$fixed["exposuresbp:exp_age_corr"] + ref_all_complete["exposuresbp:exp_age_corr"] )  * l 

  blup_SMOKE_tmp = fit_rirs$coefficients$fixed["exposuresmokbin"] + ref_all_complete["exposuresmokbin"] +
    (fit_rirs$coefficients$fixed["exposuresmokbin:exp_age_corr"] + ref_all_complete["exposuresmokbin:exp_age_corr"] )  * l 

  blup_TCHOL_tmp = fit_rirs$coefficients$fixed["exposuretchol"] + fit_rirs$coefficients$fixed["statin_bin:tchol_ind"] * dati_baseline[,statin_ind] + ref_all_complete["exposuretchol"] +
    (fit_rirs$coefficients$fixed["exposuretchol:exp_age_corr"] + ref_all_complete["exposuretchol:exp_age_corr"] )  * l 

  blups_tmp = cbind( blup_BMI_tmp, blup_HDL_tmp, blup_SBP_tmp, blup_SMOKE_tmp, blup_TCHOL_tmp)
  head(blups_tmp)

  idx_cov = which( endsWith( colnames(covariates_matrix), paste0("_", l) ))  
  covariates_matrix[ , idx_cov ] = blups_tmp
}


for( l in 0:GLOB_horizon )
{
  time_origin = j + l
  time_horizon = j + l + 5
  
  #remove who died or who experienced cvd event before time_origin
  patid_out_time = data_ext[ (!is.na( death_date ) & death_age <= time_origin ) | ( !is.na( cvd_date ) & cvd_age  <= time_origin ) | end_age <= time_origin , patid ]
  
  data_ext_in_time = data_ext[ !patid %in% patid_out_time, ]
  
  patid_in_time = unique( data_ext_in_time$patid )
  n_patid = length( unique( data_ext_in_time$patid ) )
  
  pt_analysed[ l+1, 1 ] = n_patid
  pt_analysed[ l+1, 2 ] = data_ext_in_time[ exp_count == 1, length( which( death_age <= time_origin + 1 & ( cvd_age > time_origin + 1| is.na( cvd_age ) ) ) ) ]
  pt_analysed[ l+1, 3 ] = data_ext_in_time[ exp_count == 1, length( which( cvd_age <= time_origin + 1 & ( death_age > time_origin + 1| is.na( death_age ) ) ) )]
  pt_analysed[ l+1, 4 ] = data_ext_in_time[ exp_count == 1, length( which( cvd_age <= time_origin + 1 & death_age <= time_origin + 1 ) )]
  pt_analysed[ l+1, 5 ] = data_ext_in_time[ exp_count == 1, length( which( end_age <= time_origin + 1 & ( death_age > time_origin + 1| is.na( death_age ) ) & ( cvd_age > time_origin + 1| is.na( cvd_age ) ) ) ) ]
  pt_analysed[ l+1, 6 ] = data_ext_in_time[ exp_count == 1, length( which( death_age <= time_origin + 1 |
                                                                             cvd_age <= time_origin + 1 |
                                                                             end_age <= time_origin + 1 ) ) ]
  pt_analysed[ l+1, 7 ] = data_ext_in_time[ exp_count == 1, length( which( death_age <= time_horizon & ( cvd_age > time_horizon | is.na( cvd_age ) ) ) )]
  pt_analysed[ l+1, 8 ] = data_ext_in_time[ exp_count == 1, length( which( cvd_age <= time_horizon & ( death_age > time_horizon | is.na( death_age ) ) ) ) ]
  pt_analysed[ l+1, 9 ] = data_ext_in_time[ exp_count == 1, length( which( cvd_age <= time_horizon & death_age <= time_horizon ) ) ]
  pt_analysed[ l+1, 10 ] = data_ext_in_time[ exp_count == 1, length( which( end_age <= time_horizon & ( death_age > time_horizon| is.na( death_age ) ) & ( cvd_age > time_horizon| is.na( cvd_age ) ) ) ) ]
  pt_analysed[ l+1, 11 ] = data_ext_in_time[ exp_count == 1, length( which( death_age <= time_horizon |
                                                                              cvd_age <= time_horizon |
                                                                              end_age <= time_horizon ) ) ]
  
  
   
  #CREATING DATA for COX model
  
  #setting the outcomes of interest
  outcome_time_death  = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( is.na(x), min( time_horizon, end_age ) - time_origin, min( time_horizon, x ) - time_origin ) ), .SDcols = 'death_age',  by = patid ]
  outcome_time_cvd    = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( is.na(x), min( time_horizon, end_age ) - time_origin, min( time_horizon, x ) - time_origin ) ), .SDcols = 'cvd_age',    by = patid ]
  outcome_time_statin = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( is.na(x), min( time_horizon, end_age ) - time_origin, min( time_horizon, x ) - time_origin ) ), .SDcols = 'statin_age', by = patid ]
  
  outcome_status_death  = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( x <= time_horizon & x > time_origin & !is.na(x), 1, 0) ), .SDcols = 'death_age' ]
  outcome_status_cvd    = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( x <= time_horizon & x > time_origin & !is.na(x), 1, 0) ), .SDcols = 'cvd_age' ]
  outcome_status_statin = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( x <= time_horizon & x > time_origin & !is.na(x), 1, 0) ), .SDcols = 'statin_age' ]
  
  outcome_time_death_cvd   = apply( cbind( outcome_time_cvd$cvd_age, outcome_time_death$death_age), 1, min )
  outcome_status_death_cvd = rep( 0, length(outcome_time_cvd$cvd_age) )
  #Fatal cvd
  outcome_status_death_cvd[ outcome_time_cvd$cvd_age == outcome_time_death$death_age & outcome_status_cvd$cvd_age == 1 & outcome_status_death$death_age == 1  ] = 2
  #Non fatal CVD
  outcome_status_death_cvd[ outcome_time_cvd$cvd_age < outcome_time_death$death_age & outcome_status_cvd$cvd_age == 1  ] = 1
  #Death for other causes
  outcome_status_death_cvd[  outcome_status_death$death_age == 1 & outcome_status_cvd$cvd_age == 0  ] = 3
  
  #SETTING AS 0 THOSE PEOPLE THAT DO NOT BELONG TO THE LANDMARK SUBCOHORT
  dati_baseline[ !patid %in% patid_in_time, 'lm_subcohort_index' ]      = 0
  dati_baseline[ patid %in% patid_in_time, 'status_cvd_death' ]  = outcome_status_death_cvd
  dati_baseline[ patid %in% patid_in_time, 'status_death'  ]     = outcome_status_death
  dati_baseline[ patid %in% patid_in_time, 'status_cvd' ]       = outcome_status_cvd
  dati_baseline[ patid %in% patid_in_time, 'status_statin' ]    = outcome_status_statin
  dati_baseline[ patid %in% patid_in_time, 'time_cvd_death' ]   = outcome_time_death_cvd
  dati_baseline[ patid %in% patid_in_time, 'time_death' ]       = outcome_time_death$death_age
  dati_baseline[ patid %in% patid_in_time, 'time_cvd' ]         = outcome_time_cvd$cvd_age
  dati_baseline[ patid %in% patid_in_time, 'time_statin' ]      = outcome_time_statin$statin_age
  
  idx_smoke = which( endsWith(  colnames( covariates_matrix ), paste0("smoke_blup_", l ) ) )
  idx_hdl = which( endsWith(  colnames( covariates_matrix ), paste0("hdl_blup_", l ) ) )
  idx_sbp = which( endsWith(  colnames( covariates_matrix ), paste0("sbp_blup_", l ) ) )
  idx_tchol = which( endsWith(  colnames( covariates_matrix ), paste0("tchol_blup_", l ) ) )
  idx_bmi = which( endsWith(  colnames( covariates_matrix ), paste0("bmi_blup_", l ) ) )
  
  dati_baseline[, 'smoke_blup_FE4' ] = covariates_matrix[ , ..idx_smoke ] 
  dati_baseline[, 'hdl_blup_FE4' ]   = covariates_matrix[ , ..idx_hdl ]
  dati_baseline[, 'sbp_blup_FE4'   ] = covariates_matrix[ , ..idx_sbp ]
  dati_baseline[, 'tchol_blup_FE4' ] = covariates_matrix[ , ..idx_tchol ]
  dati_baseline[, 'bmi_blup_FE4'  ]  = covariates_matrix[ , ..idx_bmi ]
  
  

  
  cox_formula_5CVD = as.formula( paste("Surv( time_cvd, status_cvd )~", cox_cov ))
  model_cox_5CVD = coxph( cox_formula_5CVD, data = dati_baseline[derivation == "derivation" & lm_subcohort_index == 1,], x = T )
  #summary(model_cox_5CVD_new)
  
  # cox_formula_5CVD = as.formula( paste("Surv( time_cvd, status_cvd )~", cox_cov ))
  # model_cox_5CVD = coxph( cox_formula_5CVD, data = data_cox[derivation == "derivation",], x = T )
  # summary(model_cox_5CVD)
  # 
  
  save( model_cox_5CVD, file = paste0( path_to_save,"cox_model_", mod_type, j, "_", l,"_", gender, ".RData") )
  
  print( paste( l, "Cox model fitted" ) )
  
  
  #Needed for the validation:
  lp_der_5CVD    = predict( model_cox_5CVD, newdata = dati_baseline[derivation == "derivation" & lm_subcohort_index == 1, ], type = "lp")    
  lp_val_5CVD    = predict( model_cox_5CVD, newdata = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, ], type = "lp")    
  base_haz_5CVD  = basehaz( model_cox_5CVD, center = T )
  surv_pred_5CVD = exp( - matrix( base_haz_5CVD$hazard, ncol = 1 ) %*% matrix( exp( lp_val_5CVD ), nrow = 1) )
  cens_pred_5CVD = survfit( Surv( dati_baseline[derivation == "validation" & lm_subcohort_index == 1, time_cvd], 1 - dati_baseline[derivation == "validation" & lm_subcohort_index == 1, status_cvd] ) ~ 1)
  idx_time = which( base_haz_5CVD$time %in% cens_pred_5CVD$time )
  
  #Estimates for 5-year CVD risk prediction at landmark age (= on the whole landmark cohort) 
  lp_der_5CVD_all_pt_lm    = predict( model_cox_5CVD, newdata = dati_baseline[derivation == "derivation", ], type = "lp")    
  estimated_surv = exp( - matrix( base_haz_5CVD$hazard, ncol = 1 ) %*% matrix( exp( lp_der_5CVD_all_pt_lm ), nrow = 1) )
  dim(estimated_surv)
  
  
  
  
  #prediction errors CVD (j+l)
  pred_error_5CVD = pred_err_function( gender = gender, 
                                       lm_age = j,
                                       model_to_check = model_cox_5CVD,
                                       start_time = j+l,
                                       data_deriv = dati_baseline[derivation == "derivation" & lm_subcohort_index == 1,],
                                       data_valid = dati_baseline[derivation == "validation" & lm_subcohort_index == 1,],
                                       time_of_int = "time_cvd",
                                       status_of_int = "status_cvd",
                                       cox_formula = cox_formula_5CVD,
                                       cumulative_baseline = base_haz_5CVD,
                                       linear_pred_derivation = lp_der_5CVD,
                                       linear_pred_validation = lp_val_5CVD,
                                       survival_prob_mat = surv_pred_5CVD,
                                       cens_km_mat = cens_pred_5CVD,
                                       idx_time = idx_time, 
                                       max_time = 5, 
                                       type = "all_c")
  
  
  if( l == 0 )
  {
    pred_err_total_5CVD     = pred_error_5CVD
    time_valid_total_5CVD   = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, time_cvd ]
    status_valid_total_5CVD = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, status_cvd ]
    lp_val_total_5CVD       = lp_val_5CVD
    survival_total_5CVD     = surv_pred_5CVD[ dim(surv_pred_5CVD)[1],]
    cluster_patid_5CVD      = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, patid ]
    start_time_total_5CVD   = rep( j + l, length( lp_val_5CVD ) )
    
  }else{
    pred_err_total_5CVD     = rbind( pred_err_total_5CVD, pred_error_5CVD )
    time_valid_total_5CVD   = c( time_valid_total_5CVD, dati_baseline[derivation == "validation" & lm_subcohort_index == 1, time_cvd ] )
    status_valid_total_5CVD = c( status_valid_total_5CVD, dati_baseline[derivation == "validation" & lm_subcohort_index == 1, status_cvd ] )
    lp_val_total_5CVD       = c( lp_val_total_5CVD, lp_val_5CVD )
    survival_total_5CVD     = c( survival_total_5CVD, surv_pred_5CVD[ dim(surv_pred_5CVD)[1],] )
    cluster_patid_5CVD      = c( cluster_patid_5CVD,dati_baseline[derivation == "validation" & lm_subcohort_index == 1, patid ] )
    start_time_total_5CVD   = c( start_time_total_5CVD , rep( j + l, length( lp_val_5CVD ) ) )
  }
  
  
  
  if( l == 0 )
  {
    #CVD outcome
    cvd_age_glob_hor = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( is.na(x), min( j + GLOB_horizon, end_age ) - j, min( j + GLOB_horizon, x ) -j) ), .SDcols = 'cvd_age', by = patid ]
    cvd_status_glob_hor = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( x <= ( j + GLOB_horizon ) & x > j & !is.na(x), 1, 0) ), .SDcols = 'cvd_age' ]
    dati_baseline[ patid %in% patid_in_time, 'cvd_status_glob_hor'] = cvd_status_glob_hor# outcome_status_cvd_cumhaz
    dati_baseline[ patid %in% patid_in_time, 'cvd_age_glob_hor'] = cvd_age_glob_hor$cvd_age #outcome_time_cvd_cumhaz$cvd_age
    cox_formula_10CVD = as.formula( paste("Surv( cvd_age_glob_hor, cvd_status_glob_hor )~", cox_cov ))
    model_cox_10CVD = coxph( cox_formula_10CVD, data = dati_baseline[derivation == "derivation" & lm_subcohort_index == 1,], x = T )
    base_haz_10CVD = basehaz( model_cox_10CVD, centered = T )
    cum_haz_est[ 1:dim( base_haz_10CVD )[ 1 ], 1:2 ] = base_haz_10CVD
    beta_l0_cvd = model_cox_10CVD$coef 
    
    lp_der_10CVD    = predict( model_cox_10CVD, newdata = dati_baseline[derivation == "derivation" & lm_subcohort_index == 1,], type = "lp")    
    lp_val_10CVD    = predict( model_cox_10CVD, newdata = dati_baseline[derivation == "validation" & lm_subcohort_index == 1,], type = "lp")    
    surv_pred_10CVD = exp( - matrix( base_haz_10CVD$hazard, ncol = 1 ) %*% matrix( exp( lp_val_10CVD ), nrow = 1) )
    cens_pred_10CVD = survfit( Surv( dati_baseline[derivation == "validation" & lm_subcohort_index == 1, cvd_age_glob_hor ], 1 - dati_baseline[derivation == "validation" & lm_subcohort_index == 1, cvd_status_glob_hor ] ) ~ 1)
    idx_10CVD       = which( base_haz_10CVD$time %in% cens_pred_10CVD$time )
    
    
    #prediction errors CVD
    pred_error_10CVD = pred_err_function( gender                 = gender, 
                                          lm_age                 = j,
                                          model_to_check         = model_cox_10CVD,
                                          start_time             = j+l,
                                          data_deriv             = dati_baseline[derivation == "derivation" & lm_subcohort_index == 1,],
                                          data_valid             = dati_baseline[derivation == "validation" & lm_subcohort_index == 1,],
                                          time_of_int            = "cvd_age_glob_hor",
                                          status_of_int          = "cvd_status_glob_hor",
                                          cox_formula            = cox_formula_10CVD,
                                          cumulative_baseline    = base_haz_10CVD,
                                          linear_pred_derivation = lp_der_10CVD,
                                          linear_pred_validation = lp_val_10CVD,
                                          survival_prob_mat      = surv_pred_10CVD,
                                          cens_km_mat            = cens_pred_10CVD,
                                          idx_time               = idx_10CVD, 
                                          max_time               = 10, 
                                          type                   = "all_c" )
    
    
    
    
    
    data_cindex_10CVD = data.frame( time_10CVD   = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, cvd_age_glob_hor ],
                                    status_10CVD = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, cvd_status_glob_hor ],
                                    lp_val_10CVD = lp_val_10CVD, 
                                    survival_expected_10CVD = surv_pred_10CVD[ dim( surv_pred_10CVD )[ 1 ], ], 
                                    cluster_patid = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, patid ], 
                                    lm_age = j,
                                    gender = gender,
                                    type   = "10CVD")
    
    save( beta_l0_cvd, file = paste0( path_to_save,"beta_l0_cvd_outcome_", mod_type, j,"_", gender, ".RData") )
    save( pred_error_10CVD, file = paste0( path_to_save,"pred_error_10CVD_", mod_type, j,"_", gender, ".RData") )
    save( data_cindex_10CVD, file = paste0( path_to_save,"data_cindex_10CVD_", mod_type, j,"_", gender, ".RData") )
    
    
    #Death outcome
    death_age_glob_hor = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( is.na(x), min( j + GLOB_horizon, end_age ) - j, min( j + GLOB_horizon, x ) -j) ), .SDcols = 'death_age', by = patid ]
    death_status_glob_hor = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( x <= (j + GLOB_horizon ) & x > j & !is.na(x), 1, 0) ), .SDcols = 'death_age' ]
    dati_baseline[ patid %in% patid_in_time, 'death_status_glob_hor'] = death_status_glob_hor# outcome_status_death_cumhaz
    dati_baseline[ patid %in% patid_in_time, 'death_age_glob_hor'] = death_age_glob_hor$death_age #outcome_time_death_cumhaz$death_age
    cox_formula_10death = as.formula( paste("Surv( death_age_glob_hor, death_status_glob_hor )~", cox_cov ))
    model_cox_10death = coxph( cox_formula_10death, data = dati_baseline[derivation == "derivation" & lm_subcohort_index == 1,], x = T )
    base_haz_10death = basehaz( model_cox_10death, centered = T )
    cum_haz_est[ 1:dim(base_haz_10death)[1], 3:4 ] = base_haz_10death
    beta_l0_death = model_cox_10death$coef 
    
    lp_der_10death    = predict( model_cox_10death, newdata = dati_baseline[derivation == "derivation" & lm_subcohort_index == 1,], type = "lp")    
    lp_val_10death    = predict( model_cox_10death, newdata = dati_baseline[derivation == "validation" & lm_subcohort_index == 1,], type = "lp")    
    surv_pred_10death = exp( - matrix( base_haz_10death$hazard, ncol = 1 ) %*% matrix( exp( lp_val_10death ), nrow = 1) )
    cens_pred_10death = survfit( Surv( dati_baseline[derivation == "validation" & lm_subcohort_index == 1, death_age_glob_hor ], 1 - dati_baseline[derivation == "validation" & lm_subcohort_index == 1, death_status_glob_hor ] ) ~ 1)
    idx_10death       = which( base_haz_10death$time %in% cens_pred_10death$time )
    
    
    #prediction errors death
    pred_error_10death = pred_err_function( gender                 = gender, 
                                            lm_age                 = j,
                                            model_to_check         = model_cox_10death,
                                            start_time             = j+l,
                                            data_deriv             = dati_baseline[derivation == "derivation" & lm_subcohort_index == 1,],
                                            data_valid             = dati_baseline[derivation == "validation" & lm_subcohort_index == 1,],
                                            time_of_int            = "death_age_glob_hor",
                                            status_of_int          = "death_status_glob_hor",
                                            cox_formula            = cox_formula_10death,
                                            cumulative_baseline    = base_haz_10death,
                                            linear_pred_derivation = lp_der_10death,
                                            linear_pred_validation = lp_val_10death,
                                            survival_prob_mat      = surv_pred_10death,
                                            cens_km_mat            = cens_pred_10death,
                                            idx_time               = idx_10death, 
                                            max_time               = 10, 
                                            type                   = "all_c" )
    
    
    
    
    data_cindex_10death = data.frame( time_10death = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, death_age_glob_hor],
                                      status_10death = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, death_status_glob_hor],
                                      lp_val_10death = lp_val_10death, 
                                      survival_expected_10death = surv_pred_10death[ dim( surv_pred_10death )[ 1 ], ], 
                                      cluster_patid = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, patid ], 
                                      lm_age = j,
                                      gender = gender,
                                      type = "10Death")
    
    save( beta_l0_death, file = paste0( path_to_save,"beta_l0_death_outcome_", mod_type, j,"_", gender, ".RData") )
    save( pred_error_10death, file = paste0( path_to_save,"pred_error_10death_", mod_type, j, "_", gender, ".RData") )
    save( data_cindex_10death, file = paste0( path_to_save,"data_cindex_10death_", mod_type, j,"_", gender, ".RData") )
    
    
    #Composite outcome
    death_age_glob_hor = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( is.na(x), min( j + GLOB_horizon, end_age ) - j, min( j + GLOB_horizon, x ) - j ) ), .SDcols = 'death_age', by = patid ]
    death_status_glob_hor = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( x <= ( j + GLOB_horizon ) & x > j & !is.na(x), 1, 0) ), .SDcols = 'death_age' ]
    
    composite_age_glob_hor = apply( cbind( cvd_age_glob_hor$cvd_age, death_age_glob_hor$death_age), 1, min )
    #composite_status_glob_hor = ifelse( cvd_status_glob_hor == 1, 1, ifelse( death_status_glob_hor == 1, 2, 0 ) )
    composite_status_glob_hor = ifelse( cvd_status_glob_hor == 1 | death_status_glob_hor == 1, 1, 0 )
    
    
    dati_baseline[ patid %in% patid_in_time, 'composite_status_glob_hor']  = composite_status_glob_hor 
    dati_baseline[ patid %in% patid_in_time, 'composite_age_glob_hor']  =  composite_age_glob_hor  
    cox_formula_10composite = as.formula( paste("Surv( composite_age_glob_hor, composite_status_glob_hor )~", cox_cov ))
    model_cox_10composite = coxph( cox_formula_10composite, data = dati_baseline[derivation == "derivation" & lm_subcohort_index == 1,], x = T)
    base_haz_10composite = basehaz( model_cox_10composite, centered = T )
    cum_haz_est[ 1:dim(base_haz_10composite)[1], 5:6 ] = base_haz_10composite
    beta_l0_composite = model_cox_10composite$coef 
    
    lp_der_10composite    = predict( model_cox_10composite, newdata = dati_baseline[derivation == "derivation" & lm_subcohort_index == 1,], type = "lp")    
    lp_val_10composite    = predict( model_cox_10composite, newdata = dati_baseline[derivation == "validation" & lm_subcohort_index == 1,], type = "lp")    
    surv_pred_10composite = exp( - matrix( base_haz_10composite$hazard, ncol = 1 ) %*% matrix( exp( lp_val_10composite ), nrow = 1) )
    cens_pred_10composite = survfit( Surv( dati_baseline[derivation == "validation" & lm_subcohort_index == 1, composite_age_glob_hor ], 1 - dati_baseline[derivation == "validation" & lm_subcohort_index == 1, composite_status_glob_hor ] ) ~ 1)
    idx_10composite       = which( base_haz_10composite$time %in% cens_pred_10composite$time )
    
    
    #prediction errors death
    pred_error_10composite = pred_err_function( gender                 = gender, 
                                                lm_age                 = j,
                                                model_to_check         = model_cox_10composite,
                                                start_time             = j+l,
                                                data_deriv             = dati_baseline[derivation == "derivation" & lm_subcohort_index == 1,],
                                                data_valid             = dati_baseline[derivation == "validation" & lm_subcohort_index == 1,],
                                                time_of_int            = "composite_age_glob_hor",
                                                status_of_int          = "composite_status_glob_hor",
                                                cox_formula            = cox_formula_10composite,
                                                cumulative_baseline    = base_haz_10composite,
                                                linear_pred_derivation = lp_der_10composite,
                                                linear_pred_validation = lp_val_10composite,
                                                survival_prob_mat      = surv_pred_10composite,
                                                cens_km_mat            = cens_pred_10composite,
                                                idx_time               = idx_10composite, 
                                                max_time               = 10, 
                                                type                   = "all_c" )
    
    
    
    
    
    
    data_cindex_10composite = data.frame( time_10composite = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, composite_age_glob_hor ],
                                          status_10composite = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, composite_status_glob_hor ],
                                          lp_val_10composite = lp_val_10composite, 
                                          survival_expected_10composite = surv_pred_10composite[ dim( surv_pred_10composite )[ 1 ],], 
                                          cluster_patid = dati_baseline[derivation == "validation" & lm_subcohort_index == 1, patid ], 
                                          lm_age = j,
                                          gender = gender,
                                          type = "10composite")
    
    save( beta_l0_composite, file = paste0( path_to_save,"beta_l0_composite_outcome_", mod_type, j,"_", gender, ".RData") )
    save( pred_error_10composite, file = paste0( path_to_save,"pred_error_10composite_", mod_type, j,"_", gender, ".RData") )
    save( data_cindex_10composite, file = paste0( path_to_save,"data_cindex_10composite_", mod_type, j,"_", gender, ".RData") )
    
  }
  
  
  
  stopifnot( min( dati_baseline[ patid %in% patid_in_time & derivation == "derivation", time_cvd ] ) > 0 )
  
  #estimated_surv2 = exp( -matrix( base_haz_5CVD$hazard, ncol = 1 ) %*% matrix( exp( lp_der_5CVD ), nrow = 1) )
  
  # test_idx = estimated_surv[ dim(estimated_surv)[1], ] %in% estimated_surv2[ dim(estimated_surv2)[1], ] 
  # head(test_idx)
  # test_idx_2 = dati_baseline[ derivation == "derivation", patid ] %in% dati_baseline[ patid %in% patid_in_time & derivation == "derivation", patid ]
  # head(test_idx_2)
  # length(test_idx_2)
  # which( test_idx != test_idx_2)
  # length(which(estimated_surv[ dim(estimated_surv)[1], test_idx_2 ] != estimated_surv2[ dim(estimated_surv2)[1], ]))
  # 
  
  
  #SAME!!! estimated_surv = survfit( model_cox_surv, newdata = data_cox[derivation == "derivation",], id = patid )
  #dim(estimated_surv)
  
  five_year_risk = estimated_surv[ dim( estimated_surv )[1],  ] #selecting the 5-year CVD risk
  risk_prediction[ unique( data_ext$patid ) %in% dati_baseline[derivation == "derivation", patid], paste0("5.y.risk.",l) := five_year_risk  ] 
  risk_prediction[ unique( data_ext$patid ) %in% dati_baseline[derivation == "derivation", patid], paste0("5.y.risk.status.",l) := ifelse( five_year_risk <= 0.95, 1, 0 ) ]
  hitting_time = apply( estimated_surv, 2, function(x) which( x  <= 0.95 )[1]  )
  risk_prediction[ unique( data_ext$patid ) %in% dati_baseline[derivation == "derivation", patid], paste0("5.y.pred.time.",l):=  base_haz_5CVD[ hitting_time, 2 ] ]
  
  
  print(paste( l, "General CV done" ) )
  
}


colnames( cum_haz_est )= c( "hazard_cvd_outcome", "time_cvd_outcome", "hazard_death_outcome", "time_death_outcome", "hazard_composite_outcome", "time_composite_outcome" )


warnings()

data_cindex_5CVD = data.frame( time_5CVD = time_valid_total_5CVD,
                               status_5CVD = status_valid_total_5CVD,
                               lp_val_5CVD = lp_val_total_5CVD,
                               survival_expected_5CVD = survival_total_5CVD,
                               cluster_patid_5CVD = cluster_patid_5CVD,
                               lm_age = j,
                               start_time_5CVD = start_time_total_5CVD,
                               gender = gender,
                               type = "5CVD"
)

#prediction errors 
save( data_cindex_5CVD, file = paste0( path_to_save, "data_cindex_5CVD_", mod_type, j,"_", gender, ".RData") )
save( pred_err_total_5CVD, file = paste0( path_to_save, "pred_error_5CVD_", mod_type, j, "_", gender, ".RData") )

save( risk_prediction, file = paste0( path_to_save, "risk_pred_FE4_", mod_type, j,"_", gender, ".RData") )
save( pt_analysed, file = paste0( path_to_save,"pt_analysed_FE4_", mod_type, j,"_", gender, ".RData") )
save( cum_haz_est, file = paste0( path_to_save,"cum_haz_est_", mod_type, j,"_", gender, ".RData") )
save( covariates_matrix, file = paste0( path_to_save, "covariates_matrix_", mod_type, j,"_", gender, ".RData") )
