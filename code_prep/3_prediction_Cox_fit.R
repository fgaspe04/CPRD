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

source("../code_prep/mixoutsamp_v3.R")
source("../code_prep/pred_error_check.R")



if( LOCAL ){
  j = 60
  load.Rdata( paste0( "../../data_created/male/data_lm_",j,"_male.RData"), paste0( "data_ext") )
  path_to_save = "../../results/results_male/"
  mod_type = "RIRS_"
  load.Rdata( paste0( "../../results/results_male/test_fit_rirs",j,"_m.RData"), paste0( "fit_rirs") )
}else{
  load.Rdata( paste0( "../../data_created/", gender, "/data_lm_",j, "_", gender,".RData"), paste0( "data") )
  mod_type = "RIRS_"
  if(gender == "male"){
    load.Rdata( paste0( "../../results/results_", gender, "/fit_rirs",j,"_male.RData"), paste0( "fit_rirs") )
  }else{
    load.Rdata( paste0( "../../results/results_", gender, "/fit_rirs",j,"_female.RData"), paste0( "fit_rirs") )
  }
  path_to_save = paste0("../../results/results_", gender )
  
  
  # load.Rdata( paste0( "../../data_created/male/data_lm_",j,"_m_data.RData"), paste0( "data_ext") )
  # path_to_save = "/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/Francesca/results/results_male/"
  # mod_type = "RIRS_"
  # load.Rdata( paste0( "/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/Francesca/results/results_male/test_fit_rirs",j,"_m.RData"), paste0( "fit_rirs") )
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

pt_analysed = matrix( 0, nrow = GLOB_horizon + 1, ncol = 11 )


#computing fixed part of prediction
length( unique(data_ext$patid[ which( data_ext$exp_age_corr <= 0 ) ] ))
#SOME PEOPLE HAVE OBSERVATION ONLY IN THE FUTURE!!!! (25% at 60 yo)
length( unique(data_ext$patid[ which( data_ext$exp_age_corr <= 0 ) ] ))/ length(unique(data_ext$patid))

dati_baseline = data.table( 
  patid = unique( data_ext$patid ),
  bp_bin_base = as.numeric( ifelse( data_ext[ exp_count == 1, !is.na( bp_med_age ) & bp_med_age  <= j ], 1, 0 ) ),
  statin_bin_base = as.numeric( ifelse( data_ext[ exp_count == 1, !is.na( statin_age ) & statin_age  <= j ], 1, 0 ) ) )



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
  
  
  #PREDICTION step
  adding_data_ext = data.table( patid = rep( unique(data_ext_in_time$patid), each = N_outcomes ),
                              exp_age_corr = rep( l, N_outcomes*n_patid),
                              exp_age = rep( j+l, N_outcomes*n_patid ),
                              exposure = rep( outcomes_names,  n_patid ),
                              sbp_ind  = rep( as.numeric( outcomes_names %in% 'sbp' ), n_patid ),
                              tchol_ind  = rep( as.numeric( outcomes_names %in% 'tchol' ), n_patid ),
                              bp_bin = rep( dati_baseline$bp_bin_base[ dati_baseline$patid %in% patid_in_time ], each = N_outcomes ),
                              statin_bin = rep( dati_baseline$statin_bin_base[ dati_baseline$patid %in% patid_in_time ], each = N_outcomes ),
                              scaled_corr = rep( NA, N_outcomes*n_patid ) )
  
  data_for_prediction_complete = rbind( data_ext_in_time[ exp_age <= time_origin, c("patid","exp_age_corr","exp_age","exposure", "sbp_ind","tchol_ind", "bp_bin", "statin_bin", "scaled_corr")], #past data
                                        adding_data_ext ) #data that we are interested in
  
  setorder( data_for_prediction_complete, patid, exp_age_corr )
  
  glob_res    = mixoutsamp( fit_rirs, as.data.frame( data_for_prediction_complete ) )
  pred_res_oi = glob_res$preddata[is.na(glob_res$preddata$scaled_corr),]
  head(pred_res_oi)
  
  #in order to produce a easy computation
  pred_res_oi$random[is.na( pred_res_oi$random )] = 0
  pred_res_oi$fitted_complete = pred_res_oi$fixed + pred_res_oi$random 
  head( pred_res_oi$fitted_complete )
  head( pred_res_oi$fitted )
  head( pred_res_oi$fixed )
  
  
  stopifnot( length( which( is.na( pred_res_oi$fitted_complete ) ) ) == 0 )
  
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
  
  data_cox = data.table( patid  = patid_in_time,
                           derivation = data_ext_in_time[ exp_count == 1, "derivation" ],
                           outcome_status_death_cvd,
                           outcome_status_death,
                           outcome_status_cvd,
                           outcome_status_statin,
                           outcome_time_death_cvd,
                           outcome_time_death = outcome_time_death$death_age,
                           outcome_time_cvd = outcome_time_cvd$cvd_age,
                           outcome_time_statin = outcome_time_statin$statin_age,
                           diab_ind = data_ext_in_time[ exp_count == 1, ifelse( diab_age <= j  & !is.na(diab_age) & diab_type == 2 & !is.na(diab_type), 1, 0 ) ],
                           bp_bin = data_ext_in_time[ exp_count == 1, ifelse( bp_med_age <= j & !is.na(bp_med_age), 1, 0 ) ],
                           statin_ind = data_ext_in_time[ exp_count == 1, ifelse( statin_age <= j  & !is.na(statin_age), 1, 0 ) ],
                           townsend_20 = data_ext_in_time[ exp_count == 1, "townsend_20" ],
                           ethnia = data_ext_in_time[ exp_count == 1, "gen_ethnicity" ],
                           renal_disease = data_ext_in_time[ exp_count == 1, ifelse( renal_age <= j  & !is.na(renal_age), 1, 0 )  ],
                           atrial_fibrillation = data_ext_in_time[ exp_count == 1, ifelse( Atrial_fibrillation_age <= j  & !is.na(Atrial_fibrillation_age), 1, 0 )  ],
                           rheumatoid_arthritis = data_ext_in_time[ exp_count == 1, ifelse( rheumatoid_arthritis_age <= j  & !is.na(rheumatoid_arthritis_age), 1, 0 )  ],
                           Severe_mental_illness_ind = data_ext_in_time[ exp_count == 1, ifelse( Severe_mental_illness_age <= j  & !is.na(Severe_mental_illness_age), 1, 0 )  ],
                           Dementia_ind = data_ext_in_time[ exp_count == 1, ifelse( Dementia_age <= j  & !is.na(Dementia_age), 1, 0 )  ],
                           Migraine_ind = data_ext_in_time[ exp_count == 1, ifelse( Migraine_age <= j  & !is.na(Migraine_age), 1, 0 )  ],
                           depression_ind = data_ext_in_time[ exp_count == 1, ifelse( depression_age <= j  & !is.na(depression_age), 1, 0 )  ],
                           SLE_ind = data_ext_in_time[ exp_count == 1, ifelse( SLE_age <= j  & !is.na(SLE_age), 1, 0 )  ],
                           family_history = data_ext_in_time[ exp_count == 1, "family_history" ],
                           landmark_age = rep( j, length(unique(data_ext_in_time$patid))),
                           smoke_blup_FE4 = pred_res_oi$fitted_complete[ pred_res_oi$exposure == "smokbin"],
                           hdl_blup_FE4 = pred_res_oi$fitted_complete[ pred_res_oi$exposure == "hdl"],
                           sbp_blup_FE4 = pred_res_oi$fitted_complete[ pred_res_oi$exposure == "sbp"],
                           tchol_blup_FE4 = pred_res_oi$fitted_complete[ pred_res_oi$exposure == "tchol"],
                           bmi_blup_FE4 = pred_res_oi$fitted_complete[ pred_res_oi$exposure == "bmi"]
                           
  )
  
  
  colnames( data_cox ) = c('patid', 'derivation', 'status_cvd_death', 'status_death', 'status_cvd', 'status_statin', 'time_cvd_death', 'time_death','time_cvd',
                             'time_statin', 'diab_ind', 'bp_bin',  'statin_ind',
                             'townsend_20', 'ethnia','renal_disease','atrial_fibrillation',
                             'rheumatoid_arthritis', "Severe_mental_illness_ind", 
                             "Dementia_ind", "Migraine_ind", "depression_ind", "SLE_ind",
                             'family_history',
                             'landmark_age', 'smoke_blup_FE4', 
                             'hdl_blup_FE4', 'sbp_blup_FE4', 'tchol_blup_FE4', 'bmi_blup_FE4')
  
  
  covariates_matrix[ covariates_matrix$patid %in% data_cox$patid, (l*N_outcomes+17):((l+1)*N_outcomes+16)] = data_cox[ ,c( 'bmi_blup_FE4', 'hdl_blup_FE4', 'sbp_blup_FE4', 'smoke_blup_FE4', 'tchol_blup_FE4' ) ]
  
  
  cox_formula = as.formula( paste("Surv( time_cvd, status_cvd )~", cox_cov ))
  model_cox = cph( cox_formula, data = data_cox[derivation == "derivation",], x = T, y = T, surv = T )
  model_cox_surv = coxph( cox_formula, data = data_cox[derivation == "derivation",], x = T )
  
  save( model_cox_surv, file = paste0( path_to_save,"cox_model_", mod_type, j, "_", l,"_", gender, ".RData") )
  
  print( paste( l, "Cox model fitted" ) )
  
  #prediction errors CVD (j+l)
  pred_error_5 = pred_err_function( lm_age = j,
                                    model_to_check = model_cox,
                                    data_deriv = data_cox[derivation == "derivation",],
                                    data_valid = data_cox[derivation == "validation",],
                                    time_of_int = "time_cvd",
                                    status_of_int = "status_cvd",
                                    cox_formula = cox_formula,
                                    max_time = 5, 
                                    type = "all")
  save( pred_error_5, file = paste0( path_to_save,"pred_error_5_", mod_type, j, "_", l, "_", gender, ".RData") )
  
  
  lp_val_5CVD = predict( model_cox_surv, newdata = data_cox[derivation == "validation",], type = "lp")    
  data_temp = data_cox[derivation == "validation",]
  data_temp[["time_cvd"]] = 5
  data_temp$status_cvd = 0
  expect_5CVD = predict( model_cox_surv, newdata = data_temp, type = "expected") 
  survival_5CVD = exp( - expect_5CVD )
  rm( data_temp )
  
  
  if( l == 0 )
  {
    time_valid_total_5CVD   = data_cox[derivation == "validation", time_cvd]
    status_valid_total_5CVD = data_cox[derivation == "validation", status_cvd]
    lp_val_total_5CVD = lp_val_5CVD
    survival_total_5CVD = survival_5CVD
    cluster_patid_5CVD = data_cox[derivation == "validation",patid]
    
  }else{
    time_valid_total_5CVD   = c( time_valid_total_5CVD, data_cox[derivation == "validation", time_cvd] )
    status_valid_total_5CVD = c( status_valid_total_5CVD, data_cox[derivation == "validation", status_cvd])
    lp_val_total_5CVD = c( lp_val_total_5CVD, lp_val_5CVD )
    survival_total_5CVD = c( survival_total_5CVD, survival_5CVD )
    cluster_patid_5CVD = c(cluster_patid_5CVD,data_cox[derivation == "validation",patid])
  }
  
   
  
  if( l == 0 )
  {
    #CVD outcome
    cvd_age_glob_hor = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( is.na(x), min( j + GLOB_horizon, end_age ) - j, min( j + GLOB_horizon, x ) -j) ), .SDcols = 'cvd_age', by = patid ]
    cvd_status_glob_hor = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( x <= ( j + GLOB_horizon ) & x > j & !is.na(x), 1, 0) ), .SDcols = 'cvd_age' ]
    data_cox$cvd_status_glob_hor = cvd_status_glob_hor# outcome_status_cvd_cumhaz
    data_cox$cvd_age_glob_hor = cvd_age_glob_hor$cvd_age #outcome_time_cvd_cumhaz$cvd_age
    cox_formula_cvd = as.formula( paste("Surv( cvd_age_glob_hor, cvd_status_glob_hor )~", cox_cov ))
    model_cox_cvd = cph( cox_formula_cvd, data = data_cox[derivation == "derivation",], x = T, y = T, surv = T )
    model_cox_cvd_surv = coxph( cox_formula_cvd, data = data_cox[derivation == "derivation",], x = T )
    bh_temp_cvd = basehaz( model_cox_cvd, centered = T )
    cum_haz_est[ 1:dim(bh_temp_cvd)[1], 1:2 ] = bh_temp_cvd
    beta_l0_cvd = model_cox_cvd$coef 
    save( beta_l0_cvd, file = paste0( path_to_save,"beta_l0_cvd_outcome_", mod_type, j,"_", gender, ".RData") )
    
    #prediction errors CVD
    pred_error_cvd_10 = pred_err_function( lm_age = j,
                                           model_to_check = model_cox_cvd,
                                           data_deriv = data_cox[derivation == "derivation",],
                                           data_valid = data_cox[derivation == "validation",],
                                           time_of_int = "cvd_age_glob_hor",
                                           status_of_int = "cvd_status_glob_hor",
                                           cox_formula = cox_formula_cvd,
                                           max_time = 10,
                                           type = "all")
    save( pred_error_cvd_10, file = paste0( path_to_save,"pred_error_cvd_10_", mod_type, j,"_", gender, ".RData") )
    
    lp_val_10CVD = predict( model_cox_cvd_surv, newdata = data_cox[derivation == "validation",], type = "lp")    
    data_temp = data_cox[derivation == "validation",]
    data_temp[["cvd_age_glob_hor"]] = 10
    data_temp$cvd_status_glob_hor = 0
    expect_10CVD = predict( model_cox_cvd_surv, newdata = data_temp, type = "expected") 
    survival_10CVD = exp( - expect_10CVD )
    rm( data_temp )
    
    data_cindex_10CVD = data.frame( time_10CVD   = data_cox[derivation == "validation", cvd_age_glob_hor ],
                                    status_10CVD = data_cox[derivation == "validation", cvd_status_glob_hor ],
                                    lp_val_10CVD = lp_val_10CVD, 
                                    survival_expected_10CVD = survival_10CVD, 
                                    cluster_patid = data_cox[derivation == "validation", patid ], 
                                    lm_age = j,
                                    type = "10CVD")
    save( data_cindex_10CVD, file = paste0( path_to_save,"data_cindex_10CVD_", mod_type, j,"_", gender, ".RData") )
    
    
    #Death outcome
    death_age_glob_hor = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( is.na(x), min( j + GLOB_horizon, end_age ) - j, min( j + GLOB_horizon, x ) -j) ), .SDcols = 'death_age', by = patid ]
    death_status_glob_hor = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( x <= (j + GLOB_horizon ) & x > j & !is.na(x), 1, 0) ), .SDcols = 'death_age' ]
    data_cox$death_status_glob_hor = death_status_glob_hor# outcome_status_death_cumhaz
    data_cox$death_age_glob_hor = death_age_glob_hor$death_age #outcome_time_death_cumhaz$death_age
    cox_formula_death = as.formula( paste("Surv( death_age_glob_hor, death_status_glob_hor )~", cox_cov ))
    model_cox_death = cph( cox_formula_death, data = data_cox[derivation == "derivation",], x = T, y = T, surv = T )
    model_cox_death_surv = coxph( cox_formula_death, data = data_cox[derivation == "derivation",], x = T )
    bh_temp_death = basehaz( model_cox_death, centered = T )
    cum_haz_est[ 1:dim(bh_temp_death)[1], 3:4 ] = bh_temp_death
    beta_l0_death = model_cox_death$coef 
    save( beta_l0_death, file = paste0( path_to_save,"beta_l0_death_outcome_", mod_type, j,"_", gender, ".RData") )
    
    #prediction errors death
    pred_error_death_10 = pred_err_function( lm_age = j,
                                             model_to_check = model_cox_death,
                                             data_deriv = data_cox[derivation == "derivation",],
                                             data_valid = data_cox[derivation == "validation",],
                                             time_of_int = "death_age_glob_hor",
                                             status_of_int = "death_status_glob_hor",
                                             cox_formula = cox_formula_death,
                                             max_time = 10, 
                                             type = "all" )
    save( pred_error_death_10, file = paste0( path_to_save,"pred_error_death_10_", mod_type, j, "_", gender, ".RData") )
    
    
    lp_val_10death = predict( model_cox_death, newdata = data_cox[derivation == "validation",], type = "lp")    
    data_temp = data_cox[derivation == "validation",]
    data_temp[["death_age_glob_hor"]] = 10
    data_temp$death_status_glob_hor = 0
    expect_10death = predict( model_cox_death_surv, newdata = data_temp, type = "expected") 
    survival_10death = exp( - expect_10death )
    rm( data_temp )
    
    data_cindex_10death = data.frame( time_10death = data_cox[derivation == "validation", death_age_glob_hor],
                                      status_10death = data_cox[derivation == "validation", death_status_glob_hor],
                                      lp_val_10death = lp_val_10death, 
                                      survival_expected_10death = survival_10death, 
                                      cluster_patid = data_cox[derivation == "validation", patid ], 
                                      lm_age = j,
                                      type = "10Death")
    save( data_cindex_10death, file = paste0( path_to_save,"data_cindex_10death_", mod_type, j,"_", gender, ".RData") )
    
    
    #Composite outcome
    death_age_glob_hor = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( is.na(x), min( j + GLOB_horizon, end_age ) - j, min( j + GLOB_horizon, x ) - j ) ), .SDcols = 'death_age', by = patid ]
    death_status_glob_hor = data_ext_in_time[ exp_count == 1, lapply( .SD, function( x ) ifelse( x <= ( j + GLOB_horizon ) & x > j & !is.na(x), 1, 0) ), .SDcols = 'death_age' ]
    
    composite_age_glob_hor = apply( cbind( cvd_age_glob_hor$cvd_age, death_age_glob_hor$death_age), 1, min )
    #composite_status_glob_hor = ifelse( cvd_status_glob_hor == 1, 1, ifelse( death_status_glob_hor == 1, 2, 0 ) )
    composite_status_glob_hor = ifelse( cvd_status_glob_hor == 1 | death_status_glob_hor == 1, 1, 0 )
    
    
    data_cox$composite_status_glob_hor  = composite_status_glob_hor 
    data_cox$composite_age_glob_hor  =  composite_age_glob_hor  
    cox_formula_composite = as.formula( paste("Surv( composite_age_glob_hor, composite_status_glob_hor )~", cox_cov ))
    model_cox_composite = cph( cox_formula_composite, data = data_cox[derivation == "derivation",], x = T, y = T, surv = T )
    model_cox_composite_surv = coxph( cox_formula_composite, data = data_cox[derivation == "derivation",], x = T)
    bh_temp_composite = basehaz( model_cox_composite, centered = T )
    cum_haz_est[ 1:dim(bh_temp_composite)[1], 5:6 ] = bh_temp_composite
    beta_l0_composite = model_cox_composite$coef 
    save( beta_l0_composite, file = paste0( path_to_save,"beta_l0_composite_outcome_", mod_type, j,"_", gender, ".RData") )
    
    #prediction errors death
    pred_error_composite_10 = pred_err_function( lm_age = j,
                                                 model_to_check = model_cox_composite,
                                                 data_deriv = data_cox[derivation == "derivation",],
                                                 data_valid = data_cox[derivation == "validation",],
                                                 time_of_int = "composite_age_glob_hor",
                                                 status_of_int = "composite_status_glob_hor",
                                                 cox_formula = cox_formula_composite,
                                                 max_time = 10, 
                                                 type = "all" )
    save( pred_error_composite_10, file = paste0( path_to_save,"pred_error_composite_10_", mod_type, j,"_", gender, ".RData") )
    
    lp_val_10composite = predict( model_cox_composite, newdata = data_cox[derivation == "validation",], type = "lp")    
    data_temp = data_cox[derivation == "validation",]
    data_temp[["composite_age_glob_hor"]] = 10
    data_temp$composite_status_glob_hor = 0
    expect_10composite = predict( model_cox_composite_surv, newdata = data_temp, type = "expected") 
    survival_10composite = exp( - expect_10composite )
    rm( data_temp )
    
    
    data_cindex_10composite = data.frame( time_10composite = data_cox[derivation == "validation", composite_age_glob_hor ],
                                          status_10composite = data_cox[derivation == "validation", composite_status_glob_hor ],
                                          lp_val_10composite = lp_val_10composite, 
                                          survival_expected_10composite = survival_10composite, 
                                          cluster_patid = data_cox[derivation == "validation", patid ], 
                                          lm_age = j,
                                          type = "10composite")
    save( data_cindex_10composite, file = paste0( path_to_save,"data_cindex_10composite_", mod_type, j,"_", gender, ".RData") )
    
  }
  
  
  
  stopifnot(min(data_cox[derivation=="derivation",time_cvd])>0)  
  estimated_surv = survfit( model_cox_surv, newdata = data_cox[derivation == "derivation",], id = patid )
  dim(estimated_surv$surv)
  bh = basehaz( model_cox_surv, centered = T )
  # mean_cov = colMeans(data_cox[, .SD, .SDcols = c('diab_ind', 'bp_bin',  'hdl_blup_FE4',   'sbp_blup_FE4',  'smoke', 'tchol_blup_FE4', 'bmi_blup_FE4', 'townsend_20', 'renal_disease','atrial_fibrillation','rheumatoid_arthritis', 'family_history' ) ], na.rm = T)
  # lin_pred_pred = t( model_cox_F$coef) %*% apply( data_cox[, .SD, .SDcols = c('diab_ind', 'bp_bin',  'hdl_blup_FE4',   'sbp_blup_FE4',  'smoke', 'tchol_blup_FE4', 'bmi_blup_FE4', 'townsend_20', 'renal_disease','atrial_fibrillation','rheumatoid_arthritis', 'family_history' ) ],1, function(x) x - mean_cov  )
  # risk_over_time = exp( as.matrix( -bh[,1] ) %*% exp( lin_pred_pred ) )
  risk_over_time = estimated_surv$surv
  five_year_risk = risk_over_time[ dim( estimated_surv$surv )[1],  ]
  risk_prediction[ unique( data_ext$patid ) %in% data_cox[derivation == "derivation", patid], paste0("5.y.risk.",l) := five_year_risk  ] 
  risk_prediction[ unique( data_ext$patid ) %in% data_cox[derivation == "derivation", patid], paste0("5.y.risk.status.",l) := ifelse( five_year_risk <= 0.95, 1, 0 ) ]
  hitting_time = apply( risk_over_time, 2, function(x) which( x  <= 0.95 )[1]  )
  risk_prediction[ unique( data_ext$patid ) %in% data_cox[derivation == "derivation", patid], paste0("5.y.pred.time.",l):=  bh[ hitting_time, 2 ] ]
  
  
  print(paste( l, "General CV done" ) )
  
}

colnames( covariates_matrix ) = c( "patid", "diab_ind", "bp_med", "statin_bin",
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

colnames( cum_haz_est )= c( "hazard_cvd_outcome", "time_cvd_outcome", "hazard_death_outcome", "time_death_outcome", "hazard_composite_outcome", "time_composite_outcome" )


warnings()

data_cindex_5CVD = data.frame( time_5CVD = time_valid_total_5CVD,
                               status_5CVD = status_valid_total_5CVD,
                               lp_val_5CVD = lp_val_total_5CVD,
                               survival_expected_5CVD = survival_total_5CVD,
                               cluster_patid_5CVD = cluster_patid_5CVD,
                               lm_age = j,
                               type = "5CVD"
)

save( data_cindex_5CVD, file = paste0( path_to_save, "data_cindex_5CVD_", mod_type, j,"_", gender, ".RData") )
save( risk_prediction, file = paste0( path_to_save, "risk_pred_FE4_", mod_type, j,"_", gender, ".RData") )
save( pt_analysed, file = paste0( path_to_save,"pt_analysed_FE4_", mod_type, j,"_", gender, ".RData") )
save( cum_haz_est, file = paste0( path_to_save,"cum_haz_est_", mod_type, j,"_", gender, ".RData") )
save( covariates_matrix, file = paste0( path_to_save, "covariates_matrix_", mod_type, j,"_", gender, ".RData") )
