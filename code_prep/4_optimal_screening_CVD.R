# #print(commandArgs(trailingOnly = TRUE))

args=commandArgs(trailingOnly = TRUE)
i      = as.integer(args[1])
gender = as.character(args[2]) #g can be male or female


LOCAL = FALSE

if( LOCAL == T)
{
  i = 60
  setwd("~/decision_theory/results/results_male/")
  source( "~/decision_theory/code/code_prep/survival_integral_nb_functions.R" )
  path_to_save_plot = "~/decision_theory/results/results_male/"
  
}else{
  setwd( paste0( "/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/Francesca/results/results_", gender,"/") )

  source( "/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/Francesca/code/code_prep/survival_integral_nb_functions.R" )
  
  path_to_save = paste0( "/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/Francesca/results/results_", gender,"/")
}



getwd()
library(miceadds)
library(xtable)
library(nlme)
library(survival)
library(tidyr)
library(data.table)


load( paste0( "risk_pred_FE4_RIRS_",i,"_", gender, ".RData" ) ) 
risk_pred_now = risk_prediction[ !is.na(`5.y.risk.0`) ] #keeping only the derivation set
rm( risk_prediction )
load( paste0( "covariates_matrix_RIRS_",i,"_", gender, ".RData") )
covariates_now = covariates_matrix[ derivation == "derivation", ]
rm( covariates_matrix )
load( paste0( "cum_haz_est_RIRS_",i,"_", gender, ".RData") )
cum_haz_now = cum_haz_est
rm( cum_haz_est )


load( paste0( "beta_l0_cvd_outcome_RIRS_",i,"_", gender, ".RData") )

print("ora qui 0")


#setting parameters
# count = 1
# 
# bs = 60
# cv = 20
# cs = 20
# theta = 0.8
# 
# a_vec = seq(0.25, 10, length.out = 40)




#optimal_matrix = matrix( 0, length(lm_age), 5)
#comparing_schemes_list = list(NULL)
#comparing_next_fup_list = list(NULL)
colnames(risk_pred_now)
GLOB_horizon = 10

risk_pred_now$status_cvd_visible = risk_pred_now$cvd_ind
risk_pred_now$status_cvd_visible[ risk_pred_now$cvd_age > (i + GLOB_horizon + 5) ] = 0 # 14 = 9+5 which is the wider observed time-window.


baseline_risk = 1 - risk_pred_now$`5.y.risk.0` 
hist(baseline_risk)
baseline_class = ifelse( baseline_risk <= 0.025/2, "Low risk",
                         ifelse( baseline_risk <= 0.05/2, "Med-low risk",
                                 ifelse( baseline_risk <= 0.075/2, "Med-high risk",
                                         ifelse( baseline_risk <= 0.10/2, "High risk", "Very high risk"))))
risk_pred_now$risk_class = factor(  baseline_class, levels =  c("Very high risk", "High risk", "Med-high risk", "Med-low risk", "Low risk"))
risk_pred_now$end_corr = ifelse( is.na( risk_pred_now$cvd_age ), pmin( risk_pred_now$end_age, i + GLOB_horizon + 5 ),  risk_pred_now$cvd_age )


seq_index_status = grep("y.risk.status.", colnames(risk_pred_now))
colnames(risk_pred_now)[seq_index_status]= paste0("5.y.status.", 0:10)

seq_index_status = grep("y.status.", colnames(risk_pred_now))
seq_index_risk = grep("y.risk.", colnames(risk_pred_now))

#Detecting the first status = 1, that is the first time 
# 5y CVD risk crosses the threshold
pred_time = apply( risk_pred_now[ , ..seq_index_status ], 1, function(x) which( x == 1 )[1] )
#people that crosses the theshold
length(which(!is.na(pred_time)))
#people that never crosses the threshold
length(which(is.na(pred_time)))

pred_risk = apply( risk_pred_now[ , ..seq_index_risk ], 1, function(x) x[ which( x < 0.95 )[1] ] )

stopifnot( length(which(!is.na(pred_risk))) == length(which(!is.na(pred_time))) )

risk_pred_now$pred_status = ifelse( is.na(pred_time), 0, 1)

#RISK is not monotonic!!!
pred_risk_before = apply( risk_pred_now[ !is.na( pred_time ) & pred_time != 1 , ..seq_index_risk ], 1, function(x) x[ which( x < 0.95 )[1] - 1 ] )


# #visualization of risk
# seq_to_melt = c(1,seq_index_risk)
# risk_pred_now_melt = reshape2::melt(risk_pred_now[sample(1:90000,1000, replace = F),..seq_to_melt], id = "patid")
# risk_pred_now_melt$value = 1 - risk_pred_now_melt$value
# ggplot( risk_pred_now_melt, aes(x = variable, y = as.factor(patid), fill = value))+
#   geom_tile() +
#   labs( x= "", y ="", title = paste0("Landmark ", i)) +
#   scale_fill_gradient(low="white", high="red") +
#   theme(axis.text.y = element_blank())
#   
# 
# length(pred_risk_before)
# pred_risk_before[1:10]
# names_cross =  paste0("5.y.risk.", pred_time[!is.na(pred_time)]-1)
# idx_cross = match( names_cross, colnames(risk_pred_now) )
# complete_matr =cbind(which(!is.na(pred_time)),idx_cross )
# risk_pred_now[ complete_matr[,1], complete_matr[,2], with = F]

pred_risk_B = 3*( pred_time - 1 ) + 1
pred_risk_A = 3*( pred_time - 1 ) - 2

xB = pred_time[ !is.na( pred_time ) ] - 1 + i
xA = pred_time[ !is.na( pred_time ) ] - 2 + i

yB = pred_risk[ !is.na( pred_time ) ]
#as.numeric( risk_pred_now[ cbind( which(!is.na( pred_time )), pred_risk_B[ which( !is.na(pred_risk_B) )] ) ] )
yA = rep( 1, length( xB ) )
yA[ which( xB != i ) ]  = pred_risk_before 
#as.numeric( risk_pred_now[ cbind( which(!is.na( pred_time ) & pred_risk_A > 0 ), pred_risk_A[ which( pred_risk_A > 0 ) ] ) ] )
m = (yB-yA)/(xB-xA)
q = yB - m*xB 

real_pred_time = ( 0.95 - q )/m
length(real_pred_time)
#risk_pred_delta_time_per_pt = i + pred_time[ which( !is.na( pred_time ) ) ] - 1 + as.numeric( risk_pred_now[ cbind( which( !is.na(pred_time) ) , 3*pred_time[!is.na(pred_time)] ) ] )
risk_pred_now$pred_time = pmin( risk_pred_now$end_corr, i + GLOB_horizon)
risk_pred_now$pred_time[ which( !is.na( pred_time ) ) ] = real_pred_time
pred_time[!is.na(pred_time)][1:10]
yA[1:10]
yB[1:10]
real_pred_time[1:10]


summary(risk_pred_now$`5.y.risk.1`[ which( risk_pred_now$`5.y.risk.1` != 0 ) ])
summary(risk_pred_now$`5.y.risk.2`[ which( risk_pred_now$`5.y.risk.2` != 0 ) ])
summary(risk_pred_now$`5.y.risk.3`[ which( risk_pred_now$`5.y.risk.3` != 0 ) ])
summary(risk_pred_now$`5.y.risk.4`[ which( risk_pred_now$`5.y.risk.4` != 0 ) ])
summary(risk_pred_now$`5.y.risk.5`[ which( risk_pred_now$`5.y.risk.5` != 0 ) ])
summary(risk_pred_now$`5.y.risk.6`[ which( risk_pred_now$`5.y.risk.6` != 0 ) ])
summary(risk_pred_now$`5.y.risk.7`[ which( risk_pred_now$`5.y.risk.7` != 0 ) ])
summary(risk_pred_now$`5.y.risk.8`[ which( risk_pred_now$`5.y.risk.8` != 0 ) ])
summary(risk_pred_now$`5.y.risk.9`[ which( risk_pred_now$`5.y.risk.9` != 0 ) ])
summary(risk_pred_now$`5.y.risk.10`[ which( risk_pred_now$`5.y.risk.10` != 0 ) ])



print("ora qui 1")


#optimal screening schedule 
cox_cov_less60 = c( "diab_ind + bp_bin + hdl_blup_FE4 + sbp_blup_FE4 + smoke_blup_FE4 + tchol_blup_FE4 + bmi_blup_FE4 + townsend_20  + depression_ind + Severe_mental_illness_ind + Migraine_ind")
cox_cov_more60 = c( "diab_ind + bp_bin  + hdl_blup_FE4 + sbp_blup_FE4 + smoke_blup_FE4 + tchol_blup_FE4 + bmi_blup_FE4 + townsend_20 + renal_disease + atrial_fibrillation + rheumatoid_arthritis + depression_ind + Severe_mental_illness_ind + Migraine_ind") # Dementia_ind +

#metti 0
if( i >= 60 )
{
  cox_cov = c('diab_ind', 'bp_med', 'hdl_blup_1', 'sbp_blup_1', 'smoke_blup_1',
              'tchol_blup_1', 'bmi_blup_1', 
              'townsend_20', 'renal_disease','atrial_fibrillation',
              'rheumatoid_arthritis',   
              'depression_ind',"Severe_mental_illness_ind",
              "Migraine_ind") #'Dementia_ind'
}else{
  cox_cov = c('diab_ind', 'bp_med', 'hdl_blup_1', 'sbp_blup_1', 'smoke_blup_1',
              'tchol_blup_1', 'bmi_blup_1', 
              'townsend_20',
              'depression_ind',"Severe_mental_illness_ind",
              'Migraine_ind' ) 
}

head( cum_haz_now$hazard_cvd_outcome)
head( cum_haz_now$time_cvd_outcome )
length(which(cum_haz_now$time_composite_outcome == 0))
length(which(cum_haz_now$hazard_composite_outcome == 0))

#this reduction is led by the type of
#outcome we are interested in
cum_haz_now_red = cum_haz_now[ which(cum_haz_now$hazard_cvd_outcome > 0 ), ]
dim(cum_haz_now)
dim(cum_haz_now_red)

cov_bit = exp( beta_l0_cvd %*% t( covariates_now[, ..cox_cov] ) ) # beta vector, not matrix[,1]
S_est = matrix( 0, nrow = dim(covariates_now)[1] , ncol = length( cum_haz_now_red$hazard_cvd_outcome ) + GLOB_horizon*2  )
dim(S_est)
S_est[ , 1: length( cum_haz_now_red$hazard_cvd_outcome )] = matrix( exp( -cum_haz_now_red$hazard_cvd_outcome %*% cov_bit ), nrow = dim(covariates_now)[1] , ncol = length( cum_haz_now_red$hazard_cvd_outcome ), byrow = T )


L = i + GLOB_horizon
t_scheme_1  = seq( i, L, by = 1 ) 
t_scheme_2  = seq( i, L, by = 2 ) 
t_scheme_3  = seq( i, L, by = 3 ) 
t_scheme_4  = seq( i, L, by = 4 ) 
t_scheme_5  = seq( i, L, by = 5 ) 
t_scheme_6  = seq( i, L, by = 6 ) 
t_scheme_7  = seq( i, L, by = 7 ) 
t_scheme_8  = seq( i, L, by = 8 ) 
t_scheme_9  = seq( i, L, by = 9 ) 
t_scheme_10 = seq( i, L, by = 10 ) 


t_scheme_1  = if( tail( t_scheme_1, 1 ) == L ) t_scheme_1 else c( t_scheme_1, L ) 
t_scheme_2  = if( tail( t_scheme_2, 1 ) == L ) t_scheme_2 else c( t_scheme_2, L )
t_scheme_3  = if( tail( t_scheme_3, 1 ) == L ) t_scheme_3 else c( t_scheme_3, L )
t_scheme_4  = if( tail( t_scheme_4, 1 ) == L ) t_scheme_4 else c( t_scheme_4, L )
t_scheme_5  = if( tail( t_scheme_5, 1 ) == L ) t_scheme_5 else c( t_scheme_5, L )
t_scheme_6  = if( tail( t_scheme_6, 1 ) == L ) t_scheme_6 else c( t_scheme_6, L ) 
t_scheme_7  = if( tail( t_scheme_7, 1 ) == L ) t_scheme_7 else c( t_scheme_7, L ) 
t_scheme_8  = if( tail( t_scheme_8, 1 ) == L ) t_scheme_8 else c( t_scheme_8, L ) 
t_scheme_9  = if( tail( t_scheme_9, 1 ) == L ) t_scheme_9 else c( t_scheme_9, L ) 
t_scheme_10 = if( tail( t_scheme_10, 1 ) == L ) t_scheme_10 else c( t_scheme_10, L ) 


time_scheme = list( t_scheme_1, t_scheme_2, t_scheme_3, t_scheme_4, t_scheme_5,  
                    t_scheme_6, t_scheme_7, t_scheme_8, t_scheme_9, t_scheme_10 )

tot_scheme_eval = length( time_scheme )

comparing_schemes = matrix( 0, nrow = dim(risk_pred_now)[1], ncol = 3*tot_scheme_eval  + 1 ) 

time_info = list( rep( 0, tot_scheme_eval ) )
count = length( cum_haz_now_red$hazard_cvd_outcome ) + 1
for( t in 1:tot_scheme_eval )
{
  time_info[[t]] = t_star( t_scheme = time_scheme[[ t ]], risk_predicted = risk_pred_now )
  S_est[ ,count ] = time_info[[t]]$k_star
  S_est[ ,count + 1] = time_info[[t]]$tk
  count = count + 2
}


omega = sort( unique( c( i + cum_haz_now_red$time_cvd_outcome, i:( i + GLOB_horizon ) ) ) )#- 1 +
n_omega = length(omega)
S1 = apply( S_est, 1, function( z ) approx( x = i + cum_haz_now_red$time_cvd_outcome, y = z[1:length( cum_haz_now_red$hazard_cvd_outcome )] , xout = omega[ which( ! omega %in% (i + cum_haz_now_red$time_cvd_outcome) ) ], method = "linear",  rule = 2, f = 0 )$y ) 
S2 = matrix( 0, nrow = dim(S_est)[1] , ncol = n_omega)
S2[ ,which( ! omega %in% (i + cum_haz_now_red$time_cvd_outcome) ) ] = t( S1 )
S2[ ,which( omega %in% (i  + cum_haz_now_red$time_cvd_outcome) ) ] = S_est[ ,which( (i  + cum_haz_now_red$time_cvd_outcome) %in% omega ) ]
S2 = cbind(S2, S_est[, (length( cum_haz_now_red$time_cvd_outcome ) + 1) : dim(S_est)[2]])

#int_no_stat_all
comparing_schemes[,1] = apply( S2, 1, function(x) +sum((omega[-1] - omega[-n_omega])/2*(x[2:n_omega] + x[1:( n_omega - 1 ) ] ) ) )
summary(comparing_schemes[,1])

theta = 0.8


for( t in 1:tot_scheme_eval  )
{
  never_statins_idx = which( S2[, n_omega + 2*t ] == L )
  statins_idx = which( S2[, n_omega + 2*t ] != L )
  
  min_t_L9 = S2[, n_omega + 2*t ]
  idx_min_t_L9 = match( min_t_L9, omega )
  S_cutoff = S2[cbind(1:dim(S2)[1],idx_min_t_L9) ]
  
  int_NS_total = function( x ){
    omega_now =  omega[ omega >=  i & omega <= x[ n_omega + 2*t ] ]
    int_fin = +sum( ( omega_now[-1] - omega_now[-length(omega_now)] )/2*( x[2:length(omega_now)] + x[1:( length(omega_now) - 1 ) ] ))
    
    return(int_fin)
  }
  
  comparing_schemes[ , 3*(t - 1) + 2  ] = apply( S2, 1,  int_NS_total )
  
  int_S_total = function( x ){
    index_now = which( omega >= x[ n_omega + 2*t ] & omega <= L )
    omega_now =  omega[ index_now ]#
    
    int_fin = +sum( ( omega_now[-1] - omega_now[-length(omega_now)] )/2*( (x[index_now[2]: index_now[length(index_now)]]) + (x[ index_now[1] : index_now[(length(index_now)-1)]]) ))
    
    return(int_fin)
  }
  
  #S2_theta_final = S2_theta
  #S2_theta_final[ ,1:n_omega] = apply(S2_theta[ ,1:n_omega], 2, function(x) x/C_y )
  S2_theta_temp = S2[ statins_idx, ]
  S2_theta_temp[ ,1:n_omega] = (S2[ statins_idx, 1:n_omega ])^(theta) * (S_cutoff[ statins_idx ])^(1-theta)
  
  comparing_schemes[ never_statins_idx, 3*(t - 1) + 3  ] = 0
  comparing_schemes[ statins_idx, 3*(t - 1) + 3  ] = apply( S2_theta_temp , 1,  int_S_total )
  
  
  #Seventh element
  max_N_visit = GLOB_horizon/t +1 
  comparing_schemes[ ,3*(t - 1) + 4 ] =  time_info[[t]]$k_star #* 10/t
  comparing_schemes[ which( time_info[[t]]$tk == L ), 3*(t - 1) + 4 ] =  max_N_visit #* 10/t
  
  #saving Figure 4
  if( i == 60 & gender == "male" )
  {
    
    pdf(file = paste0( path_to_save, "example_survival_13_time_",t, ".pdf"))
    plot(omega[1:idx_min_t_L9[statins_idx[13]]], S2[statins_idx[13],1:idx_min_t_L9[statins_idx[13]]], type = "l", xlab = "Time", ylab = "Probability of not having CVD", xlim = c(60,70), ylim = c(0.8,1), lwd = 2 )
    lines(omega[idx_min_t_L9[statins_idx[13]] : n_omega],(S2[ statins_idx[13], idx_min_t_L9[statins_idx[13]]:n_omega ])^theta * (S_cutoff[ statins_idx[13] ])^(1-theta), col = 3, lwd = 2 )
    lines(omega[idx_min_t_L9[statins_idx[13]] : n_omega],(S2[ statins_idx[13], idx_min_t_L9[statins_idx[13]]:n_omega ]), col = 1, lwd = 2, lty = 2 )
    abline(v = time_info[[t]]$t_star[statins_idx[13]], col = 2, lwd = 2, lty = 2)
    abline(v = time_info[[t]]$tk[statins_idx[13]], col = 2, lwd = 2)
    dev.off()
    
    pdf(file = paste0( path_to_save, "example_survival_20_time_",t, ".pdf"))
    plot(omega[1:idx_min_t_L9[statins_idx[20]]], S2[statins_idx[20],1:idx_min_t_L9[statins_idx[20]]], type = "l", xlab = "Time", ylab = "Probability of not having CVD", xlim = c(60,70), ylim = c(0.8,1), lwd = 2 )
    lines(omega[idx_min_t_L9[statins_idx[20]] : n_omega],(S2[ statins_idx[20], idx_min_t_L9[statins_idx[20]]:n_omega ])^theta * (S_cutoff[ statins_idx[20] ])^(1-theta), col = 3, lwd = 2 )
    lines(omega[idx_min_t_L9[statins_idx[20]] : n_omega],(S2[ statins_idx[20], idx_min_t_L9[statins_idx[20]]:n_omega ]), col = 1, lwd = 2, lty = 2 )
    abline(v = time_info[[t]]$t_star[statins_idx[20]], col = 2, lwd = 2, lty = 2)
    abline(v = time_info[[t]]$tk[statins_idx[20]], col = 2, lwd = 2)
    dev.off()
    
    pdf(file = paste0( path_to_save, "example_survival_55_time_",t, ".pdf"))
    plot(omega[1:idx_min_t_L9[statins_idx[55]]], S2[statins_idx[55],1:idx_min_t_L9[statins_idx[55]]], type = "l", xlab = "Time", ylab = "Probability of not having CVD", xlim = c(60,70), ylim = c(0.8,1), lwd = 2 )
    lines(omega[idx_min_t_L9[statins_idx[55]] : n_omega],(S2[ statins_idx[55], idx_min_t_L9[statins_idx[55]]:n_omega ])^theta * (S_cutoff[ statins_idx[55] ])^(1-theta), col = 3, lwd = 2 )
    lines(omega[idx_min_t_L9[statins_idx[55]] : n_omega],(S2[ statins_idx[55], idx_min_t_L9[statins_idx[55]]:n_omega ]), col = 1, lwd = 2, lty = 2 )
    abline(v = time_info[[t]]$t_star[statins_idx[55]], col = 2, lwd = 2, lty = 2)
    abline(v = time_info[[t]]$tk[statins_idx[55]], col = 2, lwd = 2)
    dev.off()
    
  }
  
  print(t)
}




colnames( comparing_schemes ) = c( "int_no_stat_all",
                                   "scheme_1_NS", "scheme_1_S", "scheme_1_N",
                                   "scheme_2_NS", "scheme_2_S", "scheme_2_N",
                                   "scheme_3_NS", "scheme_3_S", "scheme_3_N",
                                   "scheme_4_NS", "scheme_4_S", "scheme_4_N",
                                   "scheme_5_NS", "scheme_5_S", "scheme_5_N",
                                   "scheme_6_NS", "scheme_6_S", "scheme_6_N",
                                   "scheme_7_NS", "scheme_7_S", "scheme_7_N",
                                   "scheme_8_NS", "scheme_8_S", "scheme_8_N",
                                   "scheme_9_NS", "scheme_9_S", "scheme_9_N",
                                   "scheme_10_NS", "scheme_10_S", "scheme_10_N")


comparing_schemes = as.data.frame(comparing_schemes)
comparing_schemes$patid = risk_pred_now$patid
comparing_schemes$risk_class = risk_pred_now$risk_class
comparing_schemes$pred_time = risk_pred_now$pred_time
comparing_schemes$pred_status = risk_pred_now$pred_status
comparing_schemes$statin_age = risk_pred_now$statin_age
comparing_schemes$cvd_age = risk_pred_now$cvd_age
comparing_schemes$cvd_ind = risk_pred_now$cvd_ind
comparing_schemes$death_age = risk_pred_now$death_age
comparing_schemes$death_ind = risk_pred_now$death_ind
comparing_schemes$end_age = risk_pred_now$end_age



#comparing_schemes = cbind( patid = risk_pred_now$patid, risk_class = risk_pred_now$risk_class, 
#                           pred_time = risk_pred_now$pred_time, pred_status = risk_pred_now$pred_status,
#                           comparing_schemes )


save( comparing_schemes, file = paste0( path_to_save,"comparing_schemes_RIRS_apply_", i, "_", gender, ".RData") )
