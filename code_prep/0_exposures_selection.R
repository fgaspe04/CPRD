#print(commandArgs(trailingOnly = TRUE))
LOCAL = FALSE

library(data.table)
library(nlme)
library(survival)
library(tidyverse)
library(haven)


if( LOCAL )
{
  outcomes = read_dta("../../real_data/patients_outcomes.dta")
  exposures = read_dta("../../real_data/cprd_2M_exposures.dta")
}else{
  load("../../../patients_outcomes.RData")
  load("../../../cprd_2M_exposures.RData")
}


#keep variables of interest
exp_oi = c("bmi","hdl", "tchol", "sbp", "smokbin")
exposures_red = exposures[ which( exposures$exposure %in% exp_oi ), ]

#remove the exposures data
rm(exposures)

#remove unrealistic data
real_sbp     = which( exposures_red$exposure == "sbp"   & exposures_red$original >= 60   & exposures_red$original <= 250 )
real_tchol   = which( exposures_red$exposure == "tchol" & exposures_red$original >= 1.75 & exposures_red$original <= 20 )
real_hdl     = which( exposures_red$exposure == "hdl"   & exposures_red$original >= 0.3  & exposures_red$original <= 3.1 )
real_bmi     = which( exposures_red$exposure == "bmi"   & exposures_red$original <= 80 )
real_smokbin = which( exposures_red$exposure == "smokbin" )

real_values   = sort( c( real_sbp, real_tchol, real_hdl, real_bmi, real_smokbin ) )
exposures_red = exposures_red[ real_values, ]

rm( real_bmi, real_sbp, real_smokbin, real_tchol, real_hdl, real_values )

#check on correctness CVD, DEATH and END times
# colnames(outcomes)
# outcomes[, c("cvd_date", "death_date", "end_date")]

#correct CVD post censoring
cvd_post_censoring = which( !is.na(outcomes$cvd_date) & outcomes$end_date < outcomes$cvd_date )
length( cvd_post_censoring )
outcomes$cvd_date[ cvd_post_censoring ] = NA
outcomes$cvd_ind[ cvd_post_censoring ] = 0

cvd_post_censoring = which( !is.na(outcomes$cvd_date) & outcomes$end_date < outcomes$cvd_date )
stopifnot( length( cvd_post_censoring ) == 0 )
rm( cvd_post_censoring )

#uniform the notation for indices of events:
length( which(is.na( outcomes$cvd_ind)))
length( which(is.na( outcomes$cvd_ind ) & is.na(outcomes$cvd_date)))

outcomes$cvd_ind[ which(is.na( outcomes$cvd_ind)) ] = 0
length( which(is.na( outcomes$death_ind)))


#correct Death post censoring
death_post_censoring = which( !is.na(outcomes$death_date) & outcomes$end_date < outcomes$death_date )
length( death_post_censoring )
outcomes$death_date[ death_post_censoring ] = NA
outcomes$death_ind[ death_post_censoring ] = 0

death_post_censoring = which( !is.na(outcomes$death_date) & outcomes$end_date < outcomes$death_date )
stopifnot( length( death_post_censoring ) == 0 )
rm( death_post_censoring )

#correct CVD post Death
cvd_post_death = which( !is.na(outcomes$cvd_date) & outcomes$death_date < outcomes$cvd_date )
stopifnot( length( cvd_post_death ) == 0 )
rm( cvd_post_death )


print("People that died with CVD (during or after diagnosis)")
length( which( outcomes$cvd_ind == 1 & outcomes$death_ind == 1 ) )/length(outcomes$cvd_ind)
print("People that died without a diagnosis of CVD")
length( which( outcomes$cvd_ind == 0 & outcomes$death_ind == 1 ) )/length(outcomes$cvd_ind)
print("People censored with a diagnosis of CVD")
length( which( outcomes$cvd_ind == 1 & outcomes$death_ind == 0 ) )/length(outcomes$cvd_ind)
print("People censored without a diagnosis of CVD")
length( which( outcomes$cvd_ind == 0 & outcomes$death_ind == 0 ) )/length(outcomes$cvd_ind)

#keep only the exposures before CVD date 
exposures_red = exposures_red[, c("patid", "exp_date", "exp_age", "exposure", "original", "scaled")]
exposures_red_merged = merge( exposures_red, outcomes[ ,c("patid", "d_yob", "cvd_date", "cvd_ind", "gender") ], by = "patid", all.x = TRUE )  
cvd_exp_cond = ( exposures_red_merged$cvd_date - exposures_red_merged$d_yob )/365.25 > exposures_red_merged$exp_age | is.na( exposures_red_merged$cvd_date )

exposures_red_merged = exposures_red_merged[ which(cvd_exp_cond), ]

print("here0")

#creating the correct scaled variable
exposures_red_merged$scaled_corr = 0
exposures_red_merged$scaled_corr[ exposures_red_merged$exposure == "smokbin"] = exposures_red_merged$original[ exposures_red_merged$exposure == "smokbin"] 
  
female_id = which( exposures_red_merged$gender == "Female" )
male_id = which( exposures_red_merged$gender == "Male" )

print("here1")

for( var in c("bmi","hdl", "tchol", "sbp") ) 
{
  #creating scaled variable
  ioi_female = which( exposures_red_merged$gender == "Female" & exposures_red_merged$exposure == var )
  ioi_male = which( exposures_red_merged$gender == "Male" & exposures_red_merged$exposure == var )
  
  female_mean = mean( exposures_red_merged[ ioi_female, "original" ] ) 
  male_mean   = mean( exposures_red_merged[ ioi_male, "original" ] ) 
  
  #overall_sd = sd( exposures_red_merged$original )
  female_sd = sd( exposures_red_merged[ ioi_female, "original" ] ) 
  male_sd   = sd( exposures_red_merged[ ioi_male, "original" ] ) 
  
  exposures_red_merged$scaled_corr[ ioi_female ] = (exposures_red_merged$original[ioi_female] - female_mean)/female_sd
  exposures_red_merged$scaled_corr[ ioi_male ]   = (exposures_red_merged$original[ioi_male] - male_mean)/male_sd
}

# table( exposures_red_merged$original[exposures_red_merged$exposure == "smokbin"] )
# table( exposures_red_merged$scaled[exposures_red_merged$exposure == "smokbin"] )
# table( exposures_red_merged$scaled_corr[exposures_red_merged$exposure == "smokbin"] )

# hist( exposures_red_merged$scaled[ which( exposures_red_merged$exposure == "bmi" & exposures_red_merged$gender == "Female" )] )
# hist( exposures_red_merged$scaled_corr[ which( exposures_red_merged$exposure == "bmi" & exposures_red_merged$gender == "Female" )] )
# 
# hist( exposures_red_merged$scaled[ which( exposures_red_merged$exposure == "bmi" & exposures_red_merged$gender == "Male" )] )
# hist( exposures_red_merged$scaled_corr[ which( exposures_red_merged$exposure == "bmi" & exposures_red_merged$gender == "Male" )] )
# 
# 
# hist( exposures_red_merged$scaled[ which( exposures_red_merged$exposure == "hdl" & exposures_red_merged$gender == "Female" )] )
# hist( exposures_red_merged$scaled_corr[ which( exposures_red_merged$exposure == "hdl" & exposures_red_merged$gender == "Female" )] )
# 
# hist( exposures_red_merged$scaled[ which( exposures_red_merged$exposure == "hdl" & exposures_red_merged$gender == "Male" )] )
# hist( exposures_red_merged$scaled_corr[ which( exposures_red_merged$exposure == "hdl" & exposures_red_merged$gender == "Male" )] )
# 

# hist( exposures_red_merged$scaled[ which( exposures_red_merged$exposure == "sbp" & exposures_red_merged$gender == "Female" )] )
# hist( exposures_red_merged$scaled_corr[ which( exposures_red_merged$exposure == "sbp" & exposures_red_merged$gender == "Female" )] )
# 
# hist( exposures_red_merged$scaled[ which( exposures_red_merged$exposure == "sbp" & exposures_red_merged$gender == "Male" )] )
# hist( exposures_red_merged$scaled_corr[ which( exposures_red_merged$exposure == "sbp" & exposures_red_merged$gender == "Male" )] )
# 
# table( exposures_red_merged$scaled[ which( exposures_red_merged$exposure == "smokbin" & exposures_red_merged$gender == "Female" )] )
# table( exposures_red_merged$scaled_corr[ which( exposures_red_merged$exposure == "smokbin" & exposures_red_merged$gender == "Female" )] )
# 
# table( exposures_red_merged$scaled[ which( exposures_red_merged$exposure == "smokbin" & exposures_red_merged$gender == "Male" )] )
# table( exposures_red_merged$scaled_corr[ which( exposures_red_merged$exposure == "smokbin" & exposures_red_merged$gender == "Male" )] )
# 
# hist( exposures_red_merged$scaled[ which( exposures_red_merged$exposure == "tchol" & exposures_red_merged$gender == "Female" )] )
# hist( exposures_red_merged$scaled_corr[ which( exposures_red_merged$exposure == "tchol" & exposures_red_merged$gender == "Female" )] )
# 
# hist( exposures_red_merged$scaled[ which( exposures_red_merged$exposure == "tchol" & exposures_red_merged$gender == "Male" )] )
# hist( exposures_red_merged$scaled_corr[ which( exposures_red_merged$exposure == "tchol" & exposures_red_merged$gender == "Male" )] )

save( exposures_red_merged, file =  "../../data_created/exposures_red_merged.RData" )
