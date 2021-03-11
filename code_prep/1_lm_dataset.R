#creating derivation landmark datasets 

# #print(commandArgs(trailingOnly = TRUE))

args=commandArgs(trailingOnly = TRUE)
j = as.integer(args[1])
gender = as.character(args[2])


library(data.table)
library(nlme)
library(survival)
library(tidyverse)
library(haven)

LOCAL = FALSE


if( LOCAL )
{
  j = 60
  outcomes = read_dta("../../real_data/patients_outcomes.dta")
  load("../../data_created/exposures_red_merged.RData")
  family_data = read_dta("../../real_data/familyhistory.dta")
  
  path_to_save = paste0( "../../data_created/", gender,"/")
  
}else{
  load("../../../patients_outcomes.RData")
  load("../../data_created/exposures_red_merged.RData")
  family_data = read_dta("../../familyhistory.dta")
  
  path_to_save = paste0( "../../data_created/", gender,"/")
}



#selecting correct patients 
cat(j,"\n")
#j = 80

#Descriptives
length(unique(exposures_red_merged$patid))
length(unique(exposures_red_merged$patid[exposures_red_merged$gender == "Female"]))
length(unique(exposures_red_merged$patid[exposures_red_merged$gender == "Male"]))
length(unique(outcomes$pracid))
length(unique(outcomes$region))
length(unique(outcomes$pracid[outcomes$derivation == 1]))
length(unique(outcomes$pracid[outcomes$derivation == 0]))
length(unique(outcomes$region[outcomes$derivation == 0]))
length(which(outcomes$derivation == 1))
length(which(outcomes$derivation == 0))
dist_fup = as.numeric((outcomes$end_date - outcomes$start_date)/365.25)
summary(dist_fup)

#Filter on patients


#Focus on end, start, cvd, statins
end_age_cond     = ( outcomes$end_date - outcomes$d_yob )/365.25 > j
statin_age_cond  = ( outcomes$statins_prscd - outcomes$d_yob )/365.25 > j | is.na( outcomes$statins_prscd )
start_age_cond   = ( outcomes$start_date - outcomes$d_yob )/365.25 < j #why not <= ?
cvd_age_cond     = ( outcomes$cvd_date - outcomes$d_yob )/365.25 > j | is.na( outcomes$cvd_date )

if(gender == "male"){
  gender_cond      = outcomes$gender == 1 | outcomes$gender == "Male"
}else{
  gender_cond      = outcomes$gender == 0 | outcomes$gender == "Female"
}

selected_id = end_age_cond & statin_age_cond & start_age_cond & cvd_age_cond & gender_cond
selected_pt = outcomes$patid[ selected_id ]


#check on correctness CVD, DEATH and END times
colnames(outcomes)
outcomes[, c("cvd_date", "death_date", "end_date")]

#correct CVD post censoring
cvd_post_censoring = which( !is.na(outcomes$cvd_date) & outcomes$end_date < outcomes$cvd_date )
cvd_during_censoring = which( !is.na(outcomes$cvd_date) & outcomes$end_date == outcomes$cvd_date & is.na( outcomes$death_date )   )

length( which( !is.na(outcomes$cvd_date) & outcomes$end_date == outcomes$cvd_date ) )

length( cvd_post_censoring )
outcomes$cvd_date[ cvd_post_censoring ] = NA
outcomes$cvd_ind[ cvd_post_censoring ] = 0

outcomes$end_date[ cvd_during_censoring ] = outcomes$end_date[ cvd_during_censoring ] + 1
outcomes$cvd_ind[ cvd_during_censoring ] = 1 #it was already 1

cvd_post_censoring = which( !is.na(outcomes$cvd_date) & outcomes$end_date < outcomes$cvd_date )
stopifnot( length( cvd_post_censoring ) == 0 )
rm( cvd_post_censoring )

#uniform the notation for indices of events:
length( which(is.na( outcomes$cvd_ind)))
length( which(is.na( outcomes$cvd_ind ) & is.na(outcomes$cvd_date)))

outcomes$cvd_ind[ which(is.na( outcomes$cvd_ind)) ] = 0


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


outcomes_selected = outcomes[ selected_id, ]

#Keep only the exposures of those patients that satisfy all the requirements
exposures_selected = exposures_red_merged[ exposures_red_merged$patid %in% selected_pt, ]

rm( exposures_red_merged, outcomes, selected_id, selected_pt, end_age_cond,
    cvd_age_cond, gender_cond, start_age_cond,
    statin_age_cond, derivation )

#adding index for expsoures in the future (after L_a)
exposures_selected$future = ifelse( exposures_selected$exp_age >=j, 1, 0 )
#adding age corrected (centralised)
exposures_selected$exp_age_corr = exposures_selected$exp_age - j

#Select columns to keep in outcomes dataset

outcome_columns_to_keep = c( "patid", "pracid", "region", "bp_med_prscd", "statins_prscd", "statinuser",
                             "townsend_20", "gen_ethnicity", 
                             "diab_type", "diab_type_uncertainty", "diab_ind" , "diab_date", 
                             "rheumatoid_arthritis_ind", "rheumatoid_arthritis_date",
                             "Atrial_fibrillation_ind", "Atrial_fibrillation_date",
                             "renal_ind", "renal_date",
                             "depression_ind", "depression_date",
                             "SLE_ind", "SLE_date",                  
                             "Migraine_ind", "Migraine_date",
                             "Severe_mental_illness_ind", "Severe_mental_illness_date",
                             "Dementia_ind", "Dementia_date", 
                             "death_date", "end_date", "death_ind", "derivation" ) 



#Creating the final data: we can have people with NO exposures
data_ext = merge( exposures_selected, outcomes_selected[ ,outcome_columns_to_keep], by = "patid", all = F )

dim( data_ext )

#statin_index
data_ext$statin_bin = ifelse( data_ext$statins_prscd <= data_ext$exp_date & !is.na(data_ext$statins_prscd), 1, 0 ) 

#bp index
data_ext$bp_bin = ifelse( data_ext$bp_med_prscd <= data_ext$exp_date & !is.na(data_ext$bp_med_prscd), 1, 0 ) 

#sbp measure ind
data_ext$sbp_ind = ifelse( data_ext$exposure == "sbp", 1, 0 ) 

#tchol measure ind
data_ext$tchol_ind = ifelse( data_ext$exposure == "tchol", 1, 0 ) 


print("almost creating data_table")
data_ext = data.table( data_ext )
print("created data_table")

setkey( data_ext, patid )

data_ext[ ,exp_count := 1:.N, by = patid ]
print("created exp_count")

#data_ext[ death_ind == 1, ]
#data_ext[ death_ind == 1 & death_date != end_date, lapply(.SD, function(x) length(unique(x))), .SDcols = "patid" ]
#2902

#warning decison!!!!!!!
#data_ext[ death_ind == 1 & death_date != end_date, death_ind := 0 ] #censoring death
table( data_ext$SLE_ind )

stopifnot( length( which( data_ext$SLE_ind == 1 ) ) ==  length( which( !is.na( data_ext$SLE_date ) ) ) )
stopifnot( length( which( data_ext$Migraine_ind == 1 ) ) ==  length( which( !is.na( data_ext$Migraine_date ) ) ) )
stopifnot( length( which( data_ext$depression_ind == 1 ) ) ==  length( which( !is.na( data_ext$depression_date ) ) ) )
stopifnot( length( which( data_ext$Dementia_ind == 1 ) ) ==  length( which( !is.na( data_ext$Dementia_date ) ) ) )
stopifnot( length( which( data_ext$Severe_mental_illness_ind == 1 ) ) ==  length( which( !is.na( data_ext$Severe_mental_illness_date ) ) ) )
stopifnot( length( which( data_ext$diab_ind == 1 ) ) ==  length( which( !is.na( data_ext$diab_date ) ) ) )
stopifnot( length( which( data_ext$Atrial_fibrillation_ind == 1 ) ) ==  length( which( !is.na( data_ext$Atrial_fibrillation_date ) ) ) )
stopifnot( length( which( data_ext$rheumatoid_arthritis_ind == 1 ) ) ==  length( which( !is.na( data_ext$rheumatoid_arthritis_date ) ) ) )
stopifnot( length( which( data_ext$renal_ind == 1 ) ) ==  length( which( !is.na( data_ext$renal_date ) ) ) )

stopifnot( length( which( data_ext$death_ind == 1 ) ) ==  length( which( !is.na( data_ext$death_date ) ) ) )
stopifnot( length( which( data_ext$cvd_ind == 1 ) ) ==  length( which( !is.na( data_ext$cvd_date ) ) ) )


#creating age variables
data_ext[ ,statin_age := ifelse( is.na( statins_prscd ), NA, ( statins_prscd - d_yob )/365.25 ) ]
data_ext[ ,end_age := ifelse( is.na( end_date ), NA, ( end_date - d_yob )/365.25 ) ]
data_ext[ ,cvd_age := ifelse( is.na( cvd_date ), NA, ( cvd_date - d_yob )/365.25 ) ]
data_ext[ ,death_age := ifelse( is.na( death_date ), NA, ( death_date - d_yob )/365.25 ) ]
data_ext[ ,diab_age := ifelse( is.na( diab_date ), NA, ( diab_date - d_yob )/365.25 ) ]
data_ext[ ,bp_med_age := ifelse( is.na( bp_med_prscd ), NA, ( bp_med_prscd - d_yob )/365.25 ) ]
data_ext[ ,renal_age := ifelse( is.na( renal_date ), NA, ( renal_date - d_yob )/365.25 ) ]
data_ext[ ,rheumatoid_arthritis_age := ifelse( is.na( rheumatoid_arthritis_date ), NA, ( rheumatoid_arthritis_date - d_yob )/365.25 ) ]
data_ext[ ,Atrial_fibrillation_age := ifelse( is.na( Atrial_fibrillation_date ), NA, ( Atrial_fibrillation_date - d_yob )/365.25 ) ]
data_ext[ ,Severe_mental_illness_age := ifelse( is.na( Severe_mental_illness_date ), NA, ( Severe_mental_illness_date - d_yob )/365.25 ) ]
data_ext[ ,Migraine_age := ifelse( is.na( Migraine_date ), NA, ( Migraine_date - d_yob )/365.25 ) ]
data_ext[ ,Dementia_age := ifelse( is.na( Dementia_date ), NA, ( Dementia_date - d_yob )/365.25 ) ]
data_ext[ ,depression_age := ifelse( is.na( depression_date ), NA, ( depression_date - d_yob )/365.25 ) ]
data_ext[ ,SLE_age := ifelse( is.na( SLE_date ), NA, ( SLE_date - d_yob )/365.25 ) ]


#data_ext$gen_ethnicity[ which( data_ext$gen_ethnicity == "" )] = "Unknown"
data_ext$gen_ethnicity = as.factor( data_ext$gen_ethnicity )

#imputing !!!!!!!!!!!!!!!!!!!!!!!!
data_ext$townsend_20[ which( is.na( data_ext$townsend_20 ) ) ] = mean( data_ext$townsend_20, na.rm = T )

#family hisotry 
data_ext[ ,family_history := 0 ]
CHD_patid = family_data$patid[grep( "FH: premature coronary heart disease", family_data$desc ) ]
data_ext[ patid %in% CHD_patid, family_history := 1 ]

stopifnot( length( which( data_ext$exp_age > data_ext$end_age ) ) == 0 )

save( data_ext, file =  paste0( path_to_save, "data_lm_",j, "_", gender,".RData" ) )

