Input: 

  -  patients_outcomes.dta
  -  exposures_red_merged.RData [defined at step 0]
  -  familyhistory.dta
  -  j = landmark age under analysis
  -  g = gender under analysis 

Output:

  - creating the landmarking datasets per each gender. [10x2]

Process:

  1) Filtering on people (in line with Zhe's paper):
  
    - end date (minimum among CVD, Death and Censoring) later than j;
    - statins starting date (if it happens) later than j
    - starting date (see Zhe's paper) before j
    - CVD diagnosis (if it happens) later than j
    - gender == g
    
   2) Correcting:
   
    - CVD post censoring;
    - Death post censoring;
    - CVD post Death.
    
   3) Filtering on exposures:
   
    - select only the exposures associated with the patients selected at step 1.
    
   4) Adding useful variables for the analysis:
   
    - future: binary variable (=1 if the exposure is collected after j, =0 otherwise);
    - exp_age_corr: centralised version of the age variable (helps for the computational effort);
    - statin_bin: binary variable (=1 if the exposure is collected after the statin prescrption date, =0 otherwise);
    - bp_bin: binary variable (=1 if the exposure is collected after the blood pressure medication prescrption date, =0 otherwise);
    - sbp_ind: binary variable (=1 if the exposure is realted to SBP, =0 otherwise);
    - tchol_ind: binary variable (=1 if the exposure is realted to TCHOL, =0 otherwise);
    - age of specific diagnosis/medication: statin, end, CVD, Death, diabetes diagnosis, blood pressure medication diagnosis,
     renal disease, rheumatoid arthritis, Atrial fibrillation, Severe mental illness, Migraine, Dementia, Depression, SLE (Lupus);
    - imputing townsend 20 with its mean;
    - family_history: binary variable (=1 if "FH: premature coronary heart disease" is collected in familyhistory.dta, =0 otherwise)
    
