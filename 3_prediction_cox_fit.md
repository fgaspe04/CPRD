**Input** :

  - j = landmark age 
  - gender = female or male
  - utilities: mixoutsample_v3.R; pred_error_check.R
  - landmark data at time j, for gender:  data_lm_j_gender_data.RData (result of script 0)
  - LMEM estimates at time j for gender: test_fit_rirs_j_m.RData (result of script 1)

**Output**:

  - pt_analysed: data frame with info on number of patients
  - cum_haz_est: data frame with estimated cumulative hazard for different outcomes (this is important for computing the INB)
  - covariates_matrix: data frame with all relevant covariates (baseline covariates + 11x5 predicted time-varying covariates)
 

  a. **5-year CVD diagnosis (analysis for t star)**:
  
    - model_cox_surv: Cox model fit. Outcome of interest 5-year CVD diagnosis 
    - pred_error_5: prediction error. Outcome of interest 5-year CVD diagnosis
    - data_cindex_5CVD: data frame with all elements for computing c-indices (both landmark-specific and overall). The outcome of interest is the 5-year CVD diagnosis.
    - risk_prediction: data frame with predicted 5-year CVD risk at time j + 0, j + 1,.., j + 10
 
  b. **10-year CVD diagnosis (analysis for INB)**:
    
    - beta_l0_cvd: \betas associated to Cox model for 10-year CVD diagnosis. 
    - pred_error_10CVD: prediction error. Outcome of interest 10-year CVD diagnosis
    - data_cindex_10CVD: data frame with all elements for computing c-indices (both landmark-specific and overall). The outcome of interest is the 10-year CVD.

  
**Process**: 
    
  a. **Prediction of BLUPs**
    - BLUPs of BMI, TCHOL, HDL, smoke and SBP are estimated based on the covariates estimated at the L_a (dati_baseline).
    - Random Effectes of BLUPs are collected in the variable ref_all through the function mixoutsample.
    - BLUPs are estimated at each time of interest: t \in \{L_a, L_a+1,..,L_a+10\}.
  
  b.**Fitting Cox models for 5-year CVD risk**
    - Defining landmark subcohort, with all people still alive, with no CVD event up to time t,  t \in \{L_a, L_a+1,..,L_a+10\}.
    - Selecting the risk factors that we want to include in the Cox model (time-fixed risk factors and time-varying risk factors).
    - Defining the outcome: CVD diagnosis between L_a+l and L_a+l+5, l \in \{0,1,..,10\}
    - Fitting the Cox model on the landmark subcohort (model_cox_5CVD).
   
   c.**Assessing the prediction errors of the Cox models**
    - Computing c-indices and Brier score via pred_err_function (pred_error_5CVD)
   
   d.**Fitting Cox models for 10-year CVD risk**
    - Defining landmark cohort, with all people still alive, with no CVD event up to time L_a.
    - Selecting the risk factors that we want to include in the Cox model (time-fixed risk factors and time-varying risk factors).
    - Defining the outcome: CVD diagnosis between L_a and L_a+10.
    - Fitting the Cox model on the landmark subcohort (model_cox_10CVD).
   
   e.**Assessing the prediction errors of the Cox models**
    - Computing c-indices and Brier score via pred_err_function (pred_error_10CVD)
    
   f.**Alternatve outcomes**
    - Step d. and e. are repeated for two other possible outcomes of interest: 10-year CVD or Death (composite outcome) or 10-year Death.
    - cum_haz_est contains the cumulative baseline hazard estimated according to the three possible outcomes (10-year CVD, 10-year composite, 10-year Death). In the paper only the 10-year CVD analyses are reported.

    
    
    
