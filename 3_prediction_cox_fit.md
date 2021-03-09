**Input** :

  - j = landmark age 
  - gender = female or male
  - utilities: mixoutsample_v3.R; pred_error_check.R
  - landmark data at time j, for gender:  data_lm_j_gender_data.RData (result of script 0)
  - LMEM estimates at time j for gender: test_fit_rirs_j_m.RData (result of script 1)

**Output**:

  - pt_analysed: data frame with info on number of patients
  - cum_haz_est: data frame with estimated cumulative hazard for different outcomes (this is important for computing the INB)
  - covariates_matrix: data frame with all relevant covariates (baseline covariates + 11x5 predicted time-varying covariates). 

  a. **5-year CVD diagnosis (analysis for t star)**:
  
    - model_cox_surv: Cox model fit. Outocme of interest 5-year CVD diagnosis 
    - pred_error_5: prediction error. Outcome of interest 5-year CVD diagnosis
    - data_cindex_5CVD: data frame with all elements for computing c-indices (both landmark-specific and overall). The outcome of interest is the 5-year CVD diagnosis.
    - risk_prediction: data frame with predicted 5-year CVD risk at time j + 0, j + 1,.., j + 10
 
  b. **10-year CVD diagnosis (analysis for INB)**:
    
    - beta_l0_cvd: \betas associated to Cox model for 10-year CVD diagnosis. 
    - pred_error_cvd_10: prediction error. Outcome of interest 10-year CVD diagnosis
    - data_cindex_10cvd: data frame with all elements for computing c-indices (both landmark-specific and overall). The outcome of interest is the 10-year CVD.

  
  c. **10-year Death (analysis for INB)**:
    
    - beta_l0_death: \betas associated to Cox model for 10-year death.
    - pred_error_death_10: prediction error. Outcome of interest 10-year death
    - data_cindex_10death: data frame with all elements for computing c-indices (both landmark-specific and overall). The outcome of interest is the 10-year death.
  
  d. **10-year Composite, Death or CVD (analysis for INB)**:

    - beta_l0_composite: \betas associated to Cox model for 10-year CVD diagnosis or death.
    - pred_error_composite_10:  prediction error. Outcome of interest 10-year CVD diagnosis or death
    - data_cindex_10composite: data frame with all elements for computing c-indices (both landmark-specific and overall). The outcome of interest is the 10-year CVD diagnosis or death.
  
**Process**:
  
  - 
