**Input** :

  - j = landmark age 
  - gender = female or male
  - utilities: t_star function
  - predicted 5-year CVD risk profile for each person in the landmark cohort: risk_pred_FE4_RIRS_j_gender.RData
  - covariates: time-fixed and time-varying covariates covariates_matrix_RIRS_j_gender.RData
  - cumulative baseline hazard: cum_haz_est_RIRS_j_gender.RData
  - estimated \betas from the 10-year CVD Cx model: beta_l0_cvd_outcome_RIRS_j_gender.RData
  
**Output**:

  - comparing_schemes: evaluation of all integrals and E[N] of the INB formula, for each risk-assessment strategy evaluated.
  
**Process**: 
    
  a. **Definition of the baseline risk categories**
  
    - Definition of the baseline risk categories: classification based on the 5-year CVD risk estimated at the landmark age L_a
    - baseline_risk <= 0.025/2  => "Low risk"
    - 0.025/2 < baseline_risk <= 0.05/2  => "Med-low risk"
    - 0.05/2 < baseline_risk <= 0.075/2 => "Med-high risk"
    - 0.075/2 < baseline_risk <= 0.10/2 => "High risk"
    - baseline_risk > 0.10/2 => "Very high risk"
    
  b.**Estimate of t***
  
    - Identification of the first time after which the 5-year CVD risk exceeds the 5% threshold: pred_risk_A
    - Identification of the first time after which the 5-year CVD risk exceeds the 5% threshold: pred_risk_B 
    - Identification of t*, real_pred_time, by linear interpolation between: pred_risk_A and pred_risk_B
   
  c.**Estimate of S^{NS}**
  
    - S^{NS} is S_est = exp(-cum_haz_est*exp(\beta*covariates))  [cumulative hazards and betas are related to Cox model with 10-year CVD as outcome]
    
  d.**Evaluation of the elements of the INB formula**
  
    - Evaluated risk-assessment strategies: time_scheme (list of all considered strategies)
    - Identification of the first visit after the 5-year CVD risk exceeds the 5% threshold. Function: t_star
    - Estimate of the reference case, comparing_schemes[,1]: \int_{L_a}^{L_a+10} S^{NS}(t) dt 
    - For loop for estimating \int_{L_a}^{\tau_k*} S^{NS}(t) dt;  \int_{\tau_k*}^{L_a+10} S^{S}(t) dt; E[N], for each risk-assessment strategy under evaluation. 
    
