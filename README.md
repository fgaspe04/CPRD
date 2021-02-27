# CPRD

0) 0_exposure_selection.R

  Input: 
    - outcomes
    - exposures
 
  Output: 
    - reduced version of exposures + computed scaled time-dependent variables (exposures_red_merged.RData)
  
  Process:
  
  I select the time-dependent variables of interest: 
    - Systolic Blood Pressure (SBP);
    - Tootal Cholesterol (TCHOL);
    - High density lypoprotein (HDL);
    - Body Mass Index (BMI);
    - Smoke.
   
  I filter their sensible value (this choice is in line with Zhe's paper):
    - 60 <= SBP <= 250
    - 1.75 <= TCHOL <= 20
    - 0.3 <= HDL <= 3.1
    - BMI <= 80
    
  I correct the following errors:
    - CVD post censoring
    - Death post censoring
    - CVD post Death
   
  This step is necessary, because I want to select those measurements that happen before CVD diagnosis (if it happens).
  I select only those measurements that happen before CVD diagnosis (if it happens).
  I create the correct scaled variables folowing this formula:
  \[
  X_f = (X_f - \bar{X}_f)/sd(X_f)
  \]
  
2) Selection of the cohort
3) LME models 
