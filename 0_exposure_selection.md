# 0_exposure_selection.R

  Input: 
    
    - outcomes
    - exposures
 
  Output: 
  
    - reduced version of exposures + computed scaled time-dependent variables (exposures_red_merged.RData)
  
  Process:
  
  1) Filtering the time-dependent variables of interest: 
  
    - Systolic Blood Pressure (SBP);
    - Tootal Cholesterol (TCHOL);
    - High density lypoprotein (HDL);
    - Body Mass Index (BMI);
    - Smoke.
   
  2) Filtering their sensible value (this choice is in line with Zhe's paper):
  
    - 60 <= SBP <= 250
    - 1.75 <= TCHOL <= 20
    - 0.3 <= HDL <= 3.1
    - BMI <= 80
    
  3) Correcting the following errors:
  
    - CVD post censoring
    - Death post censoring
    - CVD post Death
   
  This step is necessary, because I want to select those measurements that happen before CVD diagnosis (if it happens).
  
  I select only those measurements that happen before CVD diagnosis (if it happens).
  
  I create the correct scaled variables following this formula:
  
<img src="https://render.githubusercontent.com/render/math?math=X_f%20%3D%20%5Cfrac%7B(X_f%20-%20%5Cbar%7BX%7D_f)%7D%7Bsd(X_f)%7D%0A">

<img src="https://render.githubusercontent.com/render/math?math=X_m%20%3D%20%5Cfrac%7B(X_m%20-%20%5Cbar%7BX%7D_m)%7D%7Bsd(X_m)%7D%0A">

 REMARK: The scaled variables are computed separately for the derivation and the validation sets.
