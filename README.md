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
  
<img src="https://render.githubusercontent.com/render/math?math=X_f%20%3D%20%5Cfrac%7B(X_f%20-%20%5Cbar%7BX%7D_f)%7D%7Bsd(X_f)%7D%0A">

<img src="https://render.githubusercontent.com/render/math?math=X_m%20%3D%20%5Cfrac%7B(X_m%20-%20%5Cbar%7BX%7D_m)%7D%7Bsd(X_m)%7D%0A">


  
2) Selection of the cohort
3) LME models 
