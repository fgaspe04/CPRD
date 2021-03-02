Input: 

  - the landmarking dataset for a specific gender at a specific age (lm_data_j_g.RData).
  
Output:

  - output of a linear mixed effect model.
  
Process:

  - Fitting a Linear Mixed Effect model with Random Intercept and Random Slopes (RIRS). We consider 5 outcomes of interest: SBP, TCHOL, BMI, HDL and smoke.  The fitted model is: 
  
  <img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Bequation%7D%0A%5Clabel%7Beq%3Almem%7D%0A%5Cbegin%7Baligned%7D%5Bb%5D%0A%20%20%20%20smoke_%7Bij%7D%20%26%3D%20%5Cbeta_%7B10%7D%20%20%2B%20%5Cbeta_%7B11%7D%20age_%7Bij%7D%20%2B%20u_%7B10i%7D%20%2B%20u_%7B11i%7Dage_%7Bij%7D%20%2B%20%5Cvarepsilon_%7Bij%7D%20%20%5C%5C%20%0A%20%20%20%20HDL_%7Bij%7D%20%20%20%20%20%26%3D%20%5Cbeta_%7B20%7D%20%20%2B%20%5Cbeta_%7B21%7D%20age_%7Bij%7D%20%2B%20u_%7B20i%7D%20%2Bu_%7B21i%7Dage_%7Bij%7D%20%2B%20%5Cvarepsilon_%7Bij%7D%20%20%20%5C%5C%0A%20%20%20%20SBP_%7Bij%7D%20%20%20%20%20%26%3D%20%5Cbeta_%7B30%7D%20%20%2B%20%5Cbeta_%7B31%7D%20age_%7Bij%7D%20%2B%20%5Cbeta_%7B32%7D%20BPM_%7Bij%7D%20%2B%20u_%7B30i%7D%20%2B%20u_%7B31i%7Dage_%7Bij%7D%20%2B%20%5Cvarepsilon_%7Bij%7D%5C%5C%0A%20%20%20%20TCHOL_%7Bij%7D%20%20%20%26%3D%20%5Cbeta_%7B40%7D%20%20%2B%20%5Cbeta_%7B41%7D%20age_%7Bij%7D%20%2B%20%20%5Cbeta_%7B42%7D%20statin_%7Bij%7D%20%2B%20u_%7B40i%7D%20%2B%20u_%7B41i%7Dage_%7Bij%7D%20%2B%20%5Cvarepsilon_%7Bij%7D%20%20%5C%5C%0A%20%20%20%20BMI_%7Bij%7D%20%26%3D%20%5Cbeta_%7B50%7D%20%20%2B%20%5Cbeta_%7B51%7D%20age_%7Bij%7D%20%2B%20u_%7B50i%7D%20%2B%20u_%7B51i%7Dage_%7Bij%7D%20%2B%20%5Cvarepsilon_%7Bij%7D%20%0A%5Cend%7Baligned%7D%0A%5Cend%7Bequation%7D%0A">
    
    
    It is important to remind that:
      
      - Smoke, BMI, TCHOL, HDL and SBP are scaled per gender (not per landmark!)
      - The age is centered wrt the landmark age.
      - We include Random Intercepts: u_{k0}, k = 1:5.
      - We include Random Slopes: u_{k1}, k = 1:5. 
      - i is the subject index.
      - j is the exposure index.
      - the code runs per gender, per landmark age, ONLY on the derivation set (the validation set will be used later).
