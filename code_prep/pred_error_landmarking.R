pred_err_landmarking = function( data, time, status, max_time ){
  
  surv_out_valid_5CVD = with( data_cox_M[derivation == "validation",], Surv( time_cvd, status_cvd ) )
  lp_val_5CVD = predict( model_cox_M_surv, newdata = data_cox_M[derivation == "validation",], type = "lp")    
  data_temp = data_cox_M[derivation == "validation",]
  data_temp[["time_cvd"]] = 5
  data_temp$status_cvd = 0
  expect_5CVD = predict( model_cox_M_surv, newdata = data_temp, type = "expected") 
  survival_5CVD = exp( - expect_5CVD )
  rm( data_temp )
  
}