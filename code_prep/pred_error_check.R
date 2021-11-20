pred_err_function = function( gender,
                              lm_age,
                              start_time,
                              model_to_check,
                              data_deriv,
                              data_valid,
                              time_of_int,
                              status_of_int,
                              cox_formula,
                              cumulative_baseline,
                              linear_pred_derivation,
                              linear_pred_validation,
                              survival_prob_mat,
                              cens_km_mat,
                              idx_time,
                              max_time,
                              type ){
  
  stopifnot( type %in% c("only_c", "all", "cumulative","all_c") )
  concordance = data.frame( gender = rep( gender, 12 ),
                            lm_age = rep( j, 12 ),
                            start_time = rep( j+l, 12),
                            macro = c( rep("discriminatory",7), rep("concordance",5) ),
                            time_dep = c( rep("fixed",4), rep("cumulative",3), "incidence", "incidence", "incidence", "incidence", "cumulative" ),
                            type = c( "c_index_0", "harrell_c_index", "harrell_c_index_pec", "gonen_c_index", 
                                      "song_zhou", "chambless_diao", "uno", 
                                      "BS_pec", "BS_pe", "BS_pew", "BS_ref", "IBS"),
                            value = 0,
                            sd = 0)
  
  surv_out_deriv = with( data_deriv, Surv( get( time_of_int ), get( status_of_int ) ) )
  surv_out_valid = with( data_valid, Surv( get( time_of_int ), get( status_of_int ) ) )
  
  #Discriminatory measures: time fixed
  #c-index Harrell
  #p1 = predictSurvProb( model_to_check, newdata = data_valid, times = max_time)    
  #harrelC1 = rcorr.cens( p1, surv_out_valid )
  
  #concordance$value[1] = harrelC1[[1]]
  
  #using concordance function from Survival
  survival_pred_val_last = survival_prob_mat[ dim( survival_prob_mat )[1], ]
  conc_tmp = concordance( surv_out_valid ~ survival_pred_val_last )
  concordance$value[1] = conc_tmp$concordance
  concordance$sd[1] = sqrt(conc_tmp$var)
  
  
  conc_tmp = concordance(model_to_check, newdata = data_valid )
  concordance$value[2] = conc_tmp$concordance
  concordance$sd[2] = sqrt(conc_tmp$var)
  
  print(" Sono qua 1")
  
  #lp_der = predict( model_to_check, newdata = data_deriv, type = "lp")    
  #lp_val = predict( model_to_check, newdata = data_valid, type = "lp")    
  
  if(type %in% c("all","all_c"))
  {
    Cpec = pec::cindex( model_to_check,
                        formula = cox_formula,
                        data = data_valid, 
                        eval.times = max_time) 	
  
    concordance$value[3] = Cpec$AppCindex[[1]]
    
    #Gonen M, et al. Concordance probability and discriminatory power in proportional hazards regression. Biometrika 2005;92:965-970.
    gonenheller = GHCI( linear_pred_validation )
    concordance$value[4] = gonenheller
    
  }

  # harrellC2 = rcorrcens( formula=Surv( cvd_age_glob_hor, cvd_status_glob_hor )~(-1 * p1),
  #            data = data_cox_M[derivation == "validation",])
  # hareelC2[1,1]
  
  

  
  # ## C-statistic by Uno et al. (Stat Med 2011;30:1105-1117)
  # unocindex = UnoC(surv_out_deriv,
  #                  surv_out_valid,
  #                  lp_val, 
  #                  time = 10) #shall I specify t??
  # 
  # concordance$value[4] = unocindex
  # 

  # #Incident case/dynamic control ROC/AUC by Heagerty et al. 2005
  # heagerty_zheng = risksetROC::risksetROC(  Stime        = data_valid[ ,get(time_of_int) ],
  #                                           status       = data_valid[ ,get(status_of_int) ],
  #                                           marker       = lp_val,
  #                                           predict.time = max_time,
  #                                           plot         = T)
  # 
  # concordance$value[4] = heagerty_zheng$AUC

  
  if(type %in% c("all", "cumulative"))
  {
    ## *Incident case* or *Cumulative case*/dynamic control AUC by Song and Zhou (Biometrics 2011;67:906-16)
    res.AUC.sh = AUC.sh( Surv.rsp     = surv_out_deriv,
                       Surv.rsp.new   = surv_out_valid,
                       lp             = linear_pred_derivation,
                       lpnew          = linear_pred_validation,
                       times          = seq(0,max_time,0.1) )
  
    concordance$value[5] = res.AUC.sh$iauc
  

  
    ## *Cumulative case*/dynamic control AUC by Chambless and Diao (Stat Med 2006;25:3474-3486.)
    res.AUC.cd = AUC.cd( Surv.rsp   = surv_out_deriv,
                       Surv.rsp.new = surv_out_valid,
                       lp           = linear_pred_derivation,
                       lpnew        = linear_pred_validation,
                       times        = seq(0,max_time,0.1))
  
    concordance$value[6] = res.AUC.cd$iauc
  

    # #*Cumulative case*/dynamic control AUC by Uno et al.
    # ## (http://biostats.bepress.com/cgi/viewcontent.cgi?article=1041&context=harvardbiostat)
    res.AUC.uno = AUC.uno( Surv.rsp   = surv_out_deriv,
                         Surv.rsp.new = surv_out_valid,
                         lpnew        = lp_val,
                         times        = seq(0,max_time,0.1))
    concordance$value[7] = res.AUC.uno$iauc

  }
  

  #Brier score
  pf = pec( model_to_check,
            formula = cox_formula,
            data = data_valid,
            verbose = F,
            cens.model = "marginal",
            maxtime = max_time )
  #BS
  concordance$value[8] = tail( pf$AppErr$coxph, 1)
  
  #BS with dynpred pe
  idx_cens = which( cens_km_mat$time %in% cumulative_baseline$time[ idx_time ] )
  N_val    = dim( data_valid )[ 1 ]
  
  brier_pe = dynpred::pe( time  = data_valid[ ,get(time_of_int) ],
                                      status  = data_valid[ ,get(status_of_int) ],
                                      tsurv   = cumulative_baseline$time[ idx_time ],
                                      survmat = survival_prob_mat[ idx_time, ],
                                      tcens   = cumulative_baseline$time[ idx_time ],
                                      censmat = matrix( rep( cens_km_mat$surv[ idx_cens ], N_val ), nrow = length(idx_cens)),
                                      FUN     = "Brier",
                                      tout    = max_time)
  
  concordance$value[9] = brier_pe$Err
  
  #BS with dynpred pew
  brier_pew = dynpred::pew( time    = data_valid[ ,get(time_of_int) ],
                                        status  = data_valid[ ,get(status_of_int) ],
                                        tsurv   = cumulative_baseline$time[ idx_time ],
                                        survmat = survival_prob_mat[ idx_time, ],
                                        tcens   = cumulative_baseline$time[ idx_time ],
                                        censmat = matrix( rep( cens_km_mat$surv[ idx_cens ], N_val ), nrow = length(idx_cens)),
                                        FUN     = "Brier",
                                        width   = max_time,
                                        tout = 0)
  
  concordance$value[10] = brier_pew$Err
  
  #Brier ref
  km_surv = survfit( Surv( data_valid[ ,get(time_of_int) ], data_valid[ ,get(status_of_int) ] )  ~ 1  )

  brier_ref = dynpred::pe( time    = data_valid[ ,get(time_of_int) ],
                           status  = data_valid[ ,get(status_of_int) ],
                           tsurv   = km_surv$time,
                           survmat = matrix( rep( km_surv$surv, N_val ), nrow = length(km_surv$time)),
                           tcens   = cens_km_mat$time,
                           censmat = matrix( rep( cens_km_mat$surv, N_val ), nrow = length(cens_km_mat$time)),
                           FUN     = "Brier",
                           tout    = max_time)
  
  concordance$value[11] = brier_ref$Err
  
  
  
  #IBS
  concordance$value[12] = crps( pf, times = max_time )[[2]]
  
  print("DONEEE")
  
  return( concordance )
}
