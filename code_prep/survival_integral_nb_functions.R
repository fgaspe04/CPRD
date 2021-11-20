# surv_function = function( time, S_temp, screening ){
#   surv_function_a = matrix( 0, nrow = length( screening ), ncol = ifelse( is.null( dim(S_temp) ), length(S_temp), dim(S_temp)[1] ) )
#   
#   if( dim( surv_function_a )[ 1 ] == 1 | is.null(dim(S_temp)) ){ #| length( surv_function_a ) >= 1
#     surv_function_a = approx( x = time, S_temp , xout = screening, method = "linear",  rule = 2, f = 0 )$y 
#   }else{
#     for( i in 1:length(screening) )
#       surv_function_a[i, ] = apply( S_temp, 1, function( y ) approx( x = time, y , xout = screening[i], method = "linear",  rule = 2, f = 0 )$y )
#   }
#   return( surv_function_a ) 
# }
# 
# 
# int_t_S = function( time, S_temp, omega ){
#   n_omega = length(omega)
#   S1  = surv_function( time, S_temp, screening = omega )
#   result_int = 0
#   
#   if( is.null( dim( S1 ) ) )
#   {
#     result_int = +sum((omega[-1]-omega[-n_omega])/2*(S1[-1]+S1[-n_omega]))
#   }else{
#     result_int = apply( S1, 1, function(x) +sum((omega[-1]-omega[-n_omega])/2*(S1[-1,x]+S1[-n_omega,x])) )
#   }
#   return(result_int)
# }
# 
# int_t_S_bis = function( time, S_temp, omega ){
#   n_omega = length(omega)
#   S1  = surv_function( time, S_temp, screening = omega )
#   result_int = 0
#   
#   if( is.null( dim( S1 ) ) )
#   {
#     result_int = c( (omega[-1]+omega[-n_omega])/2*(S1[-1]-S1[-n_omega]), 0 ) 
#   }else{
#     result_int = apply( S1, 1, function(x) +sum((omega[-1]+omega[-n_omega])/2*(S1[-1,x]-S1[-n_omega,x])) )
#   }
#   return( result_int )
# }


t_star = function( t_scheme, risk_predicted )
{
  delta_time_matrix = matrix( 0, nrow = dim(risk_predicted)[1], length(t_scheme) )
  
  for( j in 1:length( t_scheme ) )
  {
    delta_time_matrix[ ,j ] = risk_predicted$pred_time - t_scheme[ j ]
    
  }
  
  k_star =  risk_predicted$pred_status * apply( delta_time_matrix, 1, function( x ) which( x <= 0 )[1] ) + ( 1 -  risk_predicted$pred_status ) * length( t_scheme )
  closest_time = t_scheme[ k_star ]
  return( data.frame( patid = risk_predicted$patid, k_star = k_star,  tk = closest_time,
                      t_star = risk_predicted$pred_time, status = risk_predicted$pred_status,
                      censoring = risk_predicted$end_corr,
                      risk_class =  risk_predicted$risk_class ) )
}

# get_nb_function = function( time, S_now, a, cv, cs, bs )
# {
#   Sa = surv_function( time, S_temp = S_now, screening = a + i - 1 )
#   int_t_Sa = -int_t_S_bis( time, S_temp = S_now, omega = seq( a + i - 1, L , length.out = 20 ) )
#   -cv/a*( L - i + 1 ) + (bs-cs)* int_t_Sa 
# }
