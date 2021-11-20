#Tables

library(scales)
library(RColorBrewer)
library(tidyr)
library(ggpubr)
library(stringr)
library(miceadds)
library(ggplot2)
library(reshape2)
library(xtable)
library(gridExtra)
library(gtable)
library(grid)
#library(plyr)
library(dplyr)
library(gg.gap)


setwd("~/decision_theory/results/")
path_to_save_plot = "~/decision_theory/img/"



######## Loading the scheme comparison
#i = 45 

lm_age = seq( 40, 80, 5 )
n_age = length( lm_age )
stacked_comparing_scheme = NULL


for( j in 1:2 )
{
  GENDER = ifelse( j == 1, "male", "female" )
  gender = ifelse( j == 1, "m", "f" )
  
  for(i in lm_age )
  {
    load.Rdata( paste0( "results_", GENDER, "/comparing_schemes_RIRS_apply_",i,"_", GENDER, ".RData"), paste0("comparing_schemes") ) #gender
    comparing_schemes$sex = GENDER
    comparing_schemes$landmark = i
    stacked_comparing_scheme = rbind(stacked_comparing_scheme, comparing_schemes)
  }
}

GLOB_horizon = 10

colnames( stacked_comparing_scheme )
table( table( stacked_comparing_scheme$patid ) )
length( unique( stacked_comparing_scheme$patid ) )



#############################################
#Computing INB function
bns = 25000 #c( 100, 500, 1000, 2000, 3000, 4000, 5000, 10000)
us = 0.997
cs = 150 #20
cv = 18.39

source("../code/code_plot/function_INB.R")

INB_result     = matrix( 0, nrow = dim(stacked_comparing_scheme), ncol = 10)
benefit_result = matrix( 0, nrow = dim(stacked_comparing_scheme), ncol = 10)
cost_result    = matrix( 0, nrow = dim(stacked_comparing_scheme), ncol = 10)
ref_result     = bns*stacked_comparing_scheme$int_no_stat_all

for( f in 1:10 )
{
  idx_NS_integral = grep(paste0(f,"_NS"), colnames(stacked_comparing_scheme))
  NS_integral = stacked_comparing_scheme[ ,idx_NS_integral]
  
  idx_S_integral = grep( paste0(f,"_S"), colnames(stacked_comparing_scheme))
  S_integral = stacked_comparing_scheme[ ,idx_S_integral]
  
  ref_integral = stacked_comparing_scheme$int_no_stat_all
  
  idx_N_visit = grep( paste0( f,"_N" ), colnames(stacked_comparing_scheme))[2]
  N_visit = stacked_comparing_scheme[ ,idx_N_visit]
  
  full_res = function_INB(bns, us, cs, cv, NS_integral, S_integral, ref_integral, N_visit)
  
  INB_result[ ,f ]     = full_res$uf
  benefit_result[ ,f ] = full_res$benefit
  cost_result[ ,f ]    = full_res$cost
  
}

INB_result = data.frame( INB_result )
colnames( INB_result ) = paste0("uf_",1:10)
benefit_result = data.frame( benefit_result )
colnames( benefit_result ) = paste0("benefit_",1:10)
cost_result = data.frame( cost_result )
colnames( cost_result ) = paste0("cost_",1:10)



#Defininingthe optimal strategy
best_of_new = max.col( INB_result, ties.method = "last")
table( best_of_new )

#Adding the results to stacked_comparing_scheme
stacked_comparing_scheme$best_of = best_of_new
stacked_comparing_scheme = cbind( stacked_comparing_scheme,
                                  INB_result,
                                  benefit_result,
                                  cost_result, 
                                  ref_result)

table( stacked_comparing_scheme$best_of)

best_of_now = table( stacked_comparing_scheme$best_of, stacked_comparing_scheme$risk_class )
best_of_now = sweep( best_of_now, 2, colSums( best_of_now ),`/` )
best_of_now = data.frame( best_of_now )
colnames(best_of_now) = c( "Frequency", "Risk_class", "Percentage")


tab_new = table( best_of_new )
#print(tab_new)

if( length( table( best_of_new ) ) != 10 ){
  names_table = names( tab_new )
  tab_corr = rep(0, 10)
  tab_corr[ as.character( 1:10 ) %in% names_table  ] = tab_new
}else{
  tab_corr = tab_new
}

#best_of = cbind( best_of, tab_corr ) 

stacked_comparing_scheme_with_factors = stacked_comparing_scheme 
stacked_comparing_scheme_with_factors$landmark = factor(stacked_comparing_scheme_with_factors$landmark, levels = lm_age )
stacked_comparing_scheme_with_factors$best_of = factor(stacked_comparing_scheme_with_factors$best_of, levels = 1:10 )

table(stacked_comparing_scheme_with_factors$best_of)

#Table 1 and Table 2 of the supplementary material
for( j in 1:2 )
{
  GENDER = ifelse( j == 1, "male", "female" )
  gender = ifelse( j == 1, "men", "women" )
  title_gender =  paste(toupper(substring(GENDER, 1,1)), substring(GENDER, 2), sep="", collapse=" ")
  
 
  data_for_table_sup = stacked_comparing_scheme_with_factors %>%
    filter(sex == GENDER) %>%
    group_by(landmark, risk_class, best_of, .drop = F) %>%
    tally() %>%
    ungroup() %>%
    group_by(landmark) %>%
    mutate( freq = n / sum(n), unique_label = paste(risk_class, best_of)) %>%
    filter( (risk_class == "Very high risk" & best_of == 10) | risk_class != "Very high risk" ) 
  
  data_for_table_sup$unique_label[ data_for_table_sup$risk_class == "Very high risk" & data_for_table_sup$best_of == 10 ] = "Very high risk -"
  
  lev_for_all = c(  "Very high risk -", 
                    "Very high risk Total",
                    "Very high risk Never_cross",
                    paste0( "High risk ", 1:10 ),
                    "High risk Total",
                    "High risk Never_cross",
                    paste0( "Med-high risk ", 1:10 ),
                    "Med-high risk Total",
                    "Med-high risk Never_cross",
                    paste0( "Med-low risk ", 1:10 ),
                    "Med-low risk Total",
                    "Med-low risk Never_cross",
                    paste0( "Low risk ", 1:10 ),
                    "Low risk Total",
                    "Low risk Never_cross",
                    "Total",
                    "Never_cross"
                )
  
  data_for_table_sup$unique_label = factor( data_for_table_sup$unique_label, 
                                            levels = lev_for_all )
  
  
  data_for_table_sup = data_for_table_sup %>%
    mutate(info = paste0(n, " (", round(freq*100,2), "%)" )) %>% 
    select(unique_label, landmark, info) %>%
    group_by(unique_label, .drop = F) %>%
    spread(key = landmark, value = info)

  
  data_never_cross_all_risk_cat = stacked_comparing_scheme %>%
    filter( sex == GENDER ) %>%
    group_by(landmark, risk_class, .drop = F) %>%
    summarize( count_no_cross = sum(1-pred_status),
               Total = n(),
               ratio_no_cross = round( count_no_cross/Total *100,2) ,  
               Never_cross = paste0(count_no_cross, " (", ratio_no_cross,"%)")) %>%
    select(landmark, risk_class, Total, Never_cross ) %>%
    gather(Total, Never_cross, -c(landmark, risk_class)) %>%
    mutate(unique_label = factor(paste(risk_class, Total),lev_for_all )) %>%
    select(landmark, unique_label, Never_cross) %>%
    spread( key = landmark, value = Never_cross )
  
  data_never_cross_total = stacked_comparing_scheme %>%
    filter( sex == GENDER ) %>%
    group_by(landmark, .drop = F) %>%
    summarize( count_no_cross = sum(1-pred_status),
               Total = n(),
               ratio_no_cross = round( count_no_cross/Total *100,2) ,  
               Never_cross = paste0(count_no_cross, " (", ratio_no_cross,"%)")) %>%
    select(landmark, Total, Never_cross ) %>%
    gather(Total, Never_cross, -c(landmark)) %>%
    mutate(unique_label = factor( Total,lev_for_all )) %>%
    select(landmark, unique_label, Never_cross) %>%
    spread( key = landmark, value = Never_cross )
    
  
  #Merging the dataset
  data_table_final = data_never_cross_all_risk_cat %>%
    bind_rows(data_for_table_sup) %>%
    bind_rows(data_never_cross_total) %>%
    arrange(unique_label) %>%
    filter(!unique_label %in% c( "Very high risk Total", "Very high risk Never_cross")) %>%
    mutate(risk_optim = as.character(unique_label)) %>%
    separate(risk_optim, c("risk_class", "best_of"),"risk") %>%
    select(-unique_label) %>%
    mutate(risk_class = replace( risk_class, best_of %in% c(paste0(" ", 2:10), " Total", " Never_cross"), "")) %>%
    mutate_at(vars(as.character(seq(40,80, by = 5))), ~ifelse(. %in% c("0 (0%)", "0 (NaN%)" ) , "-", .)) %>%
    mutate_at(vars(risk_class, best_of), ~ifelse(. %in% c( "Never_cross", " Never_cross") , "Never cross", .)) %>%
    mutate(best_of = replace(best_of, risk_class == "Total", "Total")) %>%
    mutate(best_of = replace(best_of, risk_class == "Never cross", "Never cross")) %>%
    mutate(risk_class = replace(risk_class, risk_class %in% c("Total", "Never cross"), "")) %>%
    select(risk_class, best_of, as.character(seq(40,80, by = 5)))
  
  
  
  print( xtable(data_table_final), include.rownames = F)
  
}


#Fig.2-3 paper

for( j in 1:2 )
{
  GENDER = ifelse( j == 1, "male", "female" )
  gender = ifelse( j == 1, "m", "f" )
  gender_title =  ifelse( j == 1, "men", "women" )
  gender_id = ifelse( j == 1, "B)", "A)" )
  
  for(i in lm_age)
  {
    comparing_scheme_to_plot =  stacked_comparing_scheme %>% 
      filter(landmark == i & sex == GENDER ) 
    
    #follow-up
    mean( ifelse( is.na(comparing_scheme_to_plot$cvd_age), comparing_scheme_to_plot$end_age - i, comparing_scheme_to_plot$cvd_age - i ) ) 
    sd(ifelse( is.na(comparing_scheme_to_plot$cvd_age), comparing_scheme_to_plot$end_age - i, comparing_scheme_to_plot$cvd_age - i ))  
    
    
    uf_extended = comparing_scheme_to_plot %>%
      select(patid, starts_with("uf") ) %>%
      gather(time_between_visit, value_uf, -patid) 
    
    benefit_extended = comparing_scheme_to_plot %>%
      select(patid, starts_with("benefit") ) %>%
      gather(time_between_visit, value_benefit, -patid)
    
    cost_extended = comparing_scheme_to_plot %>%
      select(patid, starts_with("cost") ) %>%
      gather(time_between_visit, value_cost, -patid)
    
    temp_data = bind_cols( uf_extended, benefit_extended, cost_extended ) %>%
      select(patid, time_between_visit, starts_with("value") ) %>%
      left_join(comparing_scheme_to_plot %>%
                  select(patid, risk_class, pred_time, pred_status, sex, landmark, best_of, ref_result) , by = c("patid", "patid"))%>%
      mutate( time_lapse = recode(time_between_visit,
                                  uf_1 = 1, uf_2 = 2, uf_3 = 3, 
                                  uf_4 = 4, uf_5 = 5, uf_6 = 6, 
                                  uf_7 = 7, uf_8 = 8, uf_9 = 9, uf_10 = 10 ))
    
    median_info_global = temp_data %>%
      filter(risk_class != "Very high risk" ) %>%
      group_by(risk_class, time_lapse)%>%
      summarise(median_inb = median(value_uf),
                median_benefit = median(value_benefit),
                median_cost = median(value_cost), 
                median_delta_benefit = median(value_benefit - ref_result),
                patid = 1)
    
    median_info_cross_tresh = temp_data %>%
      filter(risk_class != "Very high risk", pred_status == 1 ) %>%
      group_by(risk_class, time_lapse)%>%
      summarise(median_inb = median(value_uf),
                median_benefit = median(value_benefit),
                median_cost = median(value_cost), 
                median_delta_benefit = median(value_benefit - ref_result),
                patid = 1)
    
    median_info_NO_cross_tresh = temp_data %>%
      filter(risk_class != "Very high risk", pred_status == 0 ) %>%
      group_by(risk_class, time_lapse)%>%
      summarise(median_inb = median(value_uf),
                median_benefit = median(value_benefit),
                median_cost = median(value_cost), 
                median_delta_benefit = median(value_benefit - ref_result),
                patid = 1)
    
    if( i == 40 | i == 60 | i == 80 )
    {
      p_benefit_strategy = ggplot( subset( temp_data, risk_class != "Very high risk" ), aes( x = time_lapse, y = value_benefit, group = patid, col = risk_class ) ) +
        geom_line( show.legend = F ) +
        geom_point( show.legend = F ) +
        scale_x_continuous(breaks=1:10) +
        labs( x = "Time between visits [years]", y = sprintf('Benefit - strategy [\u00A3]'), title = paste0(gender_id, " Landmark age ", i, ", ", gender_title)) +
        theme_minimal() +
        facet_grid( . ~ risk_class ) +
        theme( text = element_text(size = 20) ) +
        geom_line(data = median_info_global, aes(y = median_benefit, x = time_lapse), size = 1, color = "gray50") +
        geom_line(data = median_info_cross_tresh, aes(y = median_benefit, x = time_lapse), size = 1, color = "black") +
        geom_line(data = median_info_NO_cross_tresh, aes(y = median_benefit, x = time_lapse), size = 1, color = "gray80")
      
      ggsave( file = paste0(  path_to_save_plot, "img_", GENDER, "/optimal_scheme_benefit_strategy_", i,"_", GENDER, ".pdf"), plot = p_benefit_strategy, width = 10, height = 10 )
      
      #For \Delta benefit, cost and INB, all individuals with pred_status == 0 follow the same path, so no need of representing all of them (the median is enough!)
      p_benefit_diff = ggplot( subset( temp_data, risk_class != "Very high risk" & pred_status == 1 ), aes( x = time_lapse, y = value_benefit-ref_result, group = patid, col = risk_class ) ) +
        geom_line( show.legend = F ) +
        geom_point( show.legend = F ) +
        scale_x_continuous(breaks=1:10) +
        labs( x = "Time between visits [years]", y = sprintf('\u0394 Benefit [\u00A3]'), title = paste0("Landmark age ", i, ", ", GENDER)) +
        theme_minimal() +
        facet_grid( . ~ risk_class ) +
        theme( text = element_text(size = 20) ) +
        geom_line(data = median_info_global, aes(y = median_delta_benefit, x = time_lapse), size = 1, color = "gray50") +
        geom_line(data = median_info_cross_tresh, aes(y = median_delta_benefit, x = time_lapse), size = 1, color = "black") +
        geom_line(data = median_info_NO_cross_tresh, aes(y = median_delta_benefit, x = time_lapse), size = 1, color = "gray80")
      
      ggsave( file = paste0(  path_to_save_plot, "img_", GENDER, "/optimal_scheme_diff_benefit_", i,"_", GENDER, ".pdf"), plot = p_benefit_diff, width = 10, height = 10 )
      
      p_cost_strategy = ggplot( subset( temp_data, risk_class != "Very high risk" & pred_status == 1), aes( x = time_lapse, y = value_cost, group = patid, col = risk_class ) ) +
        geom_line( show.legend = F ) +
        geom_point( show.legend = F ) +
        scale_x_continuous(breaks=1:10) +
        labs( x = "Time between visits [years]", y = sprintf('Cost - strategy [\u00A3]'), paste0(gender_id, " Landmark age ", i, ", ", gender_title)) +
        theme_minimal() +
        facet_grid( . ~ risk_class ) +
        theme( text = element_text(size = 20) )  +
        geom_line(data = median_info_global, aes(y = median_cost, x = time_lapse), size = 1, color = "gray50") +
        geom_line(data = median_info_cross_tresh, aes(y = median_cost, x = time_lapse), size = 1, color = "black") +
        geom_line(data = median_info_NO_cross_tresh, aes(y = median_cost, x = time_lapse), size = 1, color = "gray80")
      
      ggsave( file = paste0(  path_to_save_plot, "img_", GENDER, "/optimal_scheme_cost_strategy_", i,"_", GENDER, ".pdf"), plot = p_cost_strategy, width = 10, height = 10 )
    }
    
    p_INB = ggplot( subset( temp_data, risk_class != "Very high risk" & pred_status == 1), aes( x = time_lapse, y = value_uf, group = patid, col = risk_class ) ) +
      geom_line( show.legend = F ) +
      geom_point( show.legend = F ) +
      scale_x_continuous(breaks=1:10) +
      labs( x = "Time between visits [years]", y = sprintf('Incremental Net benefit [\u00A3]'), title = paste0(gender_id, " Landmark age ", i, ", ", gender_title)) +
      theme_minimal() +
      facet_grid( . ~ risk_class ) +
      theme( text = element_text(size = 20) )  +
      geom_line(data = median_info_global, aes(y = median_inb, x = time_lapse), size = 1, color = "black") #+
      #geom_line(data = median_info_cross_tresh, aes(y = median_inb, x = time_lapse), size = 1, color = "gray50") +
      #geom_line(data = median_info_NO_cross_tresh, aes(y = median_inb, x = time_lapse), size = 1, color = "gray80")
    
    
    ggsave( file = paste0(  path_to_save_plot, "img_", GENDER, "/optimal_scheme_INB_", i,"_", GENDER, ".pdf"), plot = p_INB, width = 10, height = 10 )
    
    #ANALYSING data of high risk people, recommended a 10-year CVD
    data_high_risk_10y = comparing_scheme_to_plot %>%
      filter(risk_class == "High risk" & best_of == 10 ) %>%
      select( patid, pred_status, pred_time, cvd_age, cvd_ind, end_age, death_age, death_ind)
    summary( ifelse( is.na(data_high_risk_10y$cvd_age), data_high_risk_10y$end_age - i, data_high_risk_10y$cvd_age - i ) ) 
    sd( ifelse( is.na(data_high_risk_10y$cvd_age), data_high_risk_10y$end_age - i, data_high_risk_10y$cvd_age - i )  )  
    table( data_high_risk_10y$pred_status )
    
    plot_detail_high_risk_10yscreen = comparing_scheme_to_plot %>%
      filter(risk_class == "High risk" & best_of == 10) %>%
      select( patid, pred_status, pred_time, cvd_age, cvd_ind, end_age, death_age, death_ind) %>%
      ggplot( aes( x = pred_time, fill = factor(pred_status) )) +
      geom_histogram(bins = 50, position = "identity", alpha = 0.5) +
      labs( x = "End age", fill = "Predicted \nstatins \ninitiation")
    
    ggsave( file = paste0(  path_to_save_plot, "img_", GENDER, "/plot_detail_high_risk_10yscreen_", i,"_", GENDER, ".pdf"), plot = plot_detail_high_risk_10yscreen, width = 10, height = 10 )
    
    

    table_novhr = comparing_scheme_to_plot %>%
      select( risk_class, best_of ) %>%
      filter( risk_class != "Very high risk" ) %>%  # now required with changes to dplyr::count()
      group_by( best_of ) %>%
      summarise(freq = n()/nrow(.), count = n())
    
    best_of_now = table( comparing_scheme_to_plot$best_of, comparing_scheme_to_plot$risk_class )
    # xtable( rbind( cbind( best_of_now, rowSums( best_of_now )), c(colSums( best_of_now ), "") ), digits = 0 )
    tab_now = cbind( rbind( table(  comparing_scheme_to_plot$best_of[which(comparing_scheme_to_plot$risk_class != "Very high risk")], comparing_scheme_to_plot$risk_class[which(comparing_scheme_to_plot$risk_class != "Very high risk")] ),
                            paste( table(comparing_scheme_to_plot$risk_class[which(comparing_scheme_to_plot$risk_class != "Very high risk")]), "(", round( table(comparing_scheme_to_plot$risk_class[which(comparing_scheme_to_plot$risk_class != "Very high risk")])/sum(table(comparing_scheme_to_plot$risk_class[which(comparing_scheme_to_plot$risk_class != "Very high risk")]))*100, 2),"%)" ) ),
                     c( paste( table_novhr$count, "(", round( table_novhr$freq * 100, 2 ), "%)"), "" ) )
    
    tab_now = tab_now[,-1]
    xtable(tab_now)
    
    best_of_now = sweep( best_of_now, 2, colSums( best_of_now ),`/` )
    best_of_now = data.frame( best_of_now )
    colnames(best_of_now) = c( "Frequency", "Risk_class", "Percentage")
    best_of_now$Frequency <- factor(best_of_now$Frequency, levels = 1:10)
    
    p_frac = ggplot( subset( best_of_now, Risk_class != "Very high risk"), aes(y = Percentage, x = Risk_class, fill = Frequency )) +
      geom_bar( stat = "identity", width = 0.75 ) +
      labs( x = "Risk category", fill = "Time \nbetween \nvisits[years]", y = "Proportion of the cohort", title = paste0(gender_id, " Landmark age ", i, ", ", gender_title) ) +
      scale_fill_brewer( palette = "Spectral", drop = F ) +
      #theme_bw( ) +
      theme(text = element_text(size = 20),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.background = element_rect(fill = "transparent",colour = NA),
            legend.background = element_rect(fill = "transparent",colour = NA) )
    # theme(axis.text=element_text(size=12),
    #       axis.title=element_text(size=14), 
    #       legend.text = element_text(size = 12), 
    #       legend.title=element_text(size=14),
    #       title = element_text(size=14))
    # 
    
    ggsave( file = paste0(  path_to_save_plot, "img_", GENDER, "/optimal_scheme_summary_", i,"_", GENDER, ".pdf"), plot = p_frac, width = 10, height = 10 )
    
    
    
    if( (i < 75 & gender == "m") | (i < 80 & gender == "f") )
    {
      p_tab = ggplot() + 
        theme_void() + 
        annotation_custom(tableGrob(tab_now, theme = ttheme_default(base_size = 6)),
                          xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)
      
      
      
      ggsave( file = paste0(  path_to_save_plot, "img_", GENDER,  "/optimal_table_summary_", i,"_", GENDER, ".pdf"), plot = p_tab, width = 10, height = 10 )
    }
    
  }
}

#Fig.4-5

for( j in 1:2 )
{
  GENDER = ifelse( j == 1, "male", "female" )
  gender = ifelse( j == 1, "men", "women" )
  title_gender =  paste(toupper(substring(GENDER, 1,1)), substring(GENDER, 2), sep="", collapse=" ")
  gender_id = ifelse( j == 1, "B)", "A)" )
  
  
  #we represent the y-axis as percentage
  
  plot_by_gender_freq_optimal_cvd_01 = stacked_comparing_scheme_with_factors %>%
    count(sex, landmark, risk_class, best_of, .drop = FALSE) %>%
    filter(sex == GENDER, risk_class != "Very high risk") %>%
    group_by(landmark) %>%
    mutate(freq_best_of_by_risk_class_lm = n/sum(n)) %>%
    ggplot( aes( x = best_of,  y = freq_best_of_by_risk_class_lm, fill = factor( best_of ))) +
    geom_col( position = position_dodge2(preserve = "single")) +
    labs(fill = "Optimal \ntime \nbetween \nvisits", x = "", y = "", title = paste0("Optimal risk-assessment distribution for ", gender)) +
    scale_fill_brewer( palette = "Spectral" ) +
    theme(text = element_text(size = 15),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          legend.background = element_rect(fill = "transparent",colour = NA) ) +
    facet_grid(rows = vars(risk_class), cols = vars(landmark))
  
  ggsave( file = paste0( path_to_save_plot, "img_", GENDER, "/plot_optimal_freq_cvd_01_yaxis_",GENDER,".pdf"), plot = plot_by_gender_freq_optimal_cvd_01, bg = "transparent", width = 10, height = 10 )
}


#Figure 6: sensitivity analysis
####################################################
#
#  Sensitivity analysis: bns
#
###################################################

bns_to_eval = seq( 20000, 30000, by = 500 ) #c( 100, 500, 1000, 2000, 3000, 4000, 5000, 10000)
cs = 150 #20
us= 0.997
cv = 18.39

sens_bns_optimal = matrix( 0, nrow = 10, ncol = length(bns_to_eval) )

count = 1
for( bns in bns_to_eval )
{
  
  INB_result = matrix( 0, nrow = dim(stacked_comparing_scheme), ncol = 10)
  
  for( f in 1:10 )
  {
    idx_NS_integral = grep(paste0(f,"_NS"), colnames(stacked_comparing_scheme))
    NS_integral = stacked_comparing_scheme[ ,idx_NS_integral]
    
    idx_S_integral = grep( paste0(f,"_S"), colnames(stacked_comparing_scheme))
    S_integral = stacked_comparing_scheme[ ,idx_S_integral]
    
    ref_integral = stacked_comparing_scheme$int_no_stat_all
    
    idx_N_visit = grep( paste0( f,"_N" ), colnames(stacked_comparing_scheme))[2]
    N_visit = stacked_comparing_scheme[ ,idx_N_visit]
    
    full_res = function_INB(bns, us, cs, cv, NS_integral, S_integral, ref_integral, N_visit)
    
    INB_result[ ,f ]     = full_res$uf
    
  }
  
  best_of_new = max.col( INB_result, ties.method = "last")
  bns_data_temp = data.frame( best_of = factor( best_of_new, levels = 1:10), risk_class = stacked_comparing_scheme$risk_class )
  
  #stacked_comparing_scheme$best_of = best_of_new
  
  
  sens_bns_optimal[  , count] = bns_data_temp %>%
    filter( risk_class != "Very high risk" ) %>%
    count( best_of, .drop = FALSE) %>%
    select(n) %>%
    unlist()

  count = count + 1        
}

sens_bns_optimal = data.frame(sens_bns_optimal)
apply( sens_bns_optimal, 1, mean)
apply( sens_bns_optimal, 1, sd)
sens_bns_optimal$optimal_freq = factor(1:10)
colnames(sens_bns_optimal)[1:length(bns_to_eval)] = as.character( bns_to_eval )
sens_bns_optimal_melt = reshape2::melt(sens_bns_optimal)

sens_bns_optimal_melt = sens_bns_optimal_melt %>% group_by(variable) %>% mutate(freq = value/sum(value))

sens_bns_plot = ggplot(sens_bns_optimal_melt, aes(x = variable, y = freq, col = optimal_freq, group = optimal_freq ), size = 0.8) +
  geom_point() +
  geom_line() +
  scale_color_brewer( palette = "Spectral" ) +
  labs(x = expression(lambda), y = "Proportion of the cohort", col = "Optimal\ntime \nbetween \nvisits") +
  theme_bw() +
  theme( text = element_text(size = 10),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) )

  #theme(text = element_text(size = 10),
  #      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


gg_sens_bns_plot = gg.gap( plot = sens_bns_plot,
                           ylim = c(0, 1),
                           segments = list(c(0.05, 0.6)),
                           tick_width = c(0.005, 0.1),
                           c(0.75,0,0.25) )

g_legend = function(a.gplot){
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)}

mylegend = g_legend(sens_bns_plot)

final_p_bns = grid.arrange(gg_sens_bns_plot + theme( plot.margin=unit(c(0.1,0.5,0.1,0.5), "cm")),
                           mylegend, ncol=2,widths=c(12, 3))
ggsave( file = paste0( path_to_save_plot, "/sensitivity_analysis_bns_all.pdf"), plot = final_p_bns, width = 7, height = 7 )

####################################################
#
#  Sensitivity analysis: us
#
###################################################

bns = 25000 #c( 100, 500, 1000, 2000, 3000, 4000, 5000, 10000)
us_vec = seq(0.997, 1, by = 0.0005)
cs = 150 #20
cv = 18.39

sens_us_optimal = matrix( 0, nrow = 10, ncol = length(us_vec) )

count = 1
for( us in us_vec )
{
  
  INB_result = matrix( 0, nrow = dim(stacked_comparing_scheme), ncol = 10)
  
  for( f in 1:10 )
  {
    idx_NS_integral = grep(paste0(f,"_NS"), colnames(stacked_comparing_scheme))
    NS_integral = stacked_comparing_scheme[ ,idx_NS_integral]
    
    idx_S_integral = grep( paste0(f,"_S"), colnames(stacked_comparing_scheme))
    S_integral = stacked_comparing_scheme[ ,idx_S_integral]
    
    ref_integral = stacked_comparing_scheme$int_no_stat_all
    
    idx_N_visit = grep( paste0( f,"_N" ), colnames(stacked_comparing_scheme))[2]
    N_visit = stacked_comparing_scheme[ ,idx_N_visit]
    
    full_res = function_INB(bns, us, cs, cv, NS_integral, S_integral, ref_integral, N_visit)
    
    INB_result[ ,f ]     = full_res$uf
    
  }
  
  best_of_new = max.col( INB_result, ties.method = "last")
  us_data_temp = data.frame( best_of = factor( best_of_new, levels = 1:10), risk_class = stacked_comparing_scheme$risk_class )
  
  sens_us_optimal[  , count] = us_data_temp %>%
    filter( risk_class != "Very high risk" ) %>%
    count( best_of, .drop = FALSE) %>%
    select(n) %>%
    unlist()
  
  
  count = count + 1        
}

sens_us_optimal = data.frame(sens_us_optimal)
apply( sens_us_optimal, 1, mean)
apply( sens_us_optimal, 1, sd)
sens_us_optimal$optimal_freq = factor(1:10)
colnames(sens_us_optimal)[1:length(us_vec)] = as.character( us_vec )
sens_us_optimal_melt = reshape2::melt(sens_us_optimal)

sens_us_optimal_melt = sens_us_optimal_melt %>% group_by(variable) %>% mutate(freq = value/sum(value))

sens_us_plot = ggplot(sens_us_optimal_melt, aes(x = variable, y = freq, col = optimal_freq, group = optimal_freq )) +
  geom_point() +
  geom_line() +
  scale_color_brewer( palette = "Spectral" ) +
  labs(x = expression(u[s]), y = "Proportion of the cohort", col = "Optimal\ntime \nbetween \nvisits") +
  theme_bw() + 
  theme( text = element_text(size = 10),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) )

gg_sens_us_plot = gg.gap( plot = sens_us_plot,
                          ylim = c(0, 1),
                          segments = list(c(0.05, 0.6)),
                          tick_width = c(0.005, 0.1),
                          c(0.75,0,0.25) )


mylegend = g_legend(sens_us_plot)

final_p_us = grid.arrange(gg_sens_us_plot + theme(plot.margin=unit(c(0.1,0,0.1,0.5), "cm")),
                          mylegend, ncol=2,widths=c(12, 3))
ggsave( file = paste0( path_to_save_plot, "/sensitivity_analysis_us_all.pdf"), plot = final_p_us, width = 7, height = 7 )



####################################################
#
#  Sensitivity analysis: cs
#
###################################################

bns = 25000 #c( 100, 500, 1000, 2000, 3000, 4000, 5000, 10000)
us = 0.997
cs_vec = seq(4,322, by = 25)  #20
cv = 18.39

sens_cs_optimal = matrix( 0, nrow = 10, ncol = length(cs_vec) )

count = 1
for( cs in cs_vec )
{
  INB_result = matrix( 0, nrow = dim(stacked_comparing_scheme), ncol = 10)
  
  for( f in 1:10 )
  {
    idx_NS_integral = grep(paste0(f,"_NS"), colnames(stacked_comparing_scheme))
    NS_integral = stacked_comparing_scheme[ ,idx_NS_integral]
    
    idx_S_integral = grep( paste0(f,"_S"), colnames(stacked_comparing_scheme))
    S_integral = stacked_comparing_scheme[ ,idx_S_integral]
    
    ref_integral = stacked_comparing_scheme$int_no_stat_all
    
    idx_N_visit = grep( paste0( f,"_N" ), colnames(stacked_comparing_scheme))[2]
    N_visit = stacked_comparing_scheme[ ,idx_N_visit]
    
    full_res = function_INB(bns, us, cs, cv, NS_integral, S_integral, ref_integral, N_visit)
    
    INB_result[ ,f ]     = full_res$uf
    
  }
  
  best_of_new = max.col( INB_result, ties.method = "last")
  cs_data_temp = data.frame( best_of = factor( best_of_new, levels = 1:10), risk_class = stacked_comparing_scheme$risk_class )
  
  sens_cs_optimal[  , count] = cs_data_temp %>%
    filter( risk_class != "Very high risk" ) %>%
    count( best_of, .drop = FALSE) %>%
    select(n) %>%
    unlist()
  
  count = count + 1        
}

sens_cs_optimal = data.frame(sens_cs_optimal)
apply( sens_cs_optimal, 1, mean)
apply( sens_cs_optimal, 1, sd)
sens_cs_optimal$optimal_freq = factor(1:10)
colnames(sens_cs_optimal)[1:length(cs_vec)] = as.character( cs_vec )
sens_cs_optimal_melt = reshape2::melt(sens_cs_optimal)

sens_cs_optimal_melt = sens_cs_optimal_melt %>% group_by(variable) %>% mutate(freq = value/sum(value))

sens_cs_plot = ggplot(sens_cs_optimal_melt, aes(x = variable, y = freq, col = optimal_freq, group = optimal_freq ), size = 0.8) +
  geom_point() +
  geom_line() +
  scale_color_brewer( palette = "Spectral" ) +
  labs(x = expression(c[s]), y = "Proportion of the cohort", col = "Optimal\ntime \nbetween \nvisits") +
  theme_bw() + 
  theme( text = element_text(size = 10),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) )


gg_sens_cs_plot = gg.gap( plot = sens_cs_plot,
                          ylim = c(0, 1),
                          segments = list(c(0.05, 0.6)),
                          tick_width = c(0.005, 0.1),
                          c(0.75,0,0.25) )


mylegend = g_legend(sens_cs_plot)

final_p_cs = grid.arrange(gg_sens_cs_plot + theme(text = element_text(size = 10), plot.margin=unit(c(0.1,0,0.1,0.5), "cm")),
                          mylegend, ncol=2,widths=c(12, 3))
ggsave( file = paste0( path_to_save_plot, "/sensitivity_analysis_cs_all.pdf"), plot = final_p_cs, width = 7, height = 7 )



final_p_total = grid.arrange(gg_sens_bns_plot + labs(tag = "A)") + theme(plot.margin=unit(c(0.1,0.1,0.1,0.5), "cm")), 
                             gg_sens_us_plot + labs(tag = "B)") + theme(plot.margin=unit(c(0.1,0.1,0.1,0.5), "cm")),
                             gg_sens_cs_plot + labs(tag = "C)") + theme(plot.margin=unit(c(0.1,0.1,0.1,0.5), "cm")),
                             mylegend, ncol = 4, widths = c( 5, 5, 5, 2 ) )

ggsave( file = paste0( path_to_save_plot, "/sensitivity_analysis_all.pdf"), plot = final_p_total, width = 30, height = 21, units = "cm" )


#Table 1 and 2
i = 40

comparing_scheme_to_plot_F =  stacked_comparing_scheme %>% 
      filter(landmark == i & sex == "female" ) 

comparing_scheme_to_plot_M =  stacked_comparing_scheme %>% 
  filter(landmark == i & sex == "male" ) 

table_novhr_F = comparing_scheme_to_plot_F %>%
      select( risk_class, best_of ) %>%
      filter( risk_class != "Very high risk" ) %>%  # now required with changes to dplyr::count()
      group_by( best_of ) %>%
      summarise(freq = n()/nrow(.), count = n())

table_novhr_M = comparing_scheme_to_plot_M %>%
  select( risk_class, best_of ) %>%
  filter( risk_class != "Very high risk" ) %>%  # now required with changes to dplyr::count()
  group_by( best_of ) %>%
  summarise(freq = n()/nrow(.), count = n())


best_of_now_F = table( comparing_scheme_to_plot_F$best_of, comparing_scheme_to_plot_F$risk_class )
tab_now_F = cbind( c("women", rep("",7)), 
                   c(sort(unique(comparing_scheme_to_plot_F$best_of)), "Total"), 
                   rbind( table(  comparing_scheme_to_plot_F$best_of[which(comparing_scheme_to_plot_F$risk_class != "Very high risk")], comparing_scheme_to_plot_F$risk_class[which(comparing_scheme_to_plot_F$risk_class != "Very high risk")] ),
                        paste0( table(comparing_scheme_to_plot_F$risk_class[which(comparing_scheme_to_plot_F$risk_class != "Very high risk")]), " (", round( table(comparing_scheme_to_plot_F$risk_class[which(comparing_scheme_to_plot_F$risk_class != "Very high risk")])/sum(table(comparing_scheme_to_plot_F$risk_class[which(comparing_scheme_to_plot_F$risk_class != "Very high risk")]))*100, 2),"%)" ) ),
                 c( paste0( table_novhr_F$count, " (", round( table_novhr_F$freq * 100, 2 ), "%)"), "" ) )


best_of_now_M = table( comparing_scheme_to_plot_M$best_of, comparing_scheme_to_plot_M$risk_class )
tab_now_M = cbind( c("men", rep("",8)), 
                   c(sort(unique(comparing_scheme_to_plot_M$best_of)), "Total"), 
                   rbind( table(  comparing_scheme_to_plot_M$best_of[which(comparing_scheme_to_plot_M$risk_class != "Very high risk")], comparing_scheme_to_plot_M$risk_class[which(comparing_scheme_to_plot_M$risk_class != "Very high risk")] ),
                          paste0( table(comparing_scheme_to_plot_M$risk_class[which(comparing_scheme_to_plot_M$risk_class != "Very high risk")]), " (", round( table(comparing_scheme_to_plot_M$risk_class[which(comparing_scheme_to_plot_M$risk_class != "Very high risk")])/sum(table(comparing_scheme_to_plot_M$risk_class[which(comparing_scheme_to_plot_M$risk_class != "Very high risk")]))*100, 2),"%)" ) ),
                   c( paste0( table_novhr_M$count, " (", round( table_novhr_M$freq * 100, 2 ), "%)"), "" ) )


tab_final = rbind(tab_now_F[,-3], tab_now_M[,-3])

print( xtable(tab_final), include.rownames = F)


#Fig.3-6 supplementary
rm(list=ls())

library(survival)
path_to_save = "~/decision_theory/img/"


landmark_ages = seq( 40, 80, by = 5)
cindex_10_cvd_tot =  NULL


###Computing the Brier scores

for( j in 1:2 ) 
{  
  GENDER = ifelse( j == 1, "male", "female" )
  gender = ifelse( j == 1, "m", "f" )
  for( i in landmark_ages )
  {
    load.Rdata( paste0( "results_", GENDER,"/data_cindex_10CVD_RIRS_", i, "_", GENDER, ".RData"), paste0("cindex_10_cvd_", i ) )

    #check on colnames
    cindex_10_cvd_tot = rbind( cindex_10_cvd_tot,
                               cbind( get( paste0("cindex_10_cvd_", i ) ), gender ) )
    
  }
}


#10-year CVD risk
BS_10_CVD_lm = data.frame( gender = rep( c( "male", "female" ), each = 9 ),
                           lm_age = rep( landmark_ages, 2 ),
                           value  = 0  )

for( j in 1:2 )
{ 
  gender = ifelse( j == 1, "m", "f" )
  GENDER = ifelse( j == 1, "male", "female" )
  
  count = 1
  for( i in landmark_ages )
  {
    idx_temp = which( cindex_10_cvd_tot$lm_age == i & cindex_10_cvd_tot$gender == GENDER)
    idx_temp_status1 = which(cindex_10_cvd_tot$lm_age == i & cindex_10_cvd_tot$gender == GENDER & cindex_10_cvd_tot$status_10CVD == 1 )
    
    cens_data_temp = data.frame(time_cens = cindex_10_cvd_tot$time_10CVD[ idx_temp],
                                status_cens = 1 - cindex_10_cvd_tot$status_10CVD[ idx_temp ])
    
    cens_pred_10CVD_temp = survfit( Surv( time_cens, status_cens ) ~ 1, data = cens_data_temp)
    
    idx_cens = match( cindex_10_cvd_tot$time_10CVD[ idx_temp_status1 ], cens_pred_10CVD_temp$time)
    ipcw_temp = cens_pred_10CVD_temp$surv[ idx_cens ]
    summary( ipcw_temp )
    score_brier_i = (cindex_10_cvd_tot$survival_expected_10CVD[idx_temp_status1])^2
    summary(score_brier_i)
    
    #monitoring for too low IPCW 
    idx_not_stable = which(ipcw_temp < 0.001)
    if( length( idx_not_stable ) != 0 )
    {
      score_brier_i  = score_brier_i[ -idx_not_stable ]
      ipcw_temp      = ipcw_temp[ -idx_not_stable ]
    }
    
    BS_10_CVD_lm$value[ count + 9*( j - 1 ) ] = sum( score_brier_i/(length(idx_temp)*ipcw_temp) )
    count = count + 1
  }
}

BS_10_CVD_lm$gender_small = ifelse(BS_10_CVD_lm$gender == "female", "women", "men")


brier_10CVD_no_overall = ggplot( BS_10_CVD_lm, aes(x = lm_age, y = value, col = gender_small) ) +
  geom_point(  size = 5 ) +
  labs(y = "Brier Score", x = "Landmark age", col = "Gender") +
  scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
  theme_bw() +
  theme( text = element_text( size=30), legend.title = element_blank() )

ggsave( brier_10CVD_no_overall, file = paste0( path_to_save, "brier_10CVD_no_overall.pdf" ), width = 30, height = 21, units = "cm" )


#5-year CVD risk
pe_5_cvd = NULL

for( j in 1:2 )
{ 
  gender = ifelse( j == 1, "male", "female" )
  for( i in landmark_ages )
  {
    load.Rdata( paste0( "results_", gender,"/pred_error_5CVD_RIRS_", i, "_",gender, ".RData"), paste0("pe_5_cvd_", i ) ) 
    
    pe_5_cvd = rbind( pe_5_cvd, get( paste0("pe_5_cvd_", i ) ) )
  }
}


pe_5_cvd$gender_small = ifelse(pe_5_cvd$gender == "female", "women", "men")

brier_5CVD_no_overall = ggplot( subset( pe_5_cvd, type %in% c("BS_pec")), aes( x = start_time, y = value, col = as.factor(lm_age) )) +
  geom_point(  aes( fill = factor(lm_age), colour = factor(lm_age), alpha = 0.5 ), shape = 21, size = 4 ) +
  scale_colour_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  guides( color = FALSE, alpha = FALSE ) + 
  labs(x = "Landmark age", fill = "LM", y = "Brier score") +
  facet_wrap(vars(gender_small)) +
  theme_bw() +
  theme( text = element_text( size=30), legend.title = element_blank() )

ggsave( brier_5CVD_no_overall, file = paste0( path_to_save, "brier_5CVD_no_overall.pdf" ), width = 40, height = 25, units = "cm" )


#Computing c-index

cindex_10_composite_tot =  NULL
cindex_10_death_tot = NULL
cindex_10_cvd_tot =  NULL


###Computing the overall c-index

for( j in 1:2 ) 
{  
  GENDER = ifelse( j == 1, "male", "female" )
  gender = ifelse( j == 1, "m", "f" )
  for( i in landmark_ages )
  {
    load.Rdata( paste0( "results_", GENDER,"/data_cindex_10composite_RIRS_", i, "_", GENDER, ".RData"), paste0("cindex_10_composite_", i ) )  
    load.Rdata( paste0( "results_", GENDER,"/data_cindex_10CVD_RIRS_", i, "_", GENDER, ".RData"), paste0("cindex_10_cvd_", i ) )
    load.Rdata( paste0( "results_", GENDER,"/data_cindex_10death_RIRS_", i, "_", GENDER, ".RData"), paste0("cindex_10_death_", i ) )
    
    #check on colnames
    cindex_10_composite_tot = rbind( cindex_10_composite_tot, 
                                     cbind( get( paste0("cindex_10_composite_", i ) ), gender ) )
    cindex_10_death_tot = rbind( cindex_10_death_tot,
                                 cbind( get( paste0("cindex_10_death_", i ) ), gender ) )
    cindex_10_cvd_tot = rbind( cindex_10_cvd_tot,
                               cbind( get( paste0("cindex_10_cvd_", i ) ), gender ) )
    
  }
}



#overall c-idex
#cvd
surv_info_cvd = Surv( time = cindex_10_cvd_tot$time_10CVD, event = cindex_10_cvd_tot$status_10CVD)
cindex_temp_cvd = concordance( surv_info_cvd ~ cindex_10_cvd_tot$survival_expected_10CVD, cluster = cindex_10_cvd_tot$cluster_patid )
#cindex_temp_cvd_mine = concordance( surv_info_cvd ~ I(-cindex_10_cvd_tot$lp_val_10CVD),  cluster = cindex_10_cvd_tot$cluster_patid )
overall_cindex_cvd = data.frame( type = "CVD", conc = cindex_temp_cvd$concordance, conc_sd = sqrt( cindex_temp_cvd$var ) )

#death
surv_info_death = Surv( time = cindex_10_death_tot$time_10death, event = cindex_10_death_tot$status_10death)
cindex_temp_death = concordance( surv_info_death ~ cindex_10_death_tot$survival_expected_10death, cluster = cindex_10_death_tot$cluster_patid)
#cindex_temp_death_mine = concordance( surv_info_death ~ I(-cindex_10_death_tot$lp_val_10death), cluster = cindex_10_death_tot$cluster_patid)
overall_cindex_death = data.frame( type = "Death", conc = cindex_temp_death$concordance, conc_sd = sqrt( cindex_temp_death$var ) )

#composite
surv_info_composite = Surv( time = cindex_10_composite_tot$time_10composite, event = cindex_10_composite_tot$status_10composite)
cindex_temp_composite = concordance( surv_info_composite ~ cindex_10_composite_tot$survival_expected_10composite, cluster = cindex_10_cvd_tot$cluster_patid)
#cindex_temp_composite_mine = concordance( surv_info_composite ~ I(-cindex_10_composite_tot$lp_val_10composite), cluster = cindex_10_composite_tot$cluster_patid)
overall_cindex_composite = data.frame( type = "Composite", conc =  cindex_temp_composite$concordance, conc_sd = sqrt( cindex_temp_composite$var ) )

overall_cindex_cvd_by_gender = NULL
overall_cindex_death_by_gender = NULL
overall_cindex_composite_by_gender = NULL

cindex_cvd_10_global = NULL
cindex_death_10_global = NULL
cindex_composite_10_global = NULL

for( j in 1:2 ){
  GENDER = ifelse( j == 1, "male", "female" )
  gender = ifelse( j == 1, "m", "f" )
  
  #cvd
  idx_oi_cvd = which( cindex_10_cvd_tot$gender == GENDER )
  surv_info_cvd = Surv( time = cindex_10_cvd_tot$time_10CVD[ idx_oi_cvd ], event = cindex_10_cvd_tot$status_10CVD[ idx_oi_cvd ])
  cindex_temp_cvd = concordance( surv_info_cvd ~ cindex_10_cvd_tot$survival_expected_10CVD[ idx_oi_cvd ], cluster = cindex_10_cvd_tot$cluster_patid[ idx_oi_cvd ] )
  overall_cindex_cvd_by_gender_temp = data.frame( GENDER, type = "CVD", conc = cindex_temp_cvd$concordance, conc_sd = sqrt( cindex_temp_cvd$var ) )
  overall_cindex_cvd_by_gender = rbind( overall_cindex_cvd_by_gender,  overall_cindex_cvd_by_gender_temp )
  
  #death
  idx_oi_death = which( cindex_10_death_tot$gender == GENDER )
  surv_info_death = Surv( time = cindex_10_death_tot$time_10death[ idx_oi_death ], event = cindex_10_death_tot$status_10death[ idx_oi_death ])
  cindex_temp_death = concordance( surv_info_death ~ cindex_10_death_tot$survival_expected_10death[idx_oi_death], cluster = cindex_10_death_tot$cluster_patid[idx_oi_death])
  overall_cindex_death_by_gender_temp = data.frame( GENDER, type = "Death", conc = cindex_temp_death$concordance, conc_sd = sqrt( cindex_temp_death$var ) )
  overall_cindex_death_by_gender = rbind( overall_cindex_death_by_gender,  overall_cindex_death_by_gender_temp )
  
  #composite
  idx_oi_composite = which( cindex_10_composite_tot$gender == GENDER )
  surv_info_composite = Surv( time = cindex_10_composite_tot$time_10composite[ idx_oi_composite ], event = cindex_10_composite_tot$status_10composite[ idx_oi_composite ])
  cindex_temp_composite = concordance( surv_info_composite ~ cindex_10_composite_tot$survival_expected_10composite[idx_oi_composite], cluster = cindex_10_cvd_tot$cluster_patid[idx_oi_composite])
  overall_cindex_composite_by_gender_temp = data.frame( GENDER, type = "Composite", conc =  cindex_temp_composite$concordance, conc_sd = sqrt( cindex_temp_composite$var ) )
  overall_cindex_composite_by_gender = rbind( overall_cindex_composite_by_gender,  overall_cindex_composite_by_gender_temp )
  
  
  for( i in landmark_ages )
  {
    #cvd
    idx_oi_cvd = which(cindex_10_cvd_tot$lm_age == i & cindex_10_cvd_tot$gender == GENDER )
    surv_info_cvd = Surv( time = cindex_10_cvd_tot$time_10CVD[idx_oi_cvd], event = cindex_10_cvd_tot$status_10CVD[idx_oi_cvd])
    cindex_temp_cvd = concordance( surv_info_cvd ~ cindex_10_cvd_tot$survival_expected_10CVD[idx_oi_cvd])
    #cindex_temp_cvd_mine = concordance( surv_info_cvd ~ I(-cindex_10_cvd_tot$lp_val_10CVD[idx_oi_cvd]))
    line_to_add_cvd = c( GENDER, i, cindex_temp_cvd$concordance, sqrt( cindex_temp_cvd$var ) )
    cindex_cvd_10_global = rbind( cindex_cvd_10_global, line_to_add_cvd )
    
    #death
    idx_oi_death = which(cindex_10_death_tot$lm_age == i & cindex_10_death_tot$gender == GENDER )
    surv_info_death = Surv( time = cindex_10_death_tot$time_10death[idx_oi_death], event = cindex_10_death_tot$status_10death[idx_oi_death])
    cindex_temp_death = concordance( surv_info_death ~ cindex_10_death_tot$survival_expected_10death[idx_oi_death])
    #cindex_temp_death_mine = concordance( surv_info_death ~ I(-cindex_10_death_tot$lp_val_10eath[idx_oi_death]))
    line_to_add_death = c( GENDER, i, cindex_temp_death$concordance, sqrt( cindex_temp_death$var ) )
    cindex_death_10_global = rbind( cindex_death_10_global, line_to_add_death )
    
    #composite
    idx_oi_composite = which(cindex_10_composite_tot$lm_age == i & cindex_10_composite_tot$gender == GENDER )
    surv_info_composite = Surv( time = cindex_10_composite_tot$time_10composite[ idx_oi_composite ], event = cindex_10_composite_tot$status_10composite[ idx_oi_composite ])
    cindex_temp_composite = concordance( surv_info_composite ~ cindex_10_composite_tot$survival_expected_10composite[idx_oi_composite])
    #cindex_temp_composite_mine = concordance( surv_info_composite ~ I(-cindex_10_composite_tot$lp_val_10composite[idx_oi_composite]))
    line_to_add_composite = c( GENDER, i, cindex_temp_composite$concordance, sqrt( cindex_temp_composite$var) )
    cindex_composite_10_global = rbind( cindex_composite_10_global, line_to_add_composite )
  }
}

colnames(cindex_cvd_10_global) = colnames( cindex_death_10_global ) = colnames( cindex_composite_10_global ) = c("GENDER", "lm_age", "cindex","cindex_se")


cindex_cvd_10_global = as.data.frame(cindex_cvd_10_global)
cindex_death_10_global = as.data.frame(cindex_death_10_global)
cindex_composite_10_global = as.data.frame(cindex_composite_10_global)

cindex_cvd_10_global$type = "CVD"
cindex_death_10_global$type = "Death"
cindex_composite_10_global$type = "Composite"

cindex_cvd_10_global$overall_cindex = overall_cindex_cvd$conc
cindex_cvd_10_global$overall_cindex_se = overall_cindex_cvd$conc_sd
cindex_cvd_10_global = merge( cindex_cvd_10_global, overall_cindex_cvd_by_gender[,c("GENDER", "conc", "conc_sd")], by = "GENDER" )
colnames(cindex_cvd_10_global)[8:9] = c("overall_cindex_by_gender", "overall_cindex_se_by_gender")


cindex_death_10_global$overall_cindex = overall_cindex_death$conc
cindex_death_10_global$overall_cindex_se = overall_cindex_death$conc_sd
cindex_death_10_global = merge( cindex_death_10_global, overall_cindex_death_by_gender[,c("GENDER", "conc", "conc_sd")], by = "GENDER" )
colnames(cindex_death_10_global)[8:9] = c("overall_cindex_by_gender", "overall_cindex_se_by_gender")

cindex_composite_10_global$overall_cindex = overall_cindex_composite$conc
cindex_composite_10_global$overall_cindex_se = overall_cindex_composite$conc_sd
cindex_composite_10_global = merge( cindex_composite_10_global, overall_cindex_composite_by_gender[,c("GENDER", "conc", "conc_sd")], by = "GENDER" )
colnames(cindex_composite_10_global)[8:9] = c("overall_cindex_by_gender", "overall_cindex_se_by_gender")



cindex_10_global = rbind( cindex_cvd_10_global,
                          cindex_death_10_global, 
                          cindex_composite_10_global)

str(cindex_10_global)
cindex_10_global$lm_age = as.numeric(cindex_10_global$lm_age)
cindex_10_global$cindex = as.numeric(cindex_10_global$cindex)
cindex_10_global$cindex_se = as.numeric(cindex_10_global$cindex_se)

cindex_10_global$UB_overall = cindex_10_global$overall_cindex + qnorm(0.975) * cindex_10_global$overall_cindex_se  
cindex_10_global$LB_overall = cindex_10_global$overall_cindex - qnorm(0.975) * cindex_10_global$overall_cindex_se  

cindex_10_global$gender_small = ifelse(cindex_10_global$GENDER == "female", "women", "men")

p_10CVD =  ggplot(subset( cindex_10_global, type == "CVD"), aes( x = lm_age, y = cindex, col = gender_small ) ) +
  geom_pointrange(aes(ymin=cindex-cindex_se, ymax=cindex+cindex_se), size = 0.9 ) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 1.2) +
  geom_line( aes(y = overall_cindex), color = "black", size = 1.2) +
  geom_line( aes(y = UB_overall), color = "black", linetype = "dashed", size = 1.2) +
  geom_line( aes(y = LB_overall), color = "black", linetype = "dashed", size = 1.2) +
  labs(y = "c-index", x = "Landmark age", main = i, col = "Gender") +
  scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
  theme_bw() +
  theme( text = element_text( size=30), legend.title = element_blank() )

p_10CVD

ggsave( p_10CVD, file = paste0(path_to_save,"cindex_10CVD.pdf"), width = 30, height = 21, units = "cm")

#5year CVD

cindex_5_cvd_tot = NULL

for( j in 1:2 ) 
{  
  GENDER = ifelse( j == 1, "male", "female" )
  gender = ifelse( j == 1, "m", "f" )
  for( i in landmark_ages )
  {
    load.Rdata( paste0( "results_", GENDER,"/data_cindex_5CVD_RIRS_", i, "_", GENDER, ".RData"), paste0("cindex_5_cvd_", i ) )
    
    #check on colnames
    cindex_5_cvd_tot = rbind( cindex_5_cvd_tot,
                              cbind( get( paste0("cindex_5_cvd_", i ) ), GENDER ) )
    
  }
}

dim(cindex_5_cvd_tot)
surv_info_cvd = Surv( time = cindex_5_cvd_tot$time_5CVD, event = cindex_5_cvd_tot$status_5CVD)
cindex_temp_cvd = concordance( surv_info_cvd ~ cindex_5_cvd_tot$survival_expected_5CVD, cluster = cindex_5_cvd_tot$cluster_patid_5CVD )
overall_cindex_5cvd = data.frame( type = "5CVD", conc = cindex_temp_cvd$concordance, conc_sd = sqrt( cindex_temp_cvd$var ) )

overall_cindex_cvd_5_by_gender = NULL
cindex_cvd_5_global = NULL

for( j in 1:2 ){
  GENDER = ifelse( j == 1, "male", "female" )
  gender = ifelse( j == 1, "m", "f" )
  
  #cvd
  idx_oi_cvd = which( cindex_5_cvd_tot$gender == GENDER )
  surv_info_cvd = Surv( time = cindex_5_cvd_tot$time_5CVD[ idx_oi_cvd ], event = cindex_5_cvd_tot$status_5CVD[ idx_oi_cvd ])
  cindex_temp_cvd = concordance( surv_info_cvd ~ cindex_5_cvd_tot$survival_expected_5CVD[ idx_oi_cvd ], cluster = cindex_5_cvd_tot$cluster_patid_5CVD[ idx_oi_cvd ] )
  overall_cindex_cvd_5_by_gender_temp = data.frame( gender = GENDER, conc = cindex_temp_cvd$concordance, conc_sd = sqrt( cindex_temp_cvd$var ) )
  overall_cindex_cvd_5_by_gender = rbind( overall_cindex_cvd_5_by_gender,  overall_cindex_cvd_5_by_gender_temp )
  
  for( i in landmark_ages )
  {
    #cvd
    idx_oi_cvd = which(cindex_5_cvd_tot$lm_age == i & cindex_5_cvd_tot$gender == GENDER )
    surv_info_cvd = Surv( time = cindex_5_cvd_tot$time_5CVD[ idx_oi_cvd ], event = cindex_5_cvd_tot$status_5CVD[idx_oi_cvd])
    cindex_temp_cvd = concordance( surv_info_cvd ~ cindex_5_cvd_tot$survival_expected_5CVD[idx_oi_cvd], cluster = cindex_5_cvd_tot$cluster_patid_5CVD[ idx_oi_cvd ] )
    line_to_add_cvd = c( GENDER, i, cindex_temp_cvd$concordance, sqrt( cindex_temp_cvd$var ) )
    cindex_cvd_5_global = rbind( cindex_cvd_5_global, line_to_add_cvd )
  }
}

head(cindex_cvd_5_global)
cindex_cvd_5_global = data.frame(cindex_cvd_5_global)
colnames(cindex_cvd_5_global) = c("gender","lm_age","cindex","cindex_se")
cindex_cvd_5_global$lm_age = as.numeric(cindex_cvd_5_global$lm_age)
cindex_cvd_5_global$cindex = as.numeric(cindex_cvd_5_global$cindex)
cindex_cvd_5_global$cindex_se = as.numeric(cindex_cvd_5_global$cindex_se)
cindex_cvd_5_global$overall_cindex =  overall_cindex_5cvd$conc
cindex_cvd_5_global$overall_cindex_se =  overall_cindex_5cvd$conc_sd
cindex_cvd_5_global = merge( cindex_cvd_5_global, overall_cindex_cvd_5_by_gender, by = "gender" )

cindex_cvd_5_global$UB_overall = cindex_cvd_5_global$overall_cindex + qnorm(0.975)*cindex_cvd_5_global$overall_cindex_se  
cindex_cvd_5_global$LB_overall = cindex_cvd_5_global$overall_cindex - qnorm(0.975)*cindex_cvd_5_global$overall_cindex_se  

colnames(cindex_cvd_5_global)[7:8] = c("overall_cindex_by_gender", "overall_cindex_se_by_gender")


pe_5_cvd = NULL

for( j in 1:2 )
{ 
  gender = ifelse( j == 1, "male", "female" )
  for( i in landmark_ages )
  {
    load.Rdata( paste0( "results_", gender,"/pred_error_5CVD_RIRS_", i, "_",gender, ".RData"), paste0("pe_5_cvd_", i ) ) 
    
    pe_5_cvd = rbind( pe_5_cvd, get( paste0("pe_5_cvd_", i ) ) )
    
    
  }
}

cindex_cvd_5_global$key = paste0( cindex_cvd_5_global$gender, "_", cindex_cvd_5_global$lm_age )
pe_5_cvd$key =paste0( pe_5_cvd$gender, "_", pe_5_cvd$lm_age) 

pe_5_cvd_final = pe_5_cvd %>%
  left_join( cindex_cvd_5_global, by = c( "key" = "key" ) )%>%
  dplyr::select(ends_with(".x"), "start_time", "key", "type", "value", "sd", "cindex", "cindex_se", "overall_cindex", "overall_cindex_se" ) %>%
  mutate( gender = gender.x, lm_age = lm_age.x, cindex_punctual = value, cindex_punctual = sd,
          cindex_lm = cindex, cindex_lm_se = cindex_se )

pe_5_cvd_final$gender = relevel(pe_5_cvd_final$gender, "female")

pe_5_cvd_final$gender_small = ifelse(pe_5_cvd_final$gender == "female", "women", "men")


p_5_cvd_cindex_no_overall = ggplot( subset( pe_5_cvd_final, type == "harrell_c_index"), aes( x = start_time, y = value, col = as.factor(lm_age)) ) +
  geom_pointrange(aes(ymin=value-sd, ymax=value+sd), size = 0.9,  alpha = 0.6 ) + 
  geom_hline( yintercept = 0.5, linetype = "dashed", color = "red", size = 1.2) +
  labs(main = "Discriminatory indices - CVD", x = "Landmark age", color = "LM", y = "c-index") +
  scale_colour_brewer(palette = "Paired") +
  facet_wrap(vars(gender_small)) +
  theme_bw() +
  theme( text = element_text( size=30), legend.title = element_blank() )


p_5_cvd_cindex_no_overall

ggsave(p_5_cvd_cindex_no_overall, file = paste0(path_to_save,"cindex_5CVD_no_overall.pdf"), width = 40, height = 25, units = "cm" )


