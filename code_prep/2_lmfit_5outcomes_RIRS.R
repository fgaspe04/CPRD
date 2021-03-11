# #print(commandArgs(trailingOnly = TRUE))

#LOCAL = FALSE

args=commandArgs(trailingOnly = TRUE)
j      = as.integer(args[1])
gender = as.character(args[2]) #g can be male or female

library(data.table)
library(nlme)
library(survival)
library(tidyverse)
library(haven)
library(miceadds)


load.Rdata( paste0( "../../data_created/", gender, "/data_lm_",j, "_", gender,".RData"), paste0( "data_ext") )
path_to_save = paste0("../../results/results_", gender )



cat(j,"\n")

fit_rirs = lme( scaled_corr ~ 0 + exposure + exp_age_corr:exposure + bp_bin:sbp_ind + statin_bin:tchol_ind,  
                random=~ 0 + exposure + exp_age_corr:exposure|patid,
                weights=varIdent(form=~1|exposure),
                data = data_ext[ derivation == 1, ],
                control = lmeControl( maxIter=1000,msMaxIter=1000,msVerbose=TRUE,rel.tol=1e-6,msMaxEval=1000, niterEM = 40))

save( fit_rirs, file = paste0( path_to_save, "/test_fit_rirs",j,"_", gender,".RData" ) )

