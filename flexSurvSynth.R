library(dplyr)
library(Hmisc)
library(synthpop)
library(flexsurv)
library(synthpop)

flexSurvSynth=function(input_data,time_column,event_colum,period_column,dropout_column,covariates,k,knots,n,censtime){
  selected_cols=c(time_column,event_column,covariates,period_column,dropout_column)
  print(selected_cols)
  data_to_process=input_data[selected_cols]
  colnames(data_to_process)[colnames(data_to_process) == time_column] <- "tsurv"
  colnames(data_to_process)[colnames(data_to_process) == event_column] <- "surv"
  
  formula <- as.formula(paste("Surv(tsurv,surv) ~ ", paste(covariates, collapse = " + ")))
  surv_model <- flexsurvspline(formula, data = data_to_process, scale = scale, k = k)
  
  vars=c(covariates,period_column,dropout_column)
  vars_indices=match(vars,names(input_data))


  synth.obj = syn(input_data, visit.sequence = vars_indices,
                  k = n*nrow(data_to_process),  
                  seed = 14, 
                  method = "parametric")
  print(nrow(synth.obj))
  db_syn = synth.obj$syn  #extracting the generated dataset
  db_syn = db_syn[,vars_indices]
  
  start = min(input_data$dt_ci)
  
  db_syn$start = start
  db_syn$dt_ci = db_syn$start + db_syn$time_from_start
  
  # "simulate.flexsurvreg" -> time-to-event data froma a fitted flesurvreg
  synthetic_surv = simulate(surv_model, newdata = db_syn, nsim = 1,
                            tidy = TRUE, seed = 142, 
                            censtime = censtime) 
  # using "tidy=T", the covariates are associated to each simulated time 
  synthetic_surv$time = round(synthetic_surv$time,0)
  
  #last synthetic follow-up
  synthetic_surv$tsurv = synthetic_surv$time  
  
  synthetic_surv$surv = ifelse(synthetic_surv$dropout==1,0,
                               synthetic_surv$event)
  
  return(synthetic_surv)
}
