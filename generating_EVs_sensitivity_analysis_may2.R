library(tidyverse)
library(reshape2)

## All code relevant to estimating_drug_efficacy_curves.R ##
pfizer_data=cbind.data.frame(days=c(rep(3,697),rep(3,682),rep(5,1039),rep(5,1046)),
                             hosp=c(rep(1,5),rep(0,697-5),c(rep(1,44),rep(0,682-44),
                                                            c(rep(1,8), rep(0,1039-8),c(rep(1,66), rep(0,1046-66))))),
                             treated=c(rep('treatment',697),rep('placebo',682),rep('treatment',1039),
                                       rep('placebo',1046)))

log_binomial=glm(hosp~days+treated+days*treated,data=pfizer_data,family=binomial(link='log'))

RR_period1=exp(coefficients(log_binomial)[[3]]+0*coefficients(log_binomial)[[4]]) 
RR_period2=exp(coefficients(log_binomial)[[3]]+1*coefficients(log_binomial)[[4]]) 
RR_period3=exp(coefficients(log_binomial)[[3]]+2*coefficients(log_binomial)[[4]]) 
RR_period4=exp(coefficients(log_binomial)[[3]]+3*coefficients(log_binomial)[[4]])
RR_period5=exp(coefficients(log_binomial)[[3]]+4*coefficients(log_binomial)[[4]])
RR_period6=exp(coefficients(log_binomial)[[3]]+5*coefficients(log_binomial)[[4]])
RR_period7=exp(coefficients(log_binomial)[[3]]+6*coefficients(log_binomial)[[4]])
RR_period8=exp(coefficients(log_binomial)[[3]]+7*coefficients(log_binomial)[[4]])
RR_period9=1 # assumed to be 1 thereafter

# read in LFT and PCR positive curves from MCMC simulations
LFT_curve=read_csv("LFT_curve.csv")%>%
  filter(days_since_infection%in%seq(0,30))
colnames(LFT_curve)<-c("days_since_infection","iteration","value")

# define function to generate probability weighted reductions in hospitalization & mortality for every possible
# combination of when an LFT is administered, under testing strategies at different frequencies
generate_EVs_sensitivity<-function(j_max,data,index,incubation_period,drug_coverage_prop,delay_test_to_treatment,
                       drug_efficacy_period1,drug_efficacy_period2,
                       drug_efficacy_period3,drug_efficacy_period4,drug_efficacy_period5,drug_efficacy_period6,
                       drug_efficacy_period7,drug_efficacy_period8,drug_efficacy_period9){

  list_product_all_probs<-c()
  list_product_all_probs2<-c()
  for (i in 1:index){
    for (j in rev(seq(0,j_max))){
      if (is.na(data[i+(index*j),"value"])==TRUE) {next}
      data_subset=data[seq(i,i+(index*j),index),]
      all_probs=(1-data_subset$value)
      
      if (data_subset$days_since_infection[nrow(data_subset)]<=0+ceiling(incubation_period)-delay_test_to_treatment){
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period1*drug_coverage_prop
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==1+ceiling(incubation_period)-delay_test_to_treatment) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period2*drug_coverage_prop
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==2+ceiling(incubation_period)-delay_test_to_treatment) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period3*drug_coverage_prop
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==3+ceiling(incubation_period)-delay_test_to_treatment ) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period4*drug_coverage_prop
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==4+ceiling(incubation_period)-delay_test_to_treatment ) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period5*drug_coverage_prop
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==5+ceiling(incubation_period)-delay_test_to_treatment ) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period6*drug_coverage_prop
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==6+ceiling(incubation_period)-delay_test_to_treatment ) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period7*drug_coverage_prop
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==7+ceiling(incubation_period)-delay_test_to_treatment ) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period8*drug_coverage_prop
      }
      else {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period9*drug_coverage_prop
      }
      product_all_probs=prod(all_probs)
      list_product_all_probs=c(list_product_all_probs,product_all_probs)
    }
    
    product_all_probs2=prod(1-data[seq(i,i+(index*j_max),index),"value"],na.rm=TRUE)
    list_product_all_probs2=c(list_product_all_probs2,product_all_probs2)
    
  }
   
  return(list(list1=list_product_all_probs,list2=list_product_all_probs2))
}
  
## produce results for efficacy preserved scenario 
produce_results_efficacy_preserved<-function(index,j_max,incubation_period_scenario){
  list_sums_EV<-c()
  for (iter in 1:max(LFT_curve$iteration)){
    if(incubation_period_scenario=='longer_incubation'){
      set.seed(111)
      incubation_period=rlnorm(1,meanlog=1.63,sdlog=0.5)
    }
    else{
      set.seed(111)
      incubation_period=rweibull(1,shape=1.5,scale=3.6)
    }
    LFT_curve_iteration=LFT_curve%>%
      filter(iteration==iter)
    EV_iteration=generate_EVs_sensitivity(j_max,LFT_curve_iteration,index,incubation_period=incubation_period,
                              drug_coverage_prop=1,delay_test_to_treatment=0,
                              drug_efficacy_period1=RR_period1,
                              drug_efficacy_period2=RR_period2,
                              drug_efficacy_period3=RR_period3,
                              drug_efficacy_period4=RR_period4,
                              drug_efficacy_period5=RR_period5,
                              drug_efficacy_period6=RR_period6,
                              drug_efficacy_period7=RR_period6,
                              drug_efficacy_period8=RR_period6,
                              drug_efficacy_period9=RR_period6)

    sum_EV_iteration=(1/index)*sum(EV_iteration$list1,EV_iteration$list2)
    list_sums_EV=c(list_sums_EV,sum_EV_iteration)

  }
  median = quantile(list_sums_EV, probs = 0.5)
  LQ = quantile(list_sums_EV, probs = 0.025)
  UQ = quantile(list_sums_EV, probs = 0.975)



  return(list(df=list_sums_EV,median=median,LQ=LQ,UQ=UQ))
}

EVs_every_other_day_efficacy_preserved=produce_results_efficacy_preserved(2,15,'longer_incubation')

EVs_every_three_days_efficacy_preserved=produce_results_efficacy_preserved(3,10,'longer_incubation')

EVs_every_week_efficacy_preserved=produce_results_efficacy_preserved(7,4,'longer_incubation')

EVs_every_two_weeks_efficacy_preserved=produce_results_efficacy_preserved(14,2,'longer_incubation')

all_EVs_efficacy_preserved=cbind.data.frame(every_other_day=EVs_every_other_day_efficacy_preserved$df,
                                            every_three_days=EVs_every_three_days_efficacy_preserved$df,
                                            every_week=EVs_every_week_efficacy_preserved$df,
                                            every_two_weeks=EVs_every_two_weeks_efficacy_preserved$df)
melted_all_EVs_efficacy_preserved=melt(all_EVs_efficacy_preserved,variable.name='testing_strategy',value.name='weighted_risk_reduction')
write.csv(melted_all_EVs_efficacy_preserved,"melted_all_EVs_efficacy_preserved.csv",row.names=F)


## produce results for fast decline to zero scenario

produce_results_drop_zero<-function(index,j_max,incubation_period_scenario){
  list_sums_EV<-c()
  list_prop_given_drug<-c()
  list_prop_benefited<-c()
  list_prop_given_drug_and_benefited<-c()

  for (iter in 1:max(LFT_curve$iteration)){
    if(incubation_period_scenario=='longer_incubation'){
      incubation_period=rlnorm(1,meanlog=1.63,sdlog=0.5)
    }
    else{
      incubation_period=rweibull(1,shape=1.5,scale=3.6)
    }
    LFT_curve_iteration=LFT_curve%>%
      filter(iteration==iter)
    EV_iteration=generate_EVs_sensitivity(j_max,LFT_curve_iteration,index,incubation_period=incubation_period,
                              drug_coverage_prop=1,delay_test_to_treatment=0,
                              drug_efficacy_period1=RR_period1,
                              drug_efficacy_period2=RR_period2,
                              drug_efficacy_period3=RR_period3,
                              drug_efficacy_period4=RR_period4,
                              drug_efficacy_period5=RR_period5,
                              drug_efficacy_period6=RR_period6,
                              drug_efficacy_period7=1,
                              drug_efficacy_period8=1,
                              drug_efficacy_period9=1)

    sum_EV_iteration=(1/index)*sum(EV_iteration$list1,EV_iteration$list2)
    list_sums_EV=c(list_sums_EV,sum_EV_iteration)
  }
  median = quantile(list_sums_EV, probs = 0.5)
  LQ = quantile(list_sums_EV, probs = 0.025)
  UQ = quantile(list_sums_EV, probs = 0.975)

  return(list(df=list_sums_EV,median=median,LQ=LQ,UQ=UQ))
}

EVs_every_other_day_drop_zero=produce_results_drop_zero(2,15,"longer_incubation")

EVs_every_three_days_drop_zero=produce_results_drop_zero(3,10,"longer_incubation")

EVs_every_week_drop_zero=produce_results_drop_zero(7,4,"longer_incubation")

EVs_every_two_weeks_drop_zero=produce_results_drop_zero(14,2,"longer_incubation")

all_EVs_drop_zero=cbind.data.frame(every_other_day=EVs_every_other_day_drop_zero$df,
                                            every_three_days=EVs_every_three_days_drop_zero$df,
                                            every_week=EVs_every_week_drop_zero$df,
                                            every_two_weeks=EVs_every_two_weeks_drop_zero$df)
melted_all_EVs_drop_zero=melt(all_EVs_drop_zero,variable.name='testing_strategy',value.name='weighted_risk_reduction')
write.csv(melted_all_EVs_drop_zero,"melted_all_EVs_drop_zero.csv",row.names=F)

## re-run above analyses w/ omicron-specific incubation period

EVs_every_other_day_efficacy_preserved_incubation_sensitivity=produce_results_efficacy_preserved(2,15,'shorter_incubation')

EVs_every_three_days_efficacy_preserved_incubation_sensitivity=produce_results_efficacy_preserved(3,10,'shorter_incubation')

EVs_every_week_efficacy_preserved_incubation_sensitivity=produce_results_efficacy_preserved(7,4,'shorter_incubation')

EVs_every_two_weeks_efficacy_preserved_incubation_sensitivity=produce_results_efficacy_preserved(14,2,'shorter_incubation')

all_EVs_efficacy_preserved_incubation_sensitivity=cbind.data.frame(every_other_day=EVs_every_other_day_efficacy_preserved_incubation_sensitivity$df,
                                            every_three_days=EVs_every_three_days_efficacy_preserved_incubation_sensitivity$df,
                                            every_week=EVs_every_week_efficacy_preserved_incubation_sensitivity$df,
                                            every_two_weeks=EVs_every_two_weeks_efficacy_preserved_incubation_sensitivity$df)
melted_all_EVs_efficacy_preserved_incubation_sensitivity=melt(all_EVs_efficacy_preserved_incubation_sensitivity,variable.name='testing_strategy',value.name='weighted_risk_reduction')
write.csv(melted_all_EVs_efficacy_preserved_incubation_sensitivity,"melted_all_EVs_efficacy_preserved_incubation_sensitivity.csv",row.names=F)

EVs_every_other_day_drop_zero_incubation_sensitivity=produce_results_drop_zero(2,15,"shorter_incubation")

EVs_every_three_days_drop_zero_incubation_sensitivity=produce_results_drop_zero(3,10,"shorter_incubation")

EVs_every_week_drop_zero_incubation_sensitivity=produce_results_drop_zero(7,4,"shorter_incubation")

EVs_every_two_weeks_drop_zero_incubation_sensitivity=produce_results_drop_zero(14,2,"shorter_incubation")

all_EVs_drop_zero_incubation_sensitivity=cbind.data.frame(every_other_day=EVs_every_other_day_drop_zero_incubation_sensitivity$df,
                                   every_three_days=EVs_every_three_days_drop_zero_incubation_sensitivity$df,
                                   every_week=EVs_every_week_drop_zero_incubation_sensitivity$df,
                                   every_two_weeks=EVs_every_two_weeks_drop_zero_incubation_sensitivity$df)
melted_all_EVs_drop_zero_incubation_sensitivity=melt(all_EVs_drop_zero_incubation_sensitivity,variable.name='testing_strategy',value.name='weighted_risk_reduction')
write.csv(melted_all_EVs_drop_zero_incubation_sensitivity,"melted_all_EVs_drop_zero_incubation_sensitivity.csv",row.names=F)

