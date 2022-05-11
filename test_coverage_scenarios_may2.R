library(tidyverse)
library(reshape2)

# read in LFT and PCR positive curves from MCMC simulations
LFT_curve=read_csv("LFT_curve.csv")%>%
  filter(days_since_infection%in%seq(0,30))
colnames(LFT_curve)<-c("days_since_infection","iteration","value")

# define function to generate proportion offered drug & prop benefited for every possible
# combination of when an LFT is administered, under testing strategies at different frequencies
generate_EVs_test_covg_sensitivity<-function(j_max,data,index,incubation_period,drug_coverage_prop,
                       delay_test_to_treatment){
  
  prod_probs_all_negative_list <- c()
  list_product_benefited<-c()
  list_product_benefited2<-c()
  
  for (i in 1:index){
    for (j in rev(seq(0,j_max))){
      if (is.na(data[i+(index*j),"value"])==TRUE) {next}
      data_subset=data[seq(i,i+(index*j),index),]
      all_probs=(1-data_subset$value)
      benefited=all_probs
      if (data_subset$days_since_infection[nrow(data_subset)]<=0+ceiling(incubation_period)-delay_test_to_treatment){
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==1+ceiling(incubation_period)-delay_test_to_treatment) {
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==2+ceiling(incubation_period)-delay_test_to_treatment) {
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==3+ceiling(incubation_period)-delay_test_to_treatment ) {
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==4+ceiling(incubation_period)-delay_test_to_treatment ) {
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==5+ceiling(incubation_period)-delay_test_to_treatment ) {
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==6+ceiling(incubation_period)-delay_test_to_treatment ) {
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==7+ceiling(incubation_period)-delay_test_to_treatment ) {
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else {
        benefited[length(benefited)]=0
      }
      
      product_benefited=prod(benefited)
      list_product_benefited=c(list_product_benefited,product_benefited)
    }
    data_subset2=data[seq(i,i+(index*j_max),index),]
    data_subset2=data_subset2[complete.cases(data_subset2$value),]
    probs_all_negative=1-data_subset2[,"value"]
    prod_probs_all_negative=prod(probs_all_negative)
    prod_probs_all_negative_list=c(prod_probs_all_negative_list,prod_probs_all_negative)
    
    product_benefited2=0
    list_product_benefited2=c(list_product_benefited2,product_benefited2)
    
  }
  prop_given_drug_final_value=1-mean(prod_probs_all_negative_list)
  
  return(list(prop=prop_given_drug_final_value,prop2=list_product_benefited,prop3=list_product_benefited2))
}
# discrete hi/lo test covg and frequency scenarios

produce_results_test_covg_sensitivity<-function(index,j_max,test_coverage_prop,incubation_period_scenario){

  list_sums_EV<-c()
  list_prop_given_drug<-c()
  list_prop_benefited<-c()
  list_prop_given_drug_and_benefited<-c()

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
    EV_iteration=generate_EVs_test_covg_sensitivity(j_max,LFT_curve_iteration,index,
                                                    incubation_period=incubation_period,
                              drug_coverage_prop=1,
                              delay_test_to_treatment=0)

    prop_given_drug=EV_iteration$prop*test_coverage_prop
    list_prop_given_drug=c(list_prop_given_drug,prop_given_drug)

    prop_benefited=(1/index)*sum(EV_iteration$prop2,EV_iteration$prop3)*test_coverage_prop
    list_prop_benefited=c(list_prop_benefited,prop_benefited)

  }

  median_given_drug_probs=quantile(list_prop_given_drug,probs=0.5)
  LQ_given_drug_probs=quantile(list_prop_given_drug,probs=0.025)
  UQ_given_drug_probs=quantile(list_prop_given_drug,probs=0.975)

  median_prop_benefited=quantile(list_prop_benefited,probs=0.5)
  LQ_prop_benefited=quantile(list_prop_benefited,probs=0.025)
  UQ_prop_benefited=quantile(list_prop_benefited,probs=0.975)


  return(list(list_prop_given_drug=list_prop_given_drug,
              list_prop_benefited=list_prop_benefited,
              median_given_drug_probs=median_given_drug_probs,LQ_given_drug_probs=LQ_given_drug_probs,
              UQ_given_drug_probs=UQ_given_drug_probs,
              median_prop_benefited=median_prop_benefited,
              LQ_prop_benefited=LQ_prop_benefited,UQ_prop_benefited=UQ_prop_benefited)) # median_given_drug_and_benefited=median_given_drug_and_benefited,
}

generate_benefited_onetest<-function(days_symptoms_start,
                                     data){
  i = days_symptoms_start 
  data_subset=data[i+1,]
  benefited_or_offered=data_subset$value
  return(benefited_or_offered)
}

produce_results_onetest<-function(test_coverage_prop,incubation_period_scenario){
  list_benefited_or_offered<-c()
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
    incubation_period_rounded=ceiling(incubation_period)
    benefited_or_offered=generate_benefited_onetest(days_symptoms_start=incubation_period_rounded,
                                                    LFT_curve_iteration)*test_coverage_prop
    list_benefited_or_offered=c(list_benefited_or_offered,benefited_or_offered)
  }
  median = quantile(list_benefited_or_offered, probs = 0.5)
  LQ = quantile(list_benefited_or_offered, probs = 0.025)
  UQ = quantile(list_benefited_or_offered, probs = 0.975)
  return(list(list_benefited_or_offered=list_benefited_or_offered,median=median,LQ=LQ,UQ=UQ))
}

# output_every_other_day_test_covg_low=produce_results_test_covg_sensitivity(2,15,0.25,'longer_incubation')
# write.csv(output_every_other_day_test_covg_low$list_prop_given_drug,"given_drug_every_other_day_test_covg_low.csv")
# write.csv(output_every_other_day_test_covg_low$list_prop_benefited,"benefit_every_other_day_test_covg_low.csv")
# 
# output_every_three_days_test_covg_low=produce_results_test_covg_sensitivity(3,10,0.25,'longer_incubation')
# write.csv(output_every_three_days_test_covg_low$list_prop_given_drug,"given_drug_every_three_days_test_covg_low.csv")
# write.csv(output_every_three_days_test_covg_low$list_prop_benefited,"benefit_every_three_days_test_covg_low.csv")
# 
# output_every_week_test_covg_high=produce_results_test_covg_sensitivity(7,4,0.75,'longer_incubation')
# write.csv(output_every_week_test_covg_high$list_prop_given_drug,"given_drug_every_week_test_covg_high.csv")
# write.csv(output_every_week_test_covg_high$list_prop_benefited,"benefit_every_week_test_covg_high.csv")
# 
# output_every_two_weeks_test_covg_high=produce_results_test_covg_sensitivity(14,2,0.75,'longer_incubation')
# write.csv(output_every_two_weeks_test_covg_high$list_prop_given_drug,"given_drug_every_two_weeks_test_covg_high.csv")
# write.csv(output_every_two_weeks_test_covg_high$list_prop_benefited,"benefit_every_two_weeks_test_covg_high.csv")
# 
# output_every_other_day_test_covg_low_incubation_sensitivity=produce_results_test_covg_sensitivity(2,15,0.25,'shorter_incubation')
# write.csv(output_every_other_day_test_covg_low_incubation_sensitivity$list_prop_given_drug,"given_drug_every_other_day_test_covg_low_incubation_sensitivity.csv")
# write.csv(output_every_other_day_test_covg_low_incubation_sensitivity$list_prop_benefited,"benefit_every_other_day_test_covg_low_incubation_sensitivity.csv")
# 
# output_every_three_days_test_covg_low_incubation_sensitivity=produce_results_test_covg_sensitivity(3,10,0.25,'shorter_incubation')
# write.csv(output_every_three_days_test_covg_low_incubation_sensitivity$list_prop_given_drug,"given_drug_every_three_days_test_covg_low_incubation_sensitivity.csv")
# write.csv(output_every_three_days_test_covg_low_incubation_sensitivity$list_prop_benefited,"benefit_every_three_days_test_covg_low_incubation_sensitivity.csv")
# 
# output_every_week_test_covg_high_incubation_sensitivity=produce_results_test_covg_sensitivity(7,4,0.75,'shorter_incubation')
# write.csv(output_every_week_test_covg_high_incubation_sensitivity$list_prop_given_drug,"given_drug_every_week_test_covg_high_incubation_sensitivity.csv")
# write.csv(output_every_week_test_covg_high_incubation_sensitivity$list_prop_benefited,"benefit_every_week_test_covg_high_incubation_sensitivity.csv")
# 
# output_every_two_weeks_test_covg_high_incubation_sensitivity=produce_results_test_covg_sensitivity(14,2,0.75,'shorter_incubation')
# write.csv(output_every_two_weeks_test_covg_high_incubation_sensitivity$list_prop_given_drug,"given_drug_every_two_weeks_test_covg_high_incubation_sensitivity.csv")
# write.csv(output_every_two_weeks_test_covg_high_incubation_sensitivity$list_prop_benefited,"benefit_every_two_weeks_test_covg_high_incubation_sensitivity.csv")
# 
# output_onetest_test_covg_low=produce_results_onetest(0.25,'longer_incubation')
# write.csv(output_onetest_test_covg_low$list_benefited_or_offered,"benefit_onetest_test_covg_low.csv")
# 
# output_onetest_test_covg_hi=produce_results_onetest(0.75,'longer_incubation')
# write.csv(output_onetest_test_covg_hi$list_benefited_or_offered,"benefit_onetest_test_covg_hi.csv")
# 
# output_onetest_test_covg_low_incubation_sensitivity=produce_results_onetest(0.25,'shorter_incubation')
# write.csv(output_onetest_test_covg_low_incubation_sensitivity$list_benefited_or_offered,"benefit_onetest_test_covg_low_incubation_sensitivity.csv")
# 
# output_onetest_test_covg_hi_incubation_sensitivity=produce_results_onetest(0.75,'shorter_incubation')
# write.csv(output_onetest_test_covg_hi_incubation_sensitivity$list_benefited_or_offered,"benefit_onetest_test_covg_hi_incubation_sensitivity.csv")

# 'continuous' test covg scenarios
test_covg_seq=seq(0,1,length.out=5)

# every other day testing strategy
all_props_every_other_day=data.frame(matrix(ncol=4,nrow=0,dimnames=list(NULL,c("test_covg_prop","iteration","prop_given_drug","prop_benefited"))))
for (test_covg in test_covg_seq){
  EVs_every_other_day_test_i=produce_results_test_covg_sensitivity(2,15,test_covg,'longer_incubation')
  EVs_every_other_day_all_output_i=cbind.data.frame(test_covg_prop=rep(test_covg,4000),iteration=seq(1,4000),
                                                    prop_given_drug=EVs_every_other_day_test_i$list_prop_given_drug,
                                                    prop_benefited=EVs_every_other_day_test_i$list_prop_benefited)
  all_props_every_other_day=rbind.data.frame(all_props_every_other_day,EVs_every_other_day_all_output_i)
}
all_props_every_other_day$testing_strategy="every_other_day"

# every three days testing strategy
all_props_every_three_days=data.frame(matrix(ncol=4,nrow=0,dimnames=list(NULL,c("test_covg_prop","iteration","prop_given_drug","prop_benefited"))))

for (test_covg in test_covg_seq){
  EVs_every_three_days_test_i=produce_results_test_covg_sensitivity(3,10,test_covg,'longer_incubation')
  EVs_every_three_days_all_output_i=cbind.data.frame(test_covg_prop=rep(test_covg,4000),iteration=seq(1,4000),
                                                     prop_given_drug=EVs_every_three_days_test_i$list_prop_given_drug,
                                                     prop_benefited=EVs_every_three_days_test_i$list_prop_benefited)
  all_props_every_three_days=rbind.data.frame(all_props_every_three_days,EVs_every_three_days_all_output_i)
}
all_props_every_three_days$testing_strategy="every_three_days"

# every week testing strategy
all_props_every_week=data.frame(matrix(ncol=4,nrow=0,dimnames=list(NULL,c("test_covg_prop","iteration","prop_given_drug","prop_benefited"))))

for (test_covg in test_covg_seq){
  EVs_every_week_test_i=produce_results_test_covg_sensitivity(7,4,test_covg,'longer_incubation')
  EVs_every_week_all_output_i=cbind.data.frame(test_covg_prop=rep(test_covg,4000),iteration=seq(1,4000),
                                               prop_given_drug=EVs_every_week_test_i$list_prop_given_drug,
                                               prop_benefited=EVs_every_week_test_i$list_prop_benefited)
  all_props_every_week=rbind.data.frame(all_props_every_week,EVs_every_week_all_output_i)
}
all_props_every_week$testing_strategy="every_week"

# every two weeks testing strategy
all_props_every_two_weeks=data.frame(matrix(ncol=4,nrow=0,dimnames=list(NULL,c("test_covg_prop","iteration","prop_given_drug","prop_benefited"))))

for (test_covg in test_covg_seq){
  EVs_every_two_weeks_test_i=produce_results_test_covg_sensitivity(14,2,test_covg,'longer_incubation')
  EVs_every_two_weeks_all_output_i=cbind.data.frame(test_covg_prop=rep(test_covg,4000),iteration=seq(1,4000),
                                                    prop_given_drug=EVs_every_two_weeks_test_i$list_prop_given_drug,
                                                    prop_benefited=EVs_every_two_weeks_test_i$list_prop_benefited)
  all_props_every_two_weeks=rbind.data.frame(all_props_every_two_weeks,EVs_every_two_weeks_all_output_i)
}
all_props_every_two_weeks$testing_strategy="every_two_weeks"

#one-time testing strategy
all_props_onetest=data.frame(matrix(ncol=4,nrow=0,dimnames=list(NULL,c("test_covg_prop","iteration","prop_given_drug","prop_benefited"))))

for (test_covg in test_covg_seq){
  EVs_onetest_i=produce_results_onetest(test_covg,'longer_incubation')
  EVs_onetest_all_output_i=cbind.data.frame(test_covg_prop=rep(test_covg,4000),iteration=seq(1,4000),
                                            prop_given_drug=EVs_onetest_i$list_benefited_or_offered,
                                            prop_benefited=EVs_onetest_i$list_benefited_or_offered)
  all_props_onetest=rbind.data.frame(all_props_onetest,EVs_onetest_all_output_i)
}
all_props_onetest$testing_strategy="one-time testing"


all_props_every_strategy=rbind.data.frame(all_props_every_other_day,all_props_every_three_days,
                                          all_props_every_week,all_props_every_two_weeks,all_props_onetest)
write.csv(all_props_every_strategy,"all_props_every_strategy.csv")

## repeat for omicron-specific incubation period 

# every other day testing strategy
# all_props_every_other_day_incubation_sensitivity=data.frame()
# for (test_covg in test_covg_seq){
#   EVs_every_other_day_test_i=produce_results_test_covg_sensitivity(2,15,test_covg,'shorter_incubation')
#   EVs_every_other_day_all_output_i=cbind.data.frame(test_covg_prop=rep(test_covg,400),iteration=seq(1,400),
#                                                     prop_given_drug=EVs_every_other_day_test_i$list_prop_given_drug,
#                                                     prop_benefited=EVs_every_other_day_test_i$list_prop_benefited)
#   all_props_every_other_day_incubation_sensitivity=rbind(all_props_every_other_day_incubation_sensitivity,EVs_every_other_day_all_output_i)
# }
# all_props_every_other_day_incubation_sensitivity$testing_strategy="every_other_day"
# 
# # every three days testing strategy
# all_props_every_three_days_incubation_sensitivity=data.frame()
# for (test_covg in test_covg_seq){
#   EVs_every_three_days_test_i=produce_results_test_covg_sensitivity(3,10,test_covg,'shorter_incubation')
#   EVs_every_three_days_all_output_i=cbind.data.frame(rep(test_covg,400),seq(1,400),
#                                                      EVs_every_three_days_test_i$list_prop_given_drug,
#                                                      EVs_every_three_days_test_i$list_prop_benefited)
#   all_props_every_three_days_incubation_sensitivity=rbind(all_props_every_three_days_incubation_sensitivity,EVs_every_three_days_all_output_i)
# }
# all_props_every_three_days_incubation_sensitivity$testing_strategy="every_three_days"
# 
# # every week testing strategy
# all_props_every_week_incubation_sensitivity=data.frame()
# for (test_covg in test_covg_seq){
#   EVs_every_week_test_i=produce_results_test_covg_sensitivity(7,4,test_covg,'shorter_incubation')
#   EVs_every_week_all_output_i=cbind.data.frame(rep(test_covg,400),seq(1,400),
#                                                EVs_every_week_test_i$list_prop_given_drug,
#                                                EVs_every_week_test_i$list_prop_benefited)
#   all_props_every_week_incubation_sensitivity=rbind(all_props_every_week_incubation_sensitivity,EVs_every_week_all_output_i)
# }
# all_props_every_week_incubation_sensitivity$testing_strategy="every_week"
# 
# # every two weeks testing strategy
# all_props_every_two_weeks_incubation_sensitivity=data.frame()
# for (test_covg in test_covg_seq){
#   EVs_every_two_weeks_test_i=produce_results_test_covg_sensitivity(14,2,test_covg,'shoter_incubation')
#   EVs_every_two_weeks_all_output_i=cbind.data.frame(rep(test_covg,400),seq(1,400),
#                                                     EVs_every_two_weeks_test_i$list_prop_given_drug,
#                                                     EVs_every_two_weeks_test_i$list_prop_benefited)
#   all_props_every_two_weeks_incubation_sensitivity=rbind(all_props_every_two_weeks_incubation_sensitivity,EVs_every_two_weeks_all_output_i)
# }
# all_props_every_two_weeks_incubation_sensitivity$testing_strategy="every_two_weeks"
# 
# #one-time testing strategy
# all_props_onetime_test_incubation_sensitivity=data.frame()
# for (test_covg in test_covg_seq){
#   EVs_onetest_i=produce_results_onetest(test_covg,'shorter_incubation')
#   EVs_onetest_all_output_i=cbind.data.frame(rep(test_covg,400),seq(1,400),
#                                             EVs_onetest_i$list_prop_given_drug,
#                                             EVs_onetest_weeks_test_i$list_prop_benefited)
#   all_props_onetest_incubation_sensitivity=rbind(all_props_onetest_incubation_sensitivity,EVs_onetest_all_output_i)
# }
# all_props_onetest_incubation_sensitivity$testing_strategy="one-time testing"
# 
# 
# all_props_every_strategy_incubation_sensitivity=rbind.data.frame(all_props_every_other_day_incubation_sensitivity,
#                                                                  all_props_every_three_days_incubation_sensitivity,
#                                           all_props_every_week_incubation_sensitivity,
#                                           all_props_every_two_weeks_incubation_sensitivity,all_props_onetest_incubation_sensitivity)
# write.csv(all_props_every_strategy_incubation_sensitivity,"all_props_every_strategy_incubation_sensitivity.csv")
# 

