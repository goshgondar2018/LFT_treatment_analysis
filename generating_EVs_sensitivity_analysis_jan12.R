library(tidyverse)
library(reshape2)


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
  

produce_results_efficacy_preserved<-function(index,j_max){
  list_sums_EV<-c()
  list_prop_given_drug<-c()
  list_prop_benefited<-c()
  list_prop_given_drug_and_benefited<-c()

  for (iter in 1:max(LFT_curve$iteration)){
    incubation_period=rlnorm(1,meanlog=1.63,sdlog=0.5)
    LFT_curve_iteration=LFT_curve%>%
      filter(iteration==iter)
    EV_iteration=generate_EVs_sensitivity(j_max,LFT_curve_iteration,index,incubation_period=incubation_period,
                              drug_coverage_prop=1,delay_test_to_treatment=1,
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

EVs_every_other_day_efficacy_preserved=produce_results_efficacy_preserved(2,15)
#write.csv(EVs_every_other_day_efficacy_preserved$df,"EVs_every_other_day_df_efficacy_preserved.csv")

EVs_every_three_days_efficacy_preserved=produce_results_efficacy_preserved(3,10)
#write.csv(EVs_every_three_days_efficacy_preserved$df,"EVs_every_three_days_df_efficacy_preserved.csv")

EVs_every_week_efficacy_preserved=produce_results_efficacy_preserved(7,4)
#write.csv(EVs_every_week_efficacy_preserved$df,"EVs_every_week_df_efficacy_preserved.csv")

EVs_every_two_weeks_efficacy_preserved=produce_results_efficacy_preserved(14,2)
#write.csv(EVs_every_two_weeks_efficacy_preserved$df,"EVs_every_two_weeks_df_efficacy_preserved.csv")

all_EVs_efficacy_preserved=cbind.data.frame(every_other_day=EVs_every_other_day_efficacy_preserved$df,
                                            every_three_days=EVs_every_three_days_efficacy_preserved$df,
                                            every_week=EVs_every_week_efficacy_preserved$df,
                                            every_two_weeks=EVs_every_two_weeks_efficacy_preserved$df)
melted_all_EVs_efficacy_preserved=melt(all_EVs_efficacy_preserved,variable.name='testing_strategy',value.name='weighted_risk_reduction')
write.csv(melted_all_EVs_efficacy_preserved,"melted_all_EVs_efficacy_preserved.csv",row.names=F)



#####

produce_results_drop_zero<-function(index,j_max){
  list_sums_EV<-c()
  list_prop_given_drug<-c()
  list_prop_benefited<-c()
  list_prop_given_drug_and_benefited<-c()

  for (iter in 1:max(LFT_curve$iteration)){
    incubation_period=rlnorm(1,meanlog=1.63,sdlog=0.5)
    LFT_curve_iteration=LFT_curve%>%
      filter(iteration==iter)
    EV_iteration=generate_EVs_sensitivity(j_max,LFT_curve_iteration,index,incubation_period=incubation_period,
                              drug_coverage_prop=1,delay_test_to_treatment=1,
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

EVs_every_other_day_drop_zero=produce_results_drop_zero(2,15)
#write.csv(EVs_every_other_day_drop_zero$df,"EVs_every_other_day_df_drop_zero.csv")

EVs_every_three_days_drop_zero=produce_results_drop_zero(3,10)
#write.csv(EVs_every_three_days_drop_zero$df,"EVs_every_three_days_df_drop_zero.csv")

EVs_every_week_drop_zero=produce_results_drop_zero(7,4)
#write.csv(EVs_every_week_drop_zero$df,"EVs_every_week_df_drop_zero.csv")

EVs_every_two_weeks_drop_zero=produce_results_drop_zero(14,2)
#write.csv(EVs_every_two_weeks_drop_zero$df,"EVs_every_two_weeks_df_drop_zero.csv")

all_EVs_drop_zero=cbind.data.frame(every_other_day=EVs_every_other_day_drop_zero$df,
                                            every_three_days=EVs_every_three_days_drop_zero$df,
                                            every_week=EVs_every_week_drop_zero$df,
                                            every_two_weeks=EVs_every_two_weeks_drop_zero$df)
melted_all_EVs_drop_zero=melt(all_EVs_drop_zero,variable.name='testing_strategy',value.name='weighted_risk_reduction')
write.csv(melted_all_EVs_drop_zero,"melted_all_EVs_drop_zero.csv",row.names=F)


