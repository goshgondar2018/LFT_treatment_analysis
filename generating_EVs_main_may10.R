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
generate_EVs<-function(j_max,data,index,incubation_period,drug_coverage_prop,
                       delay_test_to_treatment,
                       drug_efficacy_period1,drug_efficacy_period2,
                       drug_efficacy_period3,drug_efficacy_period4,drug_efficacy_period5,drug_efficacy_period6,
                       drug_efficacy_period7,drug_efficacy_period8,drug_efficacy_period9){
  
  list_product_all_probs<-c()
  list_product_all_probs2<-c()
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
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period1*drug_coverage_prop
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==1+ceiling(incubation_period)-delay_test_to_treatment) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period2*drug_coverage_prop
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==2+ceiling(incubation_period)-delay_test_to_treatment) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period3*drug_coverage_prop
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==3+ceiling(incubation_period)-delay_test_to_treatment ) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period4*drug_coverage_prop
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==4+ceiling(incubation_period)-delay_test_to_treatment ) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period5*drug_coverage_prop
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==5+ceiling(incubation_period)-delay_test_to_treatment ) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period6*drug_coverage_prop
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==6+ceiling(incubation_period)-delay_test_to_treatment ) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period7*drug_coverage_prop
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else if (data_subset$days_since_infection[nrow(data_subset)]==7+ceiling(incubation_period)-delay_test_to_treatment ) {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period8*drug_coverage_prop
        benefited[length(benefited)]=(1-benefited[length(benefited)])
      }
      else {
        all_probs[length(all_probs)]=(1-all_probs[length(all_probs)])*drug_efficacy_period9*drug_coverage_prop
        benefited[length(benefited)]=0
      }
      product_all_probs=prod(all_probs)
      list_product_all_probs=c(list_product_all_probs,product_all_probs)
      
      product_benefited=prod(benefited)
      list_product_benefited=c(list_product_benefited,product_benefited)
    }
    data_subset2=data[seq(i,i+(index*j_max),index),]
    data_subset2=data_subset2[complete.cases(data_subset2$value),]
    probs_all_negative=1-data_subset2[,"value"]
    prod_probs_all_negative=prod(probs_all_negative)
    prod_probs_all_negative_list=c(prod_probs_all_negative_list,prod_probs_all_negative)
    
    product_all_probs2=prod(1-data[seq(i,i+(index*j_max),index),"value"],na.rm=TRUE)
    list_product_all_probs2=c(list_product_all_probs2,product_all_probs2)
    
    product_benefited2=0
    list_product_benefited2=c(list_product_benefited2,product_benefited2)
    
  }
  prop_given_drug_final_value=1-mean(prod_probs_all_negative_list)
  
  
  return(list(list1=list_product_all_probs,list2=list_product_all_probs2,
              prop=prop_given_drug_final_value,prop2=list_product_benefited,prop3=list_product_benefited2))
}

produce_results_main<-function(index,j_max,test_coverage_prop,incubation_period_scenario){
  
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
    EV_iteration=generate_EVs(j_max,LFT_curve_iteration,index,incubation_period=incubation_period,
                              drug_coverage_prop=1,delay_test_to_treatment=0,
                              drug_efficacy_period1=RR_period1,
                              drug_efficacy_period2=RR_period2,
                              drug_efficacy_period3=RR_period3,
                              drug_efficacy_period4=RR_period4,
                              drug_efficacy_period5=RR_period5,
                              drug_efficacy_period6=RR_period6,
                              drug_efficacy_period7=RR_period7,
                              drug_efficacy_period8=RR_period8,
                              drug_efficacy_period9=RR_period9)
    
    sum_EV_iteration=(1/index)*sum(EV_iteration$list1,EV_iteration$list2)
    list_sums_EV=c(list_sums_EV,sum_EV_iteration)
    
    prop_given_drug=EV_iteration$prop*test_coverage_prop
    list_prop_given_drug=c(list_prop_given_drug,prop_given_drug)
    
    prop_benefited=(1/index)*sum(EV_iteration$prop2,EV_iteration$prop3)*test_coverage_prop
    list_prop_benefited=c(list_prop_benefited,prop_benefited)
    
    
  }
  median = quantile(list_sums_EV, probs = 0.5)
  LQ = quantile(list_sums_EV, probs = 0.025)
  UQ = quantile(list_sums_EV, probs = 0.975)
  
  median_given_drug_probs=quantile(list_prop_given_drug,probs=0.5)
  LQ_given_drug_probs=quantile(list_prop_given_drug,probs=0.025)
  UQ_given_drug_probs=quantile(list_prop_given_drug,probs=0.975)
  
  median_prop_benefited=quantile(list_prop_benefited,probs=0.5)
  LQ_prop_benefited=quantile(list_prop_benefited,probs=0.025)
  UQ_prop_benefited=quantile(list_prop_benefited,probs=0.975)
  
  
  return(list(df=list_sums_EV,median=median,LQ=LQ,UQ=UQ,
              median_given_drug_probs=median_given_drug_probs,LQ_given_drug_probs=LQ_given_drug_probs,
              UQ_given_drug_probs=UQ_given_drug_probs,
              median_prop_benefited=median_prop_benefited,
              LQ_prop_benefited=LQ_prop_benefited,UQ_prop_benefited=UQ_prop_benefited)) # median_given_drug_and_benefited=median_given_drug_and_benefited,
}

EVs_every_other_day=produce_results_main(2,15,1,'longer_incubation')

EVs_every_three_days=produce_results_main(3,10,1,'longer_incubation')

EVs_every_week=produce_results_main(7,4,1,'longer_incubation')

EVs_every_two_weeks=produce_results_main(14,2,1,'longer_incubation')

all_EVs=cbind.data.frame(every_other_day=EVs_every_other_day$df,
                         every_three_days=EVs_every_three_days$df,
                         every_week=EVs_every_week$df,
                         every_two_weeks=EVs_every_two_weeks$df)
melted_all_EVs=melt(all_EVs,variable.name='testing_strategy',value.name='weighted_risk_reduction')
write.csv(melted_all_EVs,"melted_all_EVs.csv",row.names=F)

write.csv(cbind.data.frame(matrix(c(EVs_every_other_day$median_prop_benefited,EVs_every_other_day$LQ_prop_benefited,
                                    EVs_every_other_day$UQ_prop_benefited,EVs_every_three_days$median_prop_benefited,EVs_every_three_days$LQ_prop_benefited,
                                    EVs_every_three_days$UQ_prop_benefited,EVs_every_week$median_prop_benefited,EVs_every_week$LQ_prop_benefited,
                                    EVs_every_week$UQ_prop_benefited,EVs_every_two_weeks$median_prop_benefited,EVs_every_two_weeks$LQ_prop_benefited,
                                    EVs_every_two_weeks$UQ_prop_benefited),byrow=TRUE,ncol=3,nrow=4)),
          "prop_benefited_all.csv",row.names=F)

write.csv(cbind.data.frame(matrix(c(EVs_every_other_day$median_given_drug_probs,EVs_every_other_day$LQ_given_drug_probs,
                                    EVs_every_other_day$UQ_given_drug_probs,EVs_every_three_days$median_given_drug_probs,EVs_every_three_days$LQ_given_drug_probs,
                                    EVs_every_three_days$UQ_given_drug_probs,EVs_every_week$median_given_drug_probs,EVs_every_week$LQ_given_drug_probs,
                                    EVs_every_week$UQ_given_drug_probs,EVs_every_two_weeks$median_given_drug_probs,EVs_every_two_weeks$LQ_given_drug_probs,
                                    EVs_every_two_weeks$UQ_given_drug_probs),byrow=TRUE,ncol=3,nrow=4)),
          "prop_given_drug_probs_all.csv",row.names=F)

generate_benefited_onetest<-function(days_symptoms_start,drug_efficacy_period1,drug_coverage_prop,
                                     data){
  i = days_symptoms_start 
  data_subset=data[i+1,]
  
  all_probs=(data_subset$value*drug_efficacy_period1*drug_coverage_prop)+(1-data_subset$value)*1
  benefited_or_offered=data_subset$value
  return(list(RR=all_probs,benefited_or_offered=benefited_or_offered))
}

produce_results_onetest<-function(test_coverage_prop,incubation_period_scenario){
  list_benefited_or_offered<-c()
  list_weighted_RR<-c()
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
    weighted_RR=generate_benefited_onetest(days_symptoms_start=incubation_period_rounded,
                                                     drug_efficacy_period1 = RR_period1,drug_coverage_prop = 1,
                                                     data=LFT_curve_iteration)$RR
    list_weighted_RR=c(list_weighted_RR,weighted_RR)
    benefited_or_offered=generate_benefited_onetest(days_symptoms_start=incubation_period_rounded,
                                                    drug_efficacy_period1 = RR_period1,drug_coverage_prop = 1,
                                                    data=LFT_curve_iteration)$benefited_or_offered*test_coverage_prop
    list_benefited_or_offered=c(list_benefited_or_offered,benefited_or_offered)
  }
  median = quantile(list_benefited_or_offered, probs = 0.5)
  LQ = quantile(list_benefited_or_offered, probs = 0.025)
  UQ = quantile(list_benefited_or_offered, probs = 0.975)
  
  median_weighted_RR = quantile(list_weighted_RR, probs = 0.5)
  LQ_weighted_RR = quantile(list_weighted_RR, probs = 0.025)
  UQ_weighted_RR = quantile(list_weighted_RR, probs = 0.975)
  return(list(list_benefited_or_offered=list_benefited_or_offered,median_benefited=median,LQ_benefited=LQ,UQ_benefited=UQ,
              list_weighted_RR=list_weighted_RR,median_weighted_RR=median_weighted_RR,LQ_weighted_RR=LQ_weighted_RR,UQ_weighted_RR=UQ_weighted_RR))
}

onetest_output=produce_results_onetest(1,"longer_incubation")
write.csv(onetest_output$list_weighted_RR,"onetest_output_RR.csv",row.names=F)
write.csv(onetest_output$list_benefited_or_offered,"onetest_output_benefited_or_offered.csv",row.names=F)
## re-run above analyses with omicron-specific incubation period

EVs_every_other_day_incubation_sensitivity=produce_results_main(2,15,1,"shorter_incubation")

EVs_every_three_days_incubation_sensitivity=produce_results_main(3,10,1,'shorter_incubation')

EVs_every_week_incubation_sensitivity=produce_results_main(7,4,1,'shorter_incubation')

EVs_every_two_weeks_incubation_sensitivity=produce_results_main(14,2,1,'shorter_incubation')

all_EVs_incubation_sensitivity=cbind.data.frame(every_other_day=EVs_every_other_day_incubation_sensitivity$df,
                                                every_three_days=EVs_every_three_days_incubation_sensitivity$df,
                                                every_week=EVs_every_week_incubation_sensitivity$df,
                                                every_two_weeks=EVs_every_two_weeks_incubation_sensitivity$df)
melted_all_EVs_incubation_sensitivity=melt(all_EVs_incubation_sensitivity,variable.name='testing_strategy',value.name='weighted_risk_reduction')
write.csv(melted_all_EVs_incubation_sensitivity,"melted_all_EVs_incubation_sensitivity.csv",row.names=F)


write.csv(cbind.data.frame(matrix(c(EVs_every_other_day_incubation_sensitivity$median_prop_benefited,
                                    EVs_every_other_day_incubation_sensitivity$LQ_prop_benefited,
                                    EVs_every_other_day_incubation_sensitivity$UQ_prop_benefited,
                                    EVs_every_three_days_incubation_sensitivity$median_prop_benefited,
                                    EVs_every_three_days_incubation_sensitivity$LQ_prop_benefited,
                                    EVs_every_three_days_incubation_sensitivity$UQ_prop_benefited,
                                    EVs_every_week_incubation_sensitivity$median_prop_benefited,
                                    EVs_every_week_incubation_sensitivity$LQ_prop_benefited,
                                    EVs_every_week_incubation_sensitivity$UQ_prop_benefited,
                                    EVs_every_two_weeks_incubation_sensitivity$median_prop_benefited,
                                    EVs_every_two_weeks_incubation_sensitivity$LQ_prop_benefited,
                                    EVs_every_two_weeks_incubation_sensitivity$UQ_prop_benefited),
                                  byrow=TRUE,ncol=3,nrow=4)),
          "prop_benefited_all_incubation_sensitivity.csv",row.names=F)

write.csv(cbind.data.frame(matrix(c(EVs_every_other_day_incubation_sensitivity$median_given_drug_probs,
                                    EVs_every_other_day_incubation_sensitivity$LQ_given_drug_probs,
                                    EVs_every_other_day_incubation_sensitivity$UQ_given_drug_probs,
                                    EVs_every_three_days_incubation_sensitivity$median_given_drug_probs,
                                    EVs_every_three_days_incubation_sensitivity$LQ_given_drug_probs,
                                    EVs_every_three_days_incubation_sensitivity$UQ_given_drug_probs,
                                    EVs_every_week_incubation_sensitivity$median_given_drug_probs,
                                    EVs_every_week_incubation_sensitivity$LQ_given_drug_probs,
                                    EVs_every_week_incubation_sensitivity$UQ_given_drug_probs,
                                    EVs_every_two_weeks_incubation_sensitivity$median_given_drug_probs,
                                    EVs_every_two_weeks_incubation_sensitivity$LQ_given_drug_probs,
                                    EVs_every_two_weeks_incubation_sensitivity$UQ_given_drug_probs),
                                  byrow=TRUE,ncol=3,nrow=4)),
          "prop_given_drug_probs_all_incubation_sensitivity.csv",row.names=F)

onetest_output_incubation_sensitivity=produce_results_onetest(1,"shorter_incubation")
write.csv(onetest_output_incubation_sensitivity$list_weighted_RR,"onetest_output_RR_incubation_sensitivity.csv",row.names=F)
write.csv(onetest_output_incubation_sensitivity$list_benefited_or_offered,"onetest_output_benefited_or_offered_incubation_sensitivity.csv",row.names=F)
