library(tidyverse)
library(reshape2)
library(plotly)
# read in LFT and PCR positive curves from MCMC simulations
LFT_curve_summary=read_csv("LFT_curve_summary.csv")%>%
  filter(days_since_infection%in%seq(0,30))%>%
  select(days_since_infection,median)
colnames(LFT_curve_summary)[2]<-"value"


# define function to generate probability weighted reductions in hospitalization & mortality for every possible
# combination of when an LFT is administered, under testing strategies at different frequencies
generate_EVs<-function(drug_coverage_prop,delay_test_to_treatment,j_max,data,index,incubation_period,
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
    product_all_probs2=prod(1-data[seq(i,i+(index*j_max),index),"value"]*drug_coverage_prop*1,na.rm=TRUE)
    list_product_all_probs2=c(list_product_all_probs2,product_all_probs2)
  }
  sum_EV=(1/index)*sum(list_product_all_probs,list_product_all_probs2)

  return(sum_EV)
}

# ref: https://stackoverflow.com/questions/4391592/run-a-function-with-multiple-values-for-more-than-one-arguments-that-are-not-the
all_params=expand.grid(drug_coverage_prop=seq(0,1,length.out=11),delay_test_to_treatment=seq(1,7))
defined_function<-generate_EVs
weighted_RRs_every_other_day=mapply("defined_function",drug_coverage_prop=all_params$drug_coverage_prop,
       delay_test_to_treatment=all_params$delay_test_to_treatment,
       MoreArgs = list(j_max=15,data=LFT_curve_summary,index=2,incubation_period=5,
                       drug_efficacy_period1=RR_period1,
                       drug_efficacy_period2=RR_period2,
                       drug_efficacy_period3=RR_period3,
                       drug_efficacy_period4=RR_period4,
                       drug_efficacy_period5=RR_period5,
                       drug_efficacy_period6=RR_period6,
                       drug_efficacy_period7=RR_period7,
                       drug_efficacy_period8=RR_period8,
                       drug_efficacy_period9=RR_period9))
all_params_weighted_RRs_every_other_day=cbind.data.frame(all_params,weighted_RR=weighted_RRs_every_other_day)
colnames(all_params_weighted_RRs_every_other_day)<-c('drug_coverage','delay','weighted_RR')


plot_ly(all_params_weighted_RRs_every_other_day,x=~delay,y=~drug_coverage,z=~weighted_RR,type='contour')%>%
        layout(xaxis=list(title='delay test to treatment'),
               yaxis = list(title='drug coverage proportion'))%>%
  colorbar(title = "weighted RR")


colnames(all_params_weighted_RRs_every_other_day)<-c('drug coverage','delay','weighted RR')
plot_ly(all_params_weighted_RRs_every_other_day,x=~delay,y=~`drug coverage`,z=~`weighted RR`,
        type='scatter3d',mode='markers',marker = list(color = ~'#E69F00'))
             

weighted_RRs_every_three_days=mapply("defined_function",drug_coverage_prop=all_params$drug_coverage_prop,
                                    delay_test_to_treatment=all_params$delay_test_to_treatment,
                                    MoreArgs = list(j_max=10,data=LFT_curve_summary,index=3,incubation_period=5,
                                                    drug_efficacy_period1=RR_period1,
                                                    drug_efficacy_period2=RR_period2,
                                                    drug_efficacy_period3=RR_period3,
                                                    drug_efficacy_period4=RR_period4,
                                                    drug_efficacy_period5=RR_period5,
                                                    drug_efficacy_period6=RR_period6,
                                                    drug_efficacy_period7=RR_period7,
                                                    drug_efficacy_period8=RR_period8,
                                                    drug_efficacy_period9=RR_period9))
all_params_weighted_RRs_every_three_days=cbind.data.frame(all_params,weighted_RR=weighted_RRs_every_three_days)
colnames(all_params_weighted_RRs_every_three_days)<-c('drug_coverage','delay','weighted_RR')


plot_ly(all_params_weighted_RRs_every_three_days,x=~delay,y=~drug_coverage,z=~weighted_RR,type='contour')%>%
  layout(xaxis=list(title='delay test to treatment'),
         yaxis = list(title='drug coverage proportion'))%>%
  colorbar(title = "weighted RR")

colnames(all_params_weighted_RRs_every_three_days)<-c('drug coverage','delay','weighted RR')
plot_ly(all_params_weighted_RRs_every_three_days,x=~delay,y=~`drug coverage`,z=~`weighted RR`,
        type='scatter3d',mode='markers',marker = list(color = ~'#CC79A7'))


weighted_RRs_every_week=mapply("defined_function",drug_coverage_prop=all_params$drug_coverage_prop,
                                           delay_test_to_treatment=all_params$delay_test_to_treatment,
                                           MoreArgs = list(j_max=4,data=LFT_curve_summary,index=7,incubation_period=5,
                                                           drug_efficacy_period1=RR_period1,
                                                           drug_efficacy_period2=RR_period2,
                                                           drug_efficacy_period3=RR_period3,
                                                           drug_efficacy_period4=RR_period4,
                                                           drug_efficacy_period5=RR_period5,
                                                           drug_efficacy_period6=RR_period6,
                                                           drug_efficacy_period7=RR_period7,
                                                           drug_efficacy_period8=RR_period8,
                                                           drug_efficacy_period9=RR_period9))
all_params_weighted_RRs_every_week=cbind.data.frame(all_params,weighted_RR=weighted_RRs_every_week)
colnames(all_params_weighted_RRs_every_week)<-c('drug_coverage','delay','weighted_RR')

plot_ly(all_params_weighted_RRs_every_week,x=~delay,y=~drug_coverage,z=~weighted_RR,type='contour')%>%
  layout(xaxis=list(title='delay test to treatment'),
         yaxis = list(title='drug coverage proportion'))%>%
  colorbar(title = "weighted RR")

colnames(all_params_weighted_RRs_every_week)<-c('drug coverage','delay','weighted RR')
plot_ly(all_params_weighted_RRs_every_week,x=~delay,y=~`drug coverage`,z=~`weighted RR`,
        type='scatter3d',mode='markers',marker = list(color = ~'#009E73'))



weighted_RRs_every_two_weeks=mapply("defined_function",drug_coverage_prop=all_params$drug_coverage_prop,
                                     delay_test_to_treatment=all_params$delay_test_to_treatment,
                                     MoreArgs = list(j_max=2,data=LFT_curve_summary,index=14,incubation_period=5,
                                                     drug_efficacy_period1=RR_period1,
                                                     drug_efficacy_period2=RR_period2,
                                                     drug_efficacy_period3=RR_period3,
                                                     drug_efficacy_period4=RR_period4,
                                                     drug_efficacy_period5=RR_period5,
                                                     drug_efficacy_period6=RR_period6,
                                                     drug_efficacy_period7=RR_period7,
                                                     drug_efficacy_period8=RR_period8,
                                                     drug_efficacy_period9=RR_period9))
all_params_weighted_RRs_every_two_weeks=cbind.data.frame(all_params,weighted_RR=weighted_RRs_every_two_weeks)
colnames(all_params_weighted_RRs_every_two_weeks)<-c('drug_coverage','delay','weighted_RR')


plot_ly(all_params_weighted_RRs_every_two_weeks,x=~delay,y=~drug_coverage,z=~weighted_RR,type='contour')%>%
  layout(xaxis=list(title='delay test to treatment'),
         yaxis = list(title='drug coverage proportion'))%>%
  colorbar(title = "weighted RR")

colnames(all_params_weighted_RRs_every_two_weeks)<-c('drug coverage','delay','weighted RR')
plot_ly(all_params_weighted_RRs_every_two_weeks,x=~delay,y=~`drug coverage`,z=~`weighted RR`,
        type='scatter3d',mode='markers',marker = list(color = ~'#56B4E9'))






