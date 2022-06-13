setwd("~/Desktop/Oxford/LFT_drug_analysis/all_output_may2")
library(tidyverse)
library(ggpubr)
# main scenario
melted_all_EVs=read.csv("melted_all_EVs.csv")

prop_benefited=read.csv("prop_benefited_all.csv")
colnames(prop_benefited)<-c("median","LQ","UQ")
prop_benefited$testing_strategy=c("every_other_day","every_three_days", "every_week","every_two_weeks")

prop_given_drug=read.csv("prop_given_drug_probs_all.csv")
colnames(prop_given_drug)<-c("median","LQ","UQ")
prop_given_drug$testing_strategy=c("every_other_day","every_three_days", "every_week","every_two_weeks")

## add in prop benefited/given drug under one test strategy output
onetest_output_RR=read.csv("onetest_output_RR.csv")
melted_all_EVs_2=rbind.data.frame(melted_all_EVs[which(melted_all_EVs$testing_strategy%in%c("every_other_day","every_three_days")),],
                                cbind.data.frame(testing_strategy=rep('one_time_testing',4000),
                                                 weighted_risk_reduction=onetest_output_RR$x),
                                melted_all_EVs[which(melted_all_EVs$testing_strategy%in%c("every_week","every_two_weeks")),])%>%
  mutate(testing_strategy=factor(testing_strategy,
                                 levels=c("every_other_day","every_three_days", "one_time_testing","every_week",
                                          "every_two_weeks")))

onetest_output_benefited_or_offered=read.csv("onetest_output_benefited_or_offered.csv")
median_benefited_onetest = quantile(onetest_output_benefited_or_offered$x,probs=0.5)
LQ_benefited_onetest = quantile(onetest_output_benefited_or_offered$x,probs=0.025)
UQ_benefited_onetest = quantile(onetest_output_benefited_or_offered$x,probs=0.975)

prop_given_drug_2=rbind.data.frame(prop_given_drug[which(prop_given_drug$testing_strategy%in%c("every_other_day","every_three_days")),],
                                   cbind.data.frame(median=median_benefited_onetest,
                                                    LQ=LQ_benefited_onetest,UQ=UQ_benefited_onetest,
                                                    testing_strategy='one_time_testing'),
                                   prop_given_drug[which(prop_given_drug$testing_strategy%in%c("every_week","every_two_weeks")),])%>%
  mutate(testing_strategy=factor(testing_strategy,
                                 levels=c("every_other_day","every_three_days","one_time_testing", "every_week","every_two_weeks")))


## palette from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbPalette <- c("#CC79A7", "#E69F00",  "#009E73","#56B4E9","#D55E00")#, "#F0E442", "#0072B2", "#D55E00")

## boxplots of distribution of risk reduction estimates across testing frequencies
labels=c("every other day","every three days","one-time testing","every week","every two weeks")

plot1=ggplot(melted_all_EVs_2,aes(x=fct_reorder(testing_strategy,weighted_risk_reduction),
                                weighted_risk_reduction,fill=testing_strategy))+
  geom_boxplot()+theme_classic()+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values=cbPalette,name='testing strategy',labels=labels)+
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank())+ylab("weighted RR")

plot2=ggplot(prop_given_drug_2,aes(x=fct_reorder(testing_strategy,median,.desc=TRUE),y=median,ymin=LQ,ymax=UQ))+
  geom_point(aes(color=testing_strategy),size=3)+
  scale_color_manual(values=cbPalette,name='testing strategy',labels=labels)+
  geom_errorbar()+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank())+
  ylim(c(0.0,1.0))+xlab("testing strategy")+ylab("proportion given treatment")

ggarrange(plot1,plot2,nrow=2,legend='right',labels=c('A','B'))

ggsave("~/Desktop/Oxford/LFT_drug_analysis/all_figures_may2/main_EVs_zero_delays.pdf")

## extract medians, LQs, and UQs for each strategy
melted_all_EVs_summarized<-melted_all_EVs_2%>%
  group_by(testing_strategy)%>%
  summarise(median=quantile(weighted_risk_reduction,probs=0.5),
            LQ=quantile(weighted_risk_reduction,probs=0.02),
            UQ=quantile(weighted_risk_reduction,probs=0.975))%>%
  filter(testing_strategy%in%c("every_other_day","every_three_days","every_week","every_two_weeks"))%>%
  mutate(scenario='base case')

# repeat above w/ omicron-specific incubation period sensitivity analysis output

melted_all_EVs_incubation_sensitivity=read.csv("melted_all_EVs_incubation_sensitivity.csv")

prop_benefited_incubation_sensitivity=read.csv("prop_benefited_all_incubation_sensitivity.csv")
colnames(prop_benefited_incubation_sensitivity)<-c("median","LQ","UQ")
prop_benefited_incubation_sensitivity$testing_strategy=c("every_other_day","every_three_days", "every_week","every_two_weeks")

prop_given_drug_incubation_sensitivity=read.csv("prop_given_drug_probs_all_incubation_sensitivity.csv")
colnames(prop_given_drug_incubation_sensitivity)<-c("median","LQ","UQ")
prop_given_drug_incubation_sensitivity$testing_strategy=c("every_other_day","every_three_days", "every_week","every_two_weeks")

## output under one test strategy
onetest_output_RR_incubation_sensitivity=read.csv("onetest_output_RR_incubation_sensitivity.csv")
melted_all_EVs_2_incubation_sensitivity=rbind.data.frame(melted_all_EVs_incubation_sensitivity[which(melted_all_EVs_incubation_sensitivity$testing_strategy%in%c("every_other_day","every_three_days",
                                                                                                                                                                 "every_week")),],
                                                       cbind.data.frame(testing_strategy=rep('one_time_testing',4000),
                                                                        weighted_risk_reduction=onetest_output_RR_incubation_sensitivity$x),
                                                       melted_all_EVs_incubation_sensitivity[which(melted_all_EVs_incubation_sensitivity$testing_strategy=="every_two_weeks"),])%>%
  mutate(testing_strategy=factor(testing_strategy,
                                 levels=c("every_other_day","every_three_days", "every_week","one_time_testing",
                                          "every_two_weeks")))

onetest_output_benefited_or_offered_incubation_sensitivity=read.csv("onetest_output_benefited_or_offered_incubation_sensitivity.csv")
median_benefited_onetest_incubation_sensitivity = quantile(onetest_output_benefited_or_offered_incubation_sensitivity$x,probs=0.5)
LQ_benefited_onetest_incubation_sensitivity = quantile(onetest_output_benefited_or_offered_incubation_sensitivity$x,probs=0.025)
UQ_benefited_onetest_incubation_sensitivity = quantile(onetest_output_benefited_or_offered_incubation_sensitivity$x,probs=0.975)

prop_given_drug_2_incubation_sensitivity=rbind.data.frame(prop_given_drug[which(prop_given_drug_incubation_sensitivity$testing_strategy%in%c("every_other_day","every_three_days","every_week")),],
                                                         cbind.data.frame(median=median_benefited_onetest_incubation_sensitivity,
                                                         LQ=LQ_benefited_onetest_incubation_sensitivity,UQ=UQ_benefited_onetest_incubation_sensitivity,
                                                         testing_strategy='one_time_testing'),
                                                         prop_given_drug_incubation_sensitivity[which(prop_given_drug_incubation_sensitivity$testing_strategy=="every_two_weeks"),])%>%
  mutate(testing_strategy=factor(testing_strategy,
                                 levels=c("every_other_day","every_three_days","every_week","one_time_testing","every_two_weeks")))

# palette from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

## boxplots of distribution of risk reduction estimates across testing frequencies

labels_incubation_sensitivity=c("every other day","every three days","every week","one-time testing","every two weeks")
plot1_incubation_sensitivity=ggplot(melted_all_EVs_2_incubation_sensitivity,
                                    aes(x=fct_reorder(testing_strategy,weighted_risk_reduction),
                                        weighted_risk_reduction,fill=testing_strategy))+
  geom_boxplot()+theme_classic()+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values=cbPalette,name='testing strategy',labels=labels_incubation_sensitivity)+
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank())+ylab("weighted RR")


plot2_incubation_sensitivity=ggplot(prop_given_drug_2_incubation_sensitivity,
                                    aes(x=fct_reorder(testing_strategy,median,.desc=TRUE),
                                        y=median,ymin=LQ,ymax=UQ))+
  geom_point(aes(color=testing_strategy),size=3)+
  scale_color_manual(values=cbPalette,name='testing strategy',labels=labels_incubation_sensitivity)+
  geom_errorbar()+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank())+
  ylim(c(0.0,1.0))+xlab("testing strategy")+ylab("proportion given treatment")

ggarrange(plot1_incubation_sensitivity,plot2_incubation_sensitivity,nrow=2,legend='right',labels=c('A','B'))
ggsave("~/Desktop/Oxford/LFT_drug_analysis/all_figures_may2/main_EVs_zero_delays_incubation_sensitivity.pdf")

## extract medians, LQs, and UQs for each strategy
melted_all_EVs_summarized_incubation_sensitivity<-melted_all_EVs_2_incubation_sensitivity%>%
  group_by(testing_strategy)%>%
  summarise(median=quantile(weighted_risk_reduction,probs=0.5),
            LQ=quantile(weighted_risk_reduction,probs=0.02),
            UQ=quantile(weighted_risk_reduction,probs=0.975))%>%
  mutate(scenario='base case')

## additional drug efficacy scenarios
melted_all_EVs_drop_zero=read.csv("melted_all_EVs_drop_zero.csv")%>%
  mutate(testing_strategy=factor(testing_strategy,
                                 levels=c("every_other_day","every_three_days", "every_week","every_two_weeks")))

melted_all_EVs_drop_zero_summarized<-melted_all_EVs_drop_zero%>%
  group_by(testing_strategy)%>%
  summarise(median=quantile(weighted_risk_reduction,probs=0.5),
            LQ=quantile(weighted_risk_reduction,probs=0.02),
            UQ=quantile(weighted_risk_reduction,probs=0.975))%>%
  mutate(scenario='fast decline to zero scenario')

melted_all_EVs_efficacy_preserved=read.csv("melted_all_EVs_efficacy_preserved.csv")%>%
  mutate(testing_strategy=factor(testing_strategy,
                                 levels=c("every_other_day","every_three_days", "every_week","every_two_weeks")))

melted_all_EVs_efficacy_preserved_summarized<-melted_all_EVs_efficacy_preserved%>%
  group_by(testing_strategy)%>%
  summarise(median=quantile(weighted_risk_reduction,probs=0.5),
            LQ=quantile(weighted_risk_reduction,probs=0.02),
            UQ=quantile(weighted_risk_reduction,probs=0.975))%>%
  mutate(scenario='efficacy preserved')

# combine all drug efficacy scenarios
melted_all_EVs_summarized_combined=rbind.data.frame(melted_all_EVs_summarized,
                                                    melted_all_EVs_drop_zero_summarized,
                                                    melted_all_EVs_efficacy_preserved_summarized)
melted_all_EVs_summarized_combined$scenario=replace(melted_all_EVs_summarized_combined$scenario,
                                                     which(melted_all_EVs_summarized_combined$scenario=='fast decline to zero scenario'),'fast decline to zero')

cbPalette_scenarios=c("#CC79A7", "#E69F00",  "#009E73")#, "#F0E442", "#0072B2", "#D55E00")

## boxplots of distribution of risk reduction estimates across testing frequencies
labels_scenarios=c("base case","efficacy preserved", "fast decline to zero")

efficacy_scenarios_facet_plot=ggplot(melted_all_EVs_summarized_combined,
           aes(x=scenario,y=median,ymin=LQ,ymax=UQ))+
  geom_point(aes(color=scenario),size=3)+
  scale_color_manual(values=cbPalette_scenarios,name='scenario',labels=labels_scenarios)+
  geom_errorbar()+
  theme_classic()+
  theme(legend.position='bottom')+
  ylim(c(0.0,1.0))+xlab("testing strategy")+ylab("weighted RR")+
  #scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1))+  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

efficacy_scenarios_facet_plot+
  facet_wrap(~testing_strategy,labeller=as_labeller(facet_labels),ncol=5,nrow=1)

ggsave("~/Desktop/Oxford/LFT_drug_analysis/all_figures_may2/comparing_efficacy_scenarios.pdf")

# test coverage sensitivity analysis
test_covg_sensitivity=read_csv("all_props_every_strategy.csv")
test_covg_sensitivity_summarized=test_covg_sensitivity%>%
  group_by(test_covg_prop,testing_strategy)%>%
  summarise(median_prop_benefited=quantile(prop_benefited,probs=0.5),
            LQ_prop_benefited=quantile(prop_benefited,probs=0.025),
            UQ_prop_benefited=quantile(prop_benefited,probs=0.975),
            median_prop_given_drug=quantile(prop_given_drug,probs=0.5),
            LQ_prop_given_drug=quantile(prop_given_drug,probs=0.025),
            UQ_prop_given_drug=quantile(prop_given_drug,probs=0.975))%>%
  mutate(testing_strategy=factor(testing_strategy,
                                 levels=c("every_other_day","every_three_days", "every_week","every_two_weeks","one-time testing")))

## proportion given drug facet plot

test_covg_sensitivity_summarized$test_covg_prop_factor=as.factor(test_covg_sensitivity_summarized$test_covg_prop)
facet_labels=c('every_other_day'='every other day','every_three_days'='every three days',
               'every_week'='every week','every_two_weeks'='every two weeks','one-time testing'='one-time testing')

test_covg_sensitivity_given_drug_facet_plot=ggplot(test_covg_sensitivity_summarized,
                                                   aes(x=test_covg_prop,
                                                       y=median_prop_given_drug,ymin=LQ_prop_given_drug,ymax=UQ_prop_given_drug))+
  geom_point(aes(color=testing_strategy),size=3)+
  scale_color_manual(values=cbPalette,name='testing strategy',labels=labels)+
  geom_errorbar()+
  theme_classic()+
  theme(legend.position='bottom')+
  ylim(c(0.0,1.0))+xlab("proportion tested")+ylab("proportion given treatment")+
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1))+  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

test_covg_sensitivity_given_drug_facet_plot+
  facet_wrap(~testing_strategy,labeller=as_labeller(facet_labels),ncol=5,nrow=1)
ggsave("~/Desktop/Oxford/LFT_drug_analysis/all_figures_may2/test_covg_sensitivity_given_drug_plot.pdf")

## prop benefited
test_covg_sensitivity_prop_benefited_facet_plot=ggplot(test_covg_sensitivity_summarized,
                                                       aes(x=test_covg_prop,y=median_prop_benefited,ymin=LQ_prop_benefited,ymax=UQ_prop_benefited))+
  geom_point(aes(color=testing_strategy),size=3)+
  scale_color_manual(values=cbPalette,name='testing strategy',labels=labels)+
  geom_errorbar()+
  theme_classic()+
  theme(legend.position='bottom')+
  ylim(c(0.0,1.0))+xlab("proportion tested")+ylab("proportion benefited")+
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1))+  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

test_covg_sensitivity_prop_benefited_facet_plot+facet_wrap(~testing_strategy,labeller=as_labeller(facet_labels),
                                                           ncol=5,nrow=1)
ggsave("~/Desktop/Oxford/LFT_drug_analysis/all_figures_may2/test_covg_sensitivity_benefited_plot.pdf")


##

