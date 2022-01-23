setwd("~/Desktop/Oxford/LFT_drug_analysis/cluster_output")
library(tidyverse)
library(ggpubr)
# main scenario
melted_all_EVs=read.csv("melted_all_EVs.csv")%>%
  mutate(testing_strategy=factor(testing_strategy,
                                 levels=c("every_other_day","every_three_days", "every_week","every_two_weeks")))

prop_benefited=read.csv("prop_benefited_all.csv")

#rownames(prop_benefited)<-c("every_other_day","every_three_days", "every_week","every_two_weeks")
colnames(prop_benefited)<-c("median","LQ","UQ")
prop_benefited$testing_strategy=c("every_other_day","every_three_days", "every_week","every_two_weeks")
prop_benefited=prop_benefited%>%
  mutate(testing_strategy=factor(testing_strategy,
                                 levels=c("every_other_day","every_three_days", "every_week","every_two_weeks")))

prop_given_drug=read.csv("prop_given_drug_probs_all.csv")
colnames(prop_given_drug)<-c("median","LQ","UQ")
prop_given_drug$testing_strategy=c("every_other_day","every_three_days", "every_week","every_two_weeks")
prop_given_drug=prop_given_drug%>%
  mutate(testing_strategy=factor(testing_strategy,
                                 levels=c("every_other_day","every_three_days", "every_week","every_two_weeks")))

prop_given_drug_and_benefited=read.csv("prop_given_drug_and_benefited_all.csv")
colnames(prop_given_drug_and_benefited)<-c("median","LQ","UQ")
prop_given_drug_and_benefited$testing_strategy=c("every_other_day","every_three_days", "every_week","every_two_weeks")
prop_given_drug_and_benefited=prop_given_drug_and_benefited%>%
  mutate(testing_strategy=factor(testing_strategy,
                                 levels=c("every_other_day","every_three_days", "every_week","every_two_weeks")))

# palette from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbPalette <- c("#CC79A7", "#E69F00",  "#009E73","#56B4E9")#, "#F0E442", "#0072B2", "#D55E00")

# boxplots of distribution of risk reduction estimates across testing frequencies
labels=c("every other day","every three days", "every week","every two weeks")#,"one-time PCR") #,"one-time PCR")

plot1=ggplot(melted_all_EVs,aes(x=fct_reorder(testing_strategy,weighted_risk_reduction,.desc=FALSE),
                          weighted_risk_reduction,fill=testing_strategy))+
  geom_boxplot()+theme_classic()+scale_x_discrete(labels=labels)+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values=cbPalette,name='testing strategy',labels=labels)+
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank())+ylab("weighted RR")

plot2=ggplot(prop_given_drug,aes(x=testing_strategy,y=median,ymin=LQ,ymax=UQ))+
  geom_point(aes(color=testing_strategy),size=3)+
  scale_color_manual(values=cbPalette,name='testing strategy',labels=labels)+
  geom_errorbar()+
  theme_classic()+scale_x_discrete(labels=labels)+
  ylim(c(0.0,1.0))+xlab("testing strategy")+ylab("proportion given drug")+
  theme(axis.text.x = element_blank(),axis.title.x=element_blank(),
        legend.position='none')
  
ggarrange(plot1,plot2,nrow=2)
ggsave("~/Desktop/Oxford/LFT_drug_analysis/figures/main_EVs_zero_delays.pdf")

# extract medians, LQs, and UQs for each strategy
melted_all_EVs_summarized<-melted_all_EVs%>%
  group_by(testing_strategy)%>%
  summarise(median=quantile(weighted_risk_reduction,probs=0.5),
            LQ=quantile(weighted_risk_reduction,probs=0.02),
            UQ=quantile(weighted_risk_reduction,probs=0.975))%>%
  mutate(scenario='base case')

##

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

ggplot(melted_all_EVs_summarized_combined,aes(fill=scenario,x=fct_reorder(testing_strategy,median,.desc=FALSE),
                                              y=median))+
  geom_bar(position='dodge',stat='identity')+
  theme_classic()+scale_x_discrete(labels=labels)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("testing strategy")+ylab("weighted RR")+  scale_fill_manual(values=cbPalette)

ggsave("~/Desktop/Oxford/LFT_drug_analysis/figures/comparing_efficacy_scenarios.pdf")

