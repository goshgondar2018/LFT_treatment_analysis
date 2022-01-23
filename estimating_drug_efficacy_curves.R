library(ggpubr)

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


all_RRs=cbind.data.frame(time=c(0,1,2,3,4,5,6,7,8,9,10),RR=c(RR_period1,RR_period2,RR_period3,RR_period4,
                                                              RR_period5,RR_period6,RR_period7,
                                                              RR_period8,RR_period9,RR_period9,RR_period9))

risk_period1=exp(coefficients(log_binomial)[[1]]+0*coefficients(log_binomial)[[2]]+coefficients(log_binomial)[[3]]+0*coefficients(log_binomial)[[4]])
risk_period2=exp(coefficients(log_binomial)[[1]]+1*coefficients(log_binomial)[[2]]+coefficients(log_binomial)[[3]]+1*coefficients(log_binomial)[[4]])
risk_period3=exp(coefficients(log_binomial)[[1]]+2*coefficients(log_binomial)[[2]]+coefficients(log_binomial)[[3]]+2*coefficients(log_binomial)[[4]])
risk_period4=exp(coefficients(log_binomial)[[1]]+3*coefficients(log_binomial)[[2]]+coefficients(log_binomial)[[3]]+3*coefficients(log_binomial)[[4]])
risk_period5=exp(coefficients(log_binomial)[[1]]+4*coefficients(log_binomial)[[2]]+coefficients(log_binomial)[[3]]+4*coefficients(log_binomial)[[4]])
risk_period6=exp(coefficients(log_binomial)[[1]]+5*coefficients(log_binomial)[[2]]+coefficients(log_binomial)[[3]]+5*coefficients(log_binomial)[[4]])
risk_period7=exp(coefficients(log_binomial)[[1]]+6*coefficients(log_binomial)[[2]]+coefficients(log_binomial)[[3]]+6*coefficients(log_binomial)[[4]])
risk_period8=exp(coefficients(log_binomial)[[1]]+7*coefficients(log_binomial)[[2]]+coefficients(log_binomial)[[3]]+7*coefficients(log_binomial)[[4]])
risk_period9=exp(coefficients(log_binomial)[[1]]+4*coefficients(log_binomial)[[2]])

all_risks=cbind.data.frame(time=c(0,1,2,3,4,5,6,7,8,9,10),risk_treated=c(risk_period1,risk_period2,risk_period3,
                                        risk_period4,risk_period5,risk_period6,risk_period7,risk_period8,
                                        risk_period9,risk_period9,risk_period9))
observed_RRs=cbind.data.frame(time=c(3,5),RR=c((5/697)/(44/682),(8/1039)/(66/1046)))
observed_predicted_RRs=merge(all_RRs,observed_RRs,by='time',all.x=TRUE)
fig1=ggplot(observed_predicted_RRs,aes(x=time))+geom_line(aes(y=RR.x))+geom_point(aes(y=RR.y),col='blue',size=3)+
  scale_x_continuous(breaks=seq(0,10,2))+
  theme_classic()+ylab("RR")

#fig1=ggplot(all_RRs,aes(x=time,y=RR))+geom_line()+scale_x_continuous(breaks=seq(0,10,2))+theme_classic()
observed_risks_treated=cbind.data.frame(time=c(3,5),risk_treated=c((5/697),(8/1039)))
observed_predicted_risks_treated=merge(all_risks,observed_risks_treated,by='time',all.x=TRUE)

fig2=ggplot(observed_predicted_risks_treated,aes(x=time))+geom_line(aes(y=risk_treated.x))+
  geom_point(aes(y=risk_treated.y),col='blue',size=3)+
  scale_x_continuous(breaks=seq(0,10,2))+
  theme_classic()+
  ylab("risk among the treated")

logistic=glm(hosp~days+treated+days*treated,data=pfizer_data,family=binomial(link='logit'))

OR_period1=exp(coefficients(logistic)[[3]]+0*coefficients(logistic)[[4]])
OR_period2=exp(coefficients(logistic)[[3]]+1*coefficients(logistic)[[4]])
OR_period3=exp(coefficients(logistic)[[3]]+2*coefficients(logistic)[[4]])
OR_period4=exp(coefficients(logistic)[[3]]+3*coefficients(logistic)[[4]])
OR_period5=exp(coefficients(logistic)[[3]]+4*coefficients(logistic)[[4]])
OR_period6=exp(coefficients(logistic)[[3]]+5*coefficients(logistic)[[4]])
OR_period7=exp(coefficients(logistic)[[3]]+6*coefficients(logistic)[[4]])
OR_period8=exp(coefficients(logistic)[[3]]+7*coefficients(logistic)[[4]])
OR_period9=1

all_ORs=cbind.data.frame(time=c(0,1,2,3,4,5,6,7,8,9,10),OR=c(OR_period1,OR_period2,OR_period3,OR_period4,
                                                             OR_period5,OR_period6,OR_period7,OR_period8,
                                                             OR_period9,OR_period9,OR_period9))

odds_period1=exp(coefficients(logistic)[[1]]+0*coefficients(logistic)[[2]]+coefficients(logistic)[[3]]+0*coefficients(logistic)[[4]])
odds_period2=exp(coefficients(logistic)[[1]]+1*coefficients(logistic)[[2]]+coefficients(logistic)[[3]]+1*coefficients(logistic)[[4]])
odds_period3=exp(coefficients(logistic)[[1]]+2*coefficients(logistic)[[2]]+coefficients(logistic)[[3]]+2*coefficients(logistic)[[4]])
odds_period4=exp(coefficients(logistic)[[1]]+3*coefficients(logistic)[[2]]+coefficients(logistic)[[3]]+3*coefficients(logistic)[[4]])
odds_period5=exp(coefficients(logistic)[[1]]+4*coefficients(logistic)[[2]]+coefficients(logistic)[[3]]+4*coefficients(logistic)[[4]])
odds_period6=exp(coefficients(logistic)[[1]]+5*coefficients(logistic)[[2]]+coefficients(logistic)[[3]]+5*coefficients(logistic)[[4]])
odds_period7=exp(coefficients(logistic)[[1]]+6*coefficients(logistic)[[2]]+coefficients(logistic)[[3]]+6*coefficients(logistic)[[4]])
odds_period8=exp(coefficients(logistic)[[1]]+7*coefficients(logistic)[[2]]+coefficients(logistic)[[3]]+7*coefficients(logistic)[[4]])
odds_period9=exp(coefficients(logistic)[[1]]+4*coefficients(logistic)[[2]])

all_odds=cbind.data.frame(time=c(0,1,2,3,4,5,6,7,8,9,10),
                          odds_treated=c(odds_period1,odds_period2,odds_period3,odds_period4,
                                                      odds_period5,odds_period6,odds_period7,odds_period8,
                                                      odds_period9,odds_period9,odds_period9))


observed_ORs=cbind.data.frame(time=c(3,5),OR=c(((5/697)*(1-(44/682)))/((1-(5/697))*(44/682)),
                                               ((8/1039)*(1-(66/1046)))/((1-(8/1039))*(66/1046))))
observed_predicted_ORs=merge(all_ORs,observed_ORs,by='time',all.x=TRUE)
fig3=ggplot(observed_predicted_ORs,aes(x=time))+geom_line(aes(y=OR.x))+geom_point(aes(y=OR.y),col='blue',size=3)+
  scale_x_continuous(breaks=seq(0,10,2))+
  theme_classic()+ylab("OR")

#fig3=ggplot(all_ORs,aes(x=time,y=OR))+geom_line()+scale_x_continuous(breaks=seq(0,10,2))+theme_classic()

observed_odds_treated=cbind.data.frame(time=c(3,5),odds_treated=c((5/697)/(1-(5/697)),(8/1039)/(1-(8/1039))))
observed_predicted_odds_treated=merge(all_odds,observed_odds_treated,by='time',all.x=TRUE)
fig4=ggplot(observed_predicted_odds_treated,aes(x=time))+geom_line(aes(y=odds_treated.x))+
  geom_point(aes(y=odds_treated.y),col='blue',size=3)+
  scale_x_continuous(breaks=seq(0,10,2))+
  theme_classic()+
  ylab("odds among the treated")

ggarrange(fig1,fig2,fig3,fig4)

ggsave("drug_efficacy_curves.pdf")



