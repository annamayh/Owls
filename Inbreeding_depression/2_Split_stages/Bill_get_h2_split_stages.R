library(tidyverse)
library(brms)
library(viridis)
library(patchwork)
library(ggExtra)


setwd("/Users/ahewett1/Documents")



h2_split=readRDS("Inbreeding_depression_owls/Model_outputs/2_split_stages/2.1.Bill_split_stage_h2_rank_rand.RDS")
summary(h2_split)
#plot(h2_split)


# checking the h2 
# Ring id / RingId , RingIdpe , clutch, year, month, nestbox, (not observer as this is human induced variation)
Var.table <- as_draws_df(h2_split)

h.bill.juv <- as.mcmc((Var.table$sd_RingId__gr_stageJuvenile)^2 / ((Var.table$sd_RingId__gr_stageJuvenile)^2 + (Var.table$sd_RingId_pe__gr_stageJuvenile)^2+(Var.table$sd_month__gr_stageJuvenile)^2 + (Var.table$sd_year__gr_stageJuvenile)^2 + 
                                                                     (Var.table$sd_rank__gr_stageJuvenile)^2 + (Var.table$sd_clutch_merge__gr_stageJuvenile)^2 + (Var.table$sd_nestboxID__gr_stageJuvenile)^2 + (Var.table$b_sigma_gr_stageJuvenile)^2))
summary(h.bill.juv)
summary(h.bill.juv)$statistics[1] ## mean estimate for h2
summary(h.bill.juv)$quantiles[1] ## upper CI
summary(h.bill.juv)$quantiles[5] ## lower CI 



## now h2 for adults 
h.bill.adu <- as.mcmc((Var.table$sd_RingId__gr_stageAdult)^2 / ((Var.table$sd_RingId__gr_stageAdult)^2 + (Var.table$sd_RingId_pe__gr_stageAdult)^2+(Var.table$sd_month__gr_stageAdult)^2 + (Var.table$sd_year__gr_stageAdult)^2 + 
                                                                  (Var.table$sd_rank__gr_stageAdult)^2+(Var.table$sd_clutch_merge__gr_stageAdult)^2 + (Var.table$sd_nestboxID__gr_stageAdult)^2  + (Var.table$b_sigma_gr_stageAdult)^2))

summary(h.bill.adu)$statistics[1] ## mean estimate for h2
summary(h.bill.adu)$quantiles[1] ## upper CI
summary(h.bill.adu)$quantiles[5] ## lower CI 




total.var.juv=as.mcmc((Var.table$sd_RingId__gr_stageJuvenile)^2 + (Var.table$sd_RingId_pe__gr_stageJuvenile)^2+(Var.table$sd_month__gr_stageJuvenile)^2 + (Var.table$sd_year__gr_stageJuvenile)^2 + 
                        (Var.table$sd_rank__gr_stageJuvenile)^2 + (Var.table$sd_clutch_merge__gr_stageJuvenile)^2  + (Var.table$sd_nestboxID__gr_stageJuvenile)^2 + (Var.table$b_sigma_gr_stageJuvenile)^2)




total.var.adult=as.mcmc((Var.table$sd_RingId__gr_stageAdult)^2 + (Var.table$sd_RingId_pe__gr_stageAdult)^2+(Var.table$sd_month__gr_stageAdult)^2 + (Var.table$sd_year__gr_stageAdult)^2 + 
                          (Var.table$sd_rank__gr_stageAdult)^2+(Var.table$sd_clutch_merge__gr_stageAdult)^2 + (Var.table$sd_nestboxID__gr_stageAdult)^2 + (Var.table$b_sigma_gr_stageAdult)^2)



# variance explained by other random effects 

rands=c("sd_RingId_" , "sd_RingId_pe_", "sd_year_", "sd_month_", "sd_nestboxID_" , "sd_clutch_merge_", "sd_rank_","b_sigma")

get_variance_explained=function(stage, total.var){
    
  var_exp=NULL
  
    for (i in 1:length(rands)){
      
        eff=rands[i]
        
        rand_eff=paste0(eff, "_gr_stage",stage)
        
        var_rands <- as.mcmc((Var.table[[rand_eff]])^2) / total.var
        
        mean=summary(var_rands)$statistics[1] ## mean estimate for h2
        Upr=summary(var_rands)$quantiles[1] ## upper CI
        Lwr=summary(var_rands)$quantiles[5] ## lower CI 
        
        effect=c(mean, Upr, Lwr, stage)
        var_exp=rbind(var_exp,effect)
    
      }
    
    row.names(var_exp)=rands
    
    var_exp_df=var_exp%>%
      as.data.frame()%>%
      rownames_to_column(var="effect") %>%
      rename(stage=V4, Lwr="2.5%", Upr="97.5%")%>%
      mutate(Mean=as.numeric(Mean))%>%
      mutate(Upr=as.numeric(Upr))%>%
      mutate(Lwr=as.numeric(Lwr))
    
    
    
    var_exp_df
    
    }



juve.var=get_variance_explained("Juvenile", total.var.juv)
sum(juve.var$Mean)

adult.var=get_variance_explained("Adult", total.var.adult)
sum(adult.var$Mean)


variance_exp_per_stage=rbind(juve.var, adult.var)%>%
  mutate(effect=as.factor(effect)%>%
           fct_recode(Va="sd_RingId_", 
                      pe="sd_RingId_pe_", 
                      year="sd_year_", 
                      month="sd_month_", 
                      nestbox="sd_nestboxID_", 
                      clutch="sd_clutch_merge_",
                      rank="sd_rank_", 
                      resid="b_sigma")%>%
           fct_relevel("resid", "year", "month", "nestbox", "clutch", "rank" ,"pe", "Va")
           
           )




var=ggplot(variance_exp_per_stage, aes(fill=effect, y=Mean, x=factor(stage, level = c('Juvenile', 'Adult')))) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_viridis(discrete = T, direction = -1, option = 'magma', begin=0.1, end = 0.96, 
                     labels=c(bquote (V["Residual"]), 
                              bquote (V["Birth year"]),
                              bquote (V["Birth month"]),
                              bquote (V["Nestbox ID"]),
                              bquote (V["Clutch ID"]),
                              bquote (V["Rank"]),
                              bquote (V["Permanent env"]),
                              bquote (V["a"]))) +
  theme_bw()+
  theme(text = element_text(size = 15), 
        axis.title.x = element_blank(),
        legend.text = element_text(size = 15))+
  labs(fill="Variance \ncomponents", y="Proportion of variance explained")

var


h2=variance_exp_per_stage%>%
  filter(effect=="Va")%>%
  ggplot(aes(colour=stage, y=Mean, x=factor(stage, level = c('Juvenile', 'Adult')), ymin=Lwr, ymax=Upr))+
  geom_pointrange()+
  ylim(c(0,1))+
  theme_bw()


h2

var+h2


h22=variance_exp_per_stage%>%
  filter(effect=="Va")%>%
  add_column(trait="bill")%>%
  ggplot(aes(group=trait, y=Mean, x=factor(stage, level = c('Juvenile', 'Adult')), ymin=Lwr, ymax=Upr, colour=trait))+
  scale_x_discrete(limits=c("Juvenile", "Adult"))+
  geom_line(size=1) +
  geom_point(size=4)+
  geom_errorbar(size=1,width=0.25, position=position_dodge(0.05))+
  theme_bw()+
  ylim(0,1)+
  labs(y="Heritability")+
  theme(text = element_text(size = 15), 
        axis.title.x = element_blank(),
        legend.text = element_text(size = 15))+
  labs(colour="Trait")
  
  


var+h22


