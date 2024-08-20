library(tidyverse)
library(brms)
library(viridis)
library(patchwork)
library(ggExtra)

setwd("/Users/ahewett1/Documents")

##function to extract h2 and proportion of variance explained by random effect in the model and make it in to a nice little table
get_variance_explained_by_rands <- function(model){
  
  ## for both stages in model
  stages <- c("Juvenile", "Adult")
  
  ## random effects we want to calculate the varinace explained from 
  # Note: Observer is not included here as this is human induced variation in the trait
  rands <- c("sd_RingId_" , "sd_RingId_pe_", "sd_year_", "sd_month_", "sd_nestboxID_" , "sd_clutch_merge_", "sd_rank_", "b_sigma")
  
  var_exp_list <- list()  # Store results for each stage
  
  for(stage in stages){
    Var.table <- as_draws_df(model)
    eff_stage <- paste0("_gr_stage", stage)
    all_variance_components <- paste0(rands, eff_stage)  # Creating list of all random effects to include in variance
    
    total_variance <- 0
    for(i in all_variance_components){
      total_variance <- total_variance + (Var.table[[i]])^2 ## adding var to other vars
    }
    mcmc_total_var <- as.mcmc(total_variance)  # Total variance per stage
    
    var_exp <- NULL
    
    for (i in 1:length(rands)){
      eff <- rands[i]  # For each effect in rands 
      rand_eff <- paste0(eff, "_gr_stage", stage)  # Getting variable names 
      var_rands <- as.mcmc((Var.table[[rand_eff]])^2) / mcmc_total_var  # Dividing by total variance
      
      mean <- summary(var_rands)$statistics[1]  # Mean estimate for h2
      Upr <- summary(var_rands)$quantiles[5]    # Upper CI (97.5%)
      Lwr <- summary(var_rands)$quantiles[1]    # Lower CI (2.5%)
      
      effect <- data.frame(
        effect = eff,
        Mean = mean,
        Lwr = Lwr,
        Upr = Upr,
        stage = stage
      )
      
      var_exp <- rbind(var_exp, effect)  # Combining with others
    }
    
    var_exp_list[[stage]] <- var_exp  # Append the results for each stage to var_exp_list
  }
  
  var_exp_df <- do.call(rbind, var_exp_list)  # Combine the results after the loop
  
  # Removing row names
  rownames(var_exp_df) <- NULL
  
  # Converting to a dataframe and tidying up
  variance_exp_per_stage <- var_exp_df %>%
    mutate(effect = as.factor(effect) %>%
             fct_recode(
               Va = "sd_RingId_", 
               pe = "sd_RingId_pe_", 
               year = "sd_year_", 
               month = "sd_month_", 
               nestbox = "sd_nestboxID_", 
               clutch = "sd_clutch_merge_",
               rank="sd_rank_", 
               resid = "b_sigma"
             ) %>%
             fct_relevel("resid", "year", "month", "nestbox", "clutch", "rank" ,"pe", "Va")
    ) %>%
    mutate(Mean = as.numeric(Mean),
           Upr = as.numeric(Upr),
           Lwr = as.numeric(Lwr))
  
  return(variance_exp_per_stage)
}



### function just to get Va
get_Va_only=function(model){
    Var.table <- as_draws_df(model)
    Va_only=NULL
    stages <- c("Juvenile", "Adult")
    
    for (stage in stages){
      Va_stage <- paste0("sd_RingId__gr_stage", stage)  # Getting variable names 
    
        Va=as.mcmc(Var.table[[Va_stage]])
        summary_Va <- summary(Va)
        mean<- summary_Va$statistics[1]  # Mean estimate for h2
        Upr <- summary_Va$quantiles[5]    # Upper CI (97.5%)
        Lwr <- summary_Va$quantiles[1]  ## Lower CI
        
        Va_df <- data.frame(
          Mean = mean,
          Lwr = Lwr,
          Upr = Upr,
          stage = stage)
        
        Va_only=rbind(Va_only,Va_df)
      }
    rownames(Va_only) <- NULL
    return(Va_only) # return df
}


#############################################

h2_split_bill=readRDS("Inbreeding_depression_owls/Model_outputs/2_split_stages/2.1.Bill_split_stage_h2_rank_rand.RDS")
summary(h2_split_bill)

h2_split_mass=readRDS("Inbreeding_depression_owls/Model_outputs/2_split_stages/2.2.Mass_split_stage_h2_rank_rand.RDS")
summary(h2_split_mass)

h2_split_tarsus=readRDS("Inbreeding_depression_owls/Model_outputs/2_split_stages/2.4.Tarsus_split_stage_h2_rank_rand.RDS")
summary(h2_split_tarsus)


###### get for each trait and combine to one df for wrapping 



bill=get_variance_explained_by_rands(h2_split_bill)%>%
  add_column(trait="Bill Length")

mass=get_variance_explained_by_rands(h2_split_mass)%>%
  add_column(trait="Mass") ## dont yet have rank rand 

tarsus=get_variance_explained_by_rands(h2_split_tarsus)%>%
  add_column(trait="Tarsus Length") ## dont yet have rank rand 

split_stages_all=bill%>% ## df for plotting 
  rbind(tarsus)%>%
  rbind(mass)

#### get Va for each trait ###
va_bill=get_Va_only(h2_split_bill)%>%
  add_column(trait="Bill Length")

va_mass=get_Va_only(h2_split_mass)%>%
  add_column(trait="Mass") ## 

va_tarsus=get_Va_only(h2_split_tarsus)%>%
  add_column(trait="Tarsus Length") ## 

Va_split_stages_all=va_bill%>% ## df for plotting 
  rbind(va_tarsus)%>%
  rbind(va_mass)

Va_split_stages_all[,1]

#####################################################################################################
####################################### PLOTS #######################################
####################################### ####################################### ######################
var=ggplot(split_stages_all, aes(fill=effect, y=Mean, x=factor(stage, level = c('Juvenile', 'Adult')))) + 
  facet_wrap(~trait)+
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




h2=split_stages_all%>%
  filter(effect=="Va")%>%
  ggplot(aes(group=trait, y=Mean, x=factor(stage, level = c('Juvenile', 'Adult')), ymin=Lwr, ymax=Upr, colour=trait))+
  scale_x_discrete(limits=c("Juvenile", "Adult"))+
  geom_line(size=1,  position=position_dodge(0.2)) +
  geom_point(size=4, position=position_dodge(0.2))+
  geom_errorbar(size=1,width=0.25, position=position_dodge(0.2))+
  theme_bw()+
  ylim(0,1)+
  labs(y="Heritability")+
  theme(text = element_text(size = 15), 
        axis.title.x = element_blank(),
        legend.text = element_text(size = 15))+
  scale_color_brewer(palette = "Dark2")+
  labs(colour="Trait")
  

h2


Va=Va_split_stages_all%>%
  ggplot(aes(group=trait, y=Mean, x=factor(stage, level = c('Juvenile', 'Adult')), ymin=Lwr, ymax=Upr, colour=trait))+
  scale_x_discrete(limits=c("Juvenile", "Adult"))+
  geom_line(size=1,  position=position_dodge(0.2)) +
  geom_point(size=4, position=position_dodge(0.2))+
  geom_errorbar(size=1,width=0.25, position=position_dodge(0.2))+
  theme_bw()+
  #ylim(0,1)+
  labs(y="Va")+
  theme(text = element_text(size = 15), 
        axis.title.x = element_blank(),
        legend.text = element_text(size = 15))+
  scale_color_brewer(palette = "Dark2")+
  labs(colour="Trait")


Va 

va_and_h2=Va + h2 +plot_layout(guides = "collect")
va_and_h2

fig_split=var/h2 + plot_annotation(tag_levels = 'A') +plot_layout(heights = c(5,4))

fig_split

ggsave("Inbreeding_depression_owls/Model_outputs/2_split_stages/plots/All_split_stages_wTarsus.png",
       plot=fig_split, 
       width = 9, 
       height = 8)
