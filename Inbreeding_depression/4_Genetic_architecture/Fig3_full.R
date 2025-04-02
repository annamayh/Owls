library(tidyverse)
library(patchwork)
library(ggimage)
library(grid)
library(png)
library(cowplot)
library(ggtext)
library(viridis)

get_estimates=function(input_files){
  
  gwas_output_list <- lapply(Sys.glob(paste0(bill_files)), readRDS)
  
  unlisted_gwas=do.call(rbind,gwas_output_list)%>%
    as.data.frame()%>%
    mutate(super_scaffold = as.numeric(stringr::str_split(Window, '_', simplify = TRUE)[, 5]))%>%
    arrange(super_scaffold)%>%
    group_by(super_scaffold)%>%
    mutate(id = cur_group_id())
  
  bonf=0.05/nrow(unlisted_gwas)
  
  unlisted_gwas$sig <- ifelse(unlisted_gwas$p_val < bonf, "yes","no")      
  
  
  unlisted_gwas
  
}



get_estimates_and_plot=function(input_files, trait){
  
  gwas_output_list <- lapply(Sys.glob(paste0(input_files)), readRDS)
  unlisted_gwas=do.call(rbind,gwas_output_list)%>%
    as.data.frame()%>%
    mutate(super_scaffold = as.numeric(stringr::str_split(Window, '_', simplify = TRUE)[, 5]))%>%
    arrange(super_scaffold)%>%
    group_by(super_scaffold)%>%
    mutate(id = cur_group_id())
  
  mean_est=mean(unlisted_gwas$estimate)
  
  unlisted_gwas$Order <- 1:nrow(unlisted_gwas)
  
  axis.set <- unlisted_gwas %>% 
    group_by(super_scaffold)%>% 
    summarise(center = (max(Order) + min(Order)) / 2)
  
  colour2="grey"
  colour1="black"
  sig="red"
  
  bonf=0.05/nrow(unlisted_gwas)
  
  unlisted_gwas$alpha_group <- ifelse(unlisted_gwas$p_val <bonf, 1, 0.2)
  unlisted_gwas$color_group <- ifelse(unlisted_gwas$p_val < bonf, "sig", as.factor(unlisted_gwas$id %% 2))      
  unlisted_gwas$size_group <- ifelse(unlisted_gwas$p_val <bonf, 3, 1)  
  
  
  ## estimate plots
  brms_estimates=ggplot(unlisted_gwas, aes(x=Order, y=estimate, 
                                           col=color_group, 
                                           alpha=alpha_group,
                                           size=size_group,
                                           group=as.factor(id)))+
    
    geom_point()+
    geom_hline(yintercept=mean_est, lty=2) + 
    geom_hline(yintercept=0, lty=1) +  
    theme_bw()+
    scale_x_continuous(label = axis.set$super_scaffold, breaks = axis.set$center, limits = c(0,nrow(unlisted_gwas)), expand = c(0, 0),)+
    theme(text = element_text(size = 15),
          legend.position = "none",
          axis.text.x = element_text(angle = 40, hjust=1, size=8))+
    scale_color_manual(values=c(colour1, colour2, sig),guide=FALSE)+
    scale_alpha_identity() +
    scale_size_identity() +  # Use size_group values literally
    labs(x="Super scaffold", y='Estimate', title = paste0(trait)
    )
  
  brms_estimates
  
  
}


setwd("/Users/ahewett1/Documents")
bill_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.1_bill_brmsGWAS_FINAL/*.RDS.RDS"
mass_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.2_mass_brmsGWAS_FINAL/*.RDS.RDS"
tarsus_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.4_tarsus_brmsGWAS_FINAL/*.RDS.RDS"


bill_gen_arch=get_estimates(bill_files)
mass_gen_arch=get_estimates(mass_files)
tarsus_gen_arch=get_estimates(tarsus_files)

## summary numbers 
mean(bill_gen_arch$estimate)
mean(mass_gen_arch$estimate)
mean(tarsus_gen_arch$estimate)


colours = viridis(n = 3, end = 0.8)

bill_est = get_estimates_and_plot(bill_files, 'Bill length')+
  labs(tag = 'A')

bill_denisty=ggplot(bill_gen_arch, aes(x=estimate)) +
  geom_density(fill=colours[3], alpha=0.8)+
  theme_classic()+
  geom_vline(xintercept = 0)+
  labs(x='Estimate', y='Density')+
  theme(axis.title.y=element_blank())+
  coord_flip()



mass_est = get_estimates_and_plot(mass_files, 'Mass')+
  labs(tag = 'B')


mass_denisty=ggplot(mass_gen_arch, aes(x=estimate)) +
  geom_density(fill=colours[2], alpha=0.8)+
  theme_classic()+
  geom_vline(xintercept = 0)+
  labs(x='Estimate', y='Density')+
  theme(axis.title.y=element_blank())+
  coord_flip()



tarsus_est = get_estimates_and_plot(tarsus_files, 'Tarsus length')+
  labs(tag = 'C')

tarsus_density=ggplot(tarsus_gen_arch, aes(x=estimate)) +
  geom_density(fill=colours[1], alpha=0.65)+
  theme_classic()+
  geom_vline(xintercept = 0)+
  labs(x='Estimate', y='Density')+
  theme(axis.title.y=element_blank())+
  coord_flip()

est_plot=bill_est+bill_denisty+
  mass_est+mass_denisty+
  tarsus_est+tarsus_density+
  plot_layout(axes = 'collect', axis_titles ='collect', ncol = 2, nrow = 3, widths = c(3, 1))
  

est_plot



ggsave(est_plot,
       file = "Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4_brms_IDGWAS_ests.png",
       width = 11,
       height = 9)



