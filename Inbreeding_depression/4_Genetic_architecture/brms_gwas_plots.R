library(tidyverse)
library(patchwork)
library(ggimage)
library(grid)
library(png)
library(cowplot)
library(ggtext)

#chicken=read.table()


get_estimate_plot=function(input_files){
  
      gwas_output_list <- lapply(Sys.glob(paste0(input_files)), readRDS)
      
      # good_sc=read.table("Barn_owl_general/ss_goodQ.txt")%>%
      #   mutate(good_super_scaffolds = as.numeric(stringr::str_split(V1, '_', simplify = TRUE)[, 2]))
      
      unlisted_gwas=do.call(rbind,gwas_output_list)%>%
        as.data.frame()%>%
        mutate(super_scaffold = as.numeric(stringr::str_split(Window, '_', simplify = TRUE)[, 5]))%>%
        #filter(super_scaffold%in%good_sc$good_super_scaffolds)%>%
        arrange(super_scaffold)%>%
        group_by(super_scaffold)%>%
        mutate(id = cur_group_id())
      
      mean_est=mean(unlisted_gwas$estimate)
      
      unlisted_gwas$Order <- 1:nrow(unlisted_gwas)
      
      axis.set <- unlisted_gwas %>% 
        group_by(super_scaffold)%>% 
        summarise(center = (max(Order) + min(Order)) / 2)
      
     # unlisted_gwas$alpha_group <- ifelse(unlisted_gwas$sig == "no", 0.5, 1)
      ## set alternating colours
      #colour1=RColorBrewer::brewer.pal(1, 'Dark2')[3]
      #colour2=RColorBrewer::brewer.pal(1, 'Dark2')[2]
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
        #scale_shape_manual(values = c(13, 16))+
        scale_alpha_identity() +
        scale_size_identity() +  # Use size_group values literally
        labs(x="Super scaffold", y='Estimate'#, title = paste0(trait)
             )
      
      brms_estimates
      

}


get_gwas_plot=function(input_files, trait){
  
  gwas_output_list <- lapply(Sys.glob(paste0(input_files)), readRDS)
  
  # good_sc=read.table("Barn_owl_general/ss_goodQ.txt")%>%
  #   mutate(good_super_scaffolds = as.numeric(stringr::str_split(V1, '_', simplify = TRUE)[, 2]))
  # 
  unlisted_gwas=do.call(rbind,gwas_output_list)%>%
    as.data.frame()%>%
    mutate(super_scaffold = as.numeric(stringr::str_split(Window, '_', simplify = TRUE)[, 5]))%>%
    #filter(super_scaffold%in%good_sc$good_super_scaffolds)%>%
    arrange(super_scaffold)%>%
    group_by(super_scaffold)%>%
    mutate(id = cur_group_id())
  
  mean_est=mean(unlisted_gwas$estimate)
  
  unlisted_gwas$Order <- 1:nrow(unlisted_gwas)
  
  axis.set <- unlisted_gwas %>% 
    group_by(super_scaffold)%>% 
    summarise(center = (max(Order) + min(Order)) / 2)
  
  # unlisted_gwas$alpha_group <- ifelse(unlisted_gwas$sig == "no", 0.5, 1)
  ## set alternating colours
  colour1=RColorBrewer::brewer.pal(1, 'Dark2')[3]
  colour2="dodgerblue2"
  
  sig="red"
  
  bonf=0.05/nrow(unlisted_gwas)
  
  unlisted_gwas$alpha_group <- ifelse(unlisted_gwas$p_val <bonf, 1, 0.1)
  unlisted_gwas$color_group <- ifelse(unlisted_gwas$p_val < bonf, "sig", as.factor(unlisted_gwas$id %% 2))      
  

### and gwas like plot
brms_gwas=ggplot(unlisted_gwas, aes(x=Order, y=-log10(p_val), 
                                    col=color_group,
                                    alpha=alpha_group,
                                    group=as.factor(id)))+
  
  geom_point()+
  geom_hline(yintercept=-log10(bonf), lty=2) + 
  theme_bw()+
  scale_x_continuous(label = axis.set$super_scaffold, breaks = axis.set$center, limits = c(0,nrow(unlisted_gwas)), expand = c(0, 0),)+
  theme(text = element_text(size = 15),
        legend.position = "none",
        axis.text.x = element_text(angle = 40, hjust=1, size=8))+
  scale_color_manual(values=c(colour1, colour2, sig),guide=FALSE)+
  scale_alpha_identity() +
  scale_size_identity() +  # Use size_group values literally
  labs(x="Super scaffold", y='-Log10(p-value)'#, title = paste0(trait)
       )



brms_gwas

}


#input_files=tarsus_files


setwd("/Users/ahewett1/Documents")

bill_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.1_bill_brmsGWAS_FINAL/*.RDS.RDS"
mass_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.2_mass_brmsGWAS_FINAL/*.RDS.RDS"
tarsus_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.4_tarsus_brmsGWAS_FINAL/*.RDS.RDS"



bill_gwas=get_gwas_plot(bill_files, 'Bill Length')

mass_gwas=get_gwas_plot(mass_files, 'Mass')

tarsus_gwas=get_gwas_plot(tarsus_files, 'Tarsus Length')


gwas_plot=bill_gwas/mass_gwas/tarsus_gwas+
  plot_annotation(tag_levels = 'A')+
  plot_layout(axes = 'collect', axis_titles ='collect')

gwas_plot


## estimates plot
## bill 

bill_est <- get_estimate_plot(bill_files)+
  labs(title = "<img src = 'Inbreeding_depression_owls/Model_outputs/4_Gen_arch/plots/beak.png' height = 30>", 
       #subtitle = "Bill length"
       )+
  theme(
    plot.title = ggtext::element_markdown(hjust = 0.5))
  
 

mass_est <- get_estimate_plot(mass_files)+
  labs(title = "<img src = 'Inbreeding_depression_owls/Model_outputs/4_Gen_arch/plots/weighing_scales.png' height = 25>", 
       #subtitle = "Mass"
       )+
  theme(
    plot.title = ggtext::element_markdown(hjust = 0.5))

#<span style='font-size: 24pt'>Mass</span>

## tarsus

tarsus_est <- get_estimate_plot(tarsus_files)+
  labs(title = "<img src = 'Inbreeding_depression_owls/Model_outputs/4_Gen_arch/plots/tarsus_image.png' height = 25>",
       #subtitle = "Tarsus length"
       )+
  theme(
    plot.title = ggtext::element_markdown(hjust = 0.5))


est_plot=bill_est/mass_est/tarsus_est+
  plot_annotation(tag_levels = 'A')+
  plot_layout(axes = 'collect', axis_titles ='collect')

est_plot
# ggsave(gwas_plot,
#        file = "Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4_brms_IDGWAS_gwas.png",
#        width = 8,
#        height = 9)
# 
# 
# 
ggsave(est_plot,
       file = "Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4_brms_IDGWAS_ests3.png",
       width = 8,
       height = 9)


