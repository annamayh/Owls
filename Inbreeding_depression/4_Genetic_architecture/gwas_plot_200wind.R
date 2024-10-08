library(tidyverse)
library(patchwork)

get_gwas_plot=function(input_files, trait){


    gwas_output_list <- lapply(Sys.glob(paste0(input_files)), readRDS)
     
    # To filter out the bad quality scaffolds later 
    good_sc=read.table("Barn_owl_general/ss_goodQ.txt")%>%
      mutate(good_super_scaffolds = as.numeric(stringr::str_split(V1, '_', simplify = TRUE)[, 2]))
    
    unlisted_gwas=do.call(rbind,gwas_output_list)%>%
      as.data.frame()%>%
      rename(p_val='Pr(>|z|)')%>%
      mutate(
        log_pval = -log10(as.numeric(p_val)),
        super_scaffold = as.numeric(stringr::str_split(window, '_', simplify = TRUE)[, 5]))%>%
      filter(super_scaffold%in%good_sc$good_super_scaffolds)%>%
      arrange(super_scaffold)%>%
      group_by(super_scaffold)%>%
      mutate(id = cur_group_id())
    
    bonf=0.05/nrow(unlisted_gwas)
    
    unlisted_gwas$Order <- 1:nrow(unlisted_gwas)
    
    axis.set <- unlisted_gwas %>% 
      group_by(super_scaffold)%>% 
      summarise(center = (max(Order) + min(Order)) / 2)
    
    
    plot=ggplot(unlisted_gwas, aes(Order, log_pval, col=as.factor(id%%2)), group=as.factor(id))+
      scale_colour_manual(values = c("firebrick3","grey68")) +
      geom_point()+
      geom_hline(yintercept=-log10(bonf), linetype=2)+
      theme_bw()+
      scale_x_continuous(label = axis.set$super_scaffold, breaks = axis.set$center, limits = c(0,nrow(unlisted_gwas)), expand = c(0, 0),)+
      theme(text = element_text(size = 15),
            legend.position = "none",
            axis.text.x = element_text(angle = 60, hjust=1))+
      labs(x="Super scaffold", y='-Log10(p-value)', title = paste0(trait))
    
    plot

}



setwd("/Users/ahewett1/Documents")

bill_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.1_bill_glmmTMB_GWAS_10gens/*.RDS.RDS"
mass_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.2_mass_glmmTMB_GWAS_10gens/*.RDS.RDS"
tarsus_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.4_tarsus_glmmTMB_GWAS_10gens/*.RDS.RDS"



bill_gwas=get_gwas_plot(bill_files, 'Bill Length')

mass_gwas=get_gwas_plot(mass_files, 'Mass')

tarsus_gwas=get_gwas_plot(tarsus_files, 'Tarsus Length')


gwas_plot=bill_gwas/mass_gwas/tarsus_gwas+
  plot_annotation(tag_levels = 'A')+
  plot_layout(axes = 'collect')

gwas_plot

ggsave(gwas_plot,
       file = "Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4_glmmTMB_gwas_200wind_10gens.png",
       width = 10,
       height = 8)

