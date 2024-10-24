library(tidyverse)
library(patchwork)




get_gwas_plot=function(input_files, trait){
  
      gwas_output_list <- lapply(Sys.glob(paste0(input_files)), readRDS)
      
      good_sc=read.table("Barn_owl_general/ss_goodQ.txt")%>%
        mutate(good_super_scaffolds = as.numeric(stringr::str_split(V1, '_', simplify = TRUE)[, 2]))
      
      unlisted_gwas=do.call(rbind,gwas_output_list)%>%
        as.data.frame()%>%
        rownames_to_column(var='window')%>%
        mutate(super_scaffold = as.numeric(stringr::str_split(window, '_', simplify = TRUE)[, 5]))%>%
        filter(super_scaffold%in%good_sc$good_super_scaffolds)%>%
        arrange(super_scaffold)%>%
        group_by(super_scaffold)%>%
        mutate(id = cur_group_id())%>%
        mutate(sig=case_when(
          Q2.5>0|Q97.5<0 ~ 'yes', 
          .default = 'no'
        ))
      
      
      mean_est=mean(unlisted_gwas$Estimate)
      
      unlisted_gwas$Order <- 1:nrow(unlisted_gwas)
      
      axis.set <- unlisted_gwas %>% 
        group_by(super_scaffold)%>% 
        summarise(center = (max(Order) + min(Order)) / 2)
      
      unlisted_gwas$alpha_group <- ifelse(unlisted_gwas$sig == "no", 0.18, 1)
      ## set alternating colours
      colour1=RColorBrewer::brewer.pal(1, 'Dark2')[3]
      colour2=RColorBrewer::brewer.pal(1, 'Dark2')[2]
      
      brms_gwas=ggplot(unlisted_gwas, aes(x=Order, y=Estimate, 
                                          col=as.factor(id%%2), 
                                          alpha=alpha_group,
                                          group=as.factor(id)))+
        
        geom_point()+
        geom_hline(yintercept=0, lty=2) +  
        geom_hline(yintercept=mean_est, lty=3) + 
        theme_bw()+
        scale_x_continuous(label = axis.set$super_scaffold, breaks = axis.set$center, limits = c(0,nrow(unlisted_gwas)), expand = c(0, 0),)+
        theme(text = element_text(size = 15),
              legend.position = "none",
              axis.text.x = element_text(angle = 40, hjust=1, size=8))+
        scale_color_manual(values=c(colour1, colour2),guide=FALSE)+
        scale_alpha_identity() +
        labs(x="Super scaffold", y='Estimate', title = paste0(trait))
      
      
      brms_gwas
}




setwd("/Users/ahewett1/Documents")

bill_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.1_bill_brmsGWAS/*.RDS.RDS"
mass_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.2_mass_brmsGWAS/*.RDS.RDS"
tarsus_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.4_tarsus_brmsGWAS/*.RDS.RDS"



bill_gwas=get_gwas_plot(bill_files, 'Bill Length')

mass_gwas=get_gwas_plot(mass_files, 'Mass')

tarsus_gwas=get_gwas_plot(tarsus_files, 'Tarsus Length')


gwas_plot=bill_gwas/mass_gwas/tarsus_gwas+
  plot_annotation(tag_levels = 'A')+
  plot_layout(axes = 'collect', axis_titles ='collect')

gwas_plot



ggsave(gwas_plot,
       file = "Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4_brms_IDGWAS.png",
       width = 8,
       height = 9)

