
library(patchwork)



get_gwas_plot=function(unlisted, trait){

all_windows_gwas_plot=unlisted%>%
  rename(p_val='Pr(>|z|)')%>%
  mutate(
    log_pval = -log10(as.numeric(p_val)),
    super_scaffold = as.numeric(stringr::str_split(window, '_', simplify = TRUE)[, 5]))%>%
    arrange(super_scaffold)%>%
  group_by(super_scaffold)%>%
  mutate(id = cur_group_id())

bonf=0.05/nrow(all_windows_gwas_plot)
  
all_windows_gwas_plot$Order <- 1:nrow(all_windows_gwas_plot)

axis.set <- all_windows_gwas_plot %>% 
  group_by(super_scaffold)%>% 
  summarise(center = (max(Order) + min(Order)) / 2)


plot=ggplot(all_windows_gwas_plot, aes(Order, log_pval, col=as.factor(id%%2)), group=as.factor(id))+
  scale_colour_manual(values = c("firebrick3","grey68")) +
  geom_point()+
  geom_hline(yintercept=-log10(bonf), linetype=2)+
  theme_bw()+
  scale_x_continuous(label = axis.set$super_scaffold, breaks = axis.set$center, limits = c(0,nrow(all_windows_gwas_plot)), expand = c(0, 0),)+
  theme(text = element_text(size = 15),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust=1))+
  labs(x="Super scaffold", y='-Log10(p-value)', title = paste0(trait))

plot

}



unlisted_bill=readRDS("Inbreeding_depression_owls/Model_outputs/4_Gen_arch/glmTMB_GWAS_HBD_winds7500.RDS")
unlisted_mass=readRDS("Inbreeding_depression_owls/Model_outputs/4_Gen_arch/mass_glmTMB_GWAS_HBD_winds7500.RDS")
unlisted_tarsus=readRDS("Inbreeding_depression_owls/Model_outputs/4_Gen_arch/tarsus_glmTMB_GWAS_HBD_winds7500.RDS")



bill_gwas=get_gwas_plot(unlisted_bill, 'Bill Length')

mass_gwas=get_gwas_plot(unlisted_mass, 'Mass')

tarsus_gwas=get_gwas_plot(unlisted_tarsus, 'Tarsus Length')


gwas_plot=bill_gwas/mass_gwas/tarsus_gwas+
  plot_annotation(tag_levels = 'A')+
  plot_layout(axes = 'collect')

gwas_plot

ggsave(gwas_plot,
       file = "Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4_glmmTMB_gwas_7500wind_all.png",
       width = 10,
       height = 8)
