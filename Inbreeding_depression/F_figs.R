library(tidyverse)
library(gghalves)

setwd("/Users/ahewett1/Documents")
inbr=read.table("Inbreeding_depression_owls/All3085_FuniWE_FHBD512g.txt", header = T)

  

##Fgrm and FROH from elu
inb=read.table("Inbreeding_depression_owls/All3085_FuniWE_FHBD512g.txt", header = T)%>%
  rename(RingId=INDVs)%>%
  pivot_longer(cols = c('FHBD512gen', 'FuniWE'), names_to = 'type', values_to = 'F')

head(inb)

mean(inbr$FuniWE)

(f_dist=ggplot(inb, aes(y=F, x = type, fill = type))+
  geom_half_point(side = "l", shape = 21, alpha = 0.3, width = 0.6,) +
  geom_half_violin(side = "r",  alpha = 0.8)+ 
  stat_summary(fun = mean, 
    geom = "point", 
    color = "black", 
    size = 3)+
  theme_classic()+
  coord_flip()+
  theme(legend.position = "none", 
        axis.title.y = element_blank(),
        text = element_text(size = 18))+
  labs(y='Inbreeding coefficient value')+
  scale_x_discrete(labels = c(expression('F'[ROH]), expression('F'[uniWE])))+
  scale_fill_manual(values = c('gray40', 'gray80'))+
  geom_hline(yintercept=0, lty=2)
)


ggsave(f_dist, 
       file= "Inbreeding_depression_owls/F_dists.png",
       width = 5,
       height = 4,
       bg = 'white'
)


save(f_dist, 
       file= "Inbreeding_depression_owls/F_dists.RData"
)
