

check=read.table("Inbreeding_depression_owls/All3085_FuniWE_FHBD512g.txt", header=T)

plot(check$FuniWE)

alex_check=read.table("../Downloads/3kowls_stats.txt")%>%
  rowid_to_column()

plot(alex_check$Funiw)

merged=full_join(check, alex_check, by = join_by('INDVs'=='ind.id'))%>%
  mutate(funi_norm=(Funiw/(max(Funiw)-min(Funiw))), fhbd_norm=(FHBD512gen/(max(FHBD512gen)-min(FHBD512gen))))

ggplot(merged, aes(x=FuniWE, y=Funiw))+
  geom_point()+
  labs(x='FunW Elu', y='FuniW JG')

ggplot(merged, aes(x=fhbd_norm, y=funi_norm))+
  geom_point()



cor(merged$FuniWE, merged$Funiw)

ggplot(merged, aes(y=FuniWE, x=rowid))+
  geom_point()+
  labs(x='FunW Elu', y='FuniW JG')

ggplot(merged, aes(y=FHBD512gen, x=rowid))+
  geom_point()


ggplot(merged, aes(y=Funiw, x=Coverage))+
  geom_point()+
  labs(x='Coverage', y='FuniW JG')

ggplot(merged, aes(y=FHBD512gen, x=Coverage))+
  geom_point()+
  labs(x='Coverage', y='FROH')


ggplot(merged, aes(y=Funiw, x=rowid))+
  geom_point()+
  labs(x='RingID index', y='FuniW JG')


ggplot(merged, aes(y=Funiw, x=FHBD512gen))+
  geom_point()+
  labs(x='FROH', y='FuniW JG')+
  geom_abline()



plot(merged$Funiw)



filts=merged%>%
  filter(Coverage<5)%>%
  select(INDVs)%>%
  rename(RingId = INDVs)


write.table(filts,
            file = "Inbreeding_depression_owls/pheno_df/RingIdsCov5.txt",
            row.names = F, quote = F, sep = ",",na = "NA")



inv=read.table("../Downloads/GenotypesInversionSC22.txt", header=T)

merged_inv=left_join(merged, inv, by = join_by('INDVs'=='sample'))


ggplot(merged_inv, aes(y=FHBD512gen, x=genotypeInv))+
  geom_point()



new_Fs=read.table("../Downloads/3kowls_stats_new.txt")
plot(new_Fs$Funiw)
plot(new_Fs$Funiw.wSS22)
