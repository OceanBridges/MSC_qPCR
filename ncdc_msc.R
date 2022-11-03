rm(list=ls())

# libraries ---------------------------------------------------------------
library(maps)
library(mapdata)
library(fields)
library(ncdf4)
library(raster) # package for raster manipulation
library(ggpubr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(treemapify)
library(RColorBrewer)
library(patchwork)
library(forcats)
library(ggalluvial)
library(networkD3)
library(remotes)
library(networkD3)
library(ggalluvial)
library(ggsankey)


# data --------------------------------------------------------------------
ncdc_msc<-read.table('ncdc_msc.txt', header=T,sep="\t") %>% filter(.,!is.na(lon)) ###quito NA


ncdc_msc_noucyna<-read.table('ncdc_msc_withoutucyna.txt', header=T,sep="\t") %>% filter(.,!is.na(lon)) ###quito NA
ncdc_msc_noucyna_nonans<-ncdc_msc_noucyna[complete.cases(ncdc_msc_noucyna), ]

rates<-read.table('msc_rates.txt', header=T,sep="\t") 

connectivity<-read.table('connectivity.txt', header=T,sep="\t") %>% filter(.,!is.na(nifh)) ###quito NA
connectivity2<-read.table('connectivity2.txt', header=T,sep="\t") 
connectivity3<-read.table('connectivity3.txt', header=T,sep="\t") 

# group data of each fraction by station
ncdc_msc_noucyna_nonans<-as.data.frame(aggregate(nifh ~ fraction + station + target + lon + lat, ncdc_msc_noucyna_nonans,sum))%>% filter(.,!is.na(nifh)) ###quito NA
write.table(ncdc_msc_noucyna_nonans,file="ncdc_msc_noucyna_nonans.txt")
ncdc_msc_noucyna_nonans$station<-as.factor(ncdc_msc_noucyna_nonans$station)
ncdc_msc_noucyna_nonans$fraction<-as.factor(ncdc_msc_noucyna_nonans$fraction)
ncdc_msc_noucyna_nonans$nifh_log<-log10(ncdc_msc_noucyna_nonans$nifh)

write.table(ncdc_msc_noucyna_nonans,file="ncdc_msc_noucyna_nonans.txt")

test<-as.data.frame(aggregate(x=ncdc_msc_noucyna_nonans$nifh,by=list(ncdc_msc_noucyna_nonans$station,ncdc_msc_noucyna_nonans$target),FUN=sum))
write.table(test,file="test.txt")

percentage<-read.table('percentage.txt', header=T,sep="\t") 

percentage<-gather(percentage, key = "pct_msc_fraction", value = "pct_value", pct_fs,pct_ss,pct_susp)


ncdc_msc_noucyna_nonans<-ncdc_msc_noucyna_nonans %>%
  mutate(target=fct_relevel(target,"NP9","gETSP2","903_","gamma3","gammaA","NP10","M126"))





# satellite map ---------------------------------------------------------------------
continent<-map('worldHires',plot=F,fill=T,interior=F,xlim=c(-180,180),ylim=c(-90,90))
coast<-map('worldHires',interior=F,plot=F,xlim=c(-180,180),ylim=c(-90,90))

ssh.dir <- 'global-analysis-forecast-phy-001-024_1661624425016.nc'  
search()
ls(2)
x <- nc_open(ssh.dir)
ssh <- ncvar_get(x,'zos')
ssh<-apply(ssh, c(1,2), mean) ###promedio sobre tercera dimension (tiempo)
lat <- ncvar_get(x,'latitude')
lon <- ncvar_get(x,'longitude')

# Segunda parte-extraer valores de variable (de mapa) en funci?n de coordenadas de underway
# 
# # Rasterizo el mapa
ncdc_msc_noucyna_nonans$lon[(ncdc_msc_noucyna_nonans$lon<0)]<-360-abs(ncdc_msc_noucyna_nonans$lon)
plot(ncdc_msc_noucyna_nonans$lon,ncdc_msc_noucyna_nonans$lat)
ssh.r <- raster(as.matrix(ssh),xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon))
ssh.r<- flip(t(ssh.r), direction='y') ##lo giro
plot(ssh.r)

# #extraigo valores
ssh_ext_msc<-raster::extract(ssh.r,cbind(ncdc_msc_noucyna_nonans$lon,ncdc_msc_noucyna_nonans$lat), method="bilinear",df=TRUE) 
ssh_ext_msc_pct<-raster::extract(ssh.r,cbind(percentage$lon,percentage$lat), method="bilinear",df=TRUE) 

# #concateno en la tabla de datos 
ncdc_msc_noucyna_nonans$ssh<-ssh_ext_msc$layer 
percentage$ssh<-ssh_ext_msc_pct$layer 

# plot by group-ejemplo random para probar
ggplot(ncdc_msc_noucyna_nonans, aes(x = ssh, y = nifh, color = fraction)) +
  geom_point()


# stats -------------------------------------------------------------------

#all nifh per fraction vs ssh
ncdc_msc_ssh<-ggscatter(ncdc_msc_noucyna_nonans,x="ssh",y="nifh_log", color="fraction",
          add="reg.line",#conf.int=TRUE,
          cor.coef=TRUE,cor.method="spearman",
          xlab="SSH (m)",ylab="NCD abundance (log10 nifH l-1)")
ncdc_msc_ssh + 
  scale_fill_brewer(palette = "YlGnBu") +
  facet_grid(. ~ fraction) +
  theme(strip.text.x = element_text(size = 8, face = "bold"), 
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"))




# plots -------------------------------------------------------------------

# # m126
ggplot(data=m126_msc, aes(x=station, y=nifh, fill=fraction)) +
  geom_bar(stat="identity",position="stack") 

ggplot(data=m126_msc_percentage, aes(x=station, y=percent, fill=fraction)) +
  geom_bar(stat="identity",position="stack") 

jpeg(filename="m126_percentage.jpg",unit='cm',width = 25, height = 20, res=300)
ggscatter(m126_msc_percentage,x="ssh",y="percent",
          add="reg.line",#conf.int=TRUE,
          #cor.coef=TRUE,
          cor.method="spearman",
          xlab="SSH",ylab="MSC fraction %contribution",
          color = "fraction", palette = "Accent") +
  theme(legend.text=element_text(size=30),axis.text = element_text(size=25),axis.title=element_text(size=20)) 
dev.off()

# # gammaA
ggplot(data=gammaA_msc, aes(x=station, y=nifh, fill=fraction)) +
  geom_bar(stat="identity",position="stack") 

jpeg(filename="gammaA_percentage.jpg",unit='cm',width = 25, height = 20, res=300)
ggscatter(gammaA_msc_percentage,x="ssh",y="percent",
          add="reg.line",#conf.int=TRUE,
          #cor.coef=TRUE,
          cor.method="spearman",
          xlab="SSH",ylab="MSC fraction %contribution",
          color = "fraction", palette = "Accent") +
  theme(legend.text=element_text(size=30),axis.text = element_text(size=25),axis.title=element_text(size=20)) 
dev.off()

# # all ncds
jpeg(filename="total_nifh_per_group_underway.jpg",unit='cm',width = 25, height = 20, res=300)
underway<-read.table('underway.txt', header=T,sep="\t") %>% filter(.,!is.na(lon)) ###quito NA
ggplot(underway, aes(x = station, y = nifh, fill = target)) +
  geom_col() +
  theme_classic() +
  theme(legend.position="top")+labs(x="station",y="nifH copies l-1", title="nifH copies per NCD group") + 
  scale_x_discrete(guide = guide_axis(angle = 45))+
  guides(fill=guide_legend(title="NCD group"))
dev.off()

# bubble grid by station faceting by fraction
jpeg(filename="bubble_facets.jpg",unit='cm',width = 23, height = 12, res=300)
ncdc_msc_plot<-ggplot(ncdc_msc_noucyna_nonans, aes(x = station, y = target)) + 
  geom_point(aes(size = nifh, fill = target), alpha = 2, shape = 21) + 
  scale_size_continuous(limits = c(100, 1000000), range=c(5,20), breaks = c(100,1000,10000,100000,1000000)) + 
  labs( x= "MSC fraction", y = "NCD group", size = "nifH gene copies l-1", fill = "NCD group")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 8, face = "bold", vjust = 0.3, hjust = 1, angle=90), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 8), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.text = element_text(size = 8, face ="bold", colour ="black"), 
        legend.title = element_text(size = 8, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "left") +
  guides(
    colour = guide_legend(show = FALSE) 
  ) 
ncdc_msc_plot + 
  scale_fill_brewer(palette = "YlGnBu") +
  facet_grid(. ~ fraction) +
  theme(strip.text.x = element_text(size = 8, face = "bold"), 
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"))
dev.off()

# doughnut plot % fraction per group
donut_fs<-
  ggplot(subset(ncdc_msc_noucyna_nonans, fraction %in% c("fs")), aes(x = fraction, y = nifh, fill = target)) +
  geom_col(position="fill",width=0.5) +
  scale_x_discrete(limits = c(" ", "fs")) +
  scale_fill_brewer(palette = "YlGnBu", guide = "none") +
  coord_polar("y") +
  theme_void() +
  labs(fill="NCD group\n(nifH gene copies l-1)") 

donut_ss<-ggplot(subset(ncdc_msc_noucyna_nonans, fraction %in% c("ss")), aes(x = fraction, y = nifh, fill = target)) +
  geom_col(position="fill",width=0.5) +
  scale_x_discrete(limits = c(" ", "ss")) +
  scale_fill_brewer(palette = "YlGnBu", guide = "none") +
  coord_polar("y") +
  theme_void() +
  labs(fill="NCD group\n(nifH gene copies l-1)") 

donut_susp<-ggplot(subset(ncdc_msc_noucyna_nonans, fraction %in% c("susp")), aes(x = fraction, y = nifh, fill = target)) +
  geom_col(position="fill",width=0.5) +
  scale_x_discrete(limits = c(" ", "susp")) +
  scale_fill_brewer(palette = "YlGnBu", guide = guide_legend(reverse = TRUE)) +
  coord_polar("y") +
  theme_void() +
  labs(fill="NCD group\n(nifH gene copies l-1)") 

jpeg(filename="donut.jpg",unit='cm',width = 25, height = 20, res=300)
donut_fs+donut_ss+donut_susp
dev.off()

# corrs ssh vs %contribution each target to total nifh
ggscatter(subset(percentage, pct_msc_fraction %in% c("pct_fs")),x="ssh",y="pct_value",color="target",
                        #add="reg.line",#conf.int=TRUE,
                        #cor.coef=TRUE,
                        cor.method="spearman",
                        xlab="SSH (m)",ylab="%fs contribution", palette="YlGnBu") 

# rates
jpeg(filename="n2fix_rates_msc.jpg",unit='cm',width = 12, height = 12, res=300)
ncdc_msc_rates_plot<-ggplot(rates, aes(x = station, y = fraction)) + 
  geom_point(aes(size = n2fix, fill = fraction), alpha = 2, shape = 21) + 
  scale_size_continuous(limits = c(0.1, 31), range=c(5,20), breaks = c(0.1,2,8,12,30)) + 
  labs( x= "station", y = "MSC fraction", size = "N2 fixation rate (nmol N l-1 d-1)", fill = "MSC fraction")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 8, face = "bold", vjust = 0.3, hjust = 1, angle=90), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 8), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.text = element_text(size = 8, face ="bold", colour ="black"), 
        legend.title = element_text(size = 8, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") +
  guides(
    colour = guide_legend(show = FALSE) 
  ) 
ncdc_msc_rates_plot + 
  scale_fill_brewer(palette = "Set2") 
dev.off()


# Sankey plot connectivity surface-MSC

# create nodes dataframe
# df <- connectivity %>%
#   make_long(station,fraction,target,nifh)
# 
# ggplot(as.data.frame(connectivity),
#        aes(y = nifh, axis1 = fraction, axis2 = fraction)) +
#   geom_alluvium(aes(fill = target), width = 1/12) +
#   geom_stratum(width = 1/12, fill = "black", color = "grey") +
#   geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
#   scale_fill_brewer(type = "qual", palette = "Set1")

# # aggregate by fraction
# sankey_df<-connectivity %>%
#   dplyr::group_by(fraction) %>%
#   dplyr::summarise(udw=sum(udw),susp=sum(susp),ss=sum(ss),fs=sum(fs))
# 
# # create nodes dataframe
# fractions <- unique(as.character(connectivity$fraction))
# nodes <- data.frame(node = c(1:249), 
#                     name = c(connectivity, "udw","susp","ss","fs"))
# 
# # create a links dataframe
# connectivity<-merge(connectivity,nodes,by.x="fraction",by.y="name")
# connectivity<-merge(connectivity,nodes, by.x="nifh",by.y="name")
# links<-connectivity[,c("node.x","node.y","nifh")]
# colnames(links)<-c("source","target","value")
#   
# # draw sankey network
# networkD3::sankeyNetwork(Links = links, Nodes = nodes, 
#                          Source = 'source', 
#                          Target = 'target', 
#                          Value = 'nifH', 
#                          NodeID = 'name',
#                          units = 'gene copies l-1')


# ggplot(data = connectivity,
#        aes(axis1 = fraction,   # First variable on the X-axis
#            axis2 = target, # Second variable on the X-axis
#            y = nifh)) +
#   geom_alluvium(aes(fill = target)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("fraction", "target"),
#                    expand = c(0.15, 0.05)) +
#   theme_void() 


# df <- connectivity3 %>%
#   make_long(fraction, target)

df2<-as.data.frame(aggregate(nifh ~ fraction +  target, connectivity3,sum))%>% filter(.,!is.na(nifh)) ###quito NA

udw<-subset(df2, fraction %in% c("udw"))
susp<-subset(df2, fraction %in% c("susp"))
ss<-subset(df2, fraction %in% c("ss"))
fs<-subset(df2, fraction %in% c("fs"))

bind<bind_rows(udw,susp,ss,fs)


ggplot(df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey() +
  theme_sankey(base_size = 16)




