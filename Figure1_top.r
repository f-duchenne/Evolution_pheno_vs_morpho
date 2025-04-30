###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","ggplot2","bipartite","gridExtra","cowplot","ggpubr","scales","viridis","ggthemes") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)


preci=1000
dat=data.frame(species=rep(c("plant","poll1","poll2"),each=preci),x=seq(0,365,length.out=preci),mu=rep(c(190,150,230),each=preci),sd=rep(c(30,20,20),each=preci))
dat$abund=dnorm(dat$x,dat$mu,dat$sd)

overlaps=dcast(dat,x~species,value.var="abund")
overlaps$mut1=apply(overlaps[,c("plant","poll1")],1,min)
overlaps$mut2=apply(overlaps[,c("plant","poll2")],1,min)
overlaps$comp=apply(overlaps[,c("poll2","poll1")],1,min)
sum(overlaps$comp*0.3653654)>sum(overlaps$mut2*0.3653654)*sum(overlaps$mut1*0.3653654)

pl1=ggplot()+
geom_polygon(alpha=0.4,data=dat,aes(x=x,y=abund,color=species,fill=species))+
geom_polygon(data=overlaps,aes(x=x,y=mut1),color="#E69F00",fill="#E69F00")+
geom_polygon(data=overlaps,aes(x=x,y=mut2),color="#56B4E9",fill="#56B4E9")+
geom_polygon(data=overlaps,aes(x=x,y=comp),color="#E83F3FFF",fill="#E83F3FFF")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
ggtitle("a")+ylim(c(0,NA))+
ylab("Abundance")+
xlab("Day of the year")+
scale_colour_colorblind()+scale_fill_colorblind()

preci=1000
dat=data.frame(species=rep(c("plant","poll1","poll2"),each=preci),x=seq(0,365,length.out=preci),mu=rep(c(190,180,200),each=preci),sd=rep(c(30,20,20),each=preci))
dat$abund=dnorm(dat$x,dat$mu,dat$sd)

overlaps=dcast(dat,x~species,value.var="abund")
overlaps$mut1=apply(overlaps[,c("plant","poll1")],1,min)
overlaps$mut2=apply(overlaps[,c("plant","poll2")],1,min)
overlaps$comp=apply(overlaps[,c("poll2","poll1")],1,min)
sum(overlaps$comp*0.3653654)>sum(overlaps$mut2*0.3653654)*sum(overlaps$mut1*0.3653654)

pl2=ggplot()+
geom_polygon(alpha=0.4,data=dat,aes(x=x,y=abund,color=species,fill=species))+
geom_polygon(data=overlaps,aes(x=x,y=mut1),color="#E69F00",fill="#E69F00")+
geom_polygon(data=overlaps,aes(x=x,y=mut2),color="#56B4E9",fill="#56B4E9")+
geom_polygon(data=overlaps,aes(x=x,y=comp),color="#E83F3FFF",fill="#E83F3FFF")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
ggtitle("")+ylim(c(0,NA))+
ylab("Abundance")+
xlab("Day of the year")+
scale_colour_colorblind()+scale_fill_colorblind()


blank=ggplot()+
theme_bw()+theme(axis.line = element_blank(),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),axis.title=element_blank(),axis.text=element_blank(),legend.position="none",strip.background=element_rect(fill=NA,color=NA))


f1top=plot_grid(pl1,blank,pl2,blank,ncol=2,align="hv")

save(f1top,file="data/top_figure_1.RData")