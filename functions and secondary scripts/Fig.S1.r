pkgs <- c("ggplot2","gridExtra") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])

path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)

invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

tab=data.frame(trait=rep(c("mu","sd"),each=1000),val_lat=seq(-10,10,length.out=1000))
tab$value=c(sapply(tab$val_lat[1:1000],invlogit1),sapply(tab$val_lat[1001:2000],invlogit))

pl1=ggplot(data=tab,aes(x=val_lat,y=value))+geom_line()+theme_bw()+theme(panel.grid=element_blank(),strip.background=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+xlab("Trait parameter value on latent scale")+ylab("Trait parameter value")+facet_wrap(~trait,scales="free")+ggtitle("a")

mut=function(x){invlogit1(x+0.01)-invlogit1(x)}
mut2=function(x){invlogit(x+0.01)-invlogit(x)}

tab$mut_ef=c(sapply(tab$val_lat[1:1000],mut),sapply(tab$val_lat[1001:2000],mut2))

pl2=ggplot(data=tab,aes(x=value,y=mut_ef))+geom_line()+theme_bw()+theme(panel.grid=element_blank(),strip.background=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0))+xlab("Trait parameter value")+ylab("Evolution speed")+facet_wrap(~trait,scales="free")+ggtitle("b")


grid.arrange(pl1,pl2)

png("Figures/Fig_S1.png",width=1100,height=1000,res=140)
grid.arrange(pl1,pl2)
dev.off();