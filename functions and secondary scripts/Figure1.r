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

invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

datf=fread("data/simulated/outputs_simulations/species_level_alleg.csv")
datf$type="mu"
datf$type[grep("mu2",datf$variable)]="mu2"
datf$type[grep("sd",datf$variable)]="sd"
datf$type[grep("sd2",datf$variable)]="sd2"

ncalc=10
indf=NULL
motifs=NULL

comp_vec=2
time_vec=unique(datf$time)
ess=1
competition=6
rich=20
ti=2000
######### NETWORK COEVOLVED WITH PHENO ONLY
tr="both"
bidon=subset(datf,essai==ess & time==ti & trait==tr & rich==rich & comp==competition)

#build interaction matrix:
m = matrix(NA, rich, rich)
for(i in 1:rich){
	mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
	if(tr=="both"){
		mu1_m=invlogit1(subset(bidon,type=="mu2" & species==i)$value)
		sd1_m=invlogit(subset(bidon,type=="sd2" & species==i)$value)
	}
	for(j in 1:rich){
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+rich))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+rich))$value)
		if(tr=="both"){
			mu2_m=invlogit1(subset(bidon,type=="mu2" & species==j+rich)$value)
			sd2_m=invlogit(subset(bidon,type=="sd2" & species==j+rich)$value)
		}
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		if(tr=="both"){
			f2=function(x){pmin(dnorm(x,mu1_m,sd1_m),dnorm(x,mu2_m,sd2_m))}
			m[j,i]=(sum(f(seq(0,365,0.1)))*0.1)*(sum(f2(seq(0,365,0.1)))*0.1)
		}else{
			m[j,i]=sum(f(seq(0,365,0.1)))*0.1
		}
	}  
}
m2=round(m,digits=2)
df=melt(m2)
df$Var2=paste0("H",df$Var2)
g=graph_from_data_frame(df,directed =FALSE)
V(g)$type=TRUE
V(g)$type[grep("H",V(g)$name)]=FALSE
vec=rep(NA,length(vertex_attr(g,"name")))
vec[]="forestgreen"
vec[grep("H",vertex_attr(g,"name"))]="dodgerblue3"
vertex_attr(g,"color")=vec
edge_attr(g,"weight")=df$value
E(g)$width <- (E(g)$weight^1)*10


edge_attr(g,"color")=sapply(E(g)$weight,function(x){adjustcolor(paste0("gray",round(100-100*x/max(E(g)$weight))),alpha.f =x/max(E(g)$weight))})
size=3
net=plot_grid(base2grob(~plot(g,layout = layout_as_bipartite,vertex.label=NA,vertex.size=size,label.cex=0,size2=size,curved=T)))+
theme(plot.margin = unit(c(0,0,0,0), "cm"),plot.title=element_text(size=14,face="bold"))+ggtitle("b")


############# V MOTIFS:
colo=c("forestgreen","dodgerblue3","cyan","magenta4")

preci=1000
dat=data.frame(species=rep(c("plant","poll1","poll2"),each=preci),x=seq(0,365,length.out=preci),mu=rep(c(190,150,230),each=preci),sd=rep(c(30,20,20),each=preci))
dat$abund=dnorm(dat$x,dat$mu,dat$sd)

overlaps=dcast(dat,x~species,value.var="abund")
overlaps$mut1=apply(overlaps[,c("plant","poll1")],1,min)
overlaps$mut2=apply(overlaps[,c("plant","poll2")],1,min)
overlaps$comp=apply(overlaps[,c("poll2","poll1")],1,min)
sum(overlaps$comp)>sum(overlaps$mut2)*sum(overlaps$mut1)

pl1=ggplot()+
geom_polygon(alpha=0.4,data=dat,aes(x=x,y=abund,color=species,fill=species))+
geom_polygon(data=overlaps,aes(x=x,y=mut1),color=colo[2],fill=colo[2])+
geom_polygon(data=overlaps,aes(x=x,y=mut2),color=colo[3],fill=colo[3])+
geom_polygon(data=overlaps,aes(x=x,y=comp),color=colo[4],fill=colo[4])+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
ggtitle("c")+ylim(c(0,NA))+
ylab("Abundance")+
xlab("Day of the year")+
scale_colour_manual(values=colo[1:3])+scale_fill_manual(values=colo[1:3])

preci=1000
dat=data.frame(species=rep(c("plant","poll1","poll2"),each=preci),x=seq(0,365,length.out=preci),mu=rep(c(190,180,200),each=preci),sd=rep(c(30,20,20),each=preci))
dat$abund=dnorm(dat$x,dat$mu,dat$sd)

overlaps=dcast(dat,x~species,value.var="abund")
overlaps$mut1=apply(overlaps[,c("plant","poll1")],1,min)
overlaps$mut2=apply(overlaps[,c("plant","poll2")],1,min)
overlaps$comp=apply(overlaps[,c("poll2","poll1")],1,min)
sum(overlaps$comp)>sum(overlaps$mut2)*sum(overlaps$mut1)

pl2=ggplot()+
geom_polygon(alpha=0.4,data=dat,aes(x=x,y=abund,color=species,fill=species))+
geom_polygon(data=overlaps,aes(x=x,y=mut1),color=colo[2],fill=colo[2])+
geom_polygon(data=overlaps,aes(x=x,y=mut2),color=colo[3],fill=colo[3])+
geom_polygon(data=overlaps,aes(x=x,y=comp),color=colo[4],fill=colo[4])+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
ggtitle("")+ylim(c(0,NA))+
ylab("Abundance")+
xlab("Day of the year")+
scale_colour_manual(values=colo[1:3])+scale_fill_manual(values=colo[1:3])


blank=ggplot()+
theme_bw()+theme(axis.line = element_blank(),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),axis.title=element_blank(),axis.text=element_blank(),legend.position="none",strip.background=element_rect(fill=NA,color=NA))

motifs=plot_grid(pl1,pl2,ncol=1,align="hv")

f1top=plot_grid(blank,net,motifs,blank,ncol=2,align="hv")

pdf("Figures/Fig.1top.pdf",width=8,height=6)
f1top
dev.off();

save(f1top,file="data/top_figure_1.RData")