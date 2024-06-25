###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales","viridis","ggthemes") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")

preci=1000
dat=data.frame(species=rep(c("plant","poll1","poll2"),each=preci),x=seq(0,365,length.out=preci),mu=rep(c(190,150,220),each=preci),sd=rep(c(30,20,20),each=preci))
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


top=plot_grid(pl1,blank,pl2,blank,ncol=2,align="hv")

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
indf=fread("C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info.csv")
colo=c("dodgerblue4","chartreuse3")

indf2=subset(indf,competition==10 & rho==0.01)

p1=ggplot(data=indf2,aes(y=motifn/motifntot,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
#geom_boxplot(data=indf,aes(y=motifn/motifntot,x=time,group=paste0(trait,time,rich),color=trait,linetype=as.factor(rich)),width=100,show_guide=FALSE)+
stat_smooth()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("Proportion of 'V' motifs promoting facilitation")+
xlab("Time")+ggtitle("b")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))

p2=ggplot()+
#geom_boxplot(data=indf,aes(y=comp,x=time,group=paste0(trait,time,rich),color=trait,linetype=as.factor(rich)),width=100)+
stat_smooth(data=indf2,aes(y=comp,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("Competition for mutualistic partners")+
xlab("Time")+ggtitle("c")


leg <- ggpubr::as_ggplot(cowplot::get_legend(p1))
p1=p1+theme(legend.position="none")

bottom=plot_grid(p1,p2,leg,ncol=3,align="hv",rel_widths=c(1,1,0.3))

grid.arrange(top,bottom,ncol=1,heights=c(2,1))
pdf("Fig.1.pdf",width=9,height=9)
grid.arrange(top,bottom,ncol=1,heights=c(2,1))
dev.off();

##################### FIGURE S1

indf3=subset(indf,competition==5  & rho==0.01)

p1=ggplot()+
#geom_boxplot(data=indf,aes(y=motifn/motifntot,x=time,group=paste0(trait,time,rich),color=trait,linetype=as.factor(rich)),width=100,show_guide=FALSE)+
stat_smooth(data=indf3,aes(y=motifn/motifntot,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("Proportion of 'V' motifs promoting facilitation")+
xlab("Time")+ggtitle("a")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))

p2=ggplot()+
#geom_boxplot(data=indf,aes(y=comp,x=time,group=paste0(trait,time,rich),color=trait,linetype=as.factor(rich)),width=100)+
stat_smooth(data=indf3,aes(y=comp,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("Competition for mutualistic partners")+
xlab("Time")+ggtitle("b")

leg <- ggpubr::as_ggplot(cowplot::get_legend(p1))
p1=p1+theme(legend.position="none")

bottom=plot_grid(p1,p2,leg,ncol=3,align="hv",rel_widths=c(1,1,0.3))

png("Fig.S1.png",width=1200,height=600,res=150)
bottom
dev.off();