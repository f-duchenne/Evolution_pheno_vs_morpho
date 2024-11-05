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
ggtitle("e")+ylim(c(0,NA))+
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


bottomf2=plot_grid(pl1,blank,pl2,blank,ncol=2,align="hv")

save(bottomf2,file="bottom_figure_2.RData")

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
indf=fread("C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info_symmetric.csv")
names(indf)=gsub("weighted ","w",names(indf))
names(indf)[3]="inter. evenness"
names(indf)[5]="wNODF"
names(indf)[7]="mut_strength"
colo=c("chartreuse3","dodgerblue4","gold3")

indf2=subset(indf,competition==5 & rho==0.01 & competition_feas==competition)

p1=ggplot()+
stat_smooth(data=indf2,aes(y=mut_strength,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Average mutualism strength")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
xlab("Time")+ggtitle("b")+labs(colour="Simulations with:",fill="Simulations with:",linetype=expression(paste(n[sp]," / guild")))+scale_x_continuous(breaks=c(0,1000,2000))


indf2$comp=apply(indf2[,c("comp_a","comp_p")],1,mean)
p2=ggplot()+
stat_smooth(data=indf2,aes(y=comp,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="bottom",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Average competition stength")+guides(linetype=guide_legend(title.position="top",title.hjust = 0.5,override.aes=list(fill=NA,color="black")),colour=guide_legend(title.position="top",title.hjust = 0.5))+
xlab("Time")+ggtitle("c")+labs(colour="Simulations with:",fill="Simulations with:",linetype=expression(paste(n[sp]," / guild")))+scale_x_continuous(breaks=c(0,1000,2000))

p3=ggplot()+
stat_smooth(data=indf2,aes(y=mut_strength/(competition*comp),x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("mutualism/competition")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
xlab("Time")+ggtitle("d")+labs(colour="Simulations with:",fill="Simulations with:",linetype=expression(paste(n[sp]," / guild")))+scale_x_continuous(breaks=c(0,1000,2000))


p4=ggplot(data=subset(indf2,trait!="morpho"),aes(y=motifn/motifntot,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
stat_smooth()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"),legend.position="none")+
scale_color_manual(values=colo[-2])+scale_fill_manual(values=colo[-2])+ylab('Proportion of "V+" motifs')+
xlab("Time")+ggtitle("e")+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+scale_x_continuous(breaks=c(0,1000,2000))



leg <- ggpubr::as_ggplot(cowplot::get_legend(p2))
p2=p2+theme(legend.position="none")

bottom=plot_grid(p1,p2,p3,p4,ncol=4,align="hv",rel_widths=c(1,1,1,1,0.3))

grid.arrange(top,bottom,leg,ncol=1,heights=c(1.5,1,0.3))
pdf("Fig.1.pdf",width=8,height=8)
grid.arrange(top,bottom,leg,ncol=1,heights=c(1.5,1,0.3))
dev.off();

##################### FIGURE S1

indf3=subset(indf,competition==2  & rho==0.01 & competition_feas==competition)

p1=p1=ggplot()+
stat_smooth(data=indf3,aes(y=mut_strength,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Average mutualism strength")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
xlab("Time")+ggtitle("b")+labs(colour="Simulations with:",fill="Simulations with:",linetype=expression(paste(n[sp]," / guild")))+scale_x_continuous(breaks=c(0,1000,2000))


indf3$comp=apply(indf3[,c("comp_a","comp_p")],1,mean)
p2=ggplot()+
stat_smooth(data=indf3,aes(y=comp,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="bottom",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Average competition stength")+guides(linetype=guide_legend(title.position="top",title.hjust = 0.5,override.aes=list(fill=NA,color="black")),colour=guide_legend(title.position="top",title.hjust = 0.5))+
xlab("Time")+ggtitle("c")+labs(colour="Simulations with:",fill="Simulations with:",linetype=expression(paste(n[sp]," / guild")))+scale_x_continuous(breaks=c(0,1000,2000))

p3=ggplot()+
stat_smooth(data=indf3,aes(y=mut_strength/(comp*competition),x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("mutualism/competition")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
xlab("Time")+ggtitle("d")+labs(colour="Simulations with:",fill="Simulations with:",linetype=expression(paste(n[sp]," / guild")))+scale_x_continuous(breaks=c(0,1000,2000))


p4=ggplot(data=subset(indf3,trait!="morpho"),aes(y=motifn/motifntot,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
stat_smooth()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"),legend.position="none")+
scale_color_manual(values=colo[-2])+scale_fill_manual(values=colo[-2])+ylab('Proportion of "V+" motifs')+
xlab("Time")+ggtitle("e")+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+scale_x_continuous(breaks=c(0,1000,2000))



leg <- ggpubr::as_ggplot(cowplot::get_legend(p2))
p2=p2+theme(legend.position="none")

bottom=plot_grid(p1,p2,p3,p4,ncol=4,align="hv",rel_widths=c(1,1,1,1,0.3))

grid.arrange(bottom,leg,ncol=1,heights=c(1,0.3))

png("Fig.S1.png",width=1200,height=600,res=150)
grid.arrange(bottom,leg,ncol=1,heights=c(1,0.3))
dev.off();