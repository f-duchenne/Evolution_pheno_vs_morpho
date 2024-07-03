###################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales","piecewiseSEM",
"igraph","qgraph","car","nlme","viridis") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

colo=c("dodgerblue4","chartreuse3")
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

##################################### SPECIES LEVEL ########################
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
datf=fread("species_level_data.csv")
endpo=subset(datf,time==2000 & comp==10)


datf2=subset(datf,comp==10 & time %in% plyr::round_any(seq(sqrt(0),sqrt(2000),length.out=20)^2,10)) %>% group_by(type,species,comp,essai,rich,trait) %>% mutate(delta=sqrt((value-lag(value))^2)/(time-lag(time)))

s7=ggplot(subset(datf2,time>0),aes(x=time,y=delta,color=trait))+
stat_smooth()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+scale_color_manual(values=colo)+ylab("Changes in trait parameters")+
facet_grid(cols=vars(rich),rows=vars(type),labeller = label_bquote(cols=n[sp] == .(rich)))

png("Fig.S7.png",width=1200,height=800,res=150)
s7
dev.off();

sp1=ggplot()+
geom_density(data=subset(endpo,type=="mu"),aes(x=invlogit1(value),color=trait,fill=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("density")+
xlab("Mean value of the trait")+ggtitle("a")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_grid(. ~rich, labeller = label_bquote(cols=n[sp] == .(rich)), scales = "free_y")

sp2=ggplot()+
geom_density(data=subset(endpo,type=="sd"),aes(x=invlogit(value),color=trait,fill=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("density")+
xlab("Standard deviation of the trait")+ggtitle("b")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_wrap( ~rich, labeller = label_bquote(cols=n[sp] == .(rich)), scales = "free_y")

plot_grid(sp1,sp2,ncol=1,align="hv")

png("Fig.S2.png",width=1200,height=800,res=150)
plot_grid(sp1,sp2,ncol=1,align="hv")
dev.off();

##################################### COMMUNITY LEVEL ########################
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
indf=fread("C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info.csv")
names(indf)=gsub("weighted ","",names(indf))
names(indf)[3]="inter. evenness"
names(indf)[5]="wNODF"

indf$prop_motif=(indf$motifn/indf$motifntot)
model=glmmTMB::glmmTMB(feas~(prop_motif+rich)*trait+(1|essai),family=beta_family(link="logit"),data=subset(indf,competition==10 & rho==0.2 & time==2000))

indf2=subset(indf,competition==10 & rho==0.01)
acp=PCA(indf2[,1:6],graph=FALSE)
indf2=cbind(indf2,acp$ind$coord)


indfspeed=indf2 %>% group_by(trait,competition,essai,rich) %>% mutate(speed=sqrt((Dim.1-lag(Dim.1))^2+(Dim.2-lag(Dim.2))^2+(Dim.3-lag(Dim.3))^2)/(time-lag(time)))

s8=ggplot(subset(indfspeed,time>0),aes(x=time,y=speed,color=trait))+geom_line(alpha=0.2,aes(group=paste(essai,competition,trait,rich)))+
stat_smooth()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+scale_color_manual(values=colo)+ylab("Changes in network structure")

png("Fig.S8.png",width=1200,height=600,res=150)
s8
dev.off();


pl1=ggplot(data=subset(indf2,rich==10),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,group=paste(essai,trait)))+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point()+
geom_line(show_guide =F)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("a")+coord_fixed(ratio=1)

pl2=fviz_pca_var(acp, col.var = "black", repel =TRUE)+ggtitle("b")+
theme(plot.title=element_text(size=14,face="bold",hjust = 0))

pl3=ggplot()+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point(data=subset(indf2,time==2000),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,shape=as.factor(rich)))+
geom_point(data=subset(indf2,time==0),aes(x=Dim.1,y=Dim.2,alpha=time,shape=as.factor(rich)),color="grey",alpha=0.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("c")

pl4=ggplot(data=subset(indf2,time==2000),aes(x=as.factor(rich),y=wNODF,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(N[A] == N[P]))+ylab("Nestedness (wNODF)")+xlab(expression(paste(n[sp], " / guild")))+
ggtitle("d")

pl5=ggplot(data=subset(indf2,time==2000),aes(x=as.factor(rich),y=connectance,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(N[A] == N[P]))+ylab("Connectance")+xlab(expression(paste(n[sp], " / guild")))+ggtitle("e")

top=plot_grid(pl1,pl2,align="hv",ncol=2,rel_widths=c(1.5,1))
bottom=plot_grid(pl3,pl4,pl5,ncol=3,rel_widths=c(2,1,1))

pdf("Fig.2.pdf",width=7,height=7)
grid.arrange(top,bottom,ncol=1)
dev.off();

#################### FIG. S3
indf3=subset(indf,competition==5 & rho==0.01)
indf3=cbind(indf3,predict(acp,newdata=indf3[,1:6])$coord)

pl3=ggplot()+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point(data=subset(indf3,time==2000),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,shape=as.factor(rich)))+
geom_point(data=subset(indf3,time==0),aes(x=Dim.1,y=Dim.2,alpha=time,shape=as.factor(rich)),color="grey",alpha=0.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("a")

pl4=ggplot(data=subset(indf3,time==2000),aes(x=as.factor(rich),y=wNODF,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(N[A] == N[P]))+ylab("Nestedness (wNODF)")+xlab(expression(paste(n[sp], " / guild")))+
ggtitle("b")


pl5=ggplot(data=subset(indf3,time==2000),aes(x=as.factor(rich),y=connectance,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(N[A] == N[P]))+ylab("Connectance")+xlab(expression(paste(n[sp], " / guild")))+ggtitle("c")

bottom=plot_grid(pl3,pl4,pl5,ncol=3,rel_widths=c(2,1,1))

png("Fig.S3.png",width=1200,height=600,res=150)
bottom
dev.off();

#########################################

indf4=subset(indf,competition==10)

pl6=ggplot(data=subset(indf4,rho==0.05),aes(x=time,y=feas,color=trait,group=paste(time,trait)))+geom_violin(position=position_dodge(width=0.5),width=8,scale="width",show.legend=FALSE)+
geom_boxplot(width=3,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.background = element_blank(),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_sqrt(breaks=c(0,2000))

png("Fig.S4.png",width=1200,height=600,res=150)
pl6
dev.off();

pl6b=ggplot(data=subset(indf4,rho==0.4),aes(x=time,y=feas,color=trait,group=paste(time,trait)))+geom_violin(position=position_dodge(width=0.5),width=8,scale="width",show.legend=FALSE)+
geom_boxplot(width=3,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.background = element_blank(),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_sqrt(breaks=c(0,2000))

plot_grid(pl6,pl6b,ncol=1,align="hv")

pdf("Fig.3.pdf",width=6,height=3)
pl6b
dev.off();

############################## FIGURE S5
indf4=subset(indf,competition==5)

pl6=ggplot(data=subset(indf4,rho==0.05),aes(x=time,y=feas,color=trait,group=paste(time,trait)))+geom_violin(position=position_dodge(width=0.5),width=8,scale="width",show.legend=FALSE)+
geom_boxplot(width=3,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.background = element_blank(),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_sqrt(breaks=c(0,2000))+ggtitle("a",subtitle="low interspecific competition")

pl6b=ggplot(data=subset(indf4,rho==0.4),aes(x=time,y=feas,color=trait,group=paste(time,trait)))+geom_violin(position=position_dodge(width=0.5),width=8,scale="width",show.legend=FALSE)+
geom_boxplot(width=3,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_sqrt(breaks=c(0,2000))+ggtitle("b",subtitle="high interspecific competition")

plot_grid(pl6,pl6b,ncol=1,align="hv")
png("Fig.S5.png",width=1200,height=1000,res=150)
plot_grid(pl6,pl6b,ncol=1,align="hv")
dev.off();

############################# FIGURE 4
indf5=subset(indf4,time %in% c(0,2000))

pl7=ggplot(data=subset(indf5,time==0),aes(x=rich,y=feas,color=as.factor(rho)))+
stat_summary(fun.y=mean, geom="line", size = 0.5,position=position_dodge(width=0.5))+
stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - sd(x), fun.ymax = function(x) mean(x) + sd(x), geom = "pointrange")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Number of species per guild")+facet_wrap(~trait)+scale_colour_viridis(option="rocket",discrete=TRUE,direction=-1)+
#scale_color_brewer(palette = "Reds")+
labs(colour=expression(rho))+ggtitle("a",subtitle="Initial networks")

pl8=ggplot(data=subset(indf5,time==2000),aes(x=rich,y=feas,color=as.factor(rho)))+
stat_summary(fun.y=mean, geom="line", size = 0.5,position=position_dodge(width=0.5))+
stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - sd(x), fun.ymax = function(x) mean(x) + sd(x), geom = "pointrange")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Number of species per guild")+facet_wrap(~trait)+scale_colour_viridis(option="rocket",discrete=TRUE,direction=-1)+
#scale_color_brewer(palette = "Reds")+
labs(colour=expression(rho))+ggtitle("b",subtitle="Coevolved networks")

pdf("Fig.4.pdf",width=6,height=6)
plot_grid(pl7,pl8,ncol=1,align="hv")
dev.off();

indf %>% group_by(rich,trait,time) %>% summarise(mean(feas))



indf5=subset(indf4,time %in% c(0,2000))

pl7=ggplot(data=subset(indf5,time==0),aes(x=rich,y=feas,color=trait,alpha=as.factor(rho)))+
stat_summary(fun.y=mean, geom="line", size = 0.5,position=position_dodge(width=0.5))+
stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - sd(x), fun.ymax = function(x) mean(x) + sd(x), geom = "pointrange")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Number of species per guild")+facet_wrap(~trait)+scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
#scale_color_brewer(palette = "Reds")+
labs(alpha=expression(rho))+ggtitle("a",subtitle="Initial networks")

pl8=ggplot(data=subset(indf5,time==2000),aes(x=rich,y=feas,color=trait,alpha=as.factor(rho)))+
stat_summary(fun.y=mean, geom="line", size = 0.5,position=position_dodge(width=0.5))+
stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - sd(x), fun.ymax = function(x) mean(x) + sd(x), geom = "pointrange")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Number of species per guild")+facet_wrap(~trait)+scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
#scale_color_brewer(palette = "Reds")+
labs(alpha=expression(rho))+ggtitle("b",subtitle="Coevolved networks")

plot_grid(pl7,pl8,ncol=1,align="hv")