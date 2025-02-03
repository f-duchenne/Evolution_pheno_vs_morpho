###################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales","piecewiseSEM",
"igraph","qgraph","car","nlme","viridis","ggeffects","emmeans") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

colo=c("chartreuse3","dodgerblue4","gold3")
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

##################################### SPECIES LEVEL ########################
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
datf=fread("species_level_data_symmetric.csv")
datf$type2="pheno"
datf$type2[datf$type %in% c("mu2","sd2")]="morpho"
datf$type2[datf$trait %in% c("morpho")]="morpho"
datf$type3="mu"
datf$type3[datf$type %in% c("sd","sd2")]="sd"

datf2=subset(datf,comp==5 & time %in% plyr::round_any(seq(sqrt(0),sqrt(2000),length.out=30)^2,10))
datf2$value[datf2$type %in% c("mu","mu2")]=invlogit1(datf2$value[datf2$type %in% c("mu","mu2")])
datf2$value[datf2$type %in% c("sd","sd2")]=invlogit(datf2$value[datf2$type %in% c("sd","sd2")])
datf2=datf2 %>% group_by(type,type2,type3,species,comp,essai,rich,trait) %>% mutate(delta=sqrt((value-lag(value))^2)/(time-lag(time)))

s7_a=ggplot(subset(datf2,time>0 & type3=="mu"),aes(x=time,y=delta,color=trait,linetype=type2))+
stat_smooth()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="right",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Changes in trait parameters")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
facet_grid(cols=vars(rich),labeller = label_bquote(cols=n[sp] == .(rich)))+labs(colour="Simulations with:",linetype="trait")+ggtitle("a",subtitle="Mean of the trait(s)")

s7_b=ggplot(subset(datf2,time>0 & type3=="sd"),aes(x=time,y=delta,color=trait,linetype=type2))+
stat_smooth()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="right",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Changes in trait parameters")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
facet_grid(cols=vars(rich),labeller = label_bquote(cols=n[sp] == .(rich)))+labs(colour="Simulations with:",linetype="trait")+ggtitle("a",subtitle="Standard deviation of the trait(s)")

png("Fig.S6.png",width=1300,height=1100,res=150)
plot_grid(s7_a,s7_b,ncol=1,align="hv")
dev.off();


############### FIGURE S2
endpo=subset(datf,time==2000 & comp==5)
sp1=ggplot()+
geom_density(data=subset(endpo,type3=="mu"),aes(x=invlogit1(value),color=trait,fill=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("density")+
xlab("Mean value of the trait")+ggtitle("a")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_grid(cols=vars(rich),rows=vars(type2), labeller = label_bquote(cols=n[sp] == .(rich)), scales = "free_y")+
labs(colour="Simulations with:",fill="Simulations with:")

sp2=ggplot()+
geom_density(data=subset(endpo,type3=="sd"),aes(x=invlogit(value),color=trait,fill=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("density")+
xlab("Standard deviation of the trait")+ggtitle("b")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_grid(cols=vars(rich),rows=vars(type2), labeller = label_bquote(cols=n[sp] == .(rich)), scales = "free_y")+
labs(colour="Simulations with:",fill="Simulations with:")

plot_grid(sp1,sp2,ncol=1,align="hv")

png("Fig.S2.png",width=1200,height=1000,res=150)
plot_grid(sp1,sp2,ncol=1,align="hv")
dev.off();

##################################### COMMUNITY LEVEL ########################
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
indf=fread("C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info_symmetric.csv")
names(indf)=gsub("weighted ","w",names(indf))
names(indf)[3]="inter. evenness"
names(indf)[5]="wNODF"
names(indf)[7]="mut.strength"
indf$prop_motif=(indf$motifn/indf$motifntot)

indf2=subset(indf,rho==0.01 & competition_feas==competition)
acp=PCA(indf2[,1:7],graph=FALSE)
indf2=cbind(indf2,acp$ind$coord)
indf3=subset(indf2,competition==5)
indf3$comp=apply(indf3[,c("comp_a","comp_p")],1,mean)

indfspeed=indf3 %>% group_by(trait,competition,essai,rich) %>% mutate(speed=sqrt((Dim.1-lag(Dim.1))^2+(Dim.2-lag(Dim.2))^2+(Dim.3-lag(Dim.3))^2)/(time-lag(time)))

s8=ggplot(subset(indfspeed,time>0),aes(x=time,y=speed,color=trait))+geom_line(alpha=0.2,aes(group=paste(essai,competition,trait,rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Changes in network structure")+labs(colour="Simulations with:",fill="Simulations with:")+
guides(linetype=guide_legend(override.aes=list(fill=NA)))

png("Fig.S6.png",width=1200,height=600,res=150)
s8
dev.off();


############### FIGURE 1
pl1=ggplot(data=subset(indf2,rich==30 & time %in% c(0,20,100,220,2000)),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,group=paste(essai,trait)))+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_path(show_guide=FALSE,alpha=0.3)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
geom_point(aes(alpha=time))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ggtitle("a")+labs(colour="Simulations with:",fill="Simulations with:")+
facet_wrap(~competition,ncol=1,labeller = label_bquote(rows=italic(c) == .(competition)))+scale_alpha_continuous(limits = c(0,2000), breaks = c(0,2000),range = c(0.05, 1))

pl2=fviz_pca_var(acp, col.var = "black", repel =TRUE,labelsize = 3)+ggtitle("b")+
theme(plot.title=element_text(size=14,face="bold",hjust = 0))

pl3=ggplot()+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point(data=subset(indf3,time==2000),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,fill=trait,shape=as.factor(rich)))+
#geom_point(data=subset(indf3,time==0 & trait!="both"),aes(x=Dim.1,y=Dim.2,alpha=time,shape=as.factor(rich)),color="grey",fill="grey",alpha=0.3)+
#geom_point(data=subset(indf3,time==0 & trait=="both"),aes(x=Dim.1,y=Dim.2,alpha=time,shape=as.factor(rich)),color="black",fill="grey",alpha=0.3)+
scale_shape_manual(values = c(21,22,23))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ggtitle("c")+labs(colour="Simulations with:",fill="Simulations with:")


pl4=ggplot(data=subset(indf3,time==2000),aes(x=as.factor(rich),y=wNODF,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(N[A] == N[P]))+ylab("wNODF")+xlab(expression(paste(n[sp], " / guild")))+
ggtitle("d")

pl5=ggplot(data=subset(indf3,time==2000),aes(x=as.factor(rich),y=modularity,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(N[A] == N[P]))+ylab("Modularity")+xlab(expression(paste(n[sp], " / guild")))+ggtitle("e")

top=plot_grid(pl1,align="hv",ncol=1,rel_widths=c(1))
midle=plot_grid(pl2,pl3,ncol=1)
bottom=plot_grid(pl4,pl5,ncol=2)

lay=rbind(c(1,2),c(1,2),c(3,3))

grid.arrange(top,midle,bottom,widths=c(2,1),heights=c(1,1,0.7),layout_matrix = lay)

pdf("Fig.1.pdf",width=8,height=7)
grid.arrange(top,midle,bottom,widths=c(2,1),heights=c(1,1,0.7),layout_matrix = lay)
dev.off();

indf3b=subset(indf2,competition==2)
pl3b=ggplot()+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point(data=subset(indf3b,time==2000),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,fill=trait,shape=as.factor(rich)))+
#geom_point(data=subset(indf3,time==0 & trait!="both"),aes(x=Dim.1,y=Dim.2,alpha=time,shape=as.factor(rich)),color="grey",fill="grey",alpha=0.3)+
#geom_point(data=subset(indf3,time==0 & trait=="both"),aes(x=Dim.1,y=Dim.2,alpha=time,shape=as.factor(rich)),color="black",fill="grey",alpha=0.3)+
scale_shape_manual(values = c(21,22,23))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ggtitle("c")+labs(colour="Simulations with:",fill="Simulations with:")

png("Fig.S1.png",width=1200,height=1000,res=150)
pl3b
dev.off();


############### FIGURE 2
p1=ggplot()+
stat_smooth(data=indf3,aes(y=mut.strength,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Average mutualism strength")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
xlab("Time")+ggtitle("a")+labs(colour="Simulations with:",fill="Simulations with:",linetype=expression(paste(n[sp]," / guild")))+scale_x_continuous(breaks=c(0,1000,2000))

p2=ggplot()+
stat_smooth(data=indf3,aes(y=comp,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="bottom",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Average competition stength")+guides(linetype=guide_legend(title.position="top",title.hjust = 0.5,override.aes=list(fill=NA,color="black")),colour=guide_legend(title.position="top",title.hjust = 0.5))+
xlab("Time")+ggtitle("b")+labs(colour="Simulations with:",fill="Simulations with:",linetype=expression(paste(n[sp]," / guild")))+scale_x_continuous(breaks=c(0,1000,2000))

p3=ggplot()+
stat_smooth(data=indf3,aes(y=mut.strength/(competition*comp),x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("mutualism/competition")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
xlab("Time")+ggtitle("c")+labs(colour="Simulations with:",fill="Simulations with:")+scale_x_continuous(breaks=c(0,1000,2000))

indf3$motifn[indf3$trait=="morpho"]=0
p4=ggplot(data=subset(indf3),aes(y=motifn/motifntot,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - 1.96*sd(x)/length(x), fun.ymax = function(x) mean(x) + 1.96*sd(x)/length(x), geom = "ribbon",alpha=0.1,colour=NA)+
stat_summary(fun.y=mean, geom="line", size = 0.5,linewidth=1.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab('Proportion of "V+" motifs')+
xlab("Time")+ggtitle("d")+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+scale_x_continuous(breaks=c(0,1000,2000))

leg <- ggpubr::as_ggplot(cowplot::get_legend(p2))
p2=p2+theme(legend.position="none")
top=plot_grid(p1,p2,p3,p4,align="hv",ncol=4)
load("bottom_figure_2.RData")
grid.arrange(top,leg,bottomf2,ncol=1,heights=c(1,0.3,1))

pdf("Fig.2.pdf",width=7,height=6)
grid.arrange(top,leg,bottomf2,ncol=1,heights=c(1,0.3,1))
dev.off();

######################################### FIGURE 3
indf4=subset(indf,competition==5 & competition_feas==competition)

ggplot(data=subset(indf,rho==0.4 & competition_feas==competition),aes(x=time,y=feas_structure,color=trait,group=paste(time,trait)))+geom_violin(position=position_dodge(width=0.5),width=8,scale="width",show.legend=FALSE)+
geom_boxplot(width=3,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.background = element_blank(),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Time")+facet_grid(cols=vars(rich),rows=vars(competition), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_sqrt(breaks=c(0,1000,2000))+labs(colour="Simulations with:",fill="Simulations with:")

pl6=ggplot(data=subset(indf4,rho==0.05),aes(x=time,y=feas,color=trait,group=paste(time,trait)))+geom_violin(position=position_dodge(width=0.5),width=8,scale="width",show.legend=FALSE)+
geom_boxplot(width=3,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.background = element_blank(),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_sqrt(breaks=c(0,1000,2000))+labs(colour="Simulations with:",fill="Simulations with:")


png("Fig.S3.png",width=1200,height=600,res=150)
pl6
dev.off();

pl6b=ggplot(data=subset(indf4,rho==0.4),aes(x=time,y=feas,color=trait,group=paste(time,trait)))+geom_violin(position=position_dodge(width=0.5),width=8,scale="width",show.legend=FALSE)+
geom_boxplot(width=3,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.background = element_blank(),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_sqrt(breaks=c(0,1000,2000))+labs(colour="Simulations with:",fill="Simulations with:")

plot_grid(pl6,pl6b,ncol=1,align="hv")

pdf("Fig.3.pdf",width=6,height=3)
pl6b
dev.off();

############################## FIGURE S4
indf4=subset(indf,competition==2 & competition_feas==competition)

pl6=ggplot(data=subset(indf4,rho==0.05),aes(x=time,y=feas,color=trait,group=paste(time,trait)))+geom_violin(position=position_dodge(width=0.5),width=8,scale="width",show.legend=FALSE)+
geom_boxplot(width=3,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.background = element_blank(),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_sqrt(breaks=c(0,1000,2000))+ggtitle("a",subtitle="low interspecific competition")+labs(colour="Simulations with:",fill="Simulations with:")


pl6b=ggplot(data=subset(indf4,rho==0.4),aes(x=time,y=feas,color=trait,group=paste(time,trait)))+geom_violin(position=position_dodge(width=0.5),width=8,scale="width",show.legend=FALSE)+
geom_boxplot(width=3,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_sqrt(breaks=c(0,1000,2000))+ggtitle("b",subtitle="high interspecific competition")+labs(colour="Simulations with:",fill="Simulations with:")

plot_grid(pl6,pl6b,ncol=1,align="hv")
png("Fig.S4.png",width=1200,height=1000,res=150)
plot_grid(pl6,pl6b,ncol=1,align="hv")
dev.off();

############################# FIGURE 4
indf5=subset(indf,time %in% c(0,2000) & competition_feas==competition & competition==5)
indf5$trait=factor(indf5$trait,labels=c("both traits","only morpho","only pheno"))
indf5=indf5 %>% group_by(trait,essai,rho,competition,competition_feas,time) %>% mutate(feas_ini=feas[rich==10])  
indf5$time2=factor(indf5$time,labels=c("initial (t = 0)","coevolved (t = 2000)"))
indf5$time2=factor(indf5$time2,levels=c("coevolved (t = 2000)","initial (t = 0)"))

model=glm(feas~rich*time2*trait,data=subset(indf5,rho==0.4),family=quasibinomial)

b=ggpredict(model,c("rich[10:30]","time2","trait"))
pl7a=ggplot()+
geom_line(data=b,aes(x=x,y=predicted,color=facet,linetype=group))+
geom_point(data=subset(indf5,rho==0.4),aes(x=rich,y=feas,color=trait,shape=as.factor(time2)),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(linetype="",shape="")+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Number of species per guild")+scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
#scale_color_brewer(palette = "Reds")+
labs(alpha=expression(rho),color="Simulations with:")+scale_x_continuous(breaks=c(10,20,30))+
ggtitle("a")

b2=as.data.frame(emmeans(model, pairwise ~ rich | time2 + trait)$emmeans)
b2$asymp.LCL[b2$time2=="initial (t = 0)" & b2$trait=="only morpho"]=-14.0735
b2$asymp.UCL[b2$time2=="initial (t = 0)" & b2$trait=="only morpho"]=-0.4922
b2$time=0
b2$time[b2$time2=="coevolved (t = 2000)"]=2000
b2$time=ordered(b2$time,levels=c(0,2000))

pl7b=ggplot(data=subset(b2,trait!="only morpho"),aes(x=time,y=emmean,color=trait))+
geom_pointrange(aes(ymin=asymp.LCL,ymax=asymp.UCL,shape=as.factor(time2)))+
geom_line(aes(group=trait))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(linetype="",shape="")+ylab("Slope of the\ndiversity-stability trade off")+xlab("")+scale_color_manual(values=colo[-2],guide=FALSE)+scale_fill_manual(values=colo)+
ggtitle("b")+scale_x_discrete(breaks=c(0,2000),labels=c("initial","coevolved"))

plot_grid(pl7a,pl7b,ncol=2,align="hv")

pdf("Fig.4.pdf",width=6,height=3)
plot_grid(pl7a,pl7b,ncol=2,align="hv",rel_widths=c(2,1))
dev.off();

############################# FIGURE S5

pl7=ggplot(data=subset(indf5,time==0),aes(x=rich,y=feas,color=trait,alpha=as.factor(rho)))+
stat_summary(fun.y=mean, geom="line", size = 0.5,position=position_dodge(width=0.5))+
stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), geom = "pointrange")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Number of species per guild")+facet_wrap(~trait)+scale_color_manual(values=colo,guide=FALSE)+scale_fill_manual(values=colo)+
#scale_color_brewer(palette = "Reds")+
labs(alpha=expression(rho))+ggtitle("a",subtitle="Initial networks")+scale_x_continuous(breaks=c(10,20,30))

pl8=ggplot(data=subset(indf5,time==2000),aes(x=rich,y=feas,color=trait,alpha=as.factor(rho)))+
stat_summary(fun.y=mean, geom="line", size = 0.5,position=position_dodge(width=0.5))+
stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - sd(x), fun.ymax = function(x) mean(x) + sd(x), geom = "pointrange")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Number of species per guild")+facet_wrap(~trait)+scale_color_manual(values=colo,guide=FALSE)+scale_fill_manual(values=colo)+
#scale_color_brewer(palette = "Reds")+
labs(alpha=expression(rho))+ggtitle("b",subtitle="Coevolved networks")+scale_x_continuous(breaks=c(10,20,30))


plot_grid(pl7,pl8,ncol=1,align="hv")

png("Fig_S5.png",width=1200,height=900)
plot_grid(pl7,pl8,ncol=1,align="hv")
dev.off();


ggplot(data=subset(indf5,time==2000),aes(x=rho,y=feas,color=trait,linetype=as.factor(rich)))+
stat_smooth()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Number of species per guild")+facet_wrap(~trait)+scale_color_manual(values=colo,guide=FALSE)+scale_fill_manual(values=colo)+
#scale_color_brewer(palette = "Reds")+
labs(alpha=expression(rho))+ggtitle("a",subtitle="Initial networks")+scale_x_continuous(breaks=c(10,20,30))









indf %>% group_by(rich,trait,time) %>% summarise(mean(feas))




indf5=subset(indf4,time %in% c(0,2000))

pl7=ggplot(data=subset(indf5,time==0),aes(x=rho,y=feas,color=as.factor(rich)))+
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
