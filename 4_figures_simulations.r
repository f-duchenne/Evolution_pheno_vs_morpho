###################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales",
"car","viridis","ggeffects","emmeans") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

colo=c("chartreuse3","dodgerblue4","gold3")
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)

##################################### COMMUNITY LEVEL ########################
indf=fread("data/simulated/outputs_simulations/networks_info.csv")
names(indf)=gsub("weighted ","w",names(indf))
names(indf)[3]="inter. evenness"
names(indf)[5]="wNODF"
names(indf)[7]="mut.strength"
indf$prop_motif=(indf$motifn/indf$motifntot)

acp=PCA(indf[,1:7],graph=FALSE)
indf2=cbind(indf,acp$ind$coord)
indf2$comp=apply(indf2[,c("comp_a","comp_p")],1,mean)

bidon=subset(indf2,time==0)
s1=ggplot(data=bidon,aes(x=as.factor(competition),y=mut.strength/(competition*comp),color=trait))+geom_boxplot()+
geom_hline(yintercept=1,linetype="dashed")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
ylab("mutualism/competition in initial networks")+
xlab("Importance of competition")+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+labs(colour="Initial conditions with:")

png("Figures/Fig_S1.png",width=1200,height=900,res=150)
s1
dev.off();

indf3=subset(indf2,competition %in% c(2,6))
indf3=indf3[order(indf3$time),]

indfspeed=indf3 %>% group_by(trait,competition,essai,rich) %>% mutate(speed=sqrt((Dim.1-lag(Dim.1))^2+(Dim.2-lag(Dim.2))^2+(Dim.3-lag(Dim.3))^2)/(time-lag(time)))

ggplot(subset(indfspeed,time>0),aes(x=time,y=speed,color=trait))+geom_line(alpha=0.2,aes(group=paste(essai,competition,trait,rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Changes in network structure")+labs(colour="Simulations with:",fill="Simulations with:")+
guides(linetype=guide_legend(override.aes=list(fill=NA)))+
facet_wrap(~competition,ncol=1)

############### FIGURE 2
p1=ggplot()+
stat_smooth(data=indf3,aes(y=mut.strength,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Average mutualism strength")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
xlab("Time")+ggtitle("a")+labs(colour="Simulations with:",fill="Simulations with:",linetype=expression(paste(n[sp]," / guild")))+scale_x_continuous(breaks=c(0,1000,2000))+
facet_wrap(~competition,ncol=1,labeller = label_bquote(italic(c) == .(competition)))

p2=ggplot()+
stat_smooth(data=indf3,aes(y=comp,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="bottom",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Average competition stength")+guides(linetype=guide_legend(title.position="top",title.hjust = 0.5,override.aes=list(fill=NA,color="black")),colour=guide_legend(title.position="top",title.hjust = 0.5))+
xlab("Time")+ggtitle("b")+labs(colour="Simulations with:",fill="Simulations with:",linetype=expression(paste(n[sp]," / guild")))+scale_x_continuous(breaks=c(0,1000,2000))+
facet_wrap(~competition,ncol=1,labeller = label_bquote(italic(c) == .(competition)))

p3=ggplot()+
stat_smooth(data=indf3,aes(y=mut.strength/(competition*comp),x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("mutualism/competition")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
xlab("Time")+ggtitle("c")+labs(colour="Simulations with:",fill="Simulations with:")+scale_x_continuous(breaks=c(0,1000,2000))+
facet_wrap(~competition,ncol=1,labeller = label_bquote(italic(c) == .(competition)))

indf3$motifn[indf3$trait=="morpho"]=0
p4=ggplot(data=indf3,aes(y=motifn/motifntot,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - 1.96*sd(x)/length(x), fun.ymax = function(x) mean(x) + 1.96*sd(x)/length(x), geom = "ribbon",alpha=0.1,colour=NA)+
stat_summary(fun.y=mean, geom="line", size = 0.5,linewidth=1.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab('Proportion of "V+" motifs')+
xlab("Time")+ggtitle("d")+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+scale_x_continuous(breaks=c(0,1000,2000))+facet_wrap(~competition,ncol=1,labeller = label_bquote(italic(c) == .(competition)))

leg <- as_ggplot(get_legend(p2))
p2=p2+theme(legend.position="none")
top=plot_grid(p1,p2,p3,p4,align="hv",ncol=4)
grid.arrange(top,leg,ncol=1,heights=c(1,0.3))


pdf("Figures/Fig.2.pdf",width=7,height=5)
grid.arrange(top,leg,ncol=1,heights=c(1,0.3))
dev.off();


############### FIGURE 3
pl1=ggplot(data=subset(indf2,rich==30 & time %in% c(0,20,100,220,2000) & competition %in% c(4,6)),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,group=paste(essai,trait)))+
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
geom_point(data=subset(indf2,time==2000  & competition %in% c(4)),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,fill=trait,shape=as.factor(rich)))+
#geom_point(data=subset(indf3,time==0 & trait!="both"),aes(x=Dim.1,y=Dim.2,alpha=time,shape=as.factor(rich)),color="grey",fill="grey",alpha=0.3)+
#geom_point(data=subset(indf3,time==0 & trait=="both"),aes(x=Dim.1,y=Dim.2,alpha=time,shape=as.factor(rich)),color="black",fill="grey",alpha=0.3)+
scale_shape_manual(values = c(21,22,23))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ggtitle("c")+labs(colour="Simulations with:",fill="Simulations with:")


pl4=ggplot(data=subset(indf2,time==2000 & competition %in% c(4)),aes(x=as.factor(rich),y=wNODF,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(N[A] == N[P]))+ylab("wNODF")+xlab(expression(paste(n[sp], " / guild")))+
ggtitle("d")

pl5=ggplot(data=subset(indf2,time==2000 & competition %in% c(4)),aes(x=as.factor(rich),y=modularity,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
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

pdf("Figures/Fig.3.pdf",width=8,height=7)
grid.arrange(top,midle,bottom,widths=c(2,1),heights=c(1,1,0.7),layout_matrix = lay)
dev.off();

indf3b=subset(indf2,competition==6)
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
scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+labs(colour="Simulations with:",fill="Simulations with:")

png("Figures/Fig_S5.png",width=1200,height=1000,res=150)
pl3b
dev.off();


######################################### FIGURE 4

#EVOLUTION OF FEASIBILITY

pl6=ggplot(data=subset(indf2,competition==4),aes(x=time,y=feas_structure,color=trait,group=paste(time,trait)))+geom_violin(position=position_dodge(width=0.5),width=8,scale="width",show.legend=FALSE)+
geom_boxplot(width=3,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.background = element_blank(),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_sqrt(breaks=c(0,1000,2000))+labs(colour="Simulations with:",fill="Simulations with:")+
ggtitle("a")

pl6sup=ggplot(data=subset(indf2),aes(x=time,y=feas_structure,color=trait,group=paste(time,trait)))+geom_violin(position=position_dodge(width=0.5),width=8,scale="width",show.legend=FALSE)+
geom_boxplot(width=3,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.background = element_blank(),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Time")+facet_grid(cols=vars(rich),rows=vars(competition), labeller = label_bquote(cols=n[sp] == .(rich),rows=italic(c) == .(competition)))+scale_x_sqrt(breaks=c(0,1000,2000))+labs(colour="Simulations with:",fill="Simulations with:")

png("Figures/Fig_S6.png",width=1200,height=1200,res=150)
pl6sup
dev.off();


#Diversity-stability trade-off
indf5=subset(indf2,time %in% c(0,2000) & competition %in% c(4,6))
indf5$trait=factor(indf5$trait,labels=c("both traits","only morpho","only pheno"))
indf5=indf5 %>% group_by(trait,essai,competition,time) %>% mutate(feas_ini=feas_structure[rich==10])  
indf5$time2=factor(indf5$time,labels=c("initial (t = 0)","coevolved (t = 2000)"))
indf5$time2=factor(indf5$time2,levels=c("coevolved (t = 2000)","initial (t = 0)"))

model=glm(feas_structure~rich*time2*trait*as.factor(competition),data=subset(indf5),family=quasibinomial)

b=ggpredict(model,c("rich[10:30]","time2","trait","competition"))

pl6b=ggplot()+
geom_line(data=subset(b,panel==4),aes(x=x,y=predicted,color=facet,linetype=group))+
geom_point(data=subset(indf5,competition==4),aes(x=rich,y=feas_structure,color=trait,shape=as.factor(time2)),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(linetype="",shape="")+ylab(expression(paste("Structural stability  ",(omega))))+xlab("Number of species per guild")+scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
#scale_color_brewer(palette = "Reds")+
labs(alpha=expression(rho),color="Simulations with:")+scale_x_continuous(breaks=c(10,20,30))+
ggtitle("b")

b2=as.data.frame(emmeans(model, pairwise ~ rich | time2 + trait + competition)$emmeans)
b2$time=0
b2$time[b2$time2=="coevolved (t = 2000)"]=2000
b2$time=ordered(b2$time,levels=c(0,2000))

pl6c=ggplot(data=subset(b2),aes(x=time,y=emmean,color=trait))+
geom_hline(yintercept=0,linetype="dashed")+
geom_pointrange(aes(ymin=asymp.LCL,ymax=asymp.UCL,alpha=as.factor(competition)))+
geom_line(aes(group=paste0(trait,competition),alpha=as.factor(competition)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="right",
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(alpha=expression(italic(c)))+ylab("Slope of the\ndiversity-stability trade off")+xlab("")+scale_color_manual(values=colo,guide=FALSE)+scale_fill_manual(values=colo)+
ggtitle("c")+scale_x_discrete(breaks=c(0,2000),labels=c("initial","coevolved"))+scale_alpha_discrete(range=c(0.3,1))

bottom=plot_grid(pl6b,pl6c,ncol=2,align="hv")
grid.arrange(pl6,bottom,ncol=1)

pdf("Figures/Fig.4.pdf",width=9,height=6)
grid.arrange(pl6,bottom,ncol=1)
dev.off();
