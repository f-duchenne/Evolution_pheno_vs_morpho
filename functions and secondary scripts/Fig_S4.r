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

indf3=subset(indf2,competition %in% c(4))

time_end=5000


############### FIGURE S4
p1=ggplot()+
stat_smooth(data=indf3,aes(y=mut.strength,x=time,color=trait,fill=trait,linetype=as.factor(rich),alpha=as.factor(competition)),method="loess")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Average mutualism strength")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
xlab("Time")+ggtitle("a")+labs(colour="Simulations with:",fill="Simulations with:",linetype=expression(paste(n[sp]," / guild")))+scale_x_continuous(breaks=c(0,2500,time_end))+
facet_wrap(~competition,ncol=1,labeller = label_bquote(italic(c) == .(competition)))

p2=ggplot()+
stat_smooth(data=indf3,aes(y=comp,x=time,color=trait,fill=trait,linetype=as.factor(rich)),method="loess")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="bottom",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("Average competition stength")+guides(linetype=guide_legend(title.position="top",title.hjust = 0.5,override.aes=list(fill=NA,color="black")),colour=guide_legend(title.position="top",title.hjust = 0.5))+
xlab("Time")+ggtitle("b")+labs(colour="Simulations with:",fill="Simulations with:",linetype=expression(paste(n[sp]," / guild")))+scale_x_continuous(breaks=c(0,2500,time_end))+
facet_wrap(~competition,ncol=1,labeller = label_bquote(italic(c) == .(competition)))

p3=ggplot()+
stat_smooth(data=indf3,aes(y=mut.strength/(competition*comp),x=time,color=trait,fill=trait,linetype=as.factor(rich)),method="loess")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("mutualism/competition")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
xlab("Time")+ggtitle("c")+labs(colour="Simulations with:",fill="Simulations with:")+scale_x_continuous(breaks=c(0,2500,time_end))+
facet_wrap(~competition,ncol=1,labeller = label_bquote(italic(c) == .(competition)))

indf3$motifn[indf3$trait=="morpho"]=0
p4=ggplot(data=indf3,aes(y=motifn/motifntot,x=time,color=trait,fill=trait,linetype=as.factor(rich)),method="loess")+
stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - 1.96*sd(x)/length(x), fun.ymax = function(x) mean(x) + 1.96*sd(x)/length(x), geom = "ribbon",alpha=0.1,colour=NA)+
stat_summary(fun.y=mean, geom="line", size = 0.5,linewidth=1.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab('Proportion of "V+" motifs')+
xlab("Time")+ggtitle("d")+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+scale_x_continuous(breaks=c(0,2500,time_end))+facet_wrap(~competition,ncol=1,labeller = label_bquote(italic(c) == .(competition)))

leg <- as_ggplot(get_legend(p2))
p2=p2+theme(legend.position="none")
top=plot_grid(p1,p2,p3,p4,align="hv",ncol=4)



png("Figures/Fig_S4.png",width=1300,height=700,res=140)
grid.arrange(top,leg,ncol=1,heights=c(1,0.3))
dev.off();