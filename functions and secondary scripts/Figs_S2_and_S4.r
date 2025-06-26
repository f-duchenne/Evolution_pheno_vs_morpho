###################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","ggplot2","gridExtra","cowplot","ggpubr","scales") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

colo=c("chartreuse3","dodgerblue4","gold3")
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)

##################################### SPECIES LEVEL ########################
datf=fread("data/simulated/outputs_simulations/species_level.csv")
datf$type2="pheno"
datf$type2[datf$type %in% c("mu2","sd2")]="morpho"
datf$type2[datf$trait %in% c("morpho")]="morpho"
datf$type3="mu"
datf$type3[datf$type %in% c("sd","sd2")]="sd"

datf2=subset(datf,comp==4 & time %in% plyr::round_any(seq(sqrt(0),sqrt(2000),length.out=30)^2,10))
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

png("Figures/Fig_S2.png",width=1300,height=1100,res=150)
plot_grid(s7_a,s7_b,ncol=1,align="hv")
dev.off();


############### FIGURE S5
endpo=subset(datf,time==2000)
endpo$labi=paste0(endpo$type2,"\n c = ",endpo$comp)
sp1=ggplot()+
geom_density(data=subset(endpo,type3=="mu"),aes(x=invlogit1(value),color=trait,fill=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("density")+
xlab("Mean value of the trait")+ggtitle("a")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_grid(cols=vars(rich),rows=vars(labi), labeller = label_bquote(cols=n[sp] == .(rich)), scales = "free_y")+
labs(colour="Simulations with:",fill="Simulations with:")+
scale_y_continuous(n.breaks=3)

sp2=ggplot()+
geom_density(data=subset(endpo,type3=="sd"),aes(x=invlogit(value),color=trait,fill=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("density")+
xlab("Standard deviation of the trait")+ggtitle("b")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_grid(cols=vars(rich),rows=vars(labi), labeller = label_bquote(cols=n[sp] == .(rich)), scales = "free_y")+
labs(colour="Simulations with:",fill="Simulations with:")+
scale_y_continuous(n.breaks=3)

plot_grid(sp1,sp2,ncol=1,align="hv")

png("Figures/Fig_S4.png",width=1200,height=1600,res=150)
plot_grid(sp1,sp2,ncol=1,align="hv")
dev.off();



biche=dcast(endpo,trait+time+comp+essai+species+type2+rich~type3,value.var="value")
ggplot(data=biche,aes(x=invlogit1(mu),y=invlogit(sd),color=type2))+geom_point()+facet_wrap(~trait+comp+rich)