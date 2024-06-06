###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/scripts")

source("toolbox.R")

invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

datf=NULL
i=1
rich=10
dat=fread(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/results/ueq_",i,"_",rich,".csv"))
dive=rich*2
names(dat)[1:(dive)]=gsub("x","mu_",names(dat)[1:(dive)])
names(dat)[(dive+1):(dive*2)]=gsub("x","sd_",names(dat)[(dive+1):(dive*2)])
names(dat)[(dive*2+1)]="time"

dat2=melt(dat,id.vars=c("trait","time"))
dat2$type="mu"
dat2$type[grep("sd",dat2$variable)]="sd"
dat2$species=gsub("mu_","",dat2$variable)
dat2$species=as.numeric(gsub("sd_","",dat2$species))
dat2$species[dat2$species>dive]=dat2$species[dat2$species>dive]-dive
dat2$essai=i
dat2$rich=rich

dat2$guild="poll"
dat2$guild[dat2$species>10]="plant"
ggplot(data=subset(dat2,type=="mu" & trait=="morpho"),aes(x=time,y=invlogit1(value),color=guild,group=species))+geom_line(size=1.2,alpha=0.8)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0.5),
strip.background=element_rect(fill=NA,color=NA))+ylab("trait value")+scale_colour_manual(values=c("gold3","black"))