###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/scripts")

source("toolbox.R")

invlogit=function(x){return(1+50*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+285*exp(x)/(1+exp(x)))}

dive=10*2
datf=NULL
for(i in 1:10){
dat=fread(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/ueq_",i,".csv"))

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
datf=rbind(datf,dat2)
}


indf=NULL
for(j in c(1:2)){
for(tr in c("morpho","pheno")){
for(ti in c(1,100,300,800,1500)){
bidon=subset(datf,essai==j & time==ti & trait==tr)

#build interaction matrix:
m = matrix(NA, dive/2, dive/2)
for(i in 1:(dive/2)){
	mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
	for(j in 1:(dive/2)){
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+dive/2))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+dive/2))$value)
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		m[j,i]=integrate(f,0,365)[[1]]
	}  
}
   
ind=networklevel(m,index=c("connectance","weighted connectance","weighted nestedness","interaction evenness","H2"))
ind=as.data.frame(t(ind))
ind$modularity=computeModules(m, method="Beckett")@"likelihood"
ind$trait=tr
ind$time=ti
ind$essai=i
indf=rbind(indf,ind)
}
}
}

acp=PCA(indf[,1:6])
acp$coordinates

indf=cbind(indf,acp$ind$coord)

ggplot(data=indf,aes(x=Dim.1,y=Dim.2,alpha=time,color=trait))+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point()+
geom_line()+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))


