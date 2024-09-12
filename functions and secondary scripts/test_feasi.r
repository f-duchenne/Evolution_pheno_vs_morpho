###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/scripts/functions and secondary scripts")
source("toolbox.R")
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")
datf=fread("species_level_data.csv")

ess=1

indf=NULL
motifs=NULL
nrandom=20
#for(competition in c(5,10)){
competition=10
#for (rich in c(10,20,30)){
rich=10
#for(tr in c("morpho","pheno")){
tr="morpho"
#for(ti in plyr::round_any(seq(sqrt(0),sqrt(2000),length.out=9)^2,10)){
ti=2000
bidon=subset(datf,essai==ess & time==ti & trait==tr & rich==rich & comp==competition)

#build interaction matrix:
m = matrix(NA, rich, rich)
for(i in 1:rich){
	mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
	for(j in 1:rich){
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+rich))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+rich))$value)
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		m[j,i]=sum(f(seq(0,365,0.1)))*0.1
	}  
}

struf=NULL
for(rr in 1:nrandom){
mbi=apply(m, 1L:2L, function(p) rbinom(1, 5, p))
# mbi[mbi>0]=1
# mbi=m
stru=networklevel(mbi,index=c("connectance","nestedness","interaction evenness","H2"),H2_integer=TRUE)
stru=as.data.frame(t(stru))
stru$NODF=networklevel(round(mbi,digits=5),index=c("NODF"))
stru$modularity=computeModules(mbi, method="Beckett")@"likelihood"
struf=rbind(struf,stru)
}

ind=as.data.frame(t(apply(struf,2,mean,na.rm=T)))
ind$trait=tr
ind$time=ti
ind$essai=ess
ind$rich=rich

#build competition matrix for poll:
comp_a = matrix(NA, rich, rich)
phen_a = matrix(NA, rich, rich)
for(i in 1:rich){
	mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
	for(j in 1:(rich)){
		similarity=sum(m[,i] * m[,j]) / sum(m[,i])
		mu2=invlogit1(subset(bidon,type=="mu" & species==j)$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==j)$value)
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_a[i,j]=phen
		comp_a[i,j]=ifelse(tr=="pheno",phen*similarity,similarity)
	}  
}
diag(comp_a)=1

#build competition matrix for plants:
comp_p = matrix(NA, rich, rich)
phen_p = matrix(NA, rich, rich)
for(i in 1:rich){
	mu1=invlogit1(subset(bidon,type=="mu" & species==(i+rich))$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==(i+rich))$value)
	for(j in 1:(rich)){
		similarity=sum(t(m)[,i] * t(m)[,j]) / sum(t(m)[,i])
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+rich))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+rich))$value)
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_p[i,j]=phen
		comp_p[i,j]=ifelse(tr=="pheno",phen*similarity,similarity)
	}  
}
diag(comp_p)=1

indf=NULL
for(competition in c(0.1,0.5,1,5)){
for(rho in c(0.01,0.1,0.5)){
comp_a[,]=rho
comp_p[,]=rho
diag(comp_p)=1
diag(comp_a)=1
if(tr=="morpho"){
comp_phen_a=comp_a
}else{
comp_phen_a=comp_a*phen_a
}
diag(comp_phen_a)=1
if(tr=="morpho"){
comp_phen_p=comp_p
}else{
comp_phen_p=comp_p*phen_p
}
diag(comp_phen_p)=1

A=rbind(cbind(competition*comp_phen_p,-1*m),cbind(-1*t(m),competition*comp_phen_a))
vec_stab=c()
for(it in 1:5){vec_stab=c(vec_stab,Omega(A))}
ind=data.frame(feas=mean(vec_stab),feas_se=sd(vec_stab)/sqrt(5),competition=competition,rho=rho)

for(nul in 1:20){
mnul=m
mnul[,]=sample(c(m),nrow(m)*ncol(m))
A=rbind(cbind(competition*comp_phen_p,-1*mnul),cbind(-1*t(mnul),competition*comp_phen_a))
vec_stab=c()
for(it in 1:5){vec_stab=c(vec_stab,Omega(A))}
ind$feas_nul=mean(vec_stab)
ind$feas_nul_se=sd(vec_stab)/sqrt(5)
ind$nul=nul
indf=rbind(indf,ind)
}
}
}

indf$diag_dominance=1/indf$rho

ggplot(data=subset(indf))+geom_boxplot(aes(x=as.factor(competition),y=feas_nul),alpha=0.5,fill="grey")+geom_point(aes(x=as.factor(competition),y=feas),col="red")+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="bottom",
strip.background=element_rect(fill=NA,color=NA))+xlab("Overall comeptition strength")+
ylab("feasibility (red dot is not randomized matrix)")+facet_wrap(~rho)

ggplot(data=indf+geom_boxplot(aes(x=as.factor(competition),y=feas_nul),alpha=0.5,fill="grey")+geom_point(aes(x=as.factor(competition),y=feas),col="red")+facet_wrap(~rho,scales="free_y",ncol=3)+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="bottom",
strip.background=element_rect(fill=NA,color=NA))+xlab("Overall comeptition strength")+
ylab("feasibility (red dot is not randomized matrix)")
