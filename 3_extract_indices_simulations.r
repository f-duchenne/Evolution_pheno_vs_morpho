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
for(competition in c(5,10)){
for (rich in c(10,20,30)){
for(i in 1:100){
dat=fread(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/results/ueq_",i,"_",rich,"_",competition,".csv"))
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
dat2$comp=competition
dat2$essai=i
dat2$rich=rich
datf=rbind(datf,dat2)
}
}
}

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")
fwrite(datf,"species_level_data.csv")

###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("data.table", "dplyr","bipartite") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="/home/duchenne/pheno")
source("/home/duchenne/pheno/toolbox.R")
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
ess <- as.numeric(args_contents[[1]])
print(ess)


datf=fread("/home/duchenne/pheno/species_level_data.csv")


indf=NULL
motifs=NULL
nrandom=20
for(competition in c(5,10)){
for (rich in c(10,20,30)){
for(tr in c("morpho","pheno")){
for(ti in plyr::round_any(seq(sqrt(0),sqrt(2000),length.out=9)^2,10)){
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

motifn=0
motifntot=0
for(i in 1:rich){
id_plant=which(m[,i]>=0.01)
if(length(id_plant)>1){
combis=as.data.frame(t(combn(id_plant,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_p[x[1],x[2]]})
combis$mut=apply(combis,1,function(x){m[x[1],i]*m[x[2],i]})
motifn=motifn+sum(combis$phen_overlap<combis$mut)
motifntot=motifntot+nrow(combis)
}
}

for(i in 1:rich){
id_poll=which(m[i,]>=0.01)
if(length(id_poll)>1){
combis=as.data.frame(t(combn(id_poll,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_a[x[1],x[2]]})
combis$mut=apply(combis,1,function(x){m[i,x[1]]*m[i,x[2]]})
motifn=motifn+sum(combis$phen_overlap<combis$mut)
motifntot=motifntot+nrow(combis)
}
}

ind$motifn=motifn
ind$motifntot=motifntot
ind$comp=mean(c(comp_p,comp_a))


for(rho in c(0.01,0.05,0.2)){
comp_a[,]=rho
comp_p[,]=rho
diag(comp_p)=1
diag(comp_a)=1
comp_phen_a=comp_a*phen_a
diag(comp_phen_a)=1
comp_phen_p=comp_p*phen_p
diag(comp_phen_p)=1

A=rbind(cbind(competition*comp_p,-1*m),cbind(-1*t(m),competition*comp_a))
ind$feas=Omega(A)
ind$competition=competition
ind$rho=rho
indf=rbind(indf,ind)
}

}
}
}
}

fwrite(indf,paste0("/home/duchenne/pheno/aggreg_simues/networks_info_",ess,".csv"))

####################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

indf=NULL
for(ess in 1:100){
ind=fread(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/aggreg_simues/networks_info_",ess,".csv"))
indf=rbind(indf,ind)
}
fwrite(indf,"C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info.csv")
