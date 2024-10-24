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
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")
phenp=fread("flower_pheno_empirical.csv")
phenp[phenp$sde<1]=1
phena=fread("poll_pheno_empirical.csv")
phena[phena$sde<1]=1
load("matrices_empirical_networks.RData")

sites=names(networks)
ncalc=10
weigthed.mean=function(x,y){sum(x*y,na.rm=T)/sum(y,na.rm=T)}
empf=NULL
for (site in sites){
m=networks[[site]]

mbi=round(m,digits=5)

emp=networklevel(mbi,index=c("weighted connectance","weighted nestedness","interaction evenness","H2"))
emp=as.data.frame(t(emp))
emp$NODF=networklevel(round(mbi,digits=5),index=c("weighted NODF"))
names(emp)=gsub("weighted ","",names(emp))
emp$modularity=computeModules(mbi, method="Beckett")@"likelihood"
emp$na=ncol(m)
emp$np=nrow(m)
emp$site=site

#build competition matrix for poll:
liste=colnames(m)
comp_a = matrix(NA, ncol(m), ncol(m))
colnames(comp_a)=liste
rownames(comp_a)=liste
phen_a = matrix(NA, ncol(m), ncol(m))
colnames(phen_a)=liste
rownames(phen_a)=liste
for(i in liste){
	mu1=subset(phena,Pollinator_gen_sp==i)$mu
	sd1=subset(phena,Pollinator_gen_sp==i)$sde
	for(j in liste){
		mu2=subset(phena,Pollinator_gen_sp==j)$mu
		sd2=subset(phena,Pollinator_gen_sp==j)$sde
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_a[i,j]=phen
	}  
}
diag(comp_a)=1


#build competition matrix for plants:
liste=rownames(m)
comp_p = matrix(NA, nrow(m), nrow(m))
colnames(comp_p)=liste
rownames(comp_p)=liste
phen_p = matrix(NA, nrow(m), nrow(m))
colnames(phen_p)=liste
rownames(phen_p)=liste
for(i in liste){
	mu1=subset(phenp,Plant_gen_sp==i)$mu
	sd1=subset(phenp,Plant_gen_sp==i)$sde
	for(j in liste){
		mu1=subset(phenp,Plant_gen_sp==j)$mu
		sd1=subset(phenp,Plant_gen_sp==j)$sde
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_p[i,j]=phen
	}  
}




###3 PHEN OVERLAP AMONG PARTNERS
listea=colnames(m)
listep=rownames(m)
phen_mut = matrix(NA, nrow(m), ncol(m))
colnames(phen_mut)=listea
rownames(phen_mut)=listep
for(i in listea){
	mu1=subset(phena,Pollinator_gen_sp==i)$mu
	sd1=subset(phena,Pollinator_gen_sp==i)$sde
	for(j in listep){
		mu1=subset(phenp,Plant_gen_sp==j)$mu
		sd1=subset(phenp,Plant_gen_sp==j)$sde
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_mut[j,i]=phen
	}  
}



###### MOTIF COUNT
motifn=0
motifntot=0
liste=colnames(m)
for(i in liste){
id_plant=which(m[,i]>0)
if(length(id_plant)>1){
combis=as.data.frame(t(combn(id_plant,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_p[x[1],x[2]]})
combis$mut_overlap=sqrt(apply(combis,1,function(x){phen_mut[x[1],i]})*apply(combis,1,function(x){phen_mut[x[2],i]}))
combis$mut=apply(combis,1,function(x){sqrt(m[x[1],i]*m[x[2],i])})
motifn=motifn+sum(combis$phen_overlap<combis$mut_overlap)
motifntot=motifntot+nrow(combis)
}
}

liste=rownames(m)
for(i in liste){
id_poll=which(m[i,]>0)
if(length(id_poll)>1){
combis=as.data.frame(t(combn(id_poll,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_a[x[1],x[2]]})
combis$mut_overlap=sqrt(apply(combis,1,function(x){phen_mut[i,x[1]]})*apply(combis,1,function(x){phen_mut[i,x[2]]}))
combis$mut=apply(combis,1,function(x){sqrt(m[i,x[1]]*m[i,x[2]])})
motifn=motifn+sum(combis$phen_overlap<combis$mut_overlap)
motifntot=motifntot+nrow(combis)
}
}

emp$motifn=motifn
emp$motifntot=motifntot

for(rho in c(0.01,0.05,0.2)){
for(competition in c(0.1,1,5,10)){
comp_a[,]=rho
comp_p[,]=rho
diag(comp_p)=1
diag(comp_a)=1
comp_phen_a=comp_a*phen_a
diag(comp_phen_a)=1
comp_phen_p=comp_p*phen_p
diag(comp_phen_p)=1

A1=rbind(cbind(competition*comp_phen_p,-1*m),cbind(-1*t(m),competition*comp_phen_a))
A2=rbind(cbind(competition*comp_p,-1*m),cbind(-1*t(m),competition*comp_a))

vec_stab1=c()
vec_stab2=c()
for(it in 1:ncalc){
vec_stab1=c(vec_stab1,Omega(A1))
vec_stab2=c(vec_stab2,Omega(A2))
}
emp$feas_with_pheno=mean(vec_stab1)
emp$feas_with_pheno_se=sd(vec_stab1)/sqrt(ncalc)
emp$feas_without_pheno=mean(vec_stab2)
emp$feas_without_pheno_se=sd(vec_stab2)/sqrt(ncalc)
emp$competition=competition
emp$rho=rho
empf=rbind(empf,emp)
}
}

}

fwrite(empf,"empirical_networks.csv")

############################### EXTRACT INDICES FROM SIMULATIONS DONE ON EMPIRICAL NETWORKS
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("data.table", "dplyr","bipartite") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/results_empir")
source("C:/Users/Duchenne/Documents/evolution_pheno_morpho/scripts/functions and secondary scripts/toolbox.R")
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

liste = fread("C:/Users/Duchenne/Documents/evolution_pheno_morpho/initial_empir/liste.csv")

datf=NULL
for(ess in 1:10){
for(jj in 1:nrow(liste)){
dat=fread(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/results_empir/ueq_",liste$site[jj],"_",ess,".csv"))
nbsp_a=liste$na[jj]
nbsp_p=liste$np[jj]

for(competition in c(10)){
names(dat)[1:(nbsp_a+nbsp_p)]=gsub("x","mu_",names(dat)[1:(nbsp_a+nbsp_p)])
names(dat)[(nbsp_a+nbsp_p+1):((nbsp_a+nbsp_p)*2)]=gsub("x","sd_",names(dat)[(nbsp_a+nbsp_p+1):((nbsp_a+nbsp_p)*2)])
names(dat)[(nbsp_a+nbsp_p)*2+1]="time"

dat2=melt(dat,id.vars=c("trait","time"))
dat2$type="mu"
dat2$type[grep("sd",dat2$variable)]="sd"
dat2$species=gsub("mu_","",dat2$variable)
dat2$species=as.numeric(gsub("sd_","",dat2$species))
dat2$species[dat2$species>(nbsp_a+nbsp_p)]=dat2$species[dat2$species>(nbsp_a+nbsp_p)]-(nbsp_a+nbsp_p)
dat2$comp=competition
dat2$essai=ess
dat2$site=liste$site[jj]
datf=rbind(datf,dat2)
}
}
}
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")
fwrite(datf,"species_level_simues_empir.csv")

########################################### TO ran on HPC
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("data.table", "dplyr","bipartite") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="/home/duchenne/pheno") #run on a cluster
source("/home/duchenne/pheno/toolbox.R")  #run on a cluster
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
jj <- as.numeric(args_contents[[1]])
print(jj)


datf = fread("species_level_simues_empir.csv")
liste = fread("C:/Users/Duchenne/Documents/evolution_pheno_morpho/initial_empir/liste.csv")

indf=NULL
motifs=NULL
for(ess in 1:10){
nbsp_a=liste$na[jj]
nbsp_p=liste$np[jj]
for(competition in c(10)){
for(tr in c("morpho","pheno")){
for(ti in c(0,2000)){

bidon=subset(datf,time==ti & trait==tr & essai==ess & site==liste$site[jj])

#build interaction matrix:
m = matrix(NA, nbsp_p, nbsp_a)
for(i in 1:nbsp_a){
	mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
	for(j in 1:nbsp_p){
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+nbsp_a))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+nbsp_a))$value)
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		m[j,i]=sum(f(seq(0,365,0.1)))*0.1
	}  
}

struf=NULL
for(rr in c(0.05)){
#mbi=apply(m, 1L:2L, function(p) rbinom(1, 5, p))
#mbi[mbi>0]=1
mbi=round(m,digits=5)
stru=networklevel(mbi,index=c("weighted connectance","weighted nestedness","interaction evenness","H2"),H2_integer=FALSE)
stru=as.data.frame(t(stru))
stru$NODF=networklevel(round(mbi,digits=5),index=c("weighted NODF"))
stru$modularity=computeModules(mbi, method="Beckett")@"likelihood"
struf=rbind(struf,stru)
}

ind=as.data.frame(t(apply(struf,2,mean,na.rm=T)))
ind$trait=tr
ind$time=ti
ind$essai=ess
ind$nbsp_a=nbsp_a
ind$nbsp_p=nbsp_p
ind$site=liste$site[jj]

#build competition matrix for poll:
comp_a = matrix(NA, nbsp_a, nbsp_a)
phen_a = matrix(NA, nbsp_a, nbsp_a)
for(i in 1:nbsp_a){
	mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
	for(j in 1:nbsp_a){
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
comp_p = matrix(NA, nbsp_p, nbsp_p)
phen_p = matrix(NA, nbsp_p, nbsp_p)
for(i in 1:nbsp_p){
	mu1=invlogit1(subset(bidon,type=="mu" & species==(i+nbsp_a))$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==(i+nbsp_a))$value)
	for(j in 1:(nbsp_p)){
		similarity=sum(t(m)[,i] * t(m)[,j]) / sum(t(m)[,i])
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+nbsp_a))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+nbsp_a))$value)
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_p[i,j]=phen
		comp_p[i,j]=ifelse(tr=="pheno",phen*similarity,similarity)
	}  
}
diag(comp_p)=1

motifn=0
motifntot=0
for(i in 1:nbsp_a){
id_plant=which(m[,i]>=0.01)
if(length(id_plant)>1){
combis=as.data.frame(t(combn(id_plant,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_p[x[1],x[2]]})
combis$mut=apply(combis,1,function(x){m[x[1],i]*m[x[2],i]})
motifn=motifn+sum(combis$phen_overlap<combis$mut)
motifntot=motifntot+nrow(combis)
}
}

for(i in 1:nbsp_p){
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
for(competition in c(0.01,1,5,10)){
comp_a[,]=rho
comp_p[,]=rho
diag(comp_p)=1
diag(comp_a)=1
comp_phen_a=comp_a*phen_a
diag(comp_phen_a)=1
comp_phen_p=comp_p*phen_p
diag(comp_phen_p)=1
if(tr=="morpho"){
A=rbind(cbind(competition*comp_p,-1*m),cbind(-1*t(m),competition*comp_a))
}else{
A=rbind(cbind(competition*comp_phen_p,-1*m),cbind(-1*t(m),competition*comp_phen_a))
}

vec_stab=c()
for(it in 1:5){
vec_stab=c(vec_stab,Omega(A))
}
ind$feas_with_pheno=mean(vec_stab)
ind$feas_with_pheno_se=sd(vec_stab)/sqrt(5)
ind$competition=competition
ind$rho=rho
indf=rbind(indf,ind)
}

}
}
}
}
}

fwrite(indf,"C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info_empir.csv")
#
