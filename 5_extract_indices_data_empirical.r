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

mbi=round(m,digits=3)

emp=networklevel(mbi,index=c("weighted connectance","weighted nestedness","interaction evenness","H2"))
emp=as.data.frame(t(emp))
emp$NODF=networklevel(round(mbi,digits=5),index=c("weighted NODF"))
names(emp)=gsub("weighted ","",names(emp))
emp$modularity=computeModules(mbi, method="Beckett")@"likelihood"
emp$mut.strength=mean(m)
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
for(competition in c(0.1,1,2,5)){
comp_a[,]=rho
comp_p[,]=rho
diag(comp_p)=1
diag(comp_a)=1
comp_phen_a=comp_a*phen_a
diag(comp_phen_a)=1
comp_phen_p=comp_p*phen_p
diag(comp_phen_p)=1

A1=rbind(cbind(-1*competition*comp_phen_p,m),cbind(t(m),-1*competition*comp_phen_a))
A2=rbind(cbind(-1*competition*comp_p,m),cbind(t(m),-1*competition*comp_a))

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

competition=5

datf=NULL
for(ess in 1:10){
for(tr in c("both","pheno","morpho")){
for(jj in 1:nrow(liste)){
dat=fread(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/results_empir_symmetric/ueq_",liste$site[jj],"_",tr,"_",ess,".csv"))
nbsp_a=liste$na[jj]
nbsp_p=liste$np[jj]

dive=nbsp_a+nbsp_p

if(tr=="both"){
names(dat)[1:(dive)]=paste0("mu_",1:dive)
names(dat)[(dive+1):(dive*2)]=paste0("sd_",1:dive)
names(dat)[(dive*2+1):(dive*3)]=paste0("mu2_",1:dive)
names(dat)[(dive*3+1):(dive*4)]=paste0("sd2_",1:dive)
names(dat)[(dive*4+1)]="time"
}else{
names(dat)[1:(dive)]=paste0("mu_",1:dive)
names(dat)[(dive+1):(dive*2)]=paste0("sd_",1:dive)
names(dat)[(dive*2+1)]="time"
}

dat2=melt(dat,id.vars=c("trait","time"))
dat2$type="mu"
dat2$type[grep("mu2",dat2$variable)]="mu2"
dat2$type[grep("sd",dat2$variable)]="sd"
dat2$type[grep("sd2",dat2$variable)]="sd2"
dat2$species=gsub("mu2_","",gsub("mu_","",dat2$variable))
dat2$species=as.numeric(gsub("sd2_","",gsub("sd_","",dat2$species)))
dat2$species[dat2$species>(nbsp_a+nbsp_p)]=dat2$species[dat2$species>(nbsp_a+nbsp_p)]-(nbsp_a+nbsp_p)
dat2$comp=competition
dat2$essai=ess
dat2$site=liste$site[jj]
datf=rbind(datf,dat2)
}
}
}
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")
fwrite(datf,"species_level_simues_empir_symmetric.csv")
datf_alleg=datf[datf$time %in% c(0,2000),]
fwrite(datf_alleg,"species_level_simues_empir_symmetric_alleg.csv")

########################################### TO ran on HPC
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("data.table", "dplyr","bipartite") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="/home/duchenne/pheno") #run on a cluster
source("toolbox.R")  #run on a cluster
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
jj <- as.numeric(args_contents[[1]])
print(jj)

datf = fread("/home/duchenne/pheno/species_level_simues_empir_symmetric_alleg.csv")
liste = fread("/home/duchenne/pheno/initial_empir/liste.csv")
ncalc=10
indf=NULL


motifs=NULL
for(ess in 1:10){
nbsp_a=liste$na[jj]
nbsp_p=liste$np[jj]
for(competition in c(5)){
for(tr in c("both","pheno")){
for(ti in c(0,2000)){

bidon=subset(datf,time==ti & trait==tr & essai==ess & site==liste$site[jj])

#build interaction matrix:
#build interaction matrix:
m = matrix(NA, nbsp_p, nbsp_a)
for(i in 1:nbsp_a){
	mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
	if(tr=="both"){
		mu1_m=invlogit1(subset(bidon,type=="mu2" & species==i)$value)
		sd1_m=invlogit(subset(bidon,type=="sd2" & species==i)$value)
	}
	for(j in 1:nbsp_p){
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+nbsp_a))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+nbsp_a))$value)
		if(tr=="both"){
			mu2_m=invlogit1(subset(bidon,type=="mu2" & species==j+nbsp_a)$value)
			sd2_m=invlogit(subset(bidon,type=="sd2" & species==j+nbsp_a)$value)
		}
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		if(tr=="both"){
			f2=function(x){pmin(dnorm(x,mu1_m,sd1_m),dnorm(x,mu2_m,sd2_m))}
			m[j,i]=(sum(f(seq(0,365,0.1)))*0.1)*(sum(f2(seq(0,365,0.1)))*0.1)
		}else{
			m[j,i]=sum(f(seq(0,365,0.1)))*0.1
		}
	}  
}


mbi=round(m,digits=3)
stru=networklevel(mbi,index=c("weighted connectance","weighted nestedness","interaction evenness","H2"),H2_integer=FALSE)
struf=as.data.frame(t(stru))
struf$NODF=networklevel(round(mbi,digits=5),index=c("weighted NODF"))
struf$modularity=computeModules(mbi, method="Beckett")@"likelihood"
struf$m_mean=mean(m)

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
	for(j in i:nbsp_a){
		similarity=mean(sqrt(m[,i] * m[,j]))
		mu2=invlogit1(subset(bidon,type=="mu" & species==j)$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==j)$value)
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_a[i,j]=phen
		phen_a[j,i]=phen
		comp_a[i,j]=ifelse(tr %in% c("pheno","both"),phen*similarity,similarity)
		comp_a[j,i]=ifelse(tr %in% c("pheno","both"),phen*similarity,similarity)
	}  
}
diag(comp_a)=1

#build competition matrix for plants:
comp_p = matrix(NA, nbsp_p, nbsp_p)
phen_p = matrix(NA, nbsp_p, nbsp_p)
for(i in 1:nbsp_p){
	mu1=invlogit1(subset(bidon,type=="mu" & species==(i+nbsp_a))$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==(i+nbsp_a))$value)
	for(j in i:(nbsp_p)){
		similarity=mean(sqrt(t(m)[,i] * t(m)[,j]))
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+nbsp_a))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+nbsp_a))$value)
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_p[i,j]=phen
		phen_p[j,i]=phen
		comp_p[i,j]=ifelse(tr %in% c("pheno","both"),phen*similarity,similarity)
		comp_p[j,i]=ifelse(tr %in% c("pheno","both"),phen*similarity,similarity)
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
for(competition in c(0.01,1,2,5)){
comp_a[,]=rho
comp_p[,]=rho
diag(comp_p)=1
diag(comp_a)=1
comp_phen_a=comp_a*phen_a
diag(comp_phen_a)=1
comp_phen_p=comp_p*phen_p
diag(comp_phen_p)=1
if(tr=="morpho"){
A=rbind(cbind(-1*competition*comp_p,m),cbind(t(m),-1*competition*comp_a))
}else{
A=rbind(cbind(-1*competition*comp_phen_p,m),cbind(t(m),-1*competition*comp_phen_a))
}

vec_stab=c()
for(it in 1:ncalc){
vec_stab=c(vec_stab,Omega(A))
}
ind$feas_with_pheno=mean(vec_stab)
ind$feas_with_pheno_se=sd(vec_stab)/sqrt(ncalc)
ind$competition=competition
ind$rho=rho
indf=rbind(indf,ind)
}

}
}
}
}
}


fwrite(indf,paste0("/home/duchenne/pheno/networks_info_empir_symmetric_",jj,".csv"))
#

#################################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

indf=NULL
for(jj in 1:17){
ind=fread(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/aggreg_empir_symmetric/networks_info_empir_symmetric_",jj,".csv"))

indf=rbind(indf,ind)
}

fwrite(indf,"networks_info_empir.csv")

