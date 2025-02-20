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


datf=fread("/home/duchenne/pheno/species_level_data_new_alleg.csv")
datf$type="mu"
datf$type[grep("mu2",datf$variable)]="mu2"
datf$type[grep("sd",datf$variable)]="sd"
datf$type[grep("sd2",datf$variable)]="sd2"


ncalc=10
indf=NULL
motifs=NULL

for(competition in c(2,5)){
for (rich in c(10,20,30)){
for(tr in c("both","morpho","pheno")){
for(ti in plyr::round_any(seq(sqrt(0),sqrt(3000),length.out=10)^2,10)){
bidon=subset(datf,essai==ess & time==ti & trait==tr & rich==rich & comp==competition)

#build interaction matrix:
m = matrix(NA, rich, rich)
for(i in 1:rich){
	mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
	if(tr=="both"){
		mu1_m=invlogit1(subset(bidon,type=="mu2" & species==i)$value)
		sd1_m=invlogit(subset(bidon,type=="sd2" & species==i)$value)
	}
	for(j in 1:rich){
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+rich))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+rich))$value)
		if(tr=="both"){
			mu2_m=invlogit1(subset(bidon,type=="mu2" & species==j+rich)$value)
			sd2_m=invlogit(subset(bidon,type=="sd2" & species==j+rich)$value)
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
struf=networklevel(mbi,index=c("weighted connectance","weighted nestedness","interaction evenness","H2"),H2_integer=FALSE)
struf=as.data.frame(t(struf))
struf$NODF=networklevel(round(mbi,digits=5),index=c("weighted NODF"))
struf$modularity=computeModules(mbi, method="Beckett")@"likelihood"
struf$m_mean=mean(m)

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
	for(j in i:(rich)){
		similarity=mean(sqrt(m[,i] * m[,j]))
		mu2=invlogit1(subset(bidon,type=="mu" & species==j)$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==j)$value)
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_a[i,j]=phen
		phen_a[j,i]=phen
		comp_a[i,j]=similarity
		comp_a[j,i]=similarity
	}  
}
diag(comp_a)=1
diag(phen_a)=1

#build competition matrix for plants:
comp_p = matrix(NA, rich, rich)
phen_p = matrix(NA, rich, rich)
for(i in 1:rich){
	mu1=invlogit1(subset(bidon,type=="mu" & species==(i+rich))$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==(i+rich))$value)
	for(j in i:(rich)){
		similarity=mean(sqrt(t(m)[,i] * t(m)[,j]))
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+rich))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+rich))$value)
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_p[i,j]=phen
		phen_p[j,i]=phen
		comp_p[i,j]=similarity
		comp_p[j,i]=similarity
	}  
}
diag(comp_p)=1
diag(phen_p)=1

motifn=0
motifntot=0
for(i in 1:rich){
id_plant=which(m[,i]>=0.01)
if(length(id_plant)>1){
combis=as.data.frame(t(combn(id_plant,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_p[x[1],x[2]]}) #phenological overlap (competition) among pollinators
combis$mut=apply(combis,1,function(x){m[x[1],i]*m[x[2],i]}) #facilitation among pollinators
motifn=motifn+sum(combis$phen_overlap<combis$mut) #number of motif V+
motifntot=motifntot+nrow(combis) #number of motif V
}
}

for(i in 1:rich){
id_poll=which(m[i,]>=0.01)
if(length(id_poll)>1){
combis=as.data.frame(t(combn(id_poll,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_a[x[1],x[2]]}) #phenological overlap (competition) among plant
combis$mut=apply(combis,1,function(x){m[i,x[1]]*m[i,x[2]]}) #facilitation among plants
motifn=motifn+sum(combis$phen_overlap<combis$mut) #number of motif V+
motifntot=motifntot+nrow(combis) #number of motif V
}
}

ind$motifn=motifn
ind$motifntot=motifntot


ind$competition=competition

######################### CALCULATE FEASIBILITY WITH STRCUTRED COMPETITION NETWORKS
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
A=rbind(cbind(-1*competition*comp_phen_p,m),cbind(t(m),-1*competition*comp_phen_a))
vec_stab=c()
for(it in 1:ncalc){vec_stab=c(vec_stab,Omega(A))}
ind$feas_structure=mean(vec_stab)
ind$feas_structure_se=sd(vec_stab)/sqrt(ncalc)

######################### CALCULATE FEASIBILITY WITH MEAN FIELD BUT MEAN VALUES FROM AVERAGE COMPETITION STRENGTH IN SIMULATIONS
ind$comp_a=mean(comp_phen_a[col(comp_phen_a)!=row(comp_phen_a)])
ind$comp_p=mean(comp_phen_p[col(comp_phen_p)!=row(comp_phen_p)])

comp_phen_a[,]=mean(comp_phen_a[col(comp_phen_a)!=row(comp_phen_a)])
comp_phen_p[,]=mean(comp_phen_p[col(comp_phen_p)!=row(comp_phen_p)])
diag(comp_phen_a)=1
diag(comp_phen_p)=1
A=rbind(cbind(-1*competition*comp_phen_p,m),cbind(t(m),-1*competition*comp_phen_a))
vec_stab=c()
for(it in 1:ncalc){vec_stab=c(vec_stab,Omega(A))}
ind$feas_comp=mean(vec_stab)
ind$feas_comp_se=sd(vec_stab)/sqrt(ncalc)

######################### CALCULATE FEASIBILITY WITH FEW DIFFERENT MEAN FIELDS
for(competition_feas in c(0.1,1,2,5)){
for(rho in c(0.01,0.05,0.2,0.4,0.6)){
comp_a[,]=rho #mean(c(comp_p,comp_a))
comp_p[,]=rho #mean(c(comp_p,comp_a))
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

A=rbind(cbind(-1*competition_feas*comp_phen_p,m),cbind(t(m),-1*competition_feas*comp_phen_a))
vec_stab=c()
for(it in 1:ncalc){vec_stab=c(vec_stab,Omega(A))}
ind$feas=mean(vec_stab)
ind$feas_se=sd(vec_stab)/sqrt(ncalc)
ind$competition_feas=competition_feas
ind$rho=rho #mean(c(comp_p,comp_a)) 
indf=rbind(indf,ind)
}}

}
}
}
}

fwrite(indf,paste0("/home/duchenne/pheno/aggreg_simues/networks_info_",ess,".csv"))
#
