############################### AGGREGATE SIMULATIONS DONE ON EMPIRICAL NETWORKS
#' Check for packages and if necessary install into library 
#+ message = FALSE
pkgs <- c("data.table", "dplyr") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)


path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)

source("scripts/functions and secondary scripts/toolbox.R")
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

liste = fread("data/empirical/initial_conditions_simulations/liste.csv")

lili=list.files(path="data/empirical/outputs_simulations/results_for_each_individual_simulations/")

competition=4

datf=NULL
for(ess in 1:10){
	for(tr in c("both")){
		for(jj in 1:nrow(liste)){
			if(paste0("ueq_",liste$site[jj],"_comp_",competition,"_",tr,"_",ess,"_asym.csv") %in% lili){

				dat=fread(paste0("data/empirical/outputs_simulations/results_for_each_individual_simulations/ueq_",liste$site[jj],"_comp_",competition,"_",tr,"_",ess,".csv"))
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

				dat2=melt(dat,id.vars=c("trait","time","competition"))
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
}
fwrite(datf,"data/empirical/outputs_simulations/species_level_simues_empir.csv")
datf_alleg=datf[datf$time %in% c(0,2000),]
fwrite(datf_alleg,"data/empirical/outputs_simulations/species_level_simues_empir_alleg.csv")

###########################################  EXTRACT INDICES AND EXPORT INFORMATION AT COMMUNITY LEVEL FOR EACH SITE
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

datf = fread("/home/duchenne/pheno/species_level_simues_empir_alleg_asym.csv")
liste = fread("/home/duchenne/pheno/initial_empir/liste.csv")
ncalc=10
indf=NULL



motifs=NULL
for(ess in 1:10){
	nbsp_a=liste$na[jj]
	nbsp_p=liste$np[jj]
	for(competition in c(4)){
		for(tr in c("both")){
			for(ti in c(0,2000)){
				bidon=subset(datf,time==ti & trait==tr & essai==ess & site==liste$site[jj] & comp==competition)

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
					for(j in 1:nbsp_a){
						similarity=sum(m[,i] * m[,j])/sum(m[,i])
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
        diag(phen_a)=1

				#build competition matrix for plants:
				comp_p = matrix(NA, nbsp_p, nbsp_p)
				phen_p = matrix(NA, nbsp_p, nbsp_p)
				for(i in 1:nbsp_p){
					mu1=invlogit1(subset(bidon,type=="mu" & species==(i+nbsp_a))$value)
					sd1=invlogit(subset(bidon,type=="sd" & species==(i+nbsp_a))$value)
					for(j in 1:(nbsp_p)){
						similarity=sum(t(m)[,i] * t(m)[,j])/sum(t(m)[,i])
						mu2=invlogit1(subset(bidon,type=="mu" & species==(j+nbsp_a))$value)
						sd2=invlogit(subset(bidon,type=="sd" & species==(j+nbsp_a))$value)
						f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
						phen=sum(f(seq(0,365,0.1)))*0.1
						phen_p[i,j]=phen
						comp_p[i,j]=ifelse(tr %in% c("pheno","both"),phen*similarity,similarity)
					}  
				}
				diag(comp_p)=1
        diag(phen_p)=1

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

				A=rbind(cbind(-1*competition*comp_p,m),cbind(t(m),-1*competition*comp_a))
				vec_stab=c()
				for(it in 1:ncalc){vec_stab=c(vec_stab,Omega(A))}

				ind$feas_with_pheno=mean(vec_stab)
				ind$feas_with_pheno_se=sd(vec_stab)/sqrt(ncalc)
				ind$competition=competition
				indf=rbind(indf,ind)
			}
		}
	}
}




fwrite(indf,paste0("/home/duchenne/pheno/networks_info_empir_",jj,".csv"))
#

#################################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)


indf=NULL
for(jj in 1:17){
	ind=fread(paste0("data/empirical/outputs_simulations/results_community_level_each_site/networks_info_empir_",jj,".csv"))
	indf=rbind(indf,ind)
}

fwrite(indf,"data/empirical/outputs_simulations/networks_info_empir.csv")