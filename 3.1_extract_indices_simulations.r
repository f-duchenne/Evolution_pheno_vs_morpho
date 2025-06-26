############################### AGGREGATE SIMULATIONS
###############################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)



comp_vec=c(2,4,6)
time_vec=plyr::round_any(seq(sqrt(0),sqrt(2000),length.out=10)^2,10)

##### AGGREGATE SIMULATIONS WITH TWO TRAITS
datf=NULL
for(competition in comp_vec){
	for (rich in c(10,20,30)){
		for(i in 1:100){
			dat=fread(paste0("data/simulated/outputs_simulations/results_for_each_individual_simulations/ueq_both_",i,"_",rich,"_",competition,".csv"))
			dive=rich*2
			names(dat)[1:(dive)]=paste0("mu_",1:dive)
			names(dat)[(dive+1):(dive*2)]=paste0("sd_",1:dive)
			names(dat)[(dive*2+1):(dive*3)]=paste0("mu2_",1:dive)
			names(dat)[(dive*3+1):(dive*4)]=paste0("sd2_",1:dive)
			names(dat)[(dive*4+1)]="time"

			dat2=melt(dat,id.vars=c("trait","time"))
			dat2$type="mu"
			dat2$type[grep("mu2",dat2$variable)]="mu2"
			dat2$type[grep("sd",dat2$variable)]="sd"
			dat2$type[grep("sd2",dat2$variable)]="sd2"
			dat2$species=gsub("mu2_","",gsub("mu_","",dat2$variable))
			dat2$species=as.numeric(gsub("sd2_","",gsub("sd_","",dat2$species)))
			dat2$comp=competition
			dat2$essai=i
			dat2$rich=rich
			datf=rbind(datf,dat2)
		}
	}
}
##### AGGREGATE SIMULATIONS WITH ONE TRAITS
for(competition in comp_vec){
	for (rich in c(10,20,30)){
		for(i in 1:100){
			dat=fread(paste0("data/simulated/outputs_simulations/results_for_each_individual_simulations/ueq_",i,"_",rich,"_",competition,".csv"))
			dive=rich*2
			names(dat)[1:(dive)]=paste0("mu_",1:dive)
			names(dat)[(dive+1):(dive*2)]=paste0("sd_",1:dive)
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

fwrite(datf,"data/simulated/outputs_simulations/species_level.csv")
datf_alleg=datf[datf$time %in% time_vec,]
fwrite(datf_alleg,"data/simulated/outputs_simulations/species_level_alleg.csv")

########################################### EXTRACT INDICES AND EXPORT INFORMATION AT COMMUNITY LEVEL FOR EACH REPLICATE
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

datf=fread("/home/duchenne/pheno/species_level_alleg.csv")
datf$type="mu"
datf$type[grep("mu2",datf$variable)]="mu2"
datf$type[grep("sd",datf$variable)]="sd"
datf$type[grep("sd2",datf$variable)]="sd2"


ncalc=10
indf=NULL
motifs=NULL

comp_vec=c(2,4,6)
time_vec=unique(datf$time)

for(competition in comp_vec){
	for (rich in c(10,20,30)){
		for(tr in c("pheno","morpho","both")){
			for(ti in time_vec){
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
					for(j in 1:(rich)){
						similarity=sum(m[,i] * m[,j])/((sum(m[,i])+sum(m[,j]))/2)
						mu2=invlogit1(subset(bidon,type=="mu" & species==j)$value)
						sd2=invlogit(subset(bidon,type=="sd" & species==j)$value)
						f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
						phen=sum(f(seq(0,365,0.1)))*0.1
						phen_a[i,j]=phen
						comp_a[i,j]=similarity
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
					for(j in 1:(rich)){
						similarity=sum(t(m)[,i] * t(m)[,j])/((sum(t(m)[,i])+sum(t(m)[,j]))/2)
						mu2=invlogit1(subset(bidon,type=="mu" & species==(j+rich))$value)
						sd2=invlogit(subset(bidon,type=="sd" & species==(j+rich))$value)
						f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
						phen=sum(f(seq(0,365,0.1)))*0.1
						phen_p[i,j]=phen
						comp_p[i,j]=similarity
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
        
        ind$comp_a=mean(comp_phen_a[col(comp_phen_a)!=row(comp_phen_a)])
				ind$comp_p=mean(comp_phen_p[col(comp_phen_p)!=row(comp_phen_p)])
        
				A=rbind(cbind(-1*competition*comp_phen_p,m),cbind(t(m),-1*competition*comp_phen_a))
				vec_stab=c()
				for(it in 1:ncalc){vec_stab=c(vec_stab,Omega(A))}
				ind$feas_structure=mean(vec_stab)
				ind$feas_structure_se=sd(vec_stab)/sqrt(ncalc)
        indf=rbind(indf,ind)
        
			}
		}
	}
}

fwrite(indf,paste0("/home/duchenne/pheno/aggreg_simues_symmetric/networks_info_",ess,".csv"))
#
###############################
#################### AGGREGATE NETWORKS INFORMATION IN A UNIQUE TABLE
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
for(ess in 1:100){
	ind=fread(paste0("data/simulated/outputs_simulations/results_community_level_each_replicate/networks_info_",ess,".csv"))
	indf=rbind(indf,ind)
}
fwrite(indf,"data/simulated/outputs_simulations/networks_info.csv")


