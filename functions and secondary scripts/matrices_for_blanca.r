###################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales","piecewiseSEM",
"igraph","qgraph","car","nlme","viridis") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

colo=c("black","orange")
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)

##################################### SPECIES LEVEL ########################
datf=fread("data/simulated/outputs_simulations/species_level.csv")
competition=4
riche=20

datf=subset(datf,comp==competition & rich==riche & essai==2 & time<=500)
datf$type2="pheno"
datf$type2[datf$type %in% c("mu2","sd2")]="morpho"
datf$type2[datf$trait %in% c("morpho")]="morpho"
datf$type3="mu"
datf$type3[datf$type %in% c("sd","sd2")]="sd"

time_vec = unique(datf$time)

a=1
lili=list()
for(tr in c("pheno","morpho","both")){
	for(ti in time_vec){
	bidon=subset(datf,time==ti & trait==tr)
	if(tr=="both"){
		bidon$guild=rep(rep(c("poll","plant"),each=20),4)
	}else{
		bidon$guild=rep(rep(c("poll","plant"),each=20),2)
	}
	#build interaction matrix:
	m = matrix(NA, riche, riche)
	for(i in 1:riche){
		mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
		sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
		if(tr=="both"){
			mu1_m=invlogit1(subset(bidon,type=="mu2" & species==i)$value)
			sd1_m=invlogit(subset(bidon,type=="sd2" & species==i)$value)
		}
		for(j in 1:riche){
			mu2=invlogit1(subset(bidon,type=="mu" & species==(j+riche))$value)
			sd2=invlogit(subset(bidon,type=="sd" & species==(j+riche))$value)
			if(tr=="both"){
				mu2_m=invlogit1(subset(bidon,type=="mu2" & species==j+riche)$value)
				sd2_m=invlogit(subset(bidon,type=="sd2" & species==j+riche)$value)
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

	#build competition matrix for poll:
	comp_a = matrix(NA, riche, riche)
	phen_a = matrix(NA, riche, riche)
	for(i in 1:riche){
		mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
		sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
		for(j in 1:(riche)){
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
	comp_p = matrix(NA, riche, riche)
	phen_p = matrix(NA, riche, riche)
	for(i in 1:riche){
		mu1=invlogit1(subset(bidon,type=="mu" & species==(i+riche))$value)
		sd1=invlogit(subset(bidon,type=="sd" & species==(i+riche))$value)
		for(j in 1:(riche)){
			similarity=sum(t(m)[,i] * t(m)[,j])/((sum(t(m)[,i])+sum(t(m)[,j]))/2)
			mu2=invlogit1(subset(bidon,type=="mu" & species==(j+riche))$value)
			sd2=invlogit(subset(bidon,type=="sd" & species==(j+riche))$value)
			f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
			phen=sum(f(seq(0,365,0.1)))*0.1
			phen_p[i,j]=phen
			comp_p[i,j]=similarity
		}  
	}
	diag(comp_p)=1
	diag(phen_p)=1
	
	if(tr %in% c("pheno","both")){
		comp_af=competition*comp_a*phen_a
		comp_pf=competition*comp_p*phen_p
	}else{
		comp_af=competition*comp_a
		comp_pf=competition*comp_p
	}
	
	lili[[paste0(tr)]][[paste0("time: ",ti)]][["mutualism"]]=m
	lili[[paste0(tr)]][[paste0("time: ",ti)]][["competition poll"]]=comp_af
	lili[[paste0(tr)]][[paste0("time: ",ti)]][["competition plants"]]=comp_pf
	}
}
	