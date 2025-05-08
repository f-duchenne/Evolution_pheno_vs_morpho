###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)

source("scripts/functions and secondary scripts/toolbox.R")

weigthed.mean=function(x,y){sum(x*y,na.rm=T)/sum(y,na.rm=T)}

#calculate species phenologies from bi weekly counts
abund_counts=fread("data/empirical/interactions/weekly_counts.csv")
abund_counts$jj=yday(as.Date(abund_counts$date,format="%d/%m/%Y")) #create a column for the day of the year
phen_est=abund_counts %>% group_by(species,guild) %>% summarise(mu=weigthed.mean(jj,count),sde=sqrt(Hmisc::wtd.var(jj,count)))
phen_est=subset(phen_est,!is.na(mu))
phen_est$sde[is.na(phen_est$sde)]=1
phen_est$sde[phen_est$sde<1]=1
phena=subset(phen_est,guild=="poll")
phenp=subset(phen_est,guild=="plant")

load("data/empirical/interactions/matrices_empirical_networks.RData")

sites=names(networks)
ncalc=20
weigthed.mean=function(x,y){sum(x*y,na.rm=T)/sum(y,na.rm=T)}

library(doParallel)
# Register parallel backend
cores <- detectCores() - 1  # Use one less than available cores
registerDoParallel(cores)
final=foreach (site=sites ,.combine ="rbind")%dopar%{
	source("scripts/functions and secondary scripts/toolbox.R")
	library(bipartite)
	empf=NULL
	#mutualistic interaction matrix
	m=networks[[site]]
	
	#round for calculating network structure
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

	#build competition matrix and phenological overlap matrix for poll:
	liste=colnames(m)
	comp_a = matrix(NA, ncol(m), ncol(m))
	colnames(comp_a)=liste
	rownames(comp_a)=liste
	phen_a = matrix(NA, ncol(m), ncol(m))
	colnames(phen_a)=liste
	rownames(phen_a)=liste
	for(i in liste){
		mu1=subset(phena,species==i)$mu
		sd1=subset(phena,species==i)$sde
		for(j in liste){
			similarity=sum(m[,i] * m[,j])/((sum(m[,i])+sum(m[,j]))/2)
			mu2=subset(phena,species==j)$mu
			sd2=subset(phena,species==j)$sde
			f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
			phen=sum(f(seq(0,365,0.1)))*0.1
			phen_a[i,j]=phen
			comp_a[i,j]=similarity
		}  
	}
	diag(phen_a)=1
	diag(comp_a)=1.05

	#build competition matrix and phenological overlap matrix for plants:
	liste=rownames(m)
	comp_p = matrix(NA, nrow(m), nrow(m))
	colnames(comp_p)=liste
	rownames(comp_p)=liste
	phen_p = matrix(NA, nrow(m), nrow(m))
	colnames(phen_p)=liste
	rownames(phen_p)=liste
	for(i in liste){
		mu1=subset(phenp,species==i)$mu
		sd1=subset(phenp,species==i)$sde
		for(j in liste){
			similarity=sum(t(m)[,i] * t(m)[,j])/((sum(t(m)[,i])+sum(t(m)[,j]))/2)
			mu1=subset(phenp,species==j)$mu
			sd1=subset(phenp,species==j)$sde
			f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
			phen=sum(f(seq(0,365,0.1)))*0.1
			phen_p[i,j]=phen
			comp_p[i,j]=similarity
		}  
	}
	diag(phen_p)=1
	diag(comp_p)=1.05
	
	### PHENOLOGICAL OVERLAP AMONG PARTNERS
	listea=colnames(m)
	listep=rownames(m)
	phen_mut = matrix(NA, nrow(m), ncol(m))
	colnames(phen_mut)=listea
	rownames(phen_mut)=listep
	for(i in listea){
		mu1=subset(phena,species==i)$mu
		sd1=subset(phena,species==i)$sde
		for(j in listep){
			mu1=subset(phenp,species==j)$mu
			sd1=subset(phenp,species==j)$sde
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
	comp_phen_a=comp_a*phen_a
	comp_phen_p=comp_p*phen_p
	emp$comp_a=mean(comp_phen_a[col(comp_phen_a)!=row(comp_phen_a)])
	emp$comp_p=mean(comp_phen_p[col(comp_phen_p)!=row(comp_phen_p)])

	###### STRUCTURAL STABILITY
	for(competition in c(2,4,6)){

		A1=rbind(cbind(-1*competition*comp_phen_p,m),cbind(t(m),-1*competition*comp_phen_a))
		A2=rbind(cbind(-1*competition*comp_p,m),cbind(t(m),-1*competition*comp_a))

		vec_stab1=c()
		vec_stab2=c()
		for(it in 1:ncalc){
			vec_stab1=c(vec_stab1,Omega(A1))
			vec_stab2=c(vec_stab2,Omega(A2))
		}
		emp$feas_with_pheno=mean(vec_stab1,na.rm=TRUE)
		emp$feas_with_pheno_se=sd(vec_stab1,na.rm=TRUE)/sqrt(sum(!is.na(vec_stab1)))
		emp$feas_without_pheno=mean(vec_stab2,na.rm=TRUE)
		emp$feas_without_pheno_se=sd(vec_stab2,na.rm=TRUE)/sqrt(sum(!is.na(vec_stab2)))
		emp$competition=competition
		empf=rbind(empf,emp)
	}
	return(empf)
}

fwrite(final,"data/empirical/stability_of_empirical_networks.csv")



