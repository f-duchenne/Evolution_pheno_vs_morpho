############################### AGGREGATE SIMULATIONS
###############################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","bipartite","ggplot2") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)

invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

tr="both"
final=NULL
for(essai in 1:100){
	for (rich in c(10,20,30)){
		bidon=fread(paste0("data/simulated/initial_conditions_simulations/","pops_ini_",essai,".csv"))
		bidon=subset(bidon,sp %in% c(paste0("a",1:rich),paste0("f",1:rich)))
		
		#build interaction matrix:
		m = matrix(NA, rich, rich)
		for(i in 1:rich){
			mu1=invlogit1(subset(bidon,sp==paste0("a",i))$mu_phen)
			sd1=invlogit(subset(bidon,sp==paste0("a",i))$sd_phen)
			if(tr=="both"){
				mu1_m=invlogit1(subset(bidon,sp==paste0("a",i))$mu_morpho)
				sd1_m=invlogit(subset(bidon,sp==paste0("a",i))$sd_morpho)
			}
			for(j in 1:rich){
				mu2=invlogit1(subset(bidon,sp==paste0("f",j))$mu_phen)
				sd2=invlogit(subset(bidon,sp==paste0("f",j))$sd_phen)
				if(tr=="both"){
					mu2_m=invlogit1(subset(bidon,sp==paste0("f",j))$mu_morpho)
					sd2_m=invlogit(subset(bidon,sp==paste0("f",j))$sd_morpho)
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
		#struf$modularity=computeModules(mbi, method="Beckett")@"likelihood"
		struf$m_mean=mean(m)

		ind=as.data.frame(t(apply(struf,2,mean,na.rm=T)))
		ind$trait=tr
		ind$essai=essai
		ind$rich=rich

		#build competition matrix for poll:
		comp_a = matrix(NA, rich, rich)
		phen_a = matrix(NA, rich, rich)
		for(i in 1:rich){
			mu1=invlogit1(subset(bidon,sp==paste0("a",i))$mu_phen)
			sd1=invlogit(subset(bidon,sp==paste0("a",i))$sd_phen)
			for(j in i:(rich)){
				similarity=sum(m[,i] * m[,j])/sum(m[,i])
				mu2=invlogit1(subset(bidon,sp==paste0("a",j))$mu_phen)
				sd2=invlogit(subset(bidon,sp==paste0("a",j))$sd_phen)
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
			mu1=invlogit1(subset(bidon,sp==paste0("f",i))$mu_phen)
			sd1=invlogit(subset(bidon,sp==paste0("f",i))$sd_phen)
			for(j in i:(rich)){
				similarity=sum(t(m)[,i] * t(m)[,j])/sum(t(m)[,i])
				mu2=invlogit1(subset(bidon,sp==paste0("f",j))$mu_phen)
				sd2=invlogit(subset(bidon,sp==paste0("f",j))$sd_phen)
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
		
		#COMPETITION NETWORKS
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
		for(competition in c(2,3,4,6)){
		final=rbind(final,data.frame(essai=essai,rich=rich,ratio=mean(m)/(competition*mean(ind$comp_a,ind$comp_p)),competition=competition))
		}
	}
}



ggplot(data=final,aes(x=as.factor(competition),y=ratio))+geom_boxplot()+geom_hline(yintercept=1,linetype="dashed")+facet_wrap(~rich)