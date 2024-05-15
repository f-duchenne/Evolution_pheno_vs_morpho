###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/scripts")

source("toolbox.R")

invlogit=function(x){return(1+50*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+285*exp(x)/(1+exp(x)))}

dive=10*2
datf=NULL
for(i in 1:10){
dat=fread(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/ueq_",i,".csv"))

names(dat)[1:(dive)]=gsub("x","mu_",names(dat)[1:(dive)])
names(dat)[(dive+1):(dive*2)]=gsub("x","sd_",names(dat)[(dive+1):(dive*2)])
names(dat)[(dive*2+1)]="time"

dat2=melt(dat,id.vars=c("trait","time"))
dat2$type="mu"
dat2$type[grep("sd",dat2$variable)]="sd"
dat2$species=gsub("mu_","",dat2$variable)
dat2$species=as.numeric(gsub("sd_","",dat2$species))
dat2$species[dat2$species>dive]=dat2$species[dat2$species>dive]-dive
dat2$essai=i
datf=rbind(datf,dat2)
}


indf=NULL
for(ess in c(1:5)){
for(tr in c("morpho","pheno")){
for(ti in c(0,10,100,300,800,1500)){
bidon=subset(datf,essai==ess & time==ti & trait==tr)

#build interaction matrix:
m = matrix(NA, dive/2, dive/2)
for(i in 1:(dive/2)){
	mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
	for(j in 1:(dive/2)){
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+dive/2))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+dive/2))$value)
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		m[j,i]=sum(f(seq(0,365,0.1)))*0.1
	}  
}
   
ind=networklevel(m,index=c("weighted connectance","weighted nestedness","interaction evenness","H2"))
ind=as.data.frame(t(ind))
ind$NODF=networklevel(round(m,digits=5),index=c("weighted NODF"))
names(ind)=gsub("weighted ","",names(ind))
ind$modularity=computeModules(m, method="Beckett")@"likelihood"
ind$trait=tr
ind$time=ti
ind$essai=ess
indf=rbind(indf,ind)


#build competition matrix for plants:
comp_a = matrix(NA, dive/2, dive/2)
for(i in 1:(dive/2)){
	mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
	for(j in 1:(dive/2)){
		similarity=sum(m[,i] * m[,j]) / sum(m[,i])
		if(tr=="phen"){
			mu2=invlogit1(subset(bidon,type=="mu" & species==j)$value)
			sd2=invlogit(subset(bidon,type=="sd" & species==j)$value)
			f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
			phen=sum(f(seq(0,365,0.1)))*0.1
		}else{
			phen=1
		}
		comp_a[i,j]=phen*similarity	
	}  
}
diag(comp_a)=1

#build competition matrix for plants:
comp_p = matrix(NA, dive/2, dive/2)
for(i in 1:(dive/2)){
	mu1=invlogit1(subset(bidon,type=="mu" & species==(i+dive/2))$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==(i+dive/2))$value)
	for(j in 1:(dive/2)){
		similarity=sum(t(m)[,i] * t(m)[,j]) / sum(t(m)[,i])
		if(tr=="phen"){
			mu2=invlogit1(subset(bidon,type=="mu" & species==(j+dive/2))$value)
			sd2=invlogit(subset(bidon,type=="sd" & species==(j+dive/2))$value)
			f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
			phen=sum(f(seq(0,365,0.1)))*0.1
		}else{
			phen=1
		}
		comp_p[i,j]=phen*similarity	
	}  
}
diag(comp_p)=1

A=rbind(cbind(comp_p,-1*m),cbind(-1*m,comp_a))
Omega(A)

}
}
}

acp=PCA(indf[,1:6])

indf2=cbind(indf,acp$ind$coord)

colo=c("dodgerblue3","chartreuse3")

pl1=ggplot(data=indf2,aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,group=paste(essai,trait)))+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point()+
geom_line(show_guide =F)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("a")

pl2=fviz_pca_var(acp, col.var = "black", repel =TRUE)+ggtitle("b")

pl3=ggplot(data=subset(indf2,time==1500),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,group=paste(essai,trait)))+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("c")

top=plot_grid(pl1,pl2,align="hv",ncol=2)

grid.arrange(top,pl3,ncol=1)
