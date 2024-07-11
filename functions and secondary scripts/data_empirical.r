###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")
###### FLOWER PHENOLOGY:
flowers=fread("flower_count.csv")
flowers$date=as.Date(paste(flowers$Day,flowers$Month,flowers$Year,sep="/"),format="%d/%m/%Y")
flowers$jj=yday(flowers$date)
tab=unique(flowers[,c("Site_ID")])

empf=NULL
weigthed.mean=function(x,y){sum(x*y,na.rm=T)/sum(y,na.rm=T)}
phenp=flowers %>% group_by(Plant_gen_sp) %>% summarise(mu=weigthed.mean(jj,Flower_abundance),sde=sqrt(Hmisc::wtd.var(jj,Flower_abundance)))

fwrite(phenp,"flower_pheno_empirical.csv")

######### INTERACTIONS:
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/scripts")
source("toolbox.R")
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")
datf=as.data.frame(fread("BeeFun Master_inter.csv"))
datf$date=as.Date(paste(datf$Day,datf$Month,datf$Year,sep="/"),format="%d/%m/%Y")
datf$jj=yday(datf$date)

#pheno plants:
phenp=fread("flower_pheno_empirical.csv")
phenp$sde[is.na(phenp$sde)]=1
phenp$sde[phenp$sde==0]=1

#pheno poll:
phena=datf %>% group_by(Pollinator_gen_sp) %>% summarise(mu=weigthed.mean(jj,Frequency),sde=sqrt(Hmisc::wtd.var(jj,Frequency)))
phena=subset(phena,!is.na(mu))
phena$sde[is.na(phena$sde)]=1
phena$sde[phena$sde==0]=1

sites=unique(datf$Site_ID)
sites=sites[!(sites %in% c("El_pozo","El_Pozo"))]
weigthed.mean=function(x,y){sum(x*y,na.rm=T)/sum(y,na.rm=T)}
empf=NULL
for (site in sites){

bidon=subset(datf,Site_ID==site)
bidon=subset(bidon,!is.na(Plant_genus) & Plant_genus!="" & Pollinator_genus!="")
bidon=subset(bidon, Pollinator_gen_sp %in% phena$Pollinator_gen_sp)
bidon=subset(bidon, Plant_gen_sp %in% phenp$Plant_gen_sp)

nsampl=length(unique(paste0(bidon$Round,bidon$Year)))

m2=as.data.frame(dcast(bidon,Plant_gen_sp~Pollinator_gen_sp,value.var="Frequency",fun.aggregate =sum,na.rm=T))
m2[is.na(m2)]=0
m=as.matrix(m2[,-1])
rownames(m)=m2[,1]
m=m/nsampl
m=m[apply(m,1,sum)>0,apply(m,2,sum)>0]

emp=networklevel(m,index=c("weighted connectance","weighted nestedness","interaction evenness","H2"))
emp=as.data.frame(t(emp))
emp$NODF=networklevel(round(m,digits=5),index=c("weighted NODF"))
names(emp)=gsub("weighted ","",names(emp))
emp$modularity=computeModules(m, method="Beckett")@"likelihood"
emp$na=ncol(m)
emp$np=nrow(m)
emp$nround=nsampl
emp$nyear=length(unique(bidon$Year))
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
		similarity=sum(m[,i] * m[,j]) / sum(m[,i])
		mu2=subset(phena,Pollinator_gen_sp==j)$mu
		sd2=subset(phena,Pollinator_gen_sp==j)$sde
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_a[i,j]=phen
		comp_a[i,j]=similarity
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
		similarity=sum(t(m)[,i] * t(m)[,j]) / sum(t(m)[,i])
		mu1=subset(phenp,Plant_gen_sp==j)$mu
		sd1=subset(phenp,Plant_gen_sp==j)$sde
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_p[i,j]=phen
		comp_p[i,j]=similarity
	}  
}
diag(comp_p)=1

comp_phen_a=comp_a*phen_a
diag(comp_phen_a)=1
comp_phen_p=comp_p*phen_p
diag(comp_phen_p)=1

A1=rbind(cbind(10*comp_phen_p,-1*m),cbind(-1*t(m),10*comp_phen_a))
A2=rbind(cbind(10*comp_p,-1*m),cbind(-1*t(m),10*comp_a))

emp$feas_with_pheno=Omega(A1)
emp$feas_without_pheno=Omega(A2)

motifn=0
motifntot=0
liste=colnames(m)
for(i in liste){
id_plant=which(m[,i]>0)
if(length(id_plant)>1){
combis=as.data.frame(t(combn(id_plant,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_p[x[1],x[2]]})
combis$mut=apply(combis,1,function(x){sqrt(m[x[1],i]*m[x[2],i])})
motifn=motifn+sum(combis$phen_overlap<0.01)
motifntot=motifntot+nrow(combis)
}
}

liste=rownames(m)
for(i in liste){
id_poll=which(m[i,]>0)
if(length(id_poll)>1){
combis=as.data.frame(t(combn(id_poll,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_a[x[1],x[2]]})
combis$mut=apply(combis,1,function(x){sqrt(m[i,x[1]]*m[i,x[2]])})
motifn=motifn+sum(combis$phen_overlap<0.01)
motifntot=motifntot+nrow(combis)
}
}

emp$motifn=motifn
emp$motifntot=motifntot


empf=rbind(empf,emp)

}

fwrite(empf,"empirical_networks.csv")


b=melt(empf, measure.vars=c("feas_without_pheno","feas_with_pheno"))
b$trait="morpho"
b$trait[b$variable=="feas_with_pheno"]="pheno"

colo=c("dodgerblue4","chartreuse3")

pl3=ggplot(data=empf,aes(x=feas_without_pheno,y=feas_with_pheno))+
geom_point()+
geom_abline(intercept=0,slope=1)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+ggtitle("c")














setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
indf=fread("C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info.csv")
empf=fread("empirical_networks.csv")


acp=PCA(rbind(indf[,1:6],empf[,1:6]),graph=FALSE)
indf2=cbind(indf,acp$ind$coord[1:nrow(indf),])
empf2=cbind(empf,acp$ind$coord[(nrow(indf)+1):nrow(acp$ind$coord),])




pl3=ggplot()+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point(data=subset(indf2,time==1500),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,shape=as.factor(rich)))+
geom_point(data=subset(indf2,time==0),aes(x=Dim.1,y=Dim.2,alpha=time,shape=as.factor(rich)),color="grey",alpha=0.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("c")+
geom_point(data=empf2,aes(x=Dim.1,y=Dim.2),color="black")

