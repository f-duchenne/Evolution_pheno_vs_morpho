###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/scripts")

source("toolbox.R")

invlogit=function(x){return(1+50*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+285*exp(x)/(1+exp(x)))}

for (rich in c(10,20)){
datf=NULL
for(i in 1:50){
dat=fread(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/results/ueq_",i,"_",rich,".csv"))

dive=rich*2
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
dat2$rich=rich
datf=rbind(datf,dat2)
}
}

indf=NULL
motifs=NULL
for (rich in c(10,20)){
for(ess in c(1:50)){
for(tr in c("morpho","pheno")){
for(ti in c(0,100,250,500,1000,1500)){
bidon=subset(datf,essai==ess & time==ti & trait==tr & rich==rich)

#build interaction matrix:
m = matrix(NA, rich, rich)
for(i in 1:rich){
	mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
	for(j in 1:rich){
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+rich))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+rich))$value)
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
ind$rich=rich



#build competition matrix for plants:
comp_a = matrix(NA, rich, rich)
phen_a = matrix(NA, rich, rich)
for(i in 1:rich){
	mu1=invlogit1(subset(bidon,type=="mu" & species==i)$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==i)$value)
	for(j in 1:(rich)){
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
comp_p = matrix(NA, rich, rich)
phen_p = matrix(NA, rich, rich)
for(i in 1:rich){
	mu1=invlogit1(subset(bidon,type=="mu" & species==(i+rich))$value)
	sd1=invlogit(subset(bidon,type=="sd" & species==(i+rich))$value)
	for(j in 1:(rich)){
		similarity=sum(t(m)[,i] * t(m)[,j]) / sum(t(m)[,i])
		mu2=invlogit1(subset(bidon,type=="mu" & species==(j+rich))$value)
		sd2=invlogit(subset(bidon,type=="sd" & species==(j+rich))$value)
		f=function(x){pmin(dnorm(x,mu1,sd1),dnorm(x,mu2,sd2))}
		phen=sum(f(seq(0,365,0.1)))*0.1
		phen_p[i,j]=phen
		comp_p[i,j]=ifelse(tr=="pheno",phen*similarity,similarity)
	}  
}
diag(comp_p)=1

A=rbind(cbind(10*comp_p,-1*m),cbind(-1*m,10*comp_a))
ind$feas=Omega(A)


motifn=0
motifntot=0
for(i in 1:rich){
id_plant=which(m[,i]>0.01)
if(length(id_plant)>1){
combis=as.data.frame(t(combn(id_plant,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_p[x[1],x[2]]})
combis$mut=apply(combis,1,function(x){sqrt(m[x[1],i]*m[x[2],i])})
motifn=motifn+sum(combis$phen_overlap<combis$mut)
motifntot=motifntot+nrow(combis)
}
}

for(i in 1:rich){
id_poll=which(m[i,]>0.01)
if(length(id_poll)>1){
combis=as.data.frame(t(combn(id_poll,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_a[x[1],x[2]]})
combis$mut=apply(combis,1,function(x){sqrt(m[i,x[1]]*m[i,x[2]])})
motifn=motifn+sum(combis$phen_overlap<combis$mut)
motifntot=motifntot+nrow(combis)
}
}

ind$motifn=motifn
ind$motifntot=motifntot
ind$comp=mean(c(comp_p,comp_a))


indf=rbind(indf,ind)

}
}
}
}

fwrite(indf,"C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info_10.csv")


###################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
indf=fread("C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info_10.csv")

colo=c("dodgerblue3","chartreuse3")


p1=ggplot()+
geom_boxplot(data=indf,aes(y=motifn/motifntot,x=time,group=paste0(trait,time),color=trait))+
stat_smooth(data=indf,aes(y=motifn/motifntot,x=time,color=trait,fill=trait))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("Proportion of 'V' shape motifs without overlap")+
xlab("Time")+ggtitle("a")

p2=ggplot()+
geom_boxplot(data=indf,aes(y=comp,x=time,group=paste0(trait,time),color=trait))+
stat_smooth(data=indf,aes(y=comp,x=time,color=trait,fill=trait))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("Competition for mutualistic partners")+
xlab("Time")+ggtitle("b")


grid.arrange(p1,p2,ncol=2)
pdf("Fig.1",width=9,height=5)
grid.arrange(p1,p2,ncol=2)
dev.off();


acp=PCA(indf[,1:6],plot=FALSE)
indf2=cbind(indf,acp$ind$coord)

pl1=ggplot(data=indf2,aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,group=paste(essai,trait)))+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point()+
geom_line(show_guide =F)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("a")+coord_fixed(ratio=1)

pl2=fviz_pca_var(acp, col.var = "black", repel =TRUE)+ggtitle("b")+
theme(plot.title=element_text(size=14,face="bold",hjust = 0))

pl3=ggplot()+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point(data=subset(indf2,time==1500),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait))+
geom_point(data=subset(indf2,time==0),aes(x=Dim.1,y=Dim.2,alpha=time),color="grey",alpha=0.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("c")

top=plot_grid(pl1,pl2,align="hv",ncol=2,rel_widths=c(1.5,1))

grid.arrange(top,pl3,ncol=1)


ggplot(data=indf2,aes(x=as.factor(time),y=feas,color=trait))+geom_violin(position=position_dodge(width=0.5))+geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)

