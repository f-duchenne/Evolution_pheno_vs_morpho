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

invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

datf=NULL
for (rich in c(10,20,30)){
for(i in 1:100){
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

fwrite(datf,"species_level_data.csv")


indf=NULL
motifs=NULL
for (rich in c(10,20,30)){
for(ess in c(1:100)){
for(tr in c("morpho","pheno")){
for(ti in c(0,100,250,500,750,1000,1500)){
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



#build competition matrix for poll:
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

A=rbind(cbind(10*comp_p,-1*m),cbind(-1*t(m),10*comp_a))
ind$feas=Omega(A)


motifn=0
motifntot=0
for(i in 1:rich){
id_plant=which(m[,i]>0.01)
if(length(id_plant)>1){
combis=as.data.frame(t(combn(id_plant,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_p[x[1],x[2]]})
combis$mut=apply(combis,1,function(x){sqrt(m[x[1],i]*m[x[2],i])})
motifn=motifn+sum(combis$phen_overlap<0.01)
motifntot=motifntot+nrow(combis)
}
}

for(i in 1:rich){
id_poll=which(m[i,]>0.01)
if(length(id_poll)>1){
combis=as.data.frame(t(combn(id_poll,m=2)))
combis$phen_overlap=apply(combis,1,function(x){phen_a[x[1],x[2]]})
combis$mut=apply(combis,1,function(x){sqrt(m[i,x[1]]*m[i,x[2]])})
motifn=motifn+sum(combis$phen_overlap<0.01)
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

fwrite(indf,"C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info.csv")


###################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales","piecewiseSEM","igraph","qgraph","car") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
indf=fread("C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info.csv")

colo=c("dodgerblue4","chartreuse3")

endpo=subset(datf,time==2000)


ggplot()+
geom_density(data=subset(endpo,type=="mu"),aes(x=invlogit1(value),color=trait,fill=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("Proportion of 'V' shape motifs without overlap")+
xlab("Time")+ggtitle("b")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_wrap(~rich)

ggplot()+
geom_density(data=subset(endpo,type=="sd"),aes(x=invlogit(value),color=trait,fill=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("Proportion of 'V' shape motifs without overlap")+
xlab("Time")+ggtitle("b")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_wrap(~rich)


cortab=dcast(endpo,rich+species+essai+trait~type)

ggplot()+
geom_point(data=cortab,aes(x=invlogit1(mu),y=invlogit(sd),color=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("Proportion of 'V' shape motifs without overlap")+
xlab("Time")+ggtitle("b")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_wrap(~rich)

p1=ggplot()+
#geom_boxplot(data=indf,aes(y=motifn/motifntot,x=time,group=paste0(trait,time,rich),color=trait,linetype=as.factor(rich)),width=100,show_guide=FALSE)+
stat_smooth(data=indf,aes(y=motifn/motifntot,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("Proportion of 'V' shape motifs without overlap")+
xlab("Time")+ggtitle("b")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))

p2=ggplot()+
#geom_boxplot(data=indf,aes(y=comp,x=time,group=paste0(trait,time,rich),color=trait,linetype=as.factor(rich)),width=100)+
stat_smooth(data=indf,aes(y=comp,x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("Competition for mutualistic partners")+
xlab("Time")+ggtitle("c")


blank=ggplot()+ggtitle("a")+
theme_bw()+theme(axis.line = element_blank(),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),axis.title=element_blank(),axis.text=element_blank(),legend.position="none",strip.background=element_rect(fill=NA,color=NA))


leg <- ggpubr::as_ggplot(cowplot::get_legend(p1))
p1=p1+theme(legend.position="none")

grid.arrange(blank,p1,p2,leg,ncol=4,widths=c(1,1,1,0.2))

grid.arrange(blank,p1,p2,ncol=3)
pdf("Fig.1.pdf",width=12,height=4)
grid.arrange(blank,p1,p2,leg,ncol=4,widths=c(1,1,1,0.2))
dev.off();


acp=PCA(indf[,1:6],graph=FALSE)
indf2=cbind(indf,acp$ind$coord)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
empf=fread("empirical_networks.csv")
empf$connectance=empf$connectance+0.1
empf=cbind(empf,predict(acp,newdata=empf[,1:6])$coord)

pl1=ggplot(data=subset(indf2,rich==10),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,group=paste(essai,trait)))+
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
geom_point(data=subset(indf2,time==1500),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,shape=as.factor(rich)))+
geom_point(data=subset(indf2,time==0),aes(x=Dim.1,y=Dim.2,alpha=time,shape=as.factor(rich)),color="grey",alpha=0.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("c")+
geom_point(data=empf,aes(x=Dim.1,y=Dim.2),color="black")

top=plot_grid(pl1,pl2,align="hv",ncol=2,rel_widths=c(1.5,1))

grid.arrange(top,pl3,ncol=1)


pl4=ggplot(data=subset(indf2,time==1500),aes(x=as.factor(rich),y=NODF,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(N[A] == N[P]))+ylab("Nestedness (NODF)")+xlab("Number of species per guild")+
ggtitle("d")


pl5=ggplot(data=subset(indf2,time==1500),aes(x=as.factor(rich),y=connectance,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(N[A] == N[P]))+ylab("Connectance")+xlab("Number of species per guild")+ggtitle("e")


bottom=plot_grid(pl3,pl4,pl5,ncol=3,rel_widths=c(2,1,1))

pdf("Fig.2.pdf",width=7,height=7)
grid.arrange(top,bottom,ncol=1)
dev.off();


pl6=ggplot(data=indf,aes(x=as.factor(time),y=feas,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab("Size of the feasibility domain")+xlab("Time")+facet_wrap(~rich)

pdf("Fig.3.pdf",width=5,height=3)
pl6
dev.off();

indf %>% group_by(rich,trait,time) %>% summarise(mean(feas))

#pdf("Fig.4.pdf",width=11,height=10)
par(mfrow=c(2,2))
for(ti in c(0,1500)){
for(tr in c("morpho","pheno")){
indfb=as.data.frame(subset(indf,trait==tr & time==ti))
indfb$NODFl=logit(indfb$NODF/100,adjust=0.0001)
indfb$connectancel=logit(indfb$connectance,adjust=0.0001)
indfb$modularityl=logit(indfb$modularity,adjust=0.0001)
indfb$feasl=logit(indfb$feas,adjust=0.0001)
indfb$compl=logit(indfb$comp,adjust=0.0001)
modelco=lme(connectance~rich,data=indfb,random=~1|essai,control = lmeControl(opt = "optim"))
modelnest=lme(NODF~rich+connectance,data=indfb,random=~1|essai,control = lmeControl(opt = "optim"))
print(car::vif(modelnest))
modelcomp=lme(modularity~rich+connectance,data=indfb,random=~1|essai,control = lmeControl(opt = "optim"))
print(car::vif(modelcomp))
modelstab=lme(feas~(NODF+connectance+rich+connectance+modularity+compl),data=indfb,random=~1|essai,control = lmeControl(opt = "optim"))
print(car::vif(modelstab))

obj=piecewiseSEM::psem(modelco,modelnest,modelcomp,modelstab,data=indf)
objb=summary(obj)
left=0
right=30
center=(right-left)/2
haut=20
bas=0
cex=1.3
echelle_fleche=3
cex.arrows=0.7
l=objb$coefficients[,c("Predictor","Response","Std.Estimate","P.Value")]

l$colo=ifelse(l$Std.Estimate<0,"firebrick4","dodgerblue4")
l$lty=1
l$lty=ifelse(l$P.Value>0.05,2,l$lty)
l$curv=NA
l$curv[l$Predictor=="rich" & l$Response=="feas"]=2.7
g <- graph.data.frame(l, directed=T)
g= g %>% set_edge_attr("color", value =l$colo)
vertex_attr(g, "name") 
g= g %>% set_vertex_attr("name", value =
c("Diversity",
paste0("Connectance\n(R2 = ",objb$R2$Conditional[objb$R2$Response=="connectance"],")"),
paste0("Nestedness\n(R2 = ",objb$R2$Conditional[objb$R2$Response=="NODF"],")"),
paste0("Modularity\n(R2 = ",objb$R2$Conditional[objb$R2$Response=="modularity"],")"),
paste0("Feasibility\n(R2 = ",objb$R2$Conditional[objb$R2$Response=="feas"],")")))
coord=data.frame(label=vertex_attr(g, "name"),
x=c(center,center,right,left,center),
y=c(haut,bas+11,bas+11,bas+11,bas),vsize=25)

EL=as_edgelist(g)
EL=cbind(EL,l[,3])
asi=abs(l[,3])*echelle_fleche
asi[asi<5]=5
asi[asi>=15]=15

title=paste0(ifelse(tr=="morpho","a - Morphological trait","b - Phenological trait"),", time = ",ti)

qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
border.color="white",label.cex=cex,label.scale=F,lty=l$lty,
edge.label.cex = cex.arrows*2,edge.label.position=0.60,vsize2=9,
vsize=coord$vsize,curve=l$curv,
shape="rectangle",edge.labels=T,fade=F,asize=asi,
mar=c(5,5,5,5),knot.border.color="white",curveShape=-1,title=title)
}

}
#dev.off();