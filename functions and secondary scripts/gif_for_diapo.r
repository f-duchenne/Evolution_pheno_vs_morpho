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
datf=subset(datf,comp==2 & trait=="pheno" & rich==20 & essai==2 & time<=2000)
datf$type2="pheno"
datf$type2[datf$type %in% c("mu2","sd2")]="morpho"
datf$type2[datf$trait %in% c("morpho")]="morpho"
datf$type3="mu"
datf$type3[datf$type %in% c("sd","sd2")]="sd"

time_vec = unique(datf$time)
rich=20
tr="pheno"
a=1
for(ti in time_vec){
bidon=subset(datf,time==ti)
bidon$guild=rep(rep(c("poll","plant"),each=20),2)

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

bidon2=subset(datf,time<=ti & type=="mu")
bidon2$guild=rep(rep(c("poll","plant"),each=20),a)

traits=ggplot(data=bidon2,aes(x=time,y=invlogit1(value),color=guild,group=species,linetype=guild))+geom_line()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=c("chartreuse3","dodgerblue3"))+ylab("Trait value")+
scale_x_continuous(limits=c(0,2000),breaks=c(0,2000))

bidon2=dcast(bidon,species+guild~type)
bidon3=do.call("rbind", replicate(length(0:365), bidon2, simplify = FALSE))
bidon3$x=rep(0:365,each=nrow(bidon2))
bidon3$value=dnorm(bidon3$x,invlogit1(bidon3$mu),invlogit(bidon3$sd))


traits=ggplot(data=bidon3,aes(x=x,y=value,color=guild,group=species,linetype=guild))+geom_line(size=1.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+ylab("density")+xlab("Day of the year")+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous(breaks=c(0,365))


mat=ggplot(data=melt(mbi),aes(x=Var1,y=Var2,fill=value))+geom_tile()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank())+
scale_fill_viridis(option="rocket")+coord_fixed(expand=FALSE,ratio=1)

a=a+1
png(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/gif_pheno_c_",competition,"/",ti,".png"),height=600,width=1200,res=120)
grid.arrange(traits,mat,ncol=2,widths=c(1.5,1))
dev.off();

}



############### FIGURE S2
datf=fread("species_level_data_symmetric.csv")
datf$type2="pheno"
datf$type2[datf$type %in% c("mu2","sd2")]="morpho"
datf$type2[datf$trait %in% c("morpho")]="morpho"
datf$type3="mu"
datf$type3[datf$type %in% c("sd","sd2")]="sd"

endpo=subset(datf,time==2000 & comp==5 & trait!="both")
sp1=ggplot()+
geom_density(data=subset(endpo,type3=="mu"),aes(x=invlogit1(value),color=trait,fill=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo[-1],labels=c("only morpho","only pheno"))+scale_fill_manual(values=colo[-1],labels=c("only morpho","only pheno"))+ylab("density")+
xlab("Mean value of the trait")+ggtitle("a")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_grid(cols=vars(rich),rows=vars(type2), labeller = label_bquote(cols=n[sp] == .(rich)), scales = "free_y")+
labs(colour="Simulations with:",fill="Simulations with:")

sp2=ggplot()+
geom_density(data=subset(endpo,type3=="sd"),aes(x=invlogit(value),color=trait,fill=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo[-1],labels=c("only morpho","only pheno"))+scale_fill_manual(values=colo[-1],labels=c("only morpho","only pheno"))+ylab("density")+
xlab("Standard deviation of the trait")+ggtitle("b")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_grid(cols=vars(rich),rows=vars(type2), labeller = label_bquote(cols=n[sp] == .(rich)), scales = "free_y")+
labs(colour="Simulations with:",fill="Simulations with:")

plot_grid(sp1,sp2,ncol=1,align="hv")

pdf("fig_diapo2.pdf",width=7,height=4)
sp2
dev.off();

###############################
###################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales",
"car","viridis","ggeffects","emmeans") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

colo=c("gold3","chartreuse3","dodgerblue4")
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)

indf=fread("data/simulated/outputs_simulations/networks_info.csv")
names(indf)=gsub("weighted ","w",names(indf))
names(indf)[3]="inter. evenness"
names(indf)[5]="wNODF"
names(indf)[7]="mut.strength"
indf$prop_motif=(indf$motifn/indf$motifntot)

acp=PCA(indf[,1:7],graph=FALSE)
indf2=cbind(indf,acp$ind$coord)
indf2$comp=apply(indf2[,c("comp_a","comp_p")],1,mean)

indf2=subset(indf2,rich==30 & time==2000 & competition==4)

indf2$trait=factor(indf2$trait,levels=c("pheno","both","morpho"))

pl4=ggplot(data=indf2,aes(x=trait,y=wNODF,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(N[A] == N[P]))+ylab("Nestedness")+xlab("Simulation scenario")

pdf(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/","NODF_diapo.pdf"),width=5,height=3)
pl4
dev.off();


pl3=ggplot(data=indf2,aes(x=trait,y=modularity,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(N[A] == N[P]))+ylab("Modularity")+xlab("Simulation scenario")

pdf(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/","modularity_diapo.pdf"),width=5,height=3)
pl3
dev.off();


indf3=subset(indf2,rich==30 & competition==4 & time %in% c(0,2000))
indf3$motifn[indf3$trait=="morpho"]=0
indf3$trait=factor(indf3$trait,levels=c("pheno","both","morpho"))
indf3$trait2=as.character(indf3$trait)
indf3$trait2[indf3$time==0]="initial"
indf3$trait2=factor(indf3$trait2,levels=c("pheno","both","morpho","initial"))

p4=ggplot(data=indf3,aes(y=motifn/motifntot,x=as.factor(time),color=trait2))+
geom_boxplot()+
stat_summary(fun.y=mean, geom="line", size = 0.5,linewidth=1.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"),legend.position="none")+
scale_color_manual(values=c(colo,"grey"))+scale_fill_manual(values=colo)+ylab('Proportion of phenological motifs\npromoting facilitation')+
xlab("Time")+ggtitle("")+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_wrap(~trait,ncol=3)+scale_x_discrete(labels=c("Random","Coevolved"))


pdf(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/","vplus_diapo.pdf"),width=5,height=3)
p4
dev.off();



p3=ggplot()+
stat_smooth(data=indf3,aes(y=mut.strength/(competition*comp),x=time,color=trait,fill=trait,linetype=as.factor(rich)))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+scale_fill_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+ylab("mutualism/competition")+guides(linetype=guide_legend(override.aes=list(fill=NA,color="black")))+
xlab("Time")+ggtitle("c")+labs(colour="Simulations with:",fill="Simulations with:")+scale_x_continuous(breaks=c(0,1000,2000))+
facet_wrap(~competition,ncol=1,labeller = label_bquote(italic(c) == .(competition)))





pl5=ggplot(data=subset(indf2,time==2000 & competition %in% c(4)),aes(x=as.factor(rich),y=modularity,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo,labels=c("both traits","only morpho","only pheno"))+
labs(linetype=expression(N[A] == N[P]))+ylab("Modularity")+xlab(expression(paste(n[sp], " / guild")))+ggtitle("e")