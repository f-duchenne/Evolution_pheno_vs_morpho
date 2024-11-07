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

colo=c("chartreuse3","dodgerblue4","gold3")
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

##################################### SPECIES LEVEL ########################
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
datf=fread("species_level_data_symmetric.csv")
datf=subset(datf,comp==2 & trait=="pheno" & rich==20 & essai==1 & time<=2000)
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

bidon=subset(datf,time<=ti & type=="mu")
bidon$guild=rep(rep(c("poll","plant"),each=20),a)

traits=ggplot(data=bidon,aes(x=time,y=invlogit1(value),color=guild,group=species))+geom_line()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=c("chartreuse3","dodgerblue3"))+ylab("Trait value")+
scale_x_continuous(limits=c(0,2000),breaks=c(0,2000))

mat=ggplot(data=melt(mbi),aes(x=Var1,y=Var2,fill=value))+geom_tile()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank())+
scale_fill_viridis(option="rocket")+coord_cartesian(expand=FALSE)

a=a+1
png(paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/gif/",ti,".png"),height=600,width=1000,res=120)
grid.arrange(traits,mat,ncol=2)
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