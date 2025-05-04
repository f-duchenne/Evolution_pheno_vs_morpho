pkgs <- c("data.table", "dplyr","bipartite") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)

source("scripts/functions and secondary scripts/toolbox.R")

datf=fread("data/simulated/outputs_simulations/species_level_alleg.csv")
datf$type="mu"
datf$type[grep("mu2",datf$variable)]="mu2"
datf$type[grep("sd",datf$variable)]="sd"
datf$type[grep("sd2",datf$variable)]="sd2"


ncalc=10
indf=NULL
motifs=NULL

comp_vec=2
time_vec=unique(datf$time)

ess=1
competition=6
rich=20
ti=2000

######### NETWORK COEVOLVED WITH MORPHO ONLY
tr="morpho"
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


m2=round(m,digits=3)
obj_morpho=nestedrank(m2, method = "NODF", weighted=TRUE, normalise=TRUE, return.matrix=TRUE)$nested.matrix
obj_morpho=computeModules(m)


######### NETWORK COEVOLVED WITH PHENO ONLY
tr="pheno"
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

m2=round(m,digits=3)
obj_phen=nestedrank(m2, method = "NODF", weighted=TRUE, normalise=TRUE, return.matrix=TRUE)$nested.matrix
ggplot(data=melt(obj_phen),aes(x=as.numeric(Var1),y=as.numeric(Var2),fill=value))+geom_tile()+
scale_fill_gradient(low="white",high="black")+ scale_y_reverse()


obj_phen=computeModules(m)



pl1=as.ggplot(~plotModuleWeb(obj_morpho,displayAlabels = FALSE,displayBlabels = FALSE))+ggtitle("c")+theme(plot.title=element_text(size=14,face="bold",hjust = 0))
pl2=as.ggplot(~plotModuleWeb(obj_phen,displayAlabels = FALSE,displayBlabels = FALSE))+ggtitle("b")+theme(plot.title=element_text(size=14,face="bold",hjust = 0))

grid.arrange(pl1,pl2)



######### NETWORK COEVOLVED WITH PHENO ONLY
tr="both"
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
m2=round(m,digits=2)
df=melt(m2)
df$Var2=paste0("H",df$Var2)
g=graph_from_data_frame(df,directed =FALSE)
V(g)$type=TRUE
V(g)$type[grep("H",V(g)$name)]=FALSE
vec=rep(NA,length(vertex_attr(g,"name")))
vec[]="forestgreen"
vec[grep("H",vertex_attr(g,"name"))]="dodgerblue3"
vertex_attr(g,"color")=vec
edge_attr(g,"weight")=df$value
E(g)$width <- (E(g)$weight^1)*10


edge_attr(g,"color")=sapply(E(g)$weight,function(x){adjustcolor(paste0("gray",round(100-100*x/max(E(g)$weight))),alpha.f =x/max(E(g)$weight))})
size=10
pl1=plot_grid(base2grob(~plot(g,layout = layout_as_bipartite,vertex.label=NA,vertex.size=size,label.cex=0,size2=size,curved=T)))+
theme(plot.margin = unit(c(0,0,0,0), "cm"),plot.title=element_text(size=14,face="bold"))+ggtitle("b")