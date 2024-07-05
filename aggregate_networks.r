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
datf$Site_ID[datf$Site_ID=="Convento_de _la_luz"]="Convento_de_la_luz"
datf$Site_ID[datf$Site_ID=="La_Rocina"]="La_rocina"
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

fwrite(phena,"poll_pheno_empirical.csv")

sites=unique(datf$Site_ID)
sites=sites[!(sites %in% c("El_pozo","El_Pozo"))]
weigthed.mean=function(x,y){sum(x*y,na.rm=T)/sum(y,na.rm=T)}
empf=NULL
networks=list()
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
m=m/sqrt(matrix(apply(m,1,sum),ncol=1) %*% matrix(apply(m,2,sum),nrow=1))

networks[[site]]=m
}

save(networks, file = "matrices_empirical_networks.RData")