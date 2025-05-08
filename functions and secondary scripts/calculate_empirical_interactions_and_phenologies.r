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

weigthed.mean=function(x,y){sum(x*y,na.rm=T)/sum(y,na.rm=T)}

###### LOAD FLOWER COUNT:
flowers=fread("data/empirical/interactions/flower_count.csv")
flowers=subset(flowers,!is.na(Plant_species) & Plant_species!="sp") #remove unidentified species
flowers$date=as.Date(paste(flowers$Day,flowers$Month,flowers$Year,sep="/"),format="%d/%m/%Y")
flowers$jj=yday(flowers$date)
tab=unique(flowers[,c("Site_ID")])

flo=flowers[,c("date","Plant_gen_sp","Flower_abundance","Site_ID")]
names(flo)[2:3]=c("species","count")
flo$guild="plant"

###### LOAD INTERACTIONS:
datf=as.data.frame(fread("data/empirical/interactions/BeeFun Master_inter.csv"))
datf$Site_ID[datf$Site_ID=="Convento_de _la_luz"]="Convento_de_la_luz"
datf$Site_ID[datf$Site_ID=="La_Rocina"]="La_rocina"
datf=subset(datf,!is.na(Pollinator_species) & Pollinator_species!="sp") #remove unidentified species
datf$date=as.Date(paste(datf$Day,datf$Month,datf$Year,sep="/"),format="%d/%m/%Y")
datf$jj=yday(datf$date)

poll=datf[,c("date","Pollinator_gen_sp","Frequency","Site_ID")]
names(poll)[2:3]=c("species","count")
poll$guild="poll"

pheno_tab=rbind(flo,poll)

fwrite(pheno_tab,"data/empirical/interactions/weekly_counts.csv")

abund_counts=fread("data/empirical/interactions/weekly_counts.csv")
abund_counts$jj=yday(as.Date(abund_counts$date,format="%d/%m/%Y")) #create a column 
phen_est=abund_counts %>% group_by(species,guild) %>% summarise(mu=weigthed.mean(jj,count),sde=sqrt(Hmisc::wtd.var(jj,count)))
phen_est=subset(phen_est,!is.na(mu))
phen_est$sde[is.na(phen_est$sde)]=1
phen_est$sde[phen_est$sde<1]=1
phena=subset(phen_est,guild=="poll")
phenp=subset(phen_est,guild=="plant")

#mutualistic interactions
sites=unique(datf$Site_ID)
sites=sites[!(sites %in% c("El_pozo","El_Pozo"))]
weigthed.mean=function(x,y){sum(x*y,na.rm=T)/sum(y,na.rm=T)}
empf=NULL
networks=list()
for (site in sites){
bidon=subset(datf,Site_ID==site)
bidon=subset(bidon, Pollinator_gen_sp %in% phena$species)
bidon=subset(bidon, Plant_gen_sp %in% phenp$species)

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
names(networks)=sites
save(networks, file = "data/empirical/interactions/matrices_empirical_networks.RData")