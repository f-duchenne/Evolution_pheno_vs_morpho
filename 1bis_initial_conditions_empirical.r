###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales","viridis") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")
load("matrices_empirical_networks.RData")

liste=data.frame(site=names(networks),na=sapply(networks,ncol),np=sapply(networks,nrow))

nbrep=10

for(i in 1:nbrep){
for(j in 1:nrow(liste)){
nbsp_a=liste$na[j]
nbsp_p=liste$np[j]

duree=70
#pheno
mu_phen=runif(nbsp_a+nbsp_p,-2.5,2.5)
sd_phen=runif(nbsp_a+nbsp_p,-2.5,2.5)

#morpho
mu_morpho=runif(nbsp_a+nbsp_p,-2.5,2.5)
sd_morpho=runif(nbsp_a+nbsp_p,-2.5,2.5)


Nini <- runif(nbsp_a+nbsp_p,1,1)


#####SPECIES LEVEL INFORMATIONS TO EXPORT

final=data.frame(sp=c(paste0("a",1:nbsp_a),paste0("f",1:nbsp_p)),type=c(rep("poll",nbsp_a),rep("plant",nbsp_p)),mu_phen=mu_phen,
sd_phen=sd_phen,mu_morpho=mu_morpho,sd_morpho=sd_morpho,Nini=Nini,site=liste$site[i])
final$random=i

setwd(dir=paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/initial_empir"))
fwrite(final,paste0("pops_ini_",liste$site[j],"_",i,".csv"))
}
}
fwrite(liste,paste0("liste.csv"))


